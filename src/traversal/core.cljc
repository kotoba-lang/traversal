(ns traversal.core
  "In what order to visit the elements, given a real cache.

  `layout` decides where bytes sit; this decides the sequence they are
  touched in. The two are independent — a perfect struct-of-arrays walked
  down the wrong axis still misses every line — and the second is usually the
  larger effect, because changing an order changes how much data crosses the
  slow boundary rather than how it is arranged once it is across.

  Three families, and they answer different questions:

  - **linearizations** (`:row-major`, `:column-major`, `:morton`,
    `:hilbert`) — how an n-dimensional index becomes one address. Space-
    filling curves are supposed to keep 2-D neighbours near each other in
    1-D. `neighbour-locality` measures that instead of repeating it, and the
    measurement is worth reading: on a 16×16 grid Morton's mean neighbour
    distance is *identical* to row-major's and Hilbert's is *worse*. The win
    is real but it is in block co-residency (80% against 50%), not in the
    mean — so the intuitive metric answers the wrong question.
  - **tiling** (`tile-plan`) — how large a block may be before it stops
    fitting in a named level of cache. The capacity is asked of `machine`,
    divided by the cores that share it, and multiplied by a stated occupancy
    fraction, because a tile that fills the whole cache on paper thrashes
    against everything else running.
  - **loop order** (`loop-order`) — which loop belongs innermost. This is the
    one that produces the textbook `ikj` matmul from first principles rather
    than from memory.

  Every permutation here is checked to BE a permutation (`bijection?`). A
  reordering that silently skips or repeats an element is not a slow program,
  it is a wrong one.

  Pure `.cljc`. Depends only on `kotoba-lang/machine`."
  (:require [machine.core :as m]))

(def format-id :kotoba.traversal/v1)

(def orders #{:row-major :column-major :morton :hilbert})

;; Bit budgets. cljs bitwise operators work on signed 32-bit values, so an
;; interleave that runs past bit 30 would go negative there and not on the
;; JVM. Rather than branch per runtime, the coordinate range is capped so the
;; result always fits in a positive int32 on both.
(def max-morton-2d-bits 15)   ; side <= 32768, result < 2^30
(def max-morton-3d-bits 10)   ; side <= 1024,  result < 2^30

;; ── Morton (Z-order) ─────────────────────────────────────────────────────

(defn- spread-2 [v]
  ;; Insert a zero between each of the low 15 bits.
  (let [v (bit-and v 0x7fff)
        v (bit-and (bit-or v (bit-shift-left v 8)) 0x00ff00ff)
        v (bit-and (bit-or v (bit-shift-left v 4)) 0x0f0f0f0f)
        v (bit-and (bit-or v (bit-shift-left v 2)) 0x33333333)
        v (bit-and (bit-or v (bit-shift-left v 1)) 0x55555555)]
    v))

(defn- compact-2 [v]
  (let [v (bit-and v 0x55555555)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 1)) 0x33333333)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 2)) 0x0f0f0f0f)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 4)) 0x00ff00ff)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 8)) 0x0000ffff)]
    v))

(defn- check-bits! [vs bits what]
  (let [limit (bit-shift-left 1 bits)]
    (doseq [v vs]
      (when-not (and (integer? v) (<= 0 v) (< v limit))
        (throw (ex-info (str what " coordinate out of range")
                        {:phase :traversal/range :value v :max-exclusive limit
                         :reason "cljs bitwise ops are signed 32-bit; the range is capped so
                                  both runtimes agree"}))))))

(defn morton-encode-2d
  "Interleave the bits of `x` and `y`, x in the even positions."
  [x y]
  (check-bits! [x y] max-morton-2d-bits "morton-2d")
  (bit-or (spread-2 x) (bit-shift-left (spread-2 y) 1)))

(defn morton-decode-2d [d]
  [(compact-2 d) (compact-2 (unsigned-bit-shift-right d 1))])

(defn- spread-3 [v]
  ;; Insert two zeros between each of the low 10 bits.
  (let [v (bit-and v 0x3ff)
        v (bit-and (bit-or v (bit-shift-left v 16)) 0x030000ff)
        v (bit-and (bit-or v (bit-shift-left v 8)) 0x0300f00f)
        v (bit-and (bit-or v (bit-shift-left v 4)) 0x030c30c3)
        v (bit-and (bit-or v (bit-shift-left v 2)) 0x09249249)]
    v))

(defn- compact-3 [v]
  (let [v (bit-and v 0x09249249)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 2)) 0x030c30c3)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 4)) 0x0300f00f)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 8)) 0x030000ff)
        v (bit-and (bit-or v (unsigned-bit-shift-right v 16)) 0x000003ff)]
    v))

(defn morton-encode-3d [x y z]
  (check-bits! [x y z] max-morton-3d-bits "morton-3d")
  (bit-or (spread-3 x)
          (bit-shift-left (spread-3 y) 1)
          (bit-shift-left (spread-3 z) 2)))

(defn morton-decode-3d [d]
  [(compact-3 d)
   (compact-3 (unsigned-bit-shift-right d 1))
   (compact-3 (unsigned-bit-shift-right d 2))])

;; ── Hilbert (2-D, power-of-two side) ─────────────────────────────────────
;;
;; The classic xy2d/d2xy pair. Unlike Morton, consecutive positions on a
;; Hilbert curve are always 4-adjacent cells — no long jumps at quadrant
;; boundaries — which is why it scores better on `neighbour-locality` and why
;; it is the usual choice for spatial index keys.

(defn- pow2? [n] (and (integer? n) (pos? n) (zero? (bit-and n (dec n)))))

(defn- check-side! [side]
  (when-not (pow2? side)
    (throw (ex-info "hilbert side must be a power of two"
                    {:phase :traversal/hilbert :side side})))
  (when (> side (bit-shift-left 1 max-morton-2d-bits))
    (throw (ex-info "hilbert side out of range"
                    {:phase :traversal/hilbert :side side
                     :max (bit-shift-left 1 max-morton-2d-bits)})))
  side)

(defn- hilbert-rot [n [x y] rx ry]
  (if (zero? ry)
    (let [[x y] (if (= 1 rx) [(- n 1 x) (- n 1 y)] [x y])]
      [y x])
    [x y]))

(defn hilbert-encode-2d
  "Position of `(x, y)` along a Hilbert curve of side `side` (power of two)."
  [side x y]
  (check-side! side)
  (check-bits! [x y] max-morton-2d-bits "hilbert-2d")
  (loop [s (quot side 2) x x y y d 0]
    (if (zero? s)
      d
      (let [rx (if (pos? (bit-and x s)) 1 0)
            ry (if (pos? (bit-and y s)) 1 0)
            d (+ d (* s s (bit-xor (* 3 rx) ry)))
            [x y] (hilbert-rot side [x y] rx ry)]
        (recur (quot s 2) x y d)))))

(defn hilbert-decode-2d [side d]
  (check-side! side)
  (loop [s 1 t d x 0 y 0]
    (if (>= s side)
      [x y]
      (let [rx (bit-and 1 (quot t 2))
            ry (bit-and 1 (bit-xor t rx))
            [x y] (hilbert-rot s [x y] rx ry)
            x (+ x (* s rx))
            y (+ y (* s ry))]
        (recur (* s 2) (quot t 4) x y)))))

;; ── generic linearization ────────────────────────────────────────────────

(defn linearize
  "`[x y] -> position` under `order`, for a 2-D extent `[w h]`.

  `:morton` and `:hilbert` need a square power-of-two extent; a curve does
  not cover a ragged box without leaving holes, and quietly returning
  positions with holes in them is how a \"reordering\" starts skipping data."
  [order [w h] [x y]]
  (case order
    :row-major    (+ x (* y w))
    :column-major (+ y (* x h))
    :morton       (do (when-not (= w h)
                        (throw (ex-info "morton needs a square extent"
                                        {:phase :traversal/linearize :extent [w h]})))
                      (morton-encode-2d x y))
    :hilbert      (do (when-not (and (= w h) (pow2? w))
                        (throw (ex-info "hilbert needs a square power-of-two extent"
                                        {:phase :traversal/linearize :extent [w h]})))
                      (hilbert-encode-2d w x y))))

(defn visit-sequence
  "Every coordinate of a 2-D extent, in the order `order` visits them."
  [order [w h]]
  (->> (for [y (range h) x (range w)] [x y])
       (sort-by #(linearize order [w h] %))
       vec))

(defn bijection?
  "Does `order` visit every cell of the extent exactly once?

  Checked, not assumed. A reordering that silently skips or repeats an
  element is not a slow program, it is a wrong one — and both Morton and
  Hilbert do exactly that on an extent they were not given (ragged, or not a
  power of two), which is why `linearize` refuses those extents outright."
  [order [w h]]
  (let [n (* w h)
        positions (for [y (range h) x (range w)] (linearize order [w h] [x y]))]
    (= (set positions) (set (range n)))))

;; ── measuring the locality claim ─────────────────────────────────────────

(defn neighbour-locality
  "How well an order keeps spatial neighbours together, measured three ways.

  For every 4-adjacent pair of cells: how far apart does this order put them
  (`:mean`, `:worst`), and — the one that matters — do they land in the same
  block of `block` consecutive positions (`:same-block`)?

  **Measure this rather than repeating the folklore.** Running it on a 16×16
  grid says something people do not expect: Morton's mean neighbour distance
  is *exactly* row-major's (8.5), and Hilbert's is *worse* (9.92). Space-
  filling curves do not win on the mean, because the mean is dominated by the
  few very long jumps and the curves trade many medium jumps for a few huge
  ones.

  Where they do win is `:same-block`, which is the question a cache line or a
  page actually asks: with `block` 16, both curves put 80% of neighbouring
  pairs in one block against row-major's 50%, because a 4×4 tile of the grid
  is 16 consecutive, block-aligned positions on either curve while row-major
  spreads it across four separate rows.

  So: use `:same-block` to choose an order, and keep `:mean` around as the
  reminder that the intuitive metric answers the wrong question."
  ([order extent] (neighbour-locality order extent nil))
  ([order [w h] block]
   (let [pos (fn [x y] (linearize order [w h] [x y]))
         pairs (vec (concat (for [y (range h) x (range (dec w))]
                              [(pos x y) (pos (inc x) y)])
                            (for [y (range (dec h)) x (range w)]
                              [(pos x y) (pos x (inc y))])))
         deltas (mapv (fn [[a b]] (Math/abs (- a b))) pairs)
         n (count pairs)]
     (cond-> {:order order
              :extent [w h]
              :pairs n
              :mean (if (pos? n) (double (/ (reduce + deltas) n)) 0.0)
              :worst (if (pos? n) (apply max deltas) 0)}
       block
       (assoc :block block
              :same-block
              (if (pos? n)
                (double (/ (count (filter (fn [[a b]] (= (quot a block) (quot b block)))
                                          pairs))
                           n))
                0.0))))))

;; ── tiling ───────────────────────────────────────────────────────────────

(def tile-model
  {:model/id :kotoba.traversal.tile/square-working-set-v1
   :model/rule "arrays * tile^2 * element-bytes <= capacity * occupancy"
   :model/assumes
   ["a square tile touched by `arrays` operands, the classic blocked-matmul shape"
    "capacity is the PRIVATE share of the level: total / cores sharing it"
    "occupancy < 1 because the tile shares the level with stack, code and neighbours"]
   :model/does-not-model [:associativity-conflicts :tlb-reach :prefetch :nuca-latency]})

(defn tile-plan
  "The largest square tile whose working set still fits a named cache level.

  `:arrays` is how many operands the inner block touches at once — 3 for
  `C += A*B`. `:occupancy` defaults to 0.5: a tile sized to the full capacity
  fits on paper and thrashes in practice, because the level also holds the
  stack, the code and whatever the sibling core is doing.

  The tile is rounded DOWN to a whole number of cache lines, so a tile row
  never straddles, and clamped to the extent."
  [machine {:keys [extents element-bytes arrays level kind occupancy]
            :or {arrays 3 level 2 kind :unified occupancy 0.5}}]
  (m/require-fact machine [:cpu :cache] "a cache hierarchy")
  (let [capacity (or (m/private-cache-bytes machine level kind)
                     (throw (ex-info "machine declares no cache at that level"
                                     {:phase :traversal/tile :level level :kind kind})))
        line (m/line-bytes machine)
        usable (long (* capacity occupancy))
        raw (long (Math/floor (Math/sqrt (double (/ usable (* arrays element-bytes))))))
        per-line (max 1 (quot line element-bytes))
        rounded (max per-line (* per-line (quot raw per-line)))
        tile (max 1 (apply min rounded extents))
        working-set (* arrays tile tile element-bytes)]
    {:format format-id
     :tile/size tile
     :tile/arrays arrays
     :tile/level level
     :tile/kind kind
     :tile/capacity-bytes capacity
     :tile/occupancy occupancy
     :tile/usable-bytes usable
     :tile/working-set-bytes working-set
     :tile/fits? (<= working-set usable)
     :tile/elements-per-line per-line
     :tile/tiles (mapv #(quot (+ % (dec tile)) tile) extents)
     :tile/model tile-model
     :tile/machine (:machine/id machine)}))

(defn tile-origins
  "Origins of every tile covering `extents`, in row-major tile order."
  [[w h] tile]
  (vec (for [ty (range 0 h tile) tx (range 0 w tile)] [tx ty])))

;; ── loop order ───────────────────────────────────────────────────────────

(defn- access-cost
  "Cost of a loop's stride for one array, in lines.

  Stride 0 is free — the reference is invariant in that loop and hoists.
  Stride 1 is the unit-stride ideal. Beyond one line's worth of elements the
  cost stops growing: you were already paying a whole line per access, and a
  larger stride cannot make you pay more than that."
  [stride elements-per-line]
  (cond
    (zero? stride) 0
    :else (min stride elements-per-line)))

(defn loop-order
  "Which loop belongs innermost, derived rather than remembered.

  `accesses` is `[{:array :A :strides {:i 512 :k 1 :j 0}} …]` — the element
  stride each loop induces in each array. Loops are ordered by total cost
  descending, so the cheapest ends up innermost.

  For `C[i][j] += A[i][k] * B[k][j]` at N=512 this returns `[:i :k :j]` — the
  textbook `ikj` — because `j` is the only loop under which no array strides
  badly."
  [machine loops accesses element-bytes]
  (m/require-fact machine [:cpu :cache] "a cache hierarchy")
  (let [per-line (max 1 (quot (m/line-bytes machine) element-bytes))
        scored (mapv (fn [l]
                       {:loop l
                        :cost (reduce + 0 (map #(access-cost (get-in % [:strides l] 0) per-line)
                                               accesses))
                        :strides (into {} (map (juxt :array #(get-in % [:strides l] 0)) accesses))})
                     loops)
        ;; Worst outermost, cheapest innermost. Name breaks ties so the
        ;; result is a total order and two runs agree.
        ranked (vec (sort-by (juxt (comp - :cost) (comp str :loop)) scored))]
    {:format format-id
     :loop/order (mapv :loop ranked)
     :loop/innermost (:loop (last ranked))
     :loop/ranked ranked
     :loop/elements-per-line per-line
     :loop/machine (:machine/id machine)}))
