(ns traversal.core-test
  (:require [clojure.test :refer [deftest is testing]]
            [machine.core :as m]
            [traversal.core :as t]))

(def mach
  {:format m/format-id
   :machine/id "fixture"
   :machine/provenance :measured
   :machine/source "test fixture"
   :cpu {:arch :x86-64 :cores 8
         :simd {:name :avx2 :width-bits 256}
         :cache [{:level 1 :kind :data :bytes 32768 :line-bytes 64 :ways 8 :shared-by 1}
                 {:level 2 :kind :unified :bytes 262144 :line-bytes 64 :ways 8 :shared-by 1}
                 {:level 3 :kind :unified :bytes 8388608 :line-bytes 64 :ways 16 :shared-by 8}]}
   :page {:base-bytes 4096 :huge [2097152]}})

(deftest morton-round-trips
  (doseq [x (range 0 40 7) y (range 0 40 5)]
    (is (= [x y] (t/morton-decode-2d (t/morton-encode-2d x y)))))
  (doseq [x [0 1 2 511 1023] y [0 3 1023] z [0 7 1023]]
    (is (= [x y z] (t/morton-decode-3d (t/morton-encode-3d x y z)))))
  (testing "the low bits interleave as advertised"
    (is (= 0 (t/morton-encode-2d 0 0)))
    (is (= 1 (t/morton-encode-2d 1 0)))
    (is (= 2 (t/morton-encode-2d 0 1)))
    (is (= 3 (t/morton-encode-2d 1 1)))
    (is (= 4 (t/morton-encode-2d 2 0)))))

(deftest morton-range-is-enforced-not-wrapped
  (testing "past the declared bit budget a cljs interleave would go negative;
            the range is refused instead"
    (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                 (t/morton-encode-2d 32768 0)))
    (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                 (t/morton-encode-3d 1024 0 0)))
    (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                 (t/morton-encode-2d -1 0)))))

(deftest hilbert-round-trips
  (doseq [side [2 4 8 16 32]
          d (range (* side side))]
    (let [[x y] (t/hilbert-decode-2d side d)]
      (is (= d (t/hilbert-encode-2d side x y))
          (str "side " side " d " d)))))

(deftest hilbert-steps-are-always-adjacent
  (testing "this is the property Morton lacks: no jump at a quadrant boundary"
    (doseq [side [4 8 16]]
      (let [cells (mapv #(t/hilbert-decode-2d side %) (range (* side side)))]
        (doseq [[[x1 y1] [x2 y2]] (partition 2 1 cells)]
          (is (= 1 (+ (Math/abs (- x1 x2)) (Math/abs (- y1 y2))))
              (str "side " side ": " [x1 y1] " -> " [x2 y2])))))))

(deftest hilbert-refuses-extents-it-cannot-cover
  (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
               (t/hilbert-encode-2d 6 0 0)))
  (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
               (t/linearize :hilbert [8 4] [0 0])))
  (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
               (t/linearize :morton [8 4] [0 0]))))

(deftest every-order-is-a-permutation
  (testing "a reordering that skips or repeats an element is not slow, it is wrong"
    (doseq [order t/orders]
      (is (t/bijection? order [16 16]) (str order))))
  (testing "the rectangular orders also hold on a ragged extent"
    (is (t/bijection? :row-major [7 5]))
    (is (t/bijection? :column-major [7 5]))))

(deftest visit-sequence-agrees-with-linearize
  (let [seq16 (t/visit-sequence :hilbert [16 16])]
    (is (= 256 (count seq16)))
    (is (= 256 (count (set seq16))))
    (is (= (first seq16) (t/hilbert-decode-2d 16 0)))))

(deftest the-mean-neighbour-distance-does-not-say-what-folklore-claims
  (let [rm (t/neighbour-locality :row-major [16 16])
        mo (t/neighbour-locality :morton [16 16])
        hi (t/neighbour-locality :hilbert [16 16])]
    (testing "row-major puts vertical neighbours a whole row apart"
      (is (= 8.5 (:mean rm)))
      (is (= 16 (:worst rm))))
    (testing "Morton's mean is IDENTICAL to row-major's, and Hilbert's is worse"
      (is (= 8.5 (:mean mo)))
      (is (< (:mean rm) (:mean hi))))
    (testing "because both curves trade many medium jumps for a few huge ones"
      (is (> (:worst mo) (:worst rm)))
      (is (> (:worst hi) (:worst rm))))
    (is (= 480 (:pairs rm) (:pairs mo) (:pairs hi)))))

(deftest curves-win-on-block-co-residency-which-is-what-a-cache-asks
  (let [block 16
        rm (t/neighbour-locality :row-major [16 16] block)
        mo (t/neighbour-locality :morton [16 16] block)
        hi (t/neighbour-locality :hilbert [16 16] block)]
    (testing "a 4x4 tile is 16 consecutive block-aligned positions on a curve;
              row-major spreads it over four separate rows"
      (is (= 0.5 (:same-block rm)))
      (is (= 0.8 (:same-block mo)))
      (is (= 0.8 (:same-block hi))))
    (testing "which is the opposite ranking to the mean — hence measuring it"
      (is (> (:same-block hi) (:same-block rm)))
      (is (> (:mean hi) (:mean rm))))))

(deftest tile-plan-derives-the-blocked-matmul-tile
  (let [p (t/tile-plan mach {:extents [512 512] :element-bytes 8 :arrays 3 :level 2})]
    (testing "half of the 256 KiB private L2, three 8-byte operands"
      (is (= 262144 (:tile/capacity-bytes p)))
      (is (= 131072 (:tile/usable-bytes p)))
      (is (= 72 (:tile/size p)))
      (is (= 124416 (:tile/working-set-bytes p)))
      (is (:tile/fits? p)))
    (testing "rounded down to whole cache lines: 8 doubles per line"
      (is (= 8 (:tile/elements-per-line p)))
      (is (zero? (mod (:tile/size p) 8))))
    (is (= [8 8] (:tile/tiles p)))
    (is (= :kotoba.traversal.tile/square-working-set-v1 (get-in p [:tile/model :model/id])))))

(deftest tile-plan-uses-the-private-share-not-the-shared-total
  (testing "L3 is 8 MiB across 8 cores; a tile sized to the whole thing would
            thrash against the seven siblings"
    (let [p (t/tile-plan mach {:extents [4096 4096] :element-bytes 8 :arrays 3 :level 3})]
      (is (= 1048576 (:tile/capacity-bytes p)))
      (is (:tile/fits? p))
      (testing "and the shared total would have licensed a much larger tile"
        (is (< (:tile/size p)
               (Math/sqrt (/ (* 8388608 0.5) 24))))))))

(deftest tile-plan-clamps-to-the-extent
  (let [p (t/tile-plan mach {:extents [16 16] :element-bytes 8 :arrays 3 :level 2})]
    (is (= 16 (:tile/size p)))
    (is (= [1 1] (:tile/tiles p)))))

(deftest tile-origins-cover-the-extent
  (let [origins (t/tile-origins [10 10] 4)]
    (is (= 9 (count origins)))
    (is (= [0 0] (first origins)))
    (is (= [8 8] (last origins)))))

(deftest loop-order-derives-ikj-from-first-principles
  (let [r (t/loop-order mach [:i :j :k]
                        [{:array :A :strides {:i 512 :k 1 :j 0}}
                         {:array :B :strides {:k 512 :j 1 :i 0}}
                         {:array :C :strides {:i 512 :j 1 :k 0}}]
                        8)]
    (testing "C[i][j] += A[i][k] * B[k][j] at N=512"
      (is (= [:i :k :j] (:loop/order r)))
      (is (= :j (:loop/innermost r))))
    (testing "j wins because it is the only loop no array strides badly under"
      (is (= [16 9 2] (mapv :cost (:loop/ranked r)))))
    (is (= 8 (:loop/elements-per-line r)))))

(deftest stride-cost-saturates-at-one-line
  (testing "past a line's worth of elements you already pay a full line per
            access, so a bigger stride cannot cost more"
    (let [r (t/loop-order mach [:a :b]
                          [{:array :X :strides {:a 8 :b 100000}}]
                          8)]
      (is (= [8 8] (mapv :cost (:loop/ranked r)))))))

(deftest an-unprobed-machine-refuses-to-be-planned-against
  (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
               (t/tile-plan m/unknown {:extents [16 16] :element-bytes 8})))
  (is (thrown? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
               (t/loop-order m/unknown [:i] [{:array :A :strides {:i 1}}] 8))))

;; ── tiling benefit (calibrated on an Apple M1 Max, 2026-08-03) ───────────

(def ^:private m1max-tile
  {:format m/format-id
   :machine/id "Apple M1 Max/performance"
   :machine/provenance :measured
   :machine/source "sysctl -a (Darwin)"
   :cpu {:arch :aarch64 :cores 8
         :cache [{:level 1 :kind :data :bytes 131072 :line-bytes 128 :shared-by 1}
                 {:level 2 :kind :unified :bytes 12582912 :line-bytes 128 :shared-by 4}]}
   :page {:base-bytes 16384 :huge []}})

(deftest a-jvm-matmul-cannot-show-tiling-and-the-model-says-so-first
  (testing "measured inner loop 3.77 ns per multiply-add, bandwidth 30 GB/s"
    (let [r (t/tiling-benefit m1max-tile
                              {:n 768 :tile 48 :element-bytes 8 :arrays 3
                               :loop-ns-per-op 3.77 :bandwidth-bytes-per-ns 30.0})]
      (testing "both arms are loop-bound, so blocking has nothing to move"
        (is (= :loop (get-in r [:blocked :bound-by])))
        (is (= :loop (get-in r [:unblocked :bound-by]))))
      (testing "predicted speedup is 1.00x — which is what the sweep measured"
        (is (= 1.0 (:achievable-speedup r)))
        (is (not (:worth-tiling? r))))
      (testing "and it says how fast the loop would have to get"
        (is (< 0.25 (:loop-ns-threshold r) 0.28))))))

(deftest with-a-fast-enough-kernel-tiling-matters-enormously
  (testing "a vectorised kernel near one cycle per multiply-add"
    (let [r (t/tiling-benefit m1max-tile
                              {:n 768 :tile 48 :element-bytes 8 :arrays 3
                               :loop-ns-per-op 0.05 :bandwidth-bytes-per-ns 30.0})]
      (is (:worth-tiling? r))
      (is (= :memory (get-in r [:unblocked :bound-by])))
      (testing "unblocked re-streams B once per row -- 23x the traffic -- but the
                speedup is only 5.3x, because blocking makes the arm loop-bound
                again. Blocking helps until the loop becomes the limit, and the
                model says where that is rather than promising the traffic ratio."
        (is (< 5.3 (:achievable-speedup r) 5.4))
        (is (= :loop (get-in r [:blocked :bound-by])))))))

(deftest the-threshold-is-where-the-two-terms-cross
  (let [args {:n 512 :tile 64 :element-bytes 8 :arrays 3 :bandwidth-bytes-per-ns 30.0}
        threshold (:loop-ns-threshold (t/tiling-benefit m1max-tile
                                                        (assoc args :loop-ns-per-op 1.0)))]
    (testing "just under it, tiling is worth something; just over it, nothing"
      (is (:worth-tiling? (t/tiling-benefit m1max-tile
                                            (assoc args :loop-ns-per-op (* 0.9 threshold)))))
      (is (not (:worth-tiling? (t/tiling-benefit m1max-tile
                                                 (assoc args :loop-ns-per-op (* 1.1 threshold)))))))))
