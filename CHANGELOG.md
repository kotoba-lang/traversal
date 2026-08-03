# Changelog

## 0.2.1 — 2026-08-03

`tiling-benefit` gains a held-out validation, on one side of its threshold.

A C -O2 blocked ikj matmul at n=768 runs at 0.517 ns per multiply-add — above
the 0.27 ns threshold the model reports — so the model predicts blocking buys
nothing. Measured: the best tile managed 1.062x, six percent, inside this
harness's spread. Predicted 1.00x. Held.

Two honest limits recorded with it. The loop constant came from the unblocked
arm, so the model reproducing that arm's absolute time (234 ms predicted
against 234.15 measured) is circular; only the speedup between arms was really
predicted. And the *other* side of the threshold is untested: a kernel below
0.27 ns/madd should show tiling mattering a great deal, and clang -O2 only got
to 0.517, so confirming it needs a hand-vectorised kernel nobody has written.


## 0.2.0 — 2026-08-03

`tiling-benefit` — is there anything for a tile to win?

`tile-plan` answers how big a tile *fits*. That is a different question from
whether blocking buys anything, and the two came apart when `tile-plan`'s
answer was measured: every block width from 8 to 768 landed within 1.2x of
every other on a JVM matmul.

`tiling-benefit` predicts that. Given the measured inner loop (3.77 ns per
multiply-add) and bandwidth (30 GB/s) it reports a 1.00x achievable speedup and
`:worth-tiling? false`, because both arms are loop-bound and blocking only
moves the memory term. It also reports `:loop-ns-threshold` — 0.27 ns per
multiply-add, roughly one cycle — which is how fast the kernel would have to be
for a tiling sweep to be worth running.

With a vectorised kernel at 0.05 ns/op the same call returns 5.3x, not the 23x
traffic ratio, because blocking makes the arm loop-bound again. Blocking helps
until the loop becomes the limit and the model says where that is.

19 tests, 1857 assertions.


## 0.1.0 — 2026-08-03

Initial implementation.

- Linearizations: row-major, column-major, Morton (2-D and 3-D), Hilbert
  (2-D). Round-trip and bijection tested; extents a curve cannot cover are
  refused rather than silently returning positions with holes.
- `neighbour-locality` — measures the space-filling-curve claim instead of
  repeating it, and finds that the mean neighbour distance does NOT favour the
  curves (Morton ties row-major, Hilbert is worse). The win is block
  co-residency, 80% against 50%.
- `tile-plan` — largest square tile fitting a named cache level, using the
  PRIVATE share of that level, a stated occupancy fraction, and rounding down
  to whole lines. Carries its model, including what it does not model.
- `loop-order` — derives the textbook `ikj` matmul order from strides, with
  cost saturating at one line's worth of elements.
- Bit budgets enforced (2-D 15 bits, 3-D 10) so JVM and cljs agree, rather
  than a cljs interleave silently going negative.

16 tests, 1846 assertions.
