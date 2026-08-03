# Changelog

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
