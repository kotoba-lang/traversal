# Changelog

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
