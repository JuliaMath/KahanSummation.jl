# KahanSummation.jl

This package uses the Kahan-Babuska-Neumaier (KBN) compensated summation algorithm for
computing sums and cumulative sums.
The `sum` and `cumsum` functions in Julia's Base library use
[pairwise summation](https://en.wikipedia.org/wiki/Pairwise_summation), which provides
high accuracy and good performance.
The [KBN algorithm](https://en.wikipedia.org/wiki/Kahan_summation_algorithm) significantly
reduces numerical error at the cost of performance and memory efficiency.

```@docs
KahanSummation.sum_kbn
KahanSummation.cumsum_kbn
```
