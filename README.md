# KahanSummation.jl

[![Travis](https://travis-ci.org/JuliaMath/KahanSummation.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/KahanSummation.jl)
[![Coveralls](https://coveralls.io/repos/github/JuliaMath/KahanSummation.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaMath/KahanSummation.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaMath.github.io/KahanSummation.jl/latest)

This package provides variants of `sum` and `cumsum`, called `sum_kbn` and `cumsum_kbn`
respectively, using the Kahan-Babuska-Neumaier (KBN) algorithm for additional precision.
These functions are typically slower and less memory efficient than `sum` and `cumsum`.

These functions were formerly part of Julia's Base library.
