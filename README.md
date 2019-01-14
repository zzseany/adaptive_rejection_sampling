# Adaptive Rejection Sampling
R Package: Adaptive Rejection Sampling (`ars`)

The aim of `ars` is to quickly generate sample with hard-to-evaluate density function. The algorithms of adaptive rejection sampling is introduced with detail in Gilks et al. (1992). The method we implemented is tangent approach instead of secant approach. The package is called 'ars' and you can install and use in R.

## Performance of Sampling

We have tested several log-concave density function and 'ars' performed well for all of those distribution. The [write up](https://github.com/zzseany/adaptive_rejection_sampling/blob/master/writeup/ars.pdf) shows how we implement the algorithm and what improvements we made.

## Contributors

This package is written by three people: Jonathan Morrell, Ziyang Zhou, and Vincent Zijlmans. If you found any bug, please accept our apology and email to me so we can fix it.
