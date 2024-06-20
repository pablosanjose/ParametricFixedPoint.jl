# ParametricFixedPoint

[![Build Status](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Problem

The goal of this package is to solve for $x$ in the parametric fixed-point equation $f(x, 1) = x$ starting from a solution $x_0$ of $f(x_0, 0) = x_0$

Here, $f(x, p)$ is a general function that returns an object of the same type as x. $p$ is a real parameter in the interval $[0,1]$, but $x$ can be anything that supports addition and multiplication by scalars (e.g. $x$ could be a scalar but also an array). The problem assumes that a solution $x_0$ can be found easily at $p=0$, while the solution of interest at $p=1$ is difficult to obtain.

## Interface

The `ParametricFixedPoint` package exports a single function, with the following interface
```julia
pfp(f, x0, p0 = 1; order = 2, converge_error = false, kw...)
```
This function returns a tuple of the form `(x, error, iterations)`, where `error` is the numerical error of solution `x`,  and `iterations` is the total number of evaluations of the function `f` performed. The initial parameter `p0` is used as a seed for `p` that also gets iterated by the algorithm in a way designed to flow to a stable `p=1` fixed point. The iteration uses Anderson acceleration and employs the `aasol` implementation from SIAMFANLEquations.jl. The `order` keyword corresponds to the `m` argument of `aasol`, and controls how many prior iterations are taken into account to compute the next one. If `converge_error = true` an error will be thrown if the solution is not converged. The rest of keywords `kw` are passed directly to `aasol`. Allowed keywords are described [here](https://ctkelley.github.io/SIAMFANLEquations.jl/dev/functions/aasol/).


## Algorithm

The algorithm used by `pfp` can only find stable fixed points, i.e. points `x` such that `x = f(x, 1)` (fixed point) and `|∂ₓf(x, 1)| < 1` (stable). Unstable fixed points require recasting the problem to be stable, e.g. by using the inverse of `f`.

If `p0 == 1.0` the `pfp` function will search for the fixed point `x = f(x, 1)` using the [Anderson acceleration](https://en.wikipedia.org/wiki/Anderson_acceleration) algorithm. If `p0 < 1.0` we extend the algorithm to allow adaptive increments of `pₙ` along with `xₙ`. Effectively, we transform the original iteration `xₙ₊₁ = f(xₙ, 1)` into an extended version of the form `xₙ₊₁, pₙ₊₁ = f(xₙ, |pₙ|), f´(pₙ)`, where we start from `p₀ = p0`. Here `f´(p)` is chosen so that `p=0` is an unstable fixed point and `p = ±1` are stable fixed points. Specifically `f´(p) = ifelse(abs(p) >= 1, 1.0, sign(p) * sqrt(1.0 - (abs(p)-1.0)^2))` (other choices that retain the desired fixed point structure can be passed with the keyword `f´`). The algorithm is therefore designed to that `p` will flow along `x` using the same Anderson acceleration scheme from the initial value `p0` to `p = ±1`, in adaptive steps that will ensure optimal convergence of `x`.
