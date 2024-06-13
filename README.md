# ParametricFixedPoint

[![Build Status](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Problem

The goal of this package is to solve for $x$ in the parametric fixed-point equation $f(x, 1) = x$ starting from a solution $x_0$ of $f(x_0, 0) = x_0$

Here, $f(x, p)$ is a general function that returns an object of the same type as x. $p$ is a real parameter in the interval $[0,1]$, but $x$ can be anything that supports addition and multiplication by scalars (e.g. $x$ could be a scalar but also an array). The problem assumes that a solution $x_0$ can be found easily at $p=0$, while the solution of interest at $p=1$ is difficult to obtain.

## Interface

The `ParametricFixedPoint` package exports a single function, with the following interface
```julia
pfp(f, x0; dp = 0.05, tol = 1e-8, checkseed = true, grad_norm = x -> maximum(abs, x), kw...)
```
This function returns a tuple of the form `(x, error, iterations)`, where `error` is the numerical error of solution `x` and `iterations` is the total number of evaluations of the function `f` performed.

### Keywords
Here `dp` is the initial parameter step, denoted $\Delta p_0 \in (0,1]$ in the algorithm description. `tol` is the absolute tolerance used to determine convergence of solutions. `checkseed` controls whether `x0` is checked to actually be a solution of `f(x0, p0) = x0`, i.e. whether `|f(x0,p0) - x0|<tol`. `grad_norm` is the function used to compute `|f(x,p) - x|`. It defines a distance in the space of $x$, with all the usual properties of a distance (e.g. `grad_norm(λ*x) = abs(λ)*grad_norm(x)>0` for a real `λ`). Finally, `kw` includes additional parameters that are passed to the fixed-point `afps` solver from the `FixedPoint.jl` package to solve each step.

## Algorithm

The general idea of the method is to gradually increase the real parameter $p$ from $p=0$ (where we have a solution $x_0$) to its final value $p=1$, so that we track the fixed point for $x$ quasi-continuously in $p$. The $p_n$ steps are chosen adaptively, depending on the rate of change of the solution $x$ along the way.

### Preparation phase

As a first step, parameter $p$ is increase from $p_0=0$ to $p_1=p_0+\Delta p_0$, where $\Delta p_0 \in (0,1]$ is a user-provided parameter (typically small).

With this $p_1$ we solve the problem $f(x_1, p_1) = x_1$ by iteration, using $x_0$ as initial guess. This should converge quickly, since $x_0$ is a good guess ($x_1$ is similar to $x_0$, assuming $f$ is well-behaved and $\epsilon_0$ is small).

We do this $p$ increment a second time, from $p_1$ to $p_2=p_1 + \Delta p_1$, where we use $\Delta p_1 = \Delta p_0$. We compute the solution $x_2$ of $f(x_2, p_2) = x_2$. In this step we use $2x_1-x_0$ as initial guess, which corresponds to a linear interpolation at $p_2 = p_1 + \Delta p_1 = p_0 + 2\Delta p_0$ from the two previous results.

This concludes the preparation phase, at the end of which we have $x_0$, $x_1$ and $x_2$, the converged solutions for $p_0$, $p_1$ and $p_2$. We also have the average distance between solutions, $|\Delta x| = (|x_2-x_1|+|x_1-x_0|)/2$. Here, $|x|$ is real, and is given by a user-provided definition of `norm`, a distance over the space of solutions.

### Tracking phase
We now continue increasing $p_n$ in $\Delta p_n$ steps until we reach the final $p=1$. The key of this tracking phase is to choose $\Delta p_n$ in an adaptive way, so that the change $|\Delta x_n| = |x_{n+1} - x_{n}|$ is more or less constant, and equal to the $|\Delta x|$ we have obtained in the preparation phase.

Assume we have computed the solutions $x_n$ up to $p_n$ and would like to compute the solution at $p_{n+1} = p_n + \Delta p_n$. To choose the appropriate $\Delta p_n$ we use information of the last three solutions $x_{n-2}$, $x_{n-1}$ and $x_{n}$, with which we build a quadratic interpolation
$$x(p) \approx (p-p_n)^2 a + (p-p_n) b + c$$
The constants $a, b, c$ are obtained from the conditions $x(p_n) = x_n$, $x(p_{n-1}) = x_{n-1}$ and $x(p_{n-2}) = x_{n-2}$. This yields

$$\begin{align*}
a &= \frac{\Delta p_{n-1} x_{n-2} -(\Delta p_{n-2}+\Delta p_{n-1})x_{n-1}+\Delta p_{n-2}x_n}{(\Delta p_{n-2}+\Delta p_{n-1})\Delta p_{n-2}\Delta p_{n-1}} \\
b &= \frac{\Delta p_{n-1}^2x_{n-2} - (\Delta p_{n-2}+\Delta p_{n-1})^2 x_{n-1}+\Delta p_{n-2}(\Delta p_{n-2}+2\Delta p_{n-1}) x_n}{(\Delta p_{n-2}+\Delta p_{n-1})\Delta p_{n-2}\Delta p_{n-1}} \\
c &= x_n
\end{align*}
$$

Using this interpolation, we can estimate $x(p_{n+1}) = \Delta p_n^2 a+\Delta p_n b + c \approx x_{n+1}$. We will define the size of the next parameter step $\Delta p_n$ such that the distance $|x_{n+1}-x_n| = \Delta p_n|a\Delta p_n + b|$ is equal to $|\Delta x|$. Assuming $\Delta p_n$ is small, we obtain
$$\Delta p_n \approx \frac{|\Delta x|}{\left|a\frac{|\Delta x|}{|b|} + b\right|}$$

Given the above $\Delta p_n$, the solution $x_n$ at $p_{n+1} = p_n + \Delta p_n$ will then be computed by iteration of $f(x, p_{n+1}) = x$ using the interpolation $\Delta p_n^2 a+\Delta p_n b+c$ as initial guess.

As soon as $p_{n+1} = p_n+\Delta p_n$ exceeds the final value $p=1$, i.e. when $p_{n+1} >= 1$, we reduce $\Delta p_n = 1 - p_n$ so that the final point is $p_n=1$. After this, we conclude.

## Dependencies

To solve $f(x_n, p_n) = x_n$ at each step of the tracking algorithm, we use the `FixedPoint.jl` package as a lightweight external dependency. `FixedPoint.jl` implements a a simple accelerated iteration algorithm (`afps` function). As long as the initial guess is close to the solution, it should converge relatively fast. See the `FixedPoint.jl` [homepage](https://github.com/francescoalemanno/FixedPoint.jl) for details on supported keywords.
