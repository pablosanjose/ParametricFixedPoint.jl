# ParametricFixedPoint

[![Build Status](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/ParametricFixedPoint.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Problem

The goal of this package is to solve for $x$ in the parametric fixed-point equation $f(x, p) = x$ starting from a solution $x_0$ of $f(x_0, p_0) = x_0$

While $p$ is a real parameter, $x$ can be anything that supports addition and multiplication by scalars (e.g. $x$ could be a scalar but also an array).

## Interface

The `ParametricFixedPoint` package exports a single function, with the following interface
```julia
gfp(f, x0, (p0, p); ep0 = 0.05, checkseed = false, kw...)
```

Here `kw` are parameters passed to the fixed-point `afps` solver from the `FixedPoint` package, `ep0` is the tracking parameter, denoted $\epsilon_0 \in (0,0.5]$ below, and `checkseed` controls whether `x0` is checked to be a solution of `f(x0, p0) = x0`.

If not specified, as in `gpf(f, x0, p; kw...)`, `p0` is assumed to be zero.

## Algorithm

The general idea of the method is to gradually increase the real parameter $p$ from $p_0$ (where we have a solution $x_0$) to its final value, so that we track the fixed point for $x$ quasi-continuously in $p$.

### Preparation phase

The initial increase of $p$, from $p_0$ to $p_1=p_0+\Delta p_0$, is given by $\Delta p_0 = p_1 - p_0 = \epsilon_0*(p-p_0)$, where $\epsilon_0 \in (0,0.5]$ is a user-provided parameter (typically small).

With this $p_1$ we solve the problem $f(x_1, p_1) = x_1$ by iteration, using $x_0$ as initial guess. This should converge quickly, since $x_0$ is a good guess ($p_1$ is similar to $p_0$).

We do this $p$ increment a second time, using $\Delta p_1 = \Delta p_0$, and obtain the solution $x_2$ of $f(x_2, p_2) = x_2$. In this step we use $2x_1-x_0$ as initial guess, which corresponds to a linear interpolation at $p_2 = p_1 + \Delta p_1 = p_0 + 2\Delta p_0$ from the two previous results.

This concludes the preparation phase, at the end of which we have $x_0$, $x_1$ and $x_2$, the converged solutions for $p_0$, $p_1$ and $p_2$. We also have the average distance between solutions, $|\Delta x| = (|x_2-x_1|+|x_1-x_0|)/2$. Here, $|x|$ is real, and is given by a user-provided definition of distance over the space of solutions, with the only requirement that $|\lambda x| = \lambda|x|$ for any real $\lambda$.

### Tracking phase
We now continue increasing $p_n$ in $\Delta p_n$ steps until we reach the final $p$. The key of this tracking phase is to choose $\Delta p_n$ in an adaptive way, so that the change $|\Delta x_n| = |x_{n+1} - x_{n}|$ is more or less constant, and equal to the $|\Delta x|$ we have obtained.

Assume we have computed the solutions $x_n$ up to $p_n$ and would like to compute the solution at $p_{n+1} = p_n + \Delta p_n$. To choose the appropriate $\Delta p_n$ we use information of the last three solutions $x_{n-2}$, $x_{n-1}$ and $x_{n}$, with which we build a quadratic interpolation
$$x(p) \approx (p-p_n)^2 a + (p-p_n) b + c$$
The constants $a, b, c$ are obtained from the conditions $x(p_n) = x_n$, $x(p_{n-1}) = x_{n-1}$ and $x(p_{n-2}) = x_{n-2}$. This yields
$$\begin{align*}
a &= \frac{\Delta p_{n-1} x_{n-2} -(\Delta p_{n-1}+\Delta p_{n-1})x_{n-1}+\Delta p_{n-2}x_n}{(\Delta p_{n-2}+\Delta p_{n-1})\Delta p_{n-2}\Delta p_{n-1}} \\
b &= \frac{\Delta p_{n-1}x_{n-2} - (\Delta p_{n-2}+\Delta p_{n-1})^2 x_{n-1}+\Delta p_{n-2}(\Delta p_{n-2}+2\Delta p_{n-1}) x_n}{(\Delta p_{n-2}+\Delta p_{n-1})\Delta p_{n-2}\Delta p_{n-1}} \\
c &= x_n
\end{align*}
$$
Using this interpolation, we can estimate $x(p_{n+1}) = \Delta p_n^2 a+\Delta p_n b + c \approx x_{n+1}$. We will define the size of the next parameter step $\Delta p_n$ such that the distance $|x_{n+1}-x_n| = \Delta p_n|a\Delta p_n + b|$ is equal to $|\Delta x|$. Assuming $\Delta p_n$ is small, we obtain
$$\Delta p_n \approx \frac{|\Delta x|}{\left|a\frac{|\Delta x|}{|b|} + b\right|}$$

Given the above step, the solution $x_n$ at $p_{n+1} = p_n + \Delta p_n$ will then be computed by iteration of $f(x, p_{n+1}) = x$ using the interpolation $\Delta p_n^2 a+\Delta p_n b+c$ as initial guess.

### Final step

As soon as $p_{n+1} = p_n+\Delta p_n$ exceeds the final $p$, i.e. when $p_{n+1} >= p$, we reduce $\Delta p_n = p - p_n$ so that the final point $p_n$ matches $p$.

### Iterations

To solve $f(x_n, p_n) = x_n$ at each step we use the `FixedPoint` package, which implements a a simple accelerated iteration algorithm (`afps` function). As long as the initial guess is close to the solution, this should converge relatively fast.
