module ParametricFixedPoint

using FixedPoint

function pfp(f, x0;
    dp = 0.05, tol = 1e-8, checkseed = true, grad_norm = x -> maximum(abs, x), kw...)

    iterations = 0
    kws = (; tol, grad_norm, kw...)
    xs = similar(x0)  # seed for reuse
    distance(x1, x2) = (xs .= x2 .- x1; grad_norm(xs))

    # Check input
    0 < dp <= 1 || error("The parameter step `dp` must be in the interval (0, 1]")
    if checkseed
        distance(f(x0, 0), x0) < tol ||
            error("Initial solution is not a fixed point at p = 0 at tolerance $tol")
        iterations = 1
    end

    # Special case: single step
    if dp ≈ 1
        (x, err, its) = afps(x->f(x, 1.0), x0; kws...)
        iterations += its
        return (x, err, iterations)
    end

    # Preparation phase
    p0, p1, p2 = 0.0, dp, max(1.0, 2dp)
    (x1, err, its) = afps(x->f(x, p1), x0; kws...)
    iterations += its
    xs .= 2 .* x1 .- x0     # linear extrapolation
    (x2, err, its) = afps(x->f(x, p2), xs; kws...)
    iterations += its
    # In case we're done on p2
    p2 ≈ 1.0 && return (x2, err, iterations)
    dx = 0.5 * (distance(x1, x0) + distance(x2, x1))

    # Tracking phase
    xₙ₋₂, xₙ₋₁, xₙ = x0, x1, x2
    pₙ₋₂, pₙ₋₁, pₙ = p0, p1, p2
    dpₙ₋₂, dpₙ₋₁ = p1 - p0, p2 - p1
    a, b = similar(x0), similar(x0)
    while pₙ < 1.0 && !(pₙ ≈ 1.0)
        k = inv((dpₙ₋₂ + dpₙ₋₁) * dpₙ₋₂ * dpₙ₋₁)
        @. a = k * (dpₙ₋₁ * xₙ₋₂ - (dpₙ₋₂ + dpₙ₋₁) * xₙ₋₁ + dpₙ₋₂ * xₙ)
        @. b = k * (dpₙ₋₁^2 * xₙ₋₂ - (dpₙ₋₂ + dpₙ₋₁)^2 * xₙ₋₁ + dpₙ₋₂*(dpₙ₋₂ + 2dpₙ₋₁) * xₙ)
        c = xₙ
        u = dx/grad_norm(b)
        @. xs = a * u + b
        dpₙ = dx/grad_norm(xs)
        pₙ₊₁ = max(pₙ + dpₙ, 1.0)
        dpₙ = pₙ₊₁ - pₙ
        @. xs = dpₙ^2 * a + dpₙ * b + c
        (xₙ₊₁, err, its) = afps(x->f(x, pₙ₊₁), xs; kws...)
        iterations += its
        xₙ₋₂, xₙ₋₁, xₙ = xₙ₋₁, xₙ, xₙ₊₁
        pₙ₋₂, pₙ₋₁, pₙ = pₙ₋₁, pₙ, pₙ₊₁
        dpₙ₋₂, dpₙ₋₁ = dpₙ₋₁, dpₙ
    end

    return (xₙ, err, iterations)
end

end # module
