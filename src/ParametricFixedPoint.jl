module ParametricFixedPoint

using FixedPoint
using Base: Fix2

export pfp

function pfp(f, x0;
    dp = 1.0, tol = 1e-8, checkseed = false, grad_norm = x -> maximum(abs, x), steptol = tol, kw...)

    iters = steps = 0
    kws = (; tol = steptol, grad_norm, kw...)
    xs = Ref_or_similar(x0)  # seed for reuse
    distance(x, x´) = grad_norm(muladd_broadcast!(xs, 1.0, x´, -1.0, x))

    # Check input
    0 < dp <= 1 || error("The parameter step `dp` must be in the interval (0, 1]")
    if checkseed
        distance(f(x0, 0), x0) < tol ||
            error("Initial solution is not a fixed point at p = 0 at tolerance $tol")
        iters = 1
    end

    # Special case: single step
    if dp ≈ 1
        (x, err, its) = fixedpoint(f, x0, 1.0; kws..., tol)
        iters += its; steps += 1
        return (; x = x, error = err, steps, iters)
    end

    # Preparation phase
    p0, p1, p2 = 0.0, dp, min(1.0, 2dp)
    (x1, err, its) = fixedpoint(f, x0, p1; kws...)
    iters += its; steps += 1
    muladd_broadcast!(xs, 2.0, x1, -1.0, x0)  # xs = 2x1-x0
    (x2, err, its) = fixedpoint(f, xs, p2; kws...)
    iters += its; steps += 1
    if p2 ≈ 1.0     # In case we're done on p2
        xₙ = x2
    else
        # Tracking phase
        dx = 0.5 * (distance(x1, x0) + distance(x2, x1))
        xₙ₋₂, xₙ₋₁, xₙ = x0, x1, x2
        pₙ₋₂, pₙ₋₁, pₙ = p0, p1, p2
        dpₙ₋₂, dpₙ₋₁ = p1 - p0, p2 - p1
        ar, br = Ref_or_similar(x0), Ref_or_similar(x0)
        while pₙ < 1.0 && !(pₙ ≈ 1.0)
            k = inv((dpₙ₋₂ + dpₙ₋₁) * dpₙ₋₂ * dpₙ₋₁)
            muladd_broadcast!(ar, k * dpₙ₋₁, xₙ₋₂, - k * (dpₙ₋₂ + dpₙ₋₁), xₙ₋₁, k * dpₙ₋₂, xₙ)
            muladd_broadcast!(br, k * dpₙ₋₁^2, xₙ₋₂, - k * (dpₙ₋₂ + dpₙ₋₁)^2, xₙ₋₁, k * dpₙ₋₂ * (dpₙ₋₂ + 2dpₙ₋₁), xₙ)
            a, b, c = maybe_unref(ar), maybe_unref(br), maybe_unref(xₙ)
            u = dx / grad_norm(b)
            muladd_broadcast!(xs, u, a, 1.0, b)
            dpₙ = dx/grad_norm(xs)
            pₙ₊₁ = min(pₙ + dpₙ, 1.0)
            dpₙ = pₙ₊₁ - pₙ
            muladd_broadcast!(xs, dpₙ^2, a, dpₙ, b, 1.0, c)
            (xₙ₊₁, err, its) = fixedpoint(f, xs, pₙ₊₁; kws...)
            iters += its; steps += 1
            xₙ₋₂, xₙ₋₁, xₙ = xₙ₋₁, xₙ, xₙ₊₁
            pₙ₋₂, pₙ₋₁, pₙ = pₙ₋₁, pₙ, pₙ₊₁
            dpₙ₋₂, dpₙ₋₁ = dpₙ₋₁, dpₙ
        end
    end

    # Final step
    if steptol > tol
        (xₙ, err, its) = fixedpoint(f, xₙ, 1.0; kws..., tol)
        iters += its
    end

    return (; x = xₙ, error = err, steps, iters)
end

distance!(xs, x, x´, grad_norm) = grad_norm(muladd_broadcast!(xs, 1.0, x´, -1.0, x))

Ref_or_similar(x0::Number) = Ref(x0)
Ref_or_similar(x0) = similar(x0)

maybe_unref(x::Ref) = x[]
maybe_unref(x) = x

fixedpoint(f, x0, p; kw...) = afps(Fix2(f, p), maybe_unref(x0); kw...)

muladd_broadcast!(xs::Ref, a1, a2, a3, a4) = (xs[] = a1*a2 + a3*a4)
muladd_broadcast!(xs, a1, a2, a3, a4) = (@. xs = a1*a2 + a3*a4; xs)
muladd_broadcast!(xs::Ref, a1, a2, a3, a4, a5, a6) = (xs[] = a1*a2 + a3*a4 + a5*a6)
muladd_broadcast!(xs, a1, a2, a3, a4, a5, a6) = (@. xs = a1*a2 + a3*a4 + a5*a6; xs)


end # module
