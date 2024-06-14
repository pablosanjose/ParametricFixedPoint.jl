module ParametricFixedPoint

using FixedPoint
using FixedPointAcceleration

export pfp, FP_solver, FPA_solver

function pfp(f, x0;
    dp = 1.0, norm = x -> maximum(abs, x),
    stepsolver = FP_solver(f; grad_norm = norm, tol = 1e-3),
    finalsolver = FP_solver(f; grad_norm = norm, tol = 1e-8))

    iters = steps = 0
    xs = Ref_or_similar(x0)  # reusable seed
    # norm(x´ - x) using in-place operations if required
    distance(x, x´) = norm(muladd_broadcast!(xs, 1.0, x´, -1.0, x))

    # Check input
    0 < dp <= 1 || error("The parameter step `dp` must be in the interval (0, 1]")

    # Converge seed
    (x0´, err, its) = stepsolver(x0, 0.0)
    iters += its;

    # Special case: single step
    if dp ≈ 1
        (x, err, its) = finalsolver(x0´, 1.0)
        iters += its; steps += 1
        return (; x = x, error = err, steps, iters)
    end

    # Preparation phase
    p0, p1, p2 = 0.0, dp, min(1.0, 2dp)
    (x1, err, its) = stepsolver(x0´, p1)
    iters += its; steps += 1
    muladd_broadcast!(xs, 2.0, x1, -1.0, x0´)  # xs = 2x1-x0´
    (x2, err, its) = stepsolver(maybe_unref(xs), p2)
    iters += its; steps += 1
    if p2 ≈ 1.0     # In case we're done on p2
        xₙ = x2
    else
        # Tracking phase
        dx = 0.5 * (distance(x1, x0´) + distance(x2, x1))
        xₙ₋₂, xₙ₋₁, xₙ = x0´, x1, x2
        pₙ₋₂, pₙ₋₁, pₙ = p0, p1, p2
        dpₙ₋₂, dpₙ₋₁ = p1 - p0, p2 - p1
        ar, br = Ref_or_similar(x0), Ref_or_similar(x0)
        while pₙ < 1.0 && !(pₙ ≈ 1.0)
            k = inv((dpₙ₋₂ + dpₙ₋₁) * dpₙ₋₂ * dpₙ₋₁)
            muladd_broadcast!(ar, k * dpₙ₋₁, xₙ₋₂, - k * (dpₙ₋₂ + dpₙ₋₁), xₙ₋₁, k * dpₙ₋₂, xₙ)
            muladd_broadcast!(br, k * dpₙ₋₁^2, xₙ₋₂, - k * (dpₙ₋₂ + dpₙ₋₁)^2, xₙ₋₁, k * dpₙ₋₂ * (dpₙ₋₂ + 2dpₙ₋₁), xₙ)
            a, b, c = maybe_unref(ar), maybe_unref(br), maybe_unref(xₙ)
            u = dx / norm(b)
            muladd_broadcast!(xs, u, a, 1.0, b)
            dpₙ = dx/norm(xs)
            pₙ₊₁ = min(pₙ + dpₙ, 1.0)
            dpₙ = pₙ₊₁ - pₙ
            muladd_broadcast!(xs, dpₙ^2, a, dpₙ, b, 1.0, c)
            (xₙ₊₁, err, its) = stepsolver(maybe_unref(xs), pₙ₊₁)
            iters += its; steps += 1
            xₙ₋₂, xₙ₋₁, xₙ = xₙ₋₁, xₙ, xₙ₊₁
            pₙ₋₂, pₙ₋₁, pₙ = pₙ₋₁, pₙ, pₙ₊₁
            dpₙ₋₂, dpₙ₋₁ = dpₙ₋₁, dpₙ
        end
    end

    # Final step
    if finalsolver !== stepsolver
        (xₙ, err, its) = finalsolver(maybe_unref(xₙ), 1.0)
        iters += its
    end

    return (; x = xₙ, error = err, steps, iters)
end

Ref_or_similar(x0::Number) = Ref(x0)
Ref_or_similar(x0) = similar(x0)

maybe_unref(x::Ref) = x[]
maybe_unref(x) = x

muladd_broadcast!(xs::Ref, a1, a2, a3, a4) = (xs[] = a1*a2 + a3*a4)
muladd_broadcast!(xs, a1, a2, a3, a4) = (@. xs = a1*a2 + a3*a4; xs)
muladd_broadcast!(xs::Ref, a1, a2, a3, a4, a5, a6) = (xs[] = a1*a2 + a3*a4 + a5*a6)
muladd_broadcast!(xs, a1, a2, a3, a4, a5, a6) = (@. xs = a1*a2 + a3*a4 + a5*a6; xs)

FP_solver(f; kw...) = (xₛ, pₙ) -> begin
    solution = afps(Base.Fix2(f, pₙ), xₛ; kw...)
    return (solution.x, solution.error, solution.iters)
end

# from FixedPointAcceleration.jl
FPA_solver(f; kw...) = (xₛ, pₙ) -> begin
    solution = fixed_point(x -> f.(x, pₙ), xₛ; kw...)
    return (maybe_only(solution.FixedPoint_, xₛ), solution.Convergence_, solution.Iterations_)
end

maybe_only(x, xₛ::Number) = only(x)
maybe_only(x, xₛ) = x

end # module
