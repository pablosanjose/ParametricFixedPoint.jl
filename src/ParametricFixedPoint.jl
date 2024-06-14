module ParametricFixedPoint

using FixedPoint
using FixedPointAcceleration
using SIAMFANLEquations

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
FPA_solver(f; kw...) = (xₛ, pₙ) -> FPA_solver(f, xₛ, pₙ; kw...)

function FPA_solver(f, xₛ::Number, pₙ; kw...)
    solution = fixed_point(x -> f.(x, pₙ), xₛ; kw...)
    return (only(solution.FixedPoint_), solution.Convergence_, solution.Iterations_)
end

function FPA_solver(f, xₛ::AbstractVector, pₙ; kw...)
    solution = fixed_point(x -> f(x, pₙ), xₛ; kw...)
    return (solution.FixedPoint_, solution.Convergence_, solution.Iterations_)
end

function FPA_solver(f, xₛ::AbstractArray, pₙ; kw...)
    s = size(xₛ)
    solution = fixed_point(x -> vec(f(reshape(x, s), pₙ)), vec(xₛ); kw...)
    return (copy(reshape(solution.FixedPoint_, s)), solution.Convergence_, solution.Iterations_)
end

#### SIAMFANLEquations

# fp(x, p) -> new x
function pfp(fp, x0, p0::Real, m = 3, vstore = zeros(length(x0)+1, 3m+3); kw...)
    x0p = wrap(copy(x0), p0)
    function f!(xp´, xp)
        p = pop!(xp)
        x = unwrap(xp, x0)
        x´ = fp(x, clamp(abs(p), 0.0, 1.0))
        push!(xp, p)
        copyto!(xp´, x´)
        xp´[end] = pcurve(p)
        return xp´
    end
    sol = aasol(f!, x0p, m, vstore; kw...)
    return process_solution(sol)
end

function process_solution(sol)
    if sol.errcode != 0
        println(sol.history)
        error("Couldn't converge. Error $(sol.errcode)")
    end
    p = pop!(sol.solution)
    if !(p ≈ 1.0)
        println(sol.history)
        error("Did not reach p = 1.0, converged to p = $p")
    end
    return (; x = sol.solution, error = last(sol.history), iters = length(sol.history))
end

# pcurve(p) = ifelse(abs(p) >= 1.0, 1.0, abs(p)*(2.0 - abs(p)))
# pcurve(p) = abs(p) >= 1.0 ? 1.0 : sign(p) * sqrt(1.0 - (abs(p)-1.0)^2)
pcurve(p) = abs(p) >= 1.0 ? 1.0 : sqrt(1.0 - (abs(p)-1.0)^2)

wrap(x::Real) = [x]
wrap(x::Complex) = [real(x), imag(x)]
wrap(x::AbstractVector{T}) where {T<:Real} = x
wrap(x::AbstractVector{T}) where {T<:Complex} = reinterpret(real(T), x)
wrap(x::AbstractArray{T}) where {T<:Real} = vec(x)
wrap(x::AbstractArray{T}) where {T<:Complex} = reinterpret(real(T), vec(x))
wrap(x, p) = push!(wrap(x), p)

unwrap(x, ::T) where {T<:Real} = only(x)
unwrap(x, ::T) where {T<:Complex} = x[1] + im*x[2]
unwrap(x, ::T) where {C<:Real,T<:AbstractVector{C}} = x
unwrap(x, ::T) where {C<:Complex,T<:AbstractVector{C}} = reinterpret(C, x)
unwrap(x, m::T) where {C<:Real,T<:AbstractArray{C}} = reshape(x, size(m))
unwrap(x, m::T) where {C<:Complex,T<:AbstractArray{C}} = reinterpret(C, reshape(x, size(m)))

end # module
