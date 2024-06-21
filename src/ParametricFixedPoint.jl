module ParametricFixedPoint

export pfp, pfp!, pfp_store

using SIAMFANLEquations: aasol

# Solve fixed point x = f(x, p)
pfp(f, x0, args...; kw...) = pfp!(f, copy(x0), args...; kw...)

function pfp!(f, x0, p0::Real = 1;
    order = 2,
    f´ = pcurve, converge_error = false,
    store = pfp_store(x0, p0, order), kw...)
    isextended, v0, vstore = store
    len = length(v0)
    function f!(vp´, vp)
        if isextended
            p = clamp(abs(last(vp)), 0.0, 1.0)
            x = unvec!(x0, view(vp, 1:len-1))
        else
            p = 1.0
            x = unvec!(x0, vp)
        end
        x´ = f(x, p)
        copyto!(vp´, x´)
        isextended && (vp´[len] = f´(p))
        return vp´
    end
    sol = aasol(f!, v0, order, vstore)
    return process_solution(sol, x0, isextended; converge_error)
end

function pfp_store(x0, p0, m)
    isextended = p0 ≈ 1.0
    v0 = tovec(x0)
    isextended && push!(v0, p0)
    vstore = zeros(length(v0), max(3m+3, 2m+4))
    return isextended, v0, vstore
end

function process_solution(sol, x0, isextended; converge_error = false)
    if converge_error && sol.errcode != 0
        error("Couldn't converge. Error $(sol.errcode), see `?ParametricFixedPoint.aasol` for details.")
    end
    if isextended
        p = pop!(sol.solution)
        if converge_error && !(p ≈ 1.0)
            error("Did not reach p = 1.0, converged to p = $p")
        end
    end
    return (; x = unvec!(x0, sol.solution), error = last(sol.history), iters = length(sol.history))
end

pcurve(p) = abs(p) >= 1.0 ? 1.0 : sign(p) * sqrt(1.0 - (abs(p)-1.0)^2)

tovec(x::Real) = [x]
tovec(x::Complex) = [real(x), imag(x)]
tovec(x::AbstractVector{T}) where {T<:Real} = copy(x)
tovec(x::AbstractVector{T}) where {T<:Complex} = copy(reinterpret(real(T), x))
tovec(x::AbstractArray{T}) where {T<:Real} = copy(vec(x))
tovec(x::AbstractArray{T}) where {T<:Complex} = copy(reinterpret(real(T), vec(x)))

unvec!(::T, x) where {T<:Real} = only(x)
unvec!(::T, x) where {T<:Complex} = x[1] + im*x[2]
unvec!(x´::T, x) where {C<:Real,T<:AbstractVector{C}} = x === x´ ? x´ : copy!(x´, x)
unvec!(x´::T, x) where {C<:Complex,T<:AbstractVector{C}} = copy!(x´, reinterpret(C, x))
unvec!(x´::T, x) where {C<:Real,T<:AbstractArray{C}} = copy!(x´, reshape(x, size(x´)))
unvec!(x´::T, x) where {C<:Complex,T<:AbstractArray{C}} = copy!(x´, reinterpret(C, reshape(x, size(x´))))

end # module
