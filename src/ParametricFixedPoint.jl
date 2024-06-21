module ParametricFixedPoint

export pfp, pfp!

using SIAMFANLEquations: aasol

# Solve fixed point x = f(x, p)
pfp(f, x0, args...; kw...) = pfp!(f, copy(x0), args...; kw...)

function pfp!(f, x0, p0::Real = 1; order = 2, f´ = pcurve, converge_error = false, kw...)
    isextended = p0 ≈ 1
    m = order
    v0 = tovec(x0)
    isextended && push!(v0, p0)
    len = length(v0)
    vstore = zeros(len, max(3m+3, 2m+4))
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
    sol = aasol(f!, v0, m, vstore; kw...)
    return process_solution(sol, x0, isextended; converge_error)
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


# function fixedpoint!(f!::Function, x0::AbstractVector{T};
#         order = 2,
#         atol = 1e-8,
#         maxiter = 100,
#         dist = dx -> maximum(abs, dx),
#         store = fp_store(T, length(x0), order)) where {T}
#     done = false
#     iters = 0
#     err = zero(T)
#     col = 1
#     nrows, ncols = size(store)
#     G, G´, b = store
#     while !done
#         f!(x, x0)
#         x0 .= x .- x0
#         dx = x0  # x0 is now dx
#         err = dist(dx)
#         iters += 1
#         (iters > maxiter || err < atol) && break
#         copycol!(G, dx, col)
#         copycol!(G´, dx, col)

#         col += mod1(col, ncols)
#     end
#     return (; x, error = err, iters)
# end

# fp_store(T, size...) = zeros(T, size...), zeros(T, size...), zeros(T, size[1])

# copycol!(store, x, col) = copyto!(store, 1+(col-1)*length(x), x, 1, length(x))

end # module
