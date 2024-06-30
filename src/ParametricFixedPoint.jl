module ParametricFixedPoint

export pfp, pfp!, pfp_store

using SIAMFANLEquations: aasol

# Solve fixed point x = f(x, p). f must be given in mutating form, f!(x´, x, p[, pdata])
# where pdata is an optional data object passed as a pdata keyword argument to pfp
pfp(f!, v0, args...; kw...) = pfp!(f!, copy(v0), args...; kw...)

function pfp!(f!, v0, p0::Real = 1;
        order = 2,
        f´ = pcurve, converge_error = false,
        store = pfp_store(v0, p0, order), kw...)
    isextended, vp0, vstore = store
    function fext!(vp´, vp, pdata...)
        p = isextended ? clamp(abs(pop!(vp)), 0.0, 1.0) : 1.0
        f!(v0, vp, p, pdata...)
        copyto!(vp´, v0)
        isextended && (push!(vp, p); vp´[end] = f´(p))
        return vp´
    end
    sol = aasol(fext!, vp0, order, vstore; kw...)
    return process_solution(sol, isextended; converge_error)
end

function pfp_store(v0, p0, m)
    isextended = !(p0 ≈ 1.0)
    vp0 = copy(v0)
    isextended && push!(vp0, p0)
    vstore = zeros(length(vp0), max(3m+3, 2m+4))
    return isextended, vp0, vstore
end

function process_solution(sol, isextended; converge_error = false)
    if converge_error && sol.errcode != 0
        error("Couldn't converge. Error $(sol.errcode), see `?ParametricFixedPoint.aasol` for details.")
    end
    if isextended
        p = pop!(sol.solution)
        if converge_error && !(p ≈ 1.0)
            error("Did not reach p = 1.0, converged to p = $p")
        end
    end
    return (; x = sol.solution, error = last(sol.history), iters = length(sol.history))
end

pcurve(p) = abs(p) >= 1.0 ? sign(p) : sign(p) * sqrt(1.0 - (abs(p)-1.0)^2)
# pcurve(p) = abs(p) >= 1.0 ? sign(p) : p*5/4-p^5*1/4
# pcurve(p) = abs(p) >= 1.0 ? sign(p) : 1.5*p-0.5*p^3

end # module
