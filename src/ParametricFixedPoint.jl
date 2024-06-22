module ParametricFixedPoint

export pfp, pfp!, pfp_store

using SIAMFANLEquations: aasol

# Solve fixed point x = f(x, p). f must be given in mutating form, f!(x´, x, p[, pdata])
# where pdata is an optional data object passed as a pdata keyword argument to pfp
pfp(f!, x0, args...; kw...) = pfp!(f!, copy(x0), args...; kw...)

function pfp!(f!, x0, p0::Real = 1;
    order = 2,
    f´ = pcurve, converge_error = false,
    store = pfp_store(x0, p0, order), kw...)
    x´, vstore = store
    len = length(v0)
    function fext!(vp´, vp, pdata...)
        if isextended
            p = clamp(abs(last(vp)), 0.0, 1.0)
            x = copy!(x0, view(vp, 1:len-1))
        else
            p = 1.0
            x = x0
        end
        f!(x´, x, p, pdata...)
        copyto!(vp´, x´)
        isextended && (vp´[len] = f´(p))
        return vp´
    end
    sol = aasol(fext!, v0, order, vstore; kw...)
    return process_solution(sol, x0, isextended; converge_error)
end

function pfp_store(x0, p0, m)
    isextended = p0 ≈ 1.0
    x0´ = similar(x0)
    vstore = zeros(length(v0)+isextended, max(3m+3, 2m+4))
    return x0´, vstore
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
    return (; x = sol.solution, error = last(sol.history), iters = length(sol.history))
end

pcurve(p) = abs(p) >= 1.0 ? 1.0 : sign(p) * sqrt(1.0 - (abs(p)-1.0)^2)

end # module
