
export BatchedJacobianFreeNewtonKrylovSolver

mutable struct BatchedJacobianFreeNewtonKrylovSolver{ET} <: AbstractNonlinearSolver
    # Tolerance
    ϵ::ET
    # Max number of Newton iterations
    M::Int
    # Maximum number of batched Newton methods (number of vertical columns)
    batch_size::Int
    # Linear solver for the Jacobian system
    linearsolver!
    function BatchedJacobianFreeNewtonKrylovSolver(
        Q::AT,
        linearsolver!,
        batch_size::Int;
        ϵ = √eps(eltype(AT)),
        M = min(20, eltype(Q)),
    ) where {AT}
        return new{typeof(ϵ)}(ϵ, M, batch_size, linearsolver!)
    end
end

function initialize!(
    implicitoperator!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    args...,
)
    return false, zero(eltype(Q))
end

function nonlineariteration!(
    implicitoperator!,
    Q,
    Qrhs,
    solver::BatchedJacobianFreeNewtonKrylovSolver,
    threshold,
    args...,
)
    return nothing
end
