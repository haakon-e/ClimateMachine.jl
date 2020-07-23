#### State types

export AbstractStateType,
    Prognostic,
    Auxiliary,
    Gradient,
    GradientFlux,
    GradientLaplacian,
    Hyperdiffusive,
    UpDownIntegrals,
    UpwardIntegrals,
    DownwardIntegrals

abstract type AbstractStateType end

"""
    Prognostic <: AbstractStateType

Prognostic variables in the PDE system,
which are specified by the `BalanceLaw`, and
solved for by the ODE solver.
"""
struct Prognostic <: AbstractStateType end

"""
    Auxiliary <: AbstractStateType

Auxiliary variables help serve several purposes:

 - Pre-compute and store "expensive" variables,
   for example, quantities computed in vertical
   integrals.
 - Diagnostic exports
"""
struct Auxiliary <: AbstractStateType end

"""
    Gradient <: AbstractStateType

Variables whose gradients must be computed.
"""
struct Gradient <: AbstractStateType end

"""
    GradientFlux <: AbstractStateType

Flux variables, which are functions of gradients.
"""
struct GradientFlux <: AbstractStateType end

"""
    GradientLaplacian <: AbstractStateType

Gradient-Laplacian variables.
"""
struct GradientLaplacian <: AbstractStateType end

"""
    Hyperdiffusive <: AbstractStateType

Hyper-diffusive variables
"""
struct Hyperdiffusive <: AbstractStateType end

"""
    UpDownIntegrals <: AbstractStateType

Super-type for UpwardIntegrals and DownwardIntegrals
"""
abstract type UpDownIntegrals <: AbstractStateType end

"""
    UpwardIntegrals <: UpDownIntegrals

Variables computed in upward integrals
"""
struct UpwardIntegrals <: UpDownIntegrals end

"""
    DownwardIntegrals <: UpDownIntegrals

Variables computed in downward integrals
"""
struct DownwardIntegrals <: UpDownIntegrals end
