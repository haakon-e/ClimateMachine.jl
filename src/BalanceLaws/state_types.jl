#### State types

export AbstractStateType,
    Conservative,
    Auxiliary,
    Gradient,
    GradientFlux,
    GradientLaplacian,
    Hyperdiffusive,
    UpwardIntegrals,
    DownwardIntegrals

abstract type AbstractStateType end
struct Conservative <: AbstractStateType end
struct Auxiliary <: AbstractStateType end
struct Gradient <: AbstractStateType end
struct GradientFlux <: AbstractStateType end
struct GradientLaplacian <: AbstractStateType end
struct Hyperdiffusive <: AbstractStateType end
struct UpwardIntegrals <: AbstractStateType end
struct DownwardIntegrals <: AbstractStateType end
