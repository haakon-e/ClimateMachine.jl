#### AbstractStateType

abstract type AbstractStateType end
struct Conservative <: AbstractStateType end
struct Auxiliary <: AbstractStateType end
struct Gradient <: AbstractStateType end
struct GradientFlux <: AbstractStateType end
struct GradientLaplacian <: AbstractStateType end
struct Hyperdiffusive <: AbstractStateType end
struct VerticalIntegrals <: AbstractStateType end
struct ReverseVerticalIntegrals <: AbstractStateType end
# struct HorizontalIntegrals <: AbstractStateType end
