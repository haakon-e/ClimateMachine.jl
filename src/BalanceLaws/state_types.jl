#### AbstractStateType

abstract type AbstractStateType end
struct StateConservative <: AbstractStateType end
struct StateAuxiliary <: AbstractStateType end
struct StateGradient <: AbstractStateType end
struct StateGradientFlux <: AbstractStateType end
struct GradientLaplacian <: AbstractStateType end
struct Hyperdiffusive <: AbstractStateType end
struct VerticalIntegrals <: AbstractStateType end
struct ReverseVerticalIntegrals <: AbstractStateType end
# struct HorizontalIntegrals <: AbstractStateType end

