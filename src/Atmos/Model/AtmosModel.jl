module Atmos

export AtmosModel,
    AtmosAcousticLinearModel, AtmosAcousticGravityLinearModel, RemainderModel

using CLIMAParameters
using CLIMAParameters.Planet: grav, cp_d
using CLIMAParameters.Atmos.SubgridScale: C_smag
using DocStringExtensions
using LinearAlgebra, StaticArrays
using ..ConfigTypes
using ..Orientations
import ..Orientations:
    vertical_unit_vector,
    altitude,
    latitude,
    longitude,
    projection_normal,
    gravitational_potential,
    ∇gravitational_potential,
    projection_tangential

using ..VariableTemplates
using ..Thermodynamics
using ..TemperatureProfiles
import ..Thermodynamics: internal_energy
using ..MPIStateArrays: MPIStateArray
using ..Mesh.Grids:
    VerticalDirection, HorizontalDirection, min_node_distance, EveryDirection

using ClimateMachine.BalanceLaws:
    BalanceLaw, number_state_conservative, num_integrals

import ClimateMachine.BalanceLaws:
    vars_state,
    vars_state,
    vars_state,
    vars_gradient_laplacian,
    vars_state,
    vars_hyperdiffusive,
    vars_integrals,
    vars_reverse_integrals,
    flux_first_order!,
    flux_second_order!,
    source!,
    wavespeed,
    boundary_state!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    transform_post_gradient_laplacian!,
    init_state_auxiliary!,
    init_state_conservative!,
    update_auxiliary_state!,
    indefinite_stack_integral!,
    reverse_indefinite_stack_integral!,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!,
    reverse_integral_load_auxiliary_state!,
    reverse_integral_set_auxiliary_state!

import ClimateMachine.DGMethods:
    LocalGeometry,
    lengthscale,
    resolutionmetric,
    DGModel,
    nodal_update_auxiliary_state!
import ..DGMethods.NumericalFluxes:
    boundary_state!,
    boundary_flux_second_order!,
    normal_boundary_flux_second_order!,
    NumericalFluxFirstOrder,
    NumericalFluxGradient,
    NumericalFluxSecondOrder,
    CentralNumericalFluxHigherOrder,
    CentralNumericalFluxDivergence

import ..Courant: advective_courant, nondiffusive_courant, diffusive_courant

"""
    AtmosModel <: BalanceLaw

A `BalanceLaw` for atmosphere modeling.
Users may over-ride prescribed default values for each field.

# Usage

    AtmosModel(
        param_set,
        orientation,
        ref_state,
        turbulence,
        hyperdiffusion,
        moisture,
        radiation,
        source,
        tracers,
        boundarycondition,
        init_state_conservative
    )


# Fields
$(DocStringExtensions.FIELDS)
"""
struct AtmosModel{FT, PS, O, RS, T, HD, M, P, R, S, TR, BC, IS, DC} <:
       BalanceLaw
    "Parameter Set (type to dispatch on, e.g., planet parameters. See CLIMAParameters.jl package)"
    param_set::PS
    "An orientation model"
    orientation::O
    "Reference State (For initial conditions, or for linearisation when using implicit solvers)"
    ref_state::RS
    "Turbulence Closure (Equations for dynamics of under-resolved turbulent flows)"
    turbulence::T
    "Hyperdiffusion Model (Equations for dynamics of high-order spatial wave attenuation)"
    hyperdiffusion::HD
    "Moisture Model (Equations for dynamics of moist variables)"
    moisture::M
    "Precipitation Model (Equations for dynamics of precipitating species)"
    precipitation::P
    "Radiation Model (Equations for radiative fluxes)"
    radiation::R
    "Source Terms (Problem specific source terms)"
    source::S
    "Tracer Terms (Equations for dynamics of active and passive tracers)"
    tracers::TR
    "Boundary condition specification"
    boundarycondition::BC
    "Initial Condition (Function to assign initial values of state variables)"
    init_state_conservative::IS
    "Data Configuration (Helper field for experiment configuration)"
    data_config::DC
end


"""
    function AtmosModel{FT}()
Constructor for `AtmosModel` (where `AtmosModel <: BalanceLaw`)

"""
function AtmosModel{FT}(
    ::Type{AtmosLESConfigType},
    param_set::AbstractParameterSet;
    orientation::O = FlatOrientation(),
    ref_state::RS = HydrostaticState(DecayingTemperatureProfile{FT}(param_set),),
    turbulence::T = SmagorinskyLilly{FT}(0.21),
    hyperdiffusion::HD = NoHyperDiffusion(),
    moisture::M = EquilMoist{FT}(),
    precipitation::P = NoPrecipitation(),
    radiation::R = NoRadiation(),
    source::S = (Gravity(), Coriolis(), GeostrophicForcing{FT}(7.62e-5, 0, 0)),
    tracers::TR = NoTracers(),
    boundarycondition::BC = AtmosBC(),
    init_state_conservative::IS = nothing,
    data_config::DC = nothing,
) where {FT <: AbstractFloat, O, RS, T, HD, M, P, R, S, TR, BC, IS, DC}
    @assert param_set ≠ nothing
    @assert init_state_conservative ≠ nothing

    atmos = (
        param_set,
        orientation,
        ref_state,
        turbulence,
        hyperdiffusion,
        moisture,
        precipitation,
        radiation,
        source,
        tracers,
        boundarycondition,
        init_state_conservative,
        data_config,
    )

    return AtmosModel{FT, typeof.(atmos)...}(atmos...)
end
function AtmosModel{FT}(
    ::Type{AtmosGCMConfigType},
    param_set::AbstractParameterSet;
    orientation::O = SphericalOrientation(),
    ref_state::RS = HydrostaticState(DecayingTemperatureProfile{FT}(param_set),),
    turbulence::T = SmagorinskyLilly{FT}(C_smag(param_set)),
    hyperdiffusion::HD = NoHyperDiffusion(),
    moisture::M = EquilMoist{FT}(),
    precipitation::P = NoPrecipitation(),
    radiation::R = NoRadiation(),
    source::S = (Gravity(), Coriolis()),
    tracers::TR = NoTracers(),
    boundarycondition::BC = AtmosBC(),
    init_state_conservative::IS = nothing,
    data_config::DC = nothing,
) where {FT <: AbstractFloat, O, RS, T, HD, M, P, R, S, TR, BC, IS, DC}
    @assert param_set ≠ nothing
    @assert init_state_conservative ≠ nothing
    atmos = (
        param_set,
        orientation,
        ref_state,
        turbulence,
        hyperdiffusion,
        moisture,
        precipitation,
        radiation,
        source,
        tracers,
        boundarycondition,
        init_state_conservative,
        data_config,
    )

    return AtmosModel{FT, typeof.(atmos)...}(atmos...)
end


"""
    vars_state(m::AtmosModel, ::Conservative, FT)
Conserved state variables (Prognostic Variables)
"""
function vars_state(m::AtmosModel, ::Conservative, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        turbulence::vars_state(m.turbulence, Conservative(), FT)
        hyperdiffusion::vars_state(m.hyperdiffusion, Conservative(), FT)
        moisture::vars_state(m.moisture, Conservative(), FT)
        radiation::vars_state(m.radiation, Conservative(), FT)
        tracers::vars_state(m.tracers, Conservative(), FT)
    end
end

"""
    vars_state(m::AtmosModel, ::Gradient, FT)
Pre-transform gradient variables
"""
function vars_state(m::AtmosModel, ::Gradient, FT)
    @vars begin
        u::SVector{3, FT}
        h_tot::FT
        turbulence::vars_state(m.turbulence, Gradient(), FT)
        hyperdiffusion::vars_state(m.hyperdiffusion, Gradient(), FT)
        moisture::vars_state(m.moisture, Gradient(), FT)
        tracers::vars_state(m.tracers, Gradient(), FT)
    end
end
"""
    vars_state(m::AtmosModel, ::GradientFlux, FT)
Post-transform gradient variables
"""
function vars_state(m::AtmosModel, ::GradientFlux, FT)
    @vars begin
        ∇h_tot::SVector{3, FT}
        turbulence::vars_state(m.turbulence, GradientFlux(), FT)
        hyperdiffusion::vars_state(m.hyperdiffusion, GradientFlux(), FT)
        moisture::vars_state(m.moisture, GradientFlux(), FT)
        tracers::vars_state(m.tracers, GradientFlux(), FT)
    end
end

"""
    vars_state(m::AtmosModel, ::GradientLaplacian, FT)
Pre-transform hyperdiffusive variables
"""
function vars_state(m::AtmosModel, ::GradientLaplacian, FT)
    @vars begin
        hyperdiffusion::vars_state(m.hyperdiffusion, GradientLaplacian(), FT)
    end
end

"""
    vars_state(m::AtmosModel, ::Hyperdiffusive, FT)
Post-transform hyperdiffusive variables
"""
function vars_state(m::AtmosModel, ::Hyperdiffusive, FT)
    @vars begin
        hyperdiffusion::vars_state(m.hyperdiffusion, Hyperdiffusive(), FT)
    end
end

"""
    vars_state(m::AtmosModel, ::Auxiliary, FT)
Auxiliary variables, such as vertical (stack)
integrals, coordinates, orientation information,
reference states, subcomponent auxiliary vars,
debug variables
"""
function vars_state(m::AtmosModel, ::Auxiliary, FT)
    @vars begin
        ∫dz::vars_state(m, UpwardIntegrals(), FT)
        ∫dnz::vars_state(m, DownwardIntegrals(), FT)
        coord::SVector{3, FT}
        orientation::vars_state(m.orientation, Auxiliary(), FT)
        ref_state::vars_state(m.ref_state, Auxiliary(), FT)
        turbulence::vars_state(m.turbulence, Auxiliary(), FT)
        hyperdiffusion::vars_state(m.hyperdiffusion, Auxiliary(), FT)
        moisture::vars_state(m.moisture, Auxiliary(), FT)
        tracers::vars_state(m.tracers, Auxiliary(), FT)
        radiation::vars_state(m.radiation, Auxiliary(), FT)
    end
end
"""
    vars_state(m::AtmosModel, ::UpwardIntegrals, FT)
"""
function vars_state(m::AtmosModel, ::UpwardIntegrals, FT)
    @vars begin
        radiation::vars_state(m.radiation, UpwardIntegrals(), FT)
    end
end
"""
    vars_state(m::AtmosModel, ::DownwardIntegrals, FT)
"""
function vars_state(m::AtmosModel, ::DownwardIntegrals, FT)
    @vars begin
        radiation::vars_state(m.radiation, DownwardIntegrals(), FT)
    end
end

####
#### Forward orientation methods
####
projection_normal(bl, aux, u⃗) =
    projection_normal(bl.orientation, bl.param_set, aux, u⃗)
projection_tangential(bl, aux, u⃗) =
    projection_tangential(bl.orientation, bl.param_set, aux, u⃗)
latitude(bl, aux) = latitude(bl.orientation, aux)
longitude(bl, aux) = longitude(bl.orientation, aux)
altitude(bl, aux) = altitude(bl.orientation, bl.param_set, aux)
vertical_unit_vector(bl, aux) =
    vertical_unit_vector(bl.orientation, bl.param_set, aux)
gravitational_potential(bl, aux) = gravitational_potential(bl.orientation, aux)
∇gravitational_potential(bl, aux) =
    ∇gravitational_potential(bl.orientation, aux)

include("ref_state.jl")
include("turbulence.jl")
include("hyperdiffusion.jl")
include("moisture.jl")
include("precipitation.jl")
include("radiation.jl")
include("source.jl")
include("tracers.jl")
include("boundaryconditions.jl")
include("linear.jl")
include("remainder.jl")
include("courant.jl")
include("filters.jl")

@doc """
    flux_first_order!(
        m::AtmosModel,
        flux::Grad,
        state::Vars,
        aux::Vars,
        t::Real
    )

Computes and assembles non-diffusive fluxes in the model
equations.
""" flux_first_order!
@inline function flux_first_order!(
    m::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρinv = 1 / ρ
    ρu = state.ρu
    u = ρinv * ρu

    # advective terms
    flux.ρ = ρ * u
    flux.ρu = ρ * u .* u'
    flux.ρe = u * state.ρe

    # pressure terms
    p = pressure(m, m.moisture, state, aux)
    if m.ref_state isa HydrostaticState
        flux.ρu += (p - aux.ref_state.p) * I
    else
        flux.ρu += p * I
    end
    flux.ρe += u * p
    flux_radiation!(m.radiation, m, flux, state, aux, t)
    flux_moisture!(m.moisture, m, flux, state, aux, t)
    flux_tracers!(m.tracers, m, flux, state, aux, t)
end

function compute_gradient_argument!(
    atmos::AtmosModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρinv = 1 / state.ρ
    transform.u = ρinv * state.ρu
    transform.h_tot = total_specific_enthalpy(atmos, atmos.moisture, state, aux)

    compute_gradient_argument!(atmos.moisture, transform, state, aux, t)
    compute_gradient_argument!(atmos.turbulence, transform, state, aux, t)
    compute_gradient_argument!(atmos.hyperdiffusion, transform, state, aux, t)
    compute_gradient_argument!(atmos.tracers, transform, state, aux, t)
end

function compute_gradient_flux!(
    atmos::AtmosModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.∇h_tot = ∇transform.h_tot

    # diffusion terms required for SGS turbulence computations
    compute_gradient_flux!(
        atmos.turbulence,
        atmos.orientation,
        diffusive,
        ∇transform,
        state,
        aux,
        t,
    )
    # diffusivity of moisture components
    compute_gradient_flux!(atmos.moisture, diffusive, ∇transform, state, aux, t)
    compute_gradient_flux!(atmos.tracers, diffusive, ∇transform, state, aux, t)
end

function transform_post_gradient_laplacian!(
    atmos::AtmosModel,
    hyperdiffusive::Vars,
    hypertransform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    transform_post_gradient_laplacian!(
        atmos.hyperdiffusion,
        hyperdiffusive,
        hypertransform,
        state,
        aux,
        t,
    )
end

@doc """
    flux_second_order!(
        atmos::AtmosModel,
        flux::Grad,
        state::Vars,
        diffusive::Vars,
        hyperdiffusive::Vars,
        aux::Vars,
        t::Real
    )
Diffusive fluxes in AtmosModel. Viscosity, diffusivity are calculated
in the turbulence subcomponent and accessed within the diffusive flux
function. Contributions from subcomponents are then assembled (pointwise).
""" flux_second_order!
@inline function flux_second_order!(
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    ν, D_t, τ = turbulence_tensors(atmos, state, diffusive, aux, t)
    d_h_tot = -D_t .* diffusive.∇h_tot
    flux_second_order!(atmos, flux, state, τ, d_h_tot)
    flux_second_order!(atmos.moisture, flux, state, diffusive, aux, t, D_t)
    flux_second_order!(
        atmos.hyperdiffusion,
        flux,
        state,
        diffusive,
        hyperdiffusive,
        aux,
        t,
    )
    flux_second_order!(atmos.tracers, flux, state, diffusive, aux, t, D_t)
end

#TODO: Consider whether to not pass ρ and ρu (not state), foc BCs reasons
@inline function flux_second_order!(
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    τ,
    d_h_tot,
)
    flux.ρu += τ * state.ρ
    flux.ρe += τ * state.ρu
    flux.ρe += d_h_tot * state.ρ
end

@inline function wavespeed(m::AtmosModel, nM, state::Vars, aux::Vars, t::Real)
    ρinv = 1 / state.ρ
    u = ρinv * state.ρu
    uN = abs(dot(nM, u))
    ss = soundspeed(m, m.moisture, state, aux)

    FT = typeof(state.ρ)
    ws = fill(uN + ss, MVector{number_state_conservative(m, FT), FT})
    vars_ws = Vars{vars_state(m, Conservative(), FT)}(ws)

    wavespeed_tracers!(m.tracers, vars_ws, nM, state, aux, t)

    return ws
end


function update_auxiliary_state!(
    dg::DGModel,
    m::AtmosModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    FT = eltype(Q)
    state_auxiliary = dg.state_auxiliary

    if num_integrals(m, FT) > 0
        indefinite_stack_integral!(dg, m, Q, state_auxiliary, t, elems)
        reverse_indefinite_stack_integral!(dg, m, Q, state_auxiliary, t, elems)
    end

    nodal_update_auxiliary_state!(
        atmos_nodal_update_auxiliary_state!,
        dg,
        m,
        Q,
        t,
        elems,
    )

    return true
end

function atmos_nodal_update_auxiliary_state!(
    m::AtmosModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    atmos_nodal_update_auxiliary_state!(m.moisture, m, state, aux, t)
    atmos_nodal_update_auxiliary_state!(m.radiation, m, state, aux, t)
    atmos_nodal_update_auxiliary_state!(m.turbulence, m, state, aux, t)
    atmos_nodal_update_auxiliary_state!(m.tracers, m, state, aux, t)
end

function integral_load_auxiliary_state!(
    m::AtmosModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
)
    integral_load_auxiliary_state!(m.radiation, integ, state, aux)
end

function integral_set_auxiliary_state!(m::AtmosModel, aux::Vars, integ::Vars)
    integral_set_auxiliary_state!(m.radiation, aux, integ)
end

function reverse_integral_load_auxiliary_state!(
    m::AtmosModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
)
    reverse_integral_load_auxiliary_state!(m.radiation, integ, state, aux)
end

function reverse_integral_set_auxiliary_state!(
    m::AtmosModel,
    aux::Vars,
    integ::Vars,
)
    reverse_integral_set_auxiliary_state!(m.radiation, aux, integ)
end


@doc """
    init_state_auxiliary!(
        m::AtmosModel,
        aux::Vars,
        geom::LocalGeometry
        )
Initialise auxiliary variables for each AtmosModel subcomponent.
Store Cartesian coordinate information in `aux.coord`.
""" init_state_auxiliary!
function init_state_auxiliary!(m::AtmosModel, aux::Vars, geom::LocalGeometry)
    aux.coord = geom.coord
    init_aux!(m.orientation, m.param_set, aux)
    atmos_init_aux!(m.ref_state, m, aux, geom)
    atmos_init_aux!(m.turbulence, m, aux, geom)
    atmos_init_aux!(m.hyperdiffusion, m, aux, geom)
    atmos_init_aux!(m.tracers, m, aux, geom)
end

@doc """
    source!(
        m::AtmosModel,
        source::Vars,
        state::Vars,
        diffusive::Vars,
        aux::Vars,
        t::Real,
        direction::Direction
    )
Computes (and assembles) source terms `S(Y)` in:
```
∂Y
-- = - ∇ • F + S(Y)
∂t
```
""" source!
function source!(
    m::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    atmos_source!(m.source, m, source, state, diffusive, aux, t, direction)
end

@doc """
    init_state_conservative!(
        m::AtmosModel,
        state::Vars,
        aux::Vars,
        coords,
        t,
        args...)
Initialise state variables.
`args...` provides an option to include configuration data
(current use cases include problem constants, spline-interpolants)
""" init_state_conservative!
function init_state_conservative!(
    m::AtmosModel,
    state::Vars,
    aux::Vars,
    coords,
    t,
    args...,
)
    m.init_state_conservative(m, state, aux, coords, t, args...)
end
end # module
