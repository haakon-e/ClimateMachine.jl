"""
    TurbulenceConvection

Turbulence convection models, for example
the Eddy-Diffusivity Mass-Flux model
"""
module TurbulenceConvection

using ..BalanceLaws: BalanceLaw, AbstractStateType
using ..VariableTemplates: @vars, Vars, Grad

export TurbulenceConvectionModel, NoTurbConv

export init_state_prognostic!,
    init_aux_turbconv!, turbconv_nodal_update_auxiliary_state!

import ..BalanceLaws:
    vars_state,
    eq_tends,
    prognostic_to_primitive!,
    primitive_to_prognostic!,
    precompute,
    prognostic_vars,
    init_state_prognostic!,
    init_state_auxiliary!,
    update_auxiliary_state!,
    boundary_state!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!

using ..MPIStateArrays: MPIStateArray
using ..DGMethods: DGModel, LocalGeometry

abstract type TurbulenceConvectionModel end

"""
    NoTurbConv <: TurbulenceConvectionModel

A "no model" type, which results in kernels that
pass through and do nothing.
"""
struct NoTurbConv <: TurbulenceConvectionModel end

prognostic_vars(::NoTurbConv) = ()

eq_tends(pv, ::NoTurbConv, tt) = ()
precompute(::NoTurbConv, bl, args, ts, tend_type) = NamedTuple()

vars_state(m::TurbulenceConvectionModel, ::AbstractStateType, FT) = @vars()

function init_aux_turbconv!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    aux::Vars,
    geom::LocalGeometry,
)
    return nothing
end

function update_auxiliary_state!(
    dg::DGModel,
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    return nothing
end

function turbconv_nodal_update_auxiliary_state!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    state::Vars,
    aux::Vars,
    t::Real,
)
    return nothing
end

function compute_gradient_argument!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    return nothing
end

function compute_gradient_flux!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    return nothing
end

function integral_load_auxiliary_state!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    integ::Vars,
    state::Vars,
    aux::Vars,
)
    return nothing
end

function integral_set_auxiliary_state!(
    m::TurbulenceConvectionModel,
    bl::BalanceLaw,
    aux::Vars,
    integ::Vars,
)
    return nothing
end

function init_state_prognostic!(
    m::NoTurbConv,
    bl::BalanceLaw,
    state,
    aux,
    localgeo,
    t,
) end

prognostic_to_primitive!(::NoTurbConv, _...) = nothing
primitive_to_prognostic!(::NoTurbConv, _...) = nothing

include("boundary_conditions.jl")
include("source.jl")

end
