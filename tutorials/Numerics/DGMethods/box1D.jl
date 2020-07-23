# A box advection test to visualise how different filters work

using MPI
using OrderedCollections
using Plots
using StaticArrays

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.BalanceLaws: BalanceLaw
using ClimateMachine.Mesh.Geometry: LocalGeometry
using ClimateMachine.Mesh.Filters
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils

import ClimateMachine.BalanceLaws:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    source!,
    flux_second_order!,
    flux_first_order!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    update_auxiliary_state!,
    nodal_update_auxiliary_state!,
    init_state_auxiliary!,
    init_state_conservative!,
    boundary_state!

ClimateMachine.init(; disable_gpu = true);
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "docs", "plothelpers.jl"));

Base.@kwdef struct Box1D{FT, _init_q, _amplitude, _velo} <: BalanceLaw
    param_set::AbstractParameterSet = param_set
    init_q::FT = _init_q
    amplitude::FT = _amplitude
    velo::FT = _velo
end

vars_state_auxiliary(::Box1D, FT) = @vars(z_dim::FT);
vars_state_conservative(::Box1D, FT) = @vars(q::FT);
vars_state_gradient(m::Box1D, FT) = @vars()
vars_state_gradient_flux(m::Box1D, FT) = @vars()

function init_state_auxiliary!(m::Box1D, aux::Vars, geom::LocalGeometry)
    aux.z_dim = geom.coord[3]
end;

function init_state_conservative!(
    m::Box1D,
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
)
    if aux.z_dim >= 75  && aux.z_dim <= 125
        state.q = m.init_q + m.amplitude
    else
        state.q = m.init_q
    end
end;

function update_auxiliary_state!(
    dg::DGModel,
    m::Box1D,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    return true
end;

function source!(m::Box1D, _...) end;

@inline function flux_first_order!(
    m::Box1D,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    _...,
)
    FT = eltype(state)
    @inbounds begin
        flux.q = SVector(FT(0), FT(0), state.q * m.velo)
    end
end

@inline function flux_second_order!(
    m::Box1D,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
) end

@inline function flux_second_order!(
    m::Box1D,
    flux::Grad,
    state::Vars,
    τ,
    d_h_tot,
) end

function boundary_state!(
    nf,
    m::Box1D,
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
end;

FT = Float64;

function run_box1D(
    N_poly::Int,
    init_q::FT,
    amplitude::FT,
    velo::FT,
    plot_name::String;
    tmar_filter::Bool = false,
    cutoff_filter::Bool = false,
    exp_filter::Bool = false,
    boyd_filter::Bool = false,
    cutoff_param::Int = 1,
    exp_param_1::Int = 0,
    exp_param_2::Int = 32,
    boyd_param_1::Int = 0,
    boyd_param_2::Int = 32,
)
    N_poly = N_poly;
    nelem = 128;
    zmax = FT(600);

    m = Box1D{FT, init_q, amplitude, velo}();

    driver_config = ClimateMachine.SingleStackConfiguration(
        "Box1D",
        N_poly,
        nelem,
        zmax,
        param_set,
        m,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
        boundary = ((0, 0), (0, 0), (0, 0)),
        periodicity = (true, true, true)
    );

    t0 = FT(0)
    timeend = FT(370)

    Δ = min_node_distance(driver_config.grid, VerticalDirection())
    max_vel = m.velo
    dt = Δ / max_vel

    solver_config =
        ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
    grid = solver_config.dg.grid;
    Q = solver_config.Q;
    aux = solver_config.dg.state_auxiliary;

    output_dir = @__DIR__;
    mkpath(output_dir);
    z_key = "z"
    z_label = "z"
    z = get_z(grid)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        Q,
        vars_state_conservative(m, FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        aux,
        vars_state_auxiliary(m, FT),
    )
    all_vars = OrderedDict(state_vars..., aux_vars...);

    n_outputs = 2
    all_data = Dict[Dict([k => Dict() for k in 0:n_outputs]...),]
    all_data[1] = all_vars # store initial condition at ``t=0``

    # output
    step = [1];
    output_freq = floor(Int, timeend / dt);
    @info(" ", timeend, dt, output_freq)

    cb_output = GenericCallbacks.EveryXSimulationSteps(output_freq) do
        state_vars = SingleStackUtils.get_vars_from_nodal_stack(
            grid,
            Q,
            vars_state_conservative(m, FT),
        )
        aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
            grid,
            aux,
            vars_state_auxiliary(m, FT);
        )
        all_vars = OrderedDict(state_vars..., aux_vars...)
        push!(all_data, all_vars)
        step[1] += 1
        nothing
    end;

    filter_freq = 1
    # tmar filter
    cb_tmar =
        GenericCallbacks.EveryXSimulationSteps(filter_freq) do (init = false)
            Filters.apply!(
                solver_config.Q,
                (:q,),
                solver_config.dg.grid,
                TMARFilter(),
            )
            nothing
        end
    # cutoff filter
    cb_cutoff =
        GenericCallbacks.EveryXSimulationSteps(filter_freq) do (init = false)
            Filters.apply!(
                solver_config.Q,
                (:q,),
                solver_config.dg.grid,
                CutoffFilter(solver_config.dg.grid, cutoff_param),
            )
            nothing
        end
    # exponential filter
    cb_exp =
        GenericCallbacks.EveryXSimulationSteps(filter_freq) do (init = false)
            Filters.apply!(
                solver_config.Q,
                (:q,),
                solver_config.dg.grid,
                ExponentialFilter(solver_config.dg.grid, exp_param_1, exp_param_2)
            )
            nothing
        end
    # exponential filter
    cb_boyd =
        GenericCallbacks.EveryXSimulationSteps(filter_freq) do (init = false)
            Filters.apply!(
                solver_config.Q,
                (:q,),
                solver_config.dg.grid,
                BoydVandevenFilter(solver_config.dg.grid, boyd_param_1, boyd_param_2)
            )
            nothing
        end

    user_cb_arr = [cb_output,]
    if tmar_filter
        push!(user_cb_arr, cb_tmar)
    end
    if cutoff_filter
        push!(user_cb_arr, cb_cutoff)
    end
    if exp_filter
        push!(user_cb_arr, cb_exp)
    end
    if boyd_filter
        push!(user_cb_arr, cb_boyd)
    end
    user_cb = (user_cb_arr...,)

    ClimateMachine.invoke!(solver_config; user_callbacks = (user_cb));

    @show keys(all_data[1])

    export_plot(
        z,
        all_data,
        ("q",),
        joinpath(output_dir, plot_name),
        z_label,
        horiz_layout = true,
    );
end

# run a bunch of box model simulations to see how different filters work
run_box1D(2, 0.0, 1.0, 1.0, "box_1D_2_no_filter.pdf")
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_no_filter.pdf")
run_box1D(8, 0.0, 1.0, 1.0, "box_1D_8_no_filter.pdf")

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar.pdf", tmar_filter = true)

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_cutoff_1.pdf", cutoff_filter = true, cutoff_param = 1)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_cutoff_3.pdf", cutoff_filter = true, cutoff_param = 3)

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_exp_0_32.pdf", exp_filter = true, exp_param_1 = 0, exp_param_2 = 32)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_exp_1_32.pdf", exp_filter = true, exp_param_1 = 1, exp_param_2 = 32)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_exp_1_8.pdf", exp_filter = true, exp_param_1 = 1, exp_param_2 = 8)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_exp_1_4.pdf", exp_filter = true, exp_param_1 = 1, exp_param_2 = 4)

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_boyd_0_32.pdf", boyd_filter = true, boyd_param_1 = 0, boyd_param_2 = 32)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_boyd_1_32.pdf", boyd_filter = true, boyd_param_1 = 1, boyd_param_2 = 32)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_boyd_1_8.pdf", boyd_filter = true, boyd_param_1 = 1, boyd_param_2 = 8)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_boyd_1_4.pdf", boyd_filter = true, boyd_param_1 = 1, boyd_param_2 = 4)

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar_exp_1_8.pdf", exp_filter = true, tmar_filter = true,  exp_param_1 = 1, exp_param_2 = 8)
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar_boyd_1_8.pdf", boyd_filter = true, tmar_filter = true,  boyd_param_1 = 1, boyd_param_2 = 8)
