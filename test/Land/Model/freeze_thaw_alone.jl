# Test that freeze thaw alone reproduces expected behavior: exponential behavior
# for liquid water content, ice content, and total water conserved
using MPI
using OrderedCollections
using StaticArrays
using Statistics

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
using CLIMAParameters.Planet: ρ_cloud_liq
using CLIMAParameters.Planet: ρ_cloud_ice

using ClimateMachine
using ClimateMachine.Land
using ClimateMachine.Land.SoilWaterParameterizations
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using ClimateMachine.BalanceLaws:
    BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux, vars_state

@testset "Freeze thaw alone" begin
    FT = Float64;
    
    function get_grid_spacing(N_poly::Int64, nelem_vert::Int64, zmax::FT, zmin::FT)
        

        test_m_soil = SoilModel(
            nothing,
            PrescribedWaterModel(FT;),
            PrescribedTemperatureModel(FT;),
        )

        test_m = LandModel(
            param_set,
            test_m_soil;
            source = (),
            init_state_prognostic = init_soil_water!,
        )
        test_config = ClimateMachine.SingleStackConfiguration(
            "LandModel",
            N_poly,
            nelem_vert,
            zmax,
            param_set,
            test_m;
            zmin = zmin,
            numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
        );

        Δ = min_node_distance(test_config.grid)
        return Δ
    end
    
    function init_soil_water!(land, state, aux, coordinates, time)
        FT = eltype(state)
        state.soil.water.ϑ_l = FT(land.soil.water.initialϑ_l(aux))
        state.soil.water.θ_ice = FT(land.soil.water.initialθ_ice(aux))
    end;


    ClimateMachine.init();

    N_poly = 5;
    nelem_vert = 50;
    zmax = FT(0);
    zmin = FT(-1)
    t0 = FT(0)
    timeend = FT(30)
    dt = FT(0.05)

    n_outputs = 30;
    every_x_simulation_time = ceil(Int, timeend / n_outputs);
    Δ = get_grid_spacing(N_poly, nelem_vert, zmax, zmin)
    κ = FT(1.5) # W/m/K
    cs = FT(3e6)
    τft = max(dt, cs * Δ^2.0 / κ)




    soil_param_functions =
        SoilParamFunctions(porosity = 0.75, Ksat = 0.0, S_s = 1e-3, τft = τft)
    bottom_flux = (aux, t) -> FT(0.0)
    surface_flux = (aux, t) -> FT(0.0)
    surface_state = nothing
    bottom_state = nothing
    ϑ_l0 = (aux) -> FT(0.33)
    soil_water_model = SoilWaterModel(
        FT;
        initialϑ_l = ϑ_l0,
        dirichlet_bc = Dirichlet(
            surface_state = surface_state,
            bottom_state = bottom_state,
        ),
        neumann_bc = Neumann(
            surface_flux = surface_flux,
            bottom_flux = bottom_flux,
        ),
    )

    soil_heat_model = PrescribedTemperatureModel(FT;T = (aux,t) -> FT(267.15))# needs to be a function

    m_soil = SoilModel(soil_param_functions, soil_water_model, soil_heat_model)
    sources = (FreezeThaw(),)
    m = LandModel(
        param_set,
        m_soil;
        source = sources,
        init_state_prognostic = init_soil_water!,
    )

    driver_config = ClimateMachine.SingleStackConfiguration(
        "LandModel",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m;
        zmin = zmin,
        numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
    );



    solver_config =
        ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
    mygrid = solver_config.dg.grid;
    Q = solver_config.Q;
    aux = solver_config.dg.state_auxiliary;

    all_data = Dict([k => Dict() for k in 1:n_outputs]...)

    step = [1];
    callback = GenericCallbacks.EveryXSimulationTime(
        every_x_simulation_time,
    ) do (init = false)
        t = ODESolvers.gettime(solver_config.solver)
        grads = SingleStackUtils.get_vars_from_nodal_stack(
            mygrid,
            solver_config.dg.state_gradient_flux,
            vars_state(m, GradientFlux(), FT),
        )
        state_vars = SingleStackUtils.get_vars_from_nodal_stack(
            mygrid,
            Q,
            vars_state(m, Prognostic(), FT),
        )
        aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
            mygrid,
            aux,
            vars_state(m, Auxiliary(), FT),
        )
        all_vars = OrderedDict(state_vars..., aux_vars..., grads...)
        all_vars["t"] = [t]
        all_data[step[1]] = all_vars

        step[1] += 1
        nothing
    end;

    ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

    t = ODESolvers.gettime(solver_config.solver)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        Q,
        vars_state(m, Prognostic(), FT),
    )
    grads = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        solver_config.dg.state_gradient_flux,
        vars_state(m, GradientFlux(), FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        mygrid,
        aux,
        vars_state(m, Auxiliary(), FT),
    )
    all_vars = OrderedDict(state_vars..., aux_vars..., grads...);
    all_vars["t"] = [t]
    all_data[n_outputs] = all_vars

    m_liq = [ρ_cloud_liq(param_set)*mean(all_data[k]["soil.water.ϑ_l"]) for k in 1:n_outputs]
    m_ice = [ρ_cloud_ice(param_set)*mean(all_data[k]["soil.water.θ_ice"]) for k in 1:n_outputs]
    t = [all_data[k]["t"][1] for k in 1:n_outputs]
    total_water = m_ice+m_liq
    m_liq_of_t = m_liq[1]*exp.(-1.0.*(t.-t[1])./τft)
    m_ice_of_t = -m_liq_of_t .+ (m_ice[1]+m_liq[1])

    @test mean(abs.(m_ice+m_liq .- total_water)) < 1e-9
    @test mean(abs.(m_liq .- m_liq_of_t)) < 1e-9
    @test mean(abs.(m_ice .- m_ice_of_t)) < 1e-9
end
