using OrderedCollections
using StaticArrays
using Test

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.Grids
using ClimateMachine.MPIStateArrays
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry

import ClimateMachine.DGMethods:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    init_state_auxiliary!,
    init_state_conservative!

struct EmptyBalLaw{FT, PS} <: BalanceLaw
    "Parameters"
    param_set::PS
    "Domain height"
    zmax::FT

end
EmptyBalLaw(param_set, zmax) =
    EmptyBalLaw{typeof(zmax), typeof(param_set)}(param_set, zmax)

vars_state_auxiliary(::EmptyBalLaw, FT) = @vars(x::FT, y::FT, z::FT)
vars_state_conservative(::EmptyBalLaw, FT) = @vars(ρ::FT)
vars_state_gradient(::EmptyBalLaw, FT) = @vars()
vars_state_gradient_flux(::EmptyBalLaw, FT) = @vars()

function empty_nodal_init_state_auxiliary!(
    m::EmptyBalLaw,
    aux::Vars,
    geom::LocalGeometry,
)
    aux.x = geom.coord[1]
    aux.y = geom.coord[2]
    aux.z = geom.coord[3]
end

function init_state_auxiliary!(
    m::EmptyBalLaw,
    state_auxiliary::MPIStateArray,
    grid,
)
    nodal_init_state_auxiliary!(
        m,
        empty_nodal_init_state_auxiliary!,
        state_auxiliary,
        grid,
    )
end

function init_state_conservative!(
    m::EmptyBalLaw,
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
)
    z = aux.z
    x = aux.x
    y = aux.y
    state.ρ = (1 - 4 * (z - m.zmax / 2)^2) * (2 - x - y)
end

function test_hmean(
    grid::DiscontinuousSpectralElementGrid{T, dim, N},
    Q::MPIStateArray,
    vars,
) where {T, dim, N}
    state_vars_avg = get_horizontal_mean(grid, Q, vars)
    target = target_meanprof(grid)
    @test state_vars_avg["ρ"] ≈ target
end

function test_hvar(
    grid::DiscontinuousSpectralElementGrid{T, dim, N},
    Q::MPIStateArray,
    vars,
) where {T, dim, N}
    state_vars_var = get_horizontal_variance(grid, Q, vars)
    target = target_varprof(grid)
    @test state_vars_var["ρ"] ≈ target
end

function target_meanprof(
    grid::DiscontinuousSpectralElementGrid{T, dim, N},
) where {T, dim, N}
    Nq = N + 1
    Ntot = Nq * grid.topology.stacksize
    z = Array(get_z(grid))
    target =
        SVector{Ntot, T}([1.0 - 4.0 * (z_i - z[Ntot] / 2.0)^2 for z_i in z])
    return target
end

function target_varprof(
    grid::DiscontinuousSpectralElementGrid{T, dim, N},
) where {T, dim, N}
    Nq = N + 1
    nvertelem = grid.topology.stacksize
    Ntot = Nq * nvertelem
    z = Array(get_z(grid))
    x = z[1:Nq] * nvertelem
    scaled_var = 0.0
    for i in 1:Nq
        for j in 1:Nq
            scaled_var = scaled_var + (2 - x[i] - x[j]) * (2 - x[i] - x[j])
        end
    end
    target = SVector{Ntot, Float64}([
        (1.0 - 4.0 * (z_i - z[Ntot] / 2.0)^2) *
        (1.0 - 4.0 * (z_i - z[Ntot] / 2.0)^2) *
        (scaled_var / Nq / Nq - 1) for z_i in z
    ])
    return target
end

function main()
    FT = Float64
    ClimateMachine.init()

    m = EmptyBalLaw(param_set, FT(1))

    # Prescribe polynomial order of basis functions in finite elements
    N_poly = 5
    # Specify the number of vertical elements
    nelem_vert = 20
    # Specify the domain height
    zmax = m.zmax
    # Initial and final times
    t0 = 0.0
    timeend = 1.0
    dt = 0.1
    # Establish a `ClimateMachine` single stack configuration
    driver_config = ClimateMachine.SingleStackConfiguration(
        "HstatsTest",
        N_poly,
        nelem_vert,
        zmax,
        param_set,
        m,
    )
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        ode_dt = dt,
    )

    @testset "horizontal_stats" begin
        test_hmean(
            driver_config.grid,
            solver_config.Q,
            vars_state_conservative(m, FT),
        )
        test_hvar(
            driver_config.grid,
            solver_config.Q,
            vars_state_conservative(m, FT),
        )
    end
end
main()