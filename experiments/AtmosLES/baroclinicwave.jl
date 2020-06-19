using ClimateMachine
ClimateMachine.cli()

using ClimateMachine.Atmos
using ClimateMachine.ConfigTypes
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.Diagnostics
using ClimateMachine.GenericCallbacks
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.ODESolvers
using ClimateMachine.Thermodynamics
using ClimateMachine.VariableTemplates

using CLIMAParameters
using CLIMAParameters.Planet: e_int_v0, grav, day, cp_d, cv_d, R_d, grav, planet_radius
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

import ClimateMachine.DGMethods: vars_state_conservative, vars_state_auxiliary
import ClimateMachine.Atmos: source!, atmos_source!, altitude
import ClimateMachine.Atmos: thermo_state

using Distributions
using Random
using StaticArrays
using Test
using DocStringExtensions
using LinearAlgebra
using Logging

"""
    baroclinic_instability_cube(...)
Initialisation helper for baroclinic-wave (channel flow) test case for iterative
determination of η = p/pₛ coordinate for given z-altitude. 
"""
function baroclinic_instability_cube(FT, (x,y,z), param_set)

  ### Unpack CLIMAParameters
  _planet_radius = FT(planet_radius(param_set))
  gravity        = FT(grav(param_set))
  cp             = FT(cp_d(param_set))
  R_gas          = FT(R_d(param_set))

  ### Global Variables
  up    = FT(1)                ## See paper: Perturbation peak value
  Lp    = FT(6e5)              ## Perturbation parameter (radius)
  Lp2   = Lp*Lp              
  xc    = FT(2e6)              ## Streamwise center of perturbation
  yc    = FT(25e6)             ## Spanwise center of perturbation
  gamma_lapse = FT(5/1000)    ## Γ Lapse Rate
  Ω     = FT(7.292e-5)         ## Rotation rate [rad/s]
  f0    = 2Ω/sqrt(2)           ## 
  beta0 = f0/_planet_radius    ##  
  beta0 = -zero(FT)
  b     = FT(2)
  b2    = b*b
  u0    = FT(35)
  Ly    = FT(6e6)
  T0    = FT(288)
  T_ref = T0                   
  x0    = FT(2e7)
  p00   = FT(1e5)              ## Surface pressure

  #Step 1: Get current coordinate value by unpacking nodal coordinates from aux state
  eta = FT(1e-7)
  global geo_phi = FT(1.0)
  global temp = FT(1.0)
  @show(geo_phi, temp, "Before iterations")

  (geo_phi, temp, η) = baroclinic_iterations(geo_phi, temp, η, 

  for n = 1:100 ## Begin convergence loop
    baroclinic_instability_cube_functions!(geo_phi,temp,y,eta,f0,beta0,u0,T0,gamma_lapse, gravity, R_gas, cp)
    num  = -gravity*z + geo_phi
    den  = -R_gas/(eta)*temp
    deta = num/den
    eta  = eta - deta
    if (abs(deta) <= FT(1e-10))
        break
    elseif (abs(deta) > FT(1e-10)) && n<100
        continue
    elseif (abs(deta) > FT(1e-10)) && n==100
      @error "Baroclinic Instability Initialisation: η iterations did not converge."
    end
  end ## End convergence loop
  @show(eta, "After iterations")
  eta = min(eta,FT(1))
  eta = max(eta,FT(0))
  logeta = log(eta)
  T=temp
  press = p00*eta
  theta = T *(p00/press)^(R_gas/cp)
  rho = press/(R_gas*T)
  thetaref = T_ref * (1 - gamma_lapse*z/T0)^(-gravity/(cp*gamma_lapse))
  rhoref = p00/(T0*R_gas) * (1 - gamma_lapse*z/T0)^(gravity/(R_gas*gamma_lapse) - 1)

  ### Balanced Flow
  u = -u0*(sin(π*y/Ly))^2  * logeta * exp(-logeta*logeta/b2)

  ### Perturbation of the balanced flow
  rc2 = (x-xc)^2 + (y-yc)^2
  du = up*exp(-rc2/Lp2)
    
  ### Return u, du and combine those in the initial condition assignment. 
  ### (For testing purposes we use these separately)
  ### Returned quantities allow ClimateMachine state variables to be built
  return rho, u, du, T
end 

function baroclinic_instability_cube_functions!(geo_phi, temp,y,eta,f0,beta0,u0,T0,gamma_lapse,gravity, R_gas, cp)
  FT    = eltype(y)
  b     = FT(2)
  Ly    = FT(6e6)  
  y0    = FT(Ly/2)
  b2    = b*b
  #Get Mean Temperature
  exp1  = R_gas*gamma_lapse/gravity
  Tmean = T0*eta^exp1
  phimean = T0*gravity/gamma_lapse * (FT(1) - eta^exp1)

  logeta = log(eta)
  fac1   = (f0-beta0*y0)*(y - FT(1/2)*Ly - Ly/2π * sin(2π*y/Ly))  
  fac2   = FT(1/2)*beta0*(y^2 - Ly*y/π*sin(2π*y/Ly) - 
                        FT(1/2)*(Ly/π)^2*cos(2π*y/Ly) - 
                        Ly^2/3 - FT(1/2)*(Ly/π)^2)
  fac3 = exp(-logeta*logeta/b2)
  phi_prime=FT(0.5)*u0*(fac1 + fac2)
  geo_phi = phimean + phi_prime*fac3*logeta
  temp = Tmean + phi_prime/R_gas*fac3*(2/b2*logeta*logeta - 1)
  @show(geo_phi, temp, "Inside BICF")
  return nothing
end 

function init_baroclinicwave!(bl, state, aux, (x, y, z), t)
    ### Problem float-type
    FT = eltype(state)
    param_set = bl.param_set
    gravity = FT(grav(param_set))

    ### Unpack initial conditions (solved by iterating for η)
    (rho, u, du, T) = baroclinic_instability_cube(FT, (x,y,z), param_set)
    
    ### Primitive variables
    u⃗ = SVector{3,FT}(u,0,0)
    e_kin = FT(1/2)*sum(abs2.(u⃗))
    e_pot = gravity * z

    ### Assign state variables for initial condition
    state.ρ = rho
    state.ρu = rho .* u⃗
    state.ρe = rho * total_energy(param_set, e_kin, e_pot, T)
end

function config_baroclinicwave(FT, N, resolution, xmax, ymax, zmax)

    ics = init_baroclinicwave!     # Initial conditions
    C_smag = FT(0.23)     # Smagorinsky coefficient
    
    # Assemble source components
    source = (
        Gravity(),
    )

    # Choose default IMEX solver
    ode_solver_type = ClimateMachine.IMEXSolverType()

    # Assemble model components
    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        turbulence = SmagorinskyLilly{FT}(C_smag),
        moisture = DryModel(),
        source = source,
        boundarycondition = (
            AtmosBC(),
            AtmosBC(),
        ),
        init_state_conservative = ics,
    )

    # Assemble configuration
    config = ClimateMachine.AtmosLESConfiguration(
        "Baroclinic Wave",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_baroclinicwave!,
        solver_type = ode_solver_type,
        model = model,
    )
    return config
end

function config_diagnostics(driver_config)

    boundaries = [
        FT(0.0) FT(0.0) FT(0.0)
        FT(4e7) FT(6e6) FT(30e3)
    ]
    resolution = (FT(100e3), FT(75e3), FT(1250)) 
    interpol = ClimateMachine.InterpolationConfiguration(
        driver_config,
        boundaries,
        resolution,
    )

    default_dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        "1shours",
        driver_config.name,
    )
    core_dgngrp = setup_atmos_core_diagnostics(
        AtmosLESConfigType(),
        "1shours",
        driver_config.name,
    )
    dump_dgngrp = setup_dump_state_and_aux_diagnostics(
        AtmosLESConfigType(),
        "1shours",
        driver_config.name,
        interpol = interpol,
    )
    return ClimateMachine.DiagnosticsConfiguration([
        default_dgngrp,
        core_dgngrp,
        dump_dgngrp,
    ])
end

function main()
    FT = Float32

    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δx = FT(100000)
    Δy = FT(75000)
    Δz = FT(1250)

    resolution = (Δx, Δy, Δz)

    # Prescribe domain parameters
    xmax = FT(4e7)
    ymax = FT(6e6)
    zmax = FT(30e3)

    t0 = FT(0)

    # For a full-run, please set the timeend to 3600*6 seconds
    # For the test we set this to == 30 minutes
    days = FT(86400)
    timeend = FT(15days)
    CFLmax = FT(0.90)

    driver_config = config_baroclinicwave(FT, N, resolution, xmax, ymax, zmax)
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = CFLmax,
    )
    dgn_config = config_diagnostics(driver_config)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do (init = false)
        Filters.apply!(
            solver_config.Q,
            ("moisture.ρq_tot",),
            solver_config.dg.grid,
            TMARFilter(),
        )
        nothing
    end
    
    filterorder = 64
    filter = ExponentialFilter(solver_config.dg.grid, 0, filterorder)
    cbfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            AtmosFilterPerturbations(driver_config.bl),
            solver_config.dg.grid,
            filter,
            state_auxiliary = solver_config.dg.state_auxiliary,
        )
        nothing
    end

    # State variable
    Q = solver_config.Q
    # Volume geometry information
    vgeo = driver_config.grid.vgeo
    M = vgeo[:, Grids._M, :]
    # Unpack prognostic vars
    ρ₀ = Q.ρ
    ρe₀ = Q.ρe
    # DG variable sums
    Σρ₀ = sum(ρ₀ .* M)
    Σρe₀ = sum(ρe₀ .* M)
    cb_check_cons =
        GenericCallbacks.EveryXSimulationSteps(3000) do (init = false)
            Q = solver_config.Q
            δρ = (sum(Q.ρ .* M) - Σρ₀) / Σρ₀
            δρe = (sum(Q.ρe .* M) .- Σρe₀) ./ Σρe₀
            @show (abs(δρ))
            @show (abs(δρe))
            nothing
        end

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cbtmarfilter, cb_check_cons, cbfilter),
        check_euclidean_distance = true,
    )
end

main()
