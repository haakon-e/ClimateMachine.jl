module TurbulenceStatistics

using DelimitedFiles
using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.Mesh.Grids
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.VariableTemplates
using MPI

include("kinetic_energy.jl")

export average_kinetic_energy,
    cb_kinetic_energy, average_kinetic_energy_and_dissipation_write

end
