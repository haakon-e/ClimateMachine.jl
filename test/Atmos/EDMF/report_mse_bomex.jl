using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));

if parse(Bool, get(ENV, "CLIMATEMACHINE_PLOT_EDMF_COMPARISON", "false"))
    plot_dir = joinpath(clima_dir, "output", "bomex_edmf", "pycles_comparison")
else
    plot_dir = nothing
end

include(joinpath(@__DIR__, "compute_mse.jl"))

data_file = Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r")

#! format: off
best_mse = OrderedDict()
best_mse["prog_ρ"] = 3.4917662870525452e-02
best_mse["prog_ρu_1"] = 3.0715053983126868e+03
best_mse["prog_ρu_2"] = 1.2895738234435595e-03
best_mse["prog_moisture_ρq_tot"] = 4.1331426262587470e-02
best_mse["prog_turbconv_environment_ρatke"] = 6.7642719033845833e+02
best_mse["prog_turbconv_environment_ρaθ_liq_cv"] = 8.5667226792770009e+01
best_mse["prog_turbconv_environment_ρaq_tot_cv"] = 1.6435640302009438e+02
best_mse["prog_turbconv_updraft_1_ρa"] = 8.0846149716812874e+01
best_mse["prog_turbconv_updraft_1_ρaw"] = 8.5071202823418249e-02
best_mse["prog_turbconv_updraft_1_ρaθ_liq"] = 9.0386633006301569e+00
best_mse["prog_turbconv_updraft_1_ρaq_tot"] = 1.0803376910477068e+01
#! format: on

computed_mse = compute_mse(
    solver_config.dg.grid,
    solver_config.dg.balance_law,
    time_data,
    dons_arr,
    data_file,
    "Bomex",
    best_mse,
    400,
    plot_dir,
)

@testset "BOMEX EDMF Solution Quality Assurance (QA) tests" begin
    #! format: off
    test_mse(computed_mse, best_mse, "prog_ρ")
    test_mse(computed_mse, best_mse, "prog_ρu_1")
    test_mse(computed_mse, best_mse, "prog_ρu_2")
    test_mse(computed_mse, best_mse, "prog_moisture_ρq_tot")
    test_mse(computed_mse, best_mse, "prog_turbconv_environment_ρatke")
    test_mse(computed_mse, best_mse, "prog_turbconv_environment_ρaθ_liq_cv")
    test_mse(computed_mse, best_mse, "prog_turbconv_environment_ρaq_tot_cv")
    test_mse(computed_mse, best_mse, "prog_turbconv_updraft_1_ρa")
    test_mse(computed_mse, best_mse, "prog_turbconv_updraft_1_ρaw")
    test_mse(computed_mse, best_mse, "prog_turbconv_updraft_1_ρaθ_liq")
    test_mse(computed_mse, best_mse, "prog_turbconv_updraft_1_ρaq_tot")
    #! format: on
end
