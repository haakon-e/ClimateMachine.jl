# # Filters
# In this tutorial we show the result of applying filters
# available in the CliMA codebase in a 1 dimensional box advection setup.
# See [Filters API](https://clima.github.io/ClimateMachine.jl/latest/APIs/Numerics/Meshes/Mesh/#Filters-1) for filters interface details.

using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "tutorials", "Numerics", "DGMethods", "Box1D.jl"))

# The unfiltered result of the box advection test for order 4 polynomial is:
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_no_filter.pdf")
# ![](box_1D_4_no_filter.pdf)

# Below we show results for the same box advection test but using different filters.

# `TMARFilter()`:
run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar.pdf", tmar_filter = true)
# ![](box_1D_4_tmar.pdf)

# `CutoffFilter(grid, Nc=1)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_cutoff_1.pdf",
    cutoff_filter = true,
    cutoff_param = 1,
)
# ![](box_1D_4_cutoff_1.pdf)

# `CutoffFilter(grid, Nc=3)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_cutoff_3.pdf",
    cutoff_filter = true,
    cutoff_param = 3,
)
# ![](box_1D_4_cutoff_3.pdf)

# `ExponentialFilter(grid, Nc=1, s=4)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_4.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 4,
)
# ![](box_1D_4_exp_1_4.pdf)

# `ExponentialFilter(grid, Nc=1, s=8)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_8.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 8,
)
# ![](box_1D_4_exp_1_8.pdf)

# `ExponentialFilter(grid, Nc=1, s=32)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_32.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 32,
)
# ![](box_1D_4_exp_1_32.pdf)

# `BoydVandevenFilter(grid, Nc=1, s=4)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_4.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 4,
)
# ![](box_1D_4_boyd_1_4.pdf)

# `BoydVandevenFilter(grid, Nc=1, s=8)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_8.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 8,
)
# ![](box_1D_4_boyd_1_8.pdf)

# `BoydVandevenFilter(grid, Nc=1, s=32)`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_32.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 32,
)
# ![](box_1D_4_boyd_1_32.pdf)

# `ExponentialFilter(grid, Nc=1, s=8)` and `TMARFilter()`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_tmar_exp_1_8.pdf",
    exp_filter = true,
    tmar_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 8,
)
# ![](box_1D_4_tmar_exp_1_8.pdf)

# `BoydVandevenFilter(grid, Nc=1, s=8)` and `TMARFilter()`:
run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_tmar_boyd_1_8.pdf",
    boyd_filter = true,
    tmar_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 8,
)
# ![](box_1D_4_tmar_boyd_1_8.pdf)
