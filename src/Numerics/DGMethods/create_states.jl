#### Create states

function create_state(balance_law::BalanceLaw, grid, vt::AbstractStateType)
    topology = grid.topology
    # FIXME: Remove after updating CUDA
    h_vgeo = Array(grid.vgeo)
    FT = eltype(h_vgeo)
    Np = dofs_per_element(grid)
    DA = arraytype(grid)

    weights = view(h_vgeo, :, grid.Mid, :)
    weights = reshape(weights, size(weights, 1), 1, size(weights, 2))

    V = vars_state(balance_law, vt, FT)
    ns = number_states(balance_law, vt, FT)
    state = MPIStateArray{FT, V}(
        topology.mpicomm,
        DA,
        Np,
        rank_multiplier(vt)*ns,
        length(topology.elems),
        realelems = topology.realelems,
        ghostelems = topology.ghostelems,
        vmaprecv = grid.vmaprecv,
        vmapsend = grid.vmapsend,
        nabrtorank = topology.nabrtorank,
        nabrtovmaprecv = grid.nabrtovmaprecv,
        nabrtovmapsend = grid.nabrtovmapsend,
        weights = weights,
    )
    return state
end

function init_state(state, balance_law, grid)
    topology = grid.topology
    Np = dofs_per_element(grid)

    dim = dimensionality(grid)
    polyorder = polynomialorder(grid)
    vgeo = grid.vgeo
    device = array_device(state)
    nrealelem = length(topology.realelems)
    event = Event(device)
    event = kernel_init_state_auxiliary!(device, min(Np, 1024), Np * nrealelem)(
        balance_law,
        Val(dim),
        Val(polyorder),
        state.data,
        vgeo,
        topology.realelems,
        dependencies = (event,),
    )
    event = MPIStateArrays.begin_ghost_exchange!(
        state;
        dependencies = event,
    )
    event = MPIStateArrays.end_ghost_exchange!(
        state;
        dependencies = event,
    )
    wait(device, event)

    return state
end

