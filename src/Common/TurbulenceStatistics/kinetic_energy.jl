function average_kinetic_energy_and_dissipation_writer(
    solver_config,
    driver_config,
    E_0,
    iter,
    outfile,
    nor,
)
    mpicomm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(mpicomm)
    if mpirank == 0
        E_k0 = average_kinetic_energy(solver_config, driver_config, nor)
        dE = -(E_k0 - E_0) / iter
        open(outfile, "a") do io
            writedlm(io, [E_k0 dE])
        end
        return E_k0
    end
end

function cb_kinetic_energy(
    solver_config,
    driver_config,
    iter,
    E_0,
    outfile,
    nor,
)
    cb = GenericCallbacks.EveryXSimulationTime(iter) do
        E_0 = average_kinetic_energy_and_dissipation_writer(
            solver_config,
            driver_config,
            E_0,
            iter,
            outfile,
            nor,
        )
    end
    return cb
end

function average_kinetic_energy(solver_config, driver_config, nor)
    Q = solver_config.Q
    # Volume geometry information
    vgeo = driver_config.grid.vgeo
    # Unpack prognostic vars
    u₀ = Q.ρu ./ Q.ρ
    u_0 = u₀[:, 1, :] ./ nor
    v_0 = u₀[:, 2, :] ./ nor
    w_0 = u₀[:, 3, :] ./ nor
    mm = size(u_0, 2)
    M = vgeo[:, Grids._M, 1:mm]
    SM = sum(M)
    E_k0 = 0.5 * sum((u_0 .^ 2 .+ v_0 .^ 2 .+ w_0 .^ 2) .* M) / SM
    mpicomm = MPI.COMM_WORLD
    mpirank = MPI.Comm_rank(mpicomm)
    n = MPI.Comm_size(mpicomm)
    E_k0 = MPI.Reduce(E_k0, +, 0, mpicomm)
    if mpirank == 0
        return E_k0 / n
    end
end
