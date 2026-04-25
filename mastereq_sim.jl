include("MicroSim.jl")
using Statistics
using .MicroSim
using Printf
using NPZ
using Dates

function main()
    length(ARGS) >= 4 || error("Usage: julia simulation_1sp.jl <lambda1> <Dn1> <T>")

    L = 1.0
    Nsites = parse(Int64, ARGS[4])
    h = L / Nsites

    Dn1 = parse(Float64, ARGS[2])
    Dn1 = 0.0
    Dn2 = 0.0
    Dc = 0.2
    gamma1 = 1.0
    gamma2 = 0.0
    kappa = gamma1
    lambda1 = parse(Float64, ARGS[1])
    lambda1 = 0.0
    lambda2 = 0.0
    Tfinal = parse(Float64, ARGS[3])

    rhon1 = 320
    rhon2 = 0
    rhoc = rhon1

    μ = 10.0

    Deff = Dn1 - lambda1 * (T1(rhoc, μ) - T0(rhoc, μ))
    @printf("Deff: %.6f\n", Deff)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    output_dir = @sprintf(
        "l1_%.2f_Dn1_%.2f_T_%.2f_%s",
        lambda1, Dn1, Tfinal, timestamp
    )

    # rescaling
    Dn1 /= h^2
    Dn2 /= h^2
    Dc /= h^2
    lambda1 /= h^2
    lambda2 /= h^2

    extra_rate = (rhon1 / 8)*(Nsites / 100)^2

    par = Params(
        Dn1, Dn2, Dc,
        gamma1, gamma2, kappa,
        μ, lambda1, lambda2,
        Tfinal,
        500000,     # save_rate
        false,       # save, false because we average and wanna save later
        output_dir
    )

    st = initialize_state(Nsites, rhon1, rhon2, rhoc, μ; drhoc=rhoc, drho1=0.0)
    run_sim!(st, par)
    
    mkpath(output_dir)
    npzwrite(joinpath(output_dir, "occupancies_t.npy"), history_array(st))
    npzwrite(joinpath(output_dir, "tau_t.npy"), collect(st.tau_history))

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)

end

main()