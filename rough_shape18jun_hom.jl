include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates

function main()
    # length(ARGS) == 2 || error("Usage: (here give kappa lambda and T) julia simulation_1sp.jl <Nsites> <T>")
    L = 35.00
    Nsites = 7000 * parse(Float64, ARGS[1]) / 200 # 100 150 200 250 
    h = L / Nsites

    Dn1 = 1.0
    Dn2 = 0.0
    Dc = 3.0 
    kappa = 0.125
    gamma1 = 2 * kappa
    gamma2 = 0.0
    lambda1 = 1.0
    lambda2 = 0.0
    Tfinal = 50.0
    zeta = 0

    dummy = 7000 * parse(Float64, ARGS[2]) / 200

    rhon2 = 0
    rhoc = 2
    rhon1 = Int(round(rhoc * kappa / gamma1))

    μ = 1.0 # f = mu c
    # μ = 1.0
    # μ = Dn1 * kappa / (lambda1 * gamma1)

    total_rate = (Dn1 * rhon1 + Dc * rhoc + μ * rhoc * rhon1 * lambda1) * Nsites

    Deff = Dn1 - lambda1 * (T1(rhoc, μ) - T0(rhoc, μ))
    println(T1(rhoc, μ))
    println(T0(rhoc, μ))
    @printf("Deff: %.6f\n", Deff)
    @printf("Total rate: %.6f\n", total_rate)

    ν = lambda1 * T1(rhoc, μ) / (Deff + lambda1 * T1(rhoc, μ))
    @printf("nu: %.6f\n", ν)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    output_dir = @sprintf(
        "/scratch03.local/gtucci/micro/julia/%.2f_sites_L_35.0Dc%.2fkappa%.2f_lambda%.2f_%s",
        Nsites, Dc, kappa,lambda1, timestamp,
     )

#    output_dir = @sprintf(
#       "window_test_L_12.56Nsites%.2fT%.2f_%s",
#       Nsites, Tfinal, timestamp
#    )

    # rescaling
    Dn1 /= h^2
    Dn2 /= h^2
    Dc /= h^2
    lambda1 /= h^2
    lambda2 /= h^2

    par = Params(
        Dn1, Dn2, Dc,
        gamma1, gamma2, kappa,
        μ, lambda1, lambda2,
        Tfinal,
        30000000,     # save_rate
        true,       # save
        output_dir
    )

    st = initialize_state(round(Int, Nsites), rhon1, rhon2, rhoc, μ; drho1=0, drhoc=0)
    println(sum(st.occc))
    println(sum(st.occ1))
    run_sim!(st, par)
    # dir_local = @sprintf(
    #     "cluster_data/03jun/Dc%.2f/kappa%.2f/",
    #     par.Dc/par.Dn1, par.kappa)
    # save_sim_dir(st, par, dir_local)

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)
end

main()