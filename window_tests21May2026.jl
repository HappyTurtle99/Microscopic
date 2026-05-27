include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates

function main()
    length(ARGS) == 2 || error("Usage: (here give kappa lambda and T) julia simulation_1sp.jl <Nsites> <T>")
    L = 3.14
    Nsites = parse(Int64, ARGS[1])
    h = L / Nsites
#
    Dn1 = 1.0
    Dn2 = 0.0
    Dc = 1.0
    kappa = 1.0
    gamma1 = 2 * kappa
    gamma2 = 0.0
    lambda1 = 1.0
    lambda2 = 0.0
    Tfinal = parse(Float64, ARGS[2])
    zeta = 0

    rhon2 = 0
    rhoc = 2
    rhon1 = Int(round(rhoc * kappa / gamma1))

    μ = 10.0
    # μ = 1.0
    # μ = Dn1 * kappa / (lambda1 * gamma1)

    Deff = Dn1 - lambda1 * (T1(rhoc, μ) - T0(rhoc, μ))
    @printf("Deff: %.6f\n", Deff)
    ν = lambda1 * T1(rhoc, μ) / (Deff + lambda1 * T1(rhoc, μ))
    @printf("nu: %.6f\n", ν)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    output_dir = @sprintf(
        "/scratch03.local/gtucci/micro/julia/window_t_L_12.56_Nsites_%.2f_Tfinal_%.2f_%s",
        Nsites, Tfinal, timestamp
    )

   # output_dir = @sprintf(
    #    "window_test_L_18.85_Nsites_%.2f_T_%.2f_%s",
     #   Nsites, Tfinal, timestamp
    #)

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

    st = initialize_state(Nsites, rhon1, rhon2, rhoc, μ; drho1=0, drhoc=0)
    run_sim!(st, par)

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)
end

main()
