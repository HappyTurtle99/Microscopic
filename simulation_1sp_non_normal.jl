include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates

function main()
    length(ARGS) == 3 || error("Usage: (here give kappa lambda and T) julia simulation_1sp.jl <lambda1> <Dn1> <T>")
# 11 1.0 1.0
    L = parse(Float64, ARGS[1])
    Nsites = Int64(round(150 * L))
    h = L / Nsites
#
    Dn1 = parse(Float64, ARGS[2])
    Dn2 = 0.0
    Dc = 1.0
    kappa = 0.0
    gamma1 = 2 * kappa
    gamma2 = 0.0
    zeta = 10.0
    lambda1 = parse(Float64, ARGS[3])
    lambda2 = 0.0
    Tfinal = 10

    rhon2 = 0
    rhoc = 2
    rhon1 = 1#Int(round(rhoc * kappa / gamma1))

    μ = 10.0
    # μ = 1.0
    # μ = Dn1 * kappa / (lambda1 * gamma1)

    Deff = Dn1 - lambda1 * (T1(rhoc, μ) - T0(rhoc, μ))
    @printf("Deff: %.6f\n", Deff)
    ν = lambda1 * T1(rhoc, μ) / (Deff + lambda1 * T1(rhoc, μ))
    @printf("nu: %.6f\n", ν)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    # output_dir = @sprintf(
    #     "/scratch.local/gtucci/micro/julia/l1_%.2f_Dc_%.2f_Nsites_%.2f_%s",
    #     lambda1, Dc, Nsites, timestamp
    # )

    output_dir = @sprintf(
        "nonnormal_zeta10_L_%.2f_Dn1_%.2f_lambda_%.2f_%s",
        L, Dn1, lambda1, timestamp
    )

    # rescaling
    Dn1 /= h^2
    Dn2 /= h^2
    Dc /= h^2
    lambda1 /= h^2
    lambda2 /= h^2

    par = Params(
        Dn1, Dn2, Dc,
        gamma1, gamma2, kappa, zeta,
        μ, lambda1, lambda2,
        Tfinal,
        3000000,     # save_rate
        true,       # save
        output_dir
    )

    st = initialize_state(Nsites, rhon1, rhon2, rhoc, μ; drho1=0, drhoc=0)
    run_sim!(st, par)

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)
end

main()