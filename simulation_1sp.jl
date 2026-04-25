include("MicroSim.jl")
using .MicroSim
using Printf
using Dates

function main()
    length(ARGS) >= 3 || error("Usage: julia simulation_1sp.jl <lambda1> <Dn1> <T>")

    L = 1.0
    Nsites = parse(Int64, ARGS[1])
    h = L / Nsites

    Dn1 = 0.1
    Dn2 = 0.0
    Dc = parse(Float64, ARGS[2])
    gamma1 = 0.2
    gamma2 = 0.0
    kappa = 0.1
    lambda1 = parse(Float64, ARGS[3])
    lambda2 = 0.0
    Tfinal = 0.001

    if lambda1 > 2
        Tfinal = 5
    else 
        Tfinal = 12.5
    end

    rhon1 = 150
    rhon2 = 0
    rhoc = Int(round(rhon1 * gamma1 / kappa))

    μ = Dn1 * kappa / (lambda1 * gamma1)

    Deff = Dn1 - lambda1 * (T1(rhoc, μ) - T0(rhoc, μ))
    @printf("Deff: %.6f\n", Deff)
    ν = lambda1 * T1(rhoc, μ) / (Deff + lambda1 * T1(rhoc, μ))
    @printf("nu: %.6f\n", ν)

    μ = μ * ν
    μ = 10.0

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    output_dir = @sprintf(
        "/scratch.local/gtucci/micro/julia/l1_%.2f_Dn1_%.2f_T_%.2f_%s",
        lambda1, Dn1, Tfinal, timestamp
    )

    # output_dir = @sprintf(
    #     "/Users/sarce/Microscopic/l1_%.2f_Dn1_%.2f_T_%.2f_%s",
    #     lambda1, Dn1, Tfinal, timestamp
    # )

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
        50000000,     # save_rate
        true,       # save
        output_dir
    )

    st = initialize_state(Nsites, rhon1, rhon2, rhoc, μ; drho1=rhon1*0.1, drhoc=rhoc*0.1)
    run_sim!(st, par)

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)
end

main()