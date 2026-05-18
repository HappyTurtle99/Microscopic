include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using NPZ

function main()
    length(ARGS) == 3 || error("Usage: (here give kappa lambda and T) julia simulation_1sp.jl <lambda1> <Dn1> <T>")
    L = parse(Float64, ARGS[1]) 
    Nsites = Int64(round(100 * L))
    h = L / Nsites
#
    Dn1 = parse(Float64, ARGS[2]) 
    Dn2 = 0.0
    Dc = 1.0
    kappa = 1.0
    gamma1 = kappa
    gamma2 = 0.0
    lambda1 = parse(Float64, ARGS[3])
    lambda2 = 0.0
    Tfinal = 0.01

    rhon2 = 0
    rhoc = 10
    rhon1 = Int(round(rhoc * kappa / gamma1))

    μ = 0.1
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

    # output_dir = @sprintf(
    #     "/Users/sarce/Microscopic/stepcheck",
    #     L, Dn1, lambda1, timestamp
    # )
    output_dir = "/Users/sarce/Microscopic/stepcheck"

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
        3000000,     # save_rate
        false,       # save
        output_dir
    )

    st = initialize_state(Nsites, rhon1, rhon2, rhoc, μ; drho1=rhon1, drhoc=rhoc)
    run_sim!(st, par)

    sample_weights_empirical =  zeros(Float64, 2 * Nsites)

    for i in 1:2*20*Nsites
        sample_weights_empirical[MicroSimFast.sample(st.chemo_rate_tree1)] += 1
    end
    outdir = par.output_dir

    mkpath(outdir)
    npzwrite(joinpath(outdir, "chemo_rates1_sim.npy"), st.chemo_rates1)
    npzwrite(joinpath(outdir, "chemo_rate_tree1_sim.npy"), st.chemo_rate_tree1.tree)
    npzwrite(joinpath(outdir, "sample_weights.npy"), sample_weights_empirical)
    npzwrite(joinpath(outdir, "occc.npy"), st.occc)

    println("Job finished! Number of sites: ", Nsites, " Time: ", Tfinal)
end

main()