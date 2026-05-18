include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/Users/sarce/Microscopic/L_12.56_Dn1_1.00_lambda_10.00_2026-05-15_131930"
    
    par_old = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))
    
    Tfinal = st.tau + 1

    par = Params(
        par_old.Dn1, par_old.Dn2, par_old.Dc,
        par_old.gamma1, par_old.gamma2, par_old.kappa,
        par_old.μ, par_old.lambda1, par_old.lambda2,
        Tfinal,
        par_old.save_rate,
        par_old.save,
        par_old.output_dir
    )

    rhoc = sum(st.occc) / length(st.occc)
    Deff = par.Dn1 - par.lambda1 * (T1(rhoc, par.μ) - T0(rhoc, par.μ))
    @printf("Deff: %.6f\n", Deff)
    ν = par.lambda1 * T1(rhoc, par.μ) / (Deff + par.lambda1 * T1(rhoc, par.μ))
    @printf("nu: %.6f\n", ν)

    run_sim!(st, par)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()