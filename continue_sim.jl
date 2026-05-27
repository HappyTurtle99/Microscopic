include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/scratch03.local/gtucci/micro/julia/window_t_L_12.56_Nsites_3000.00_Tfinal_150.00_2026-05-19_155535"
    
    par_old = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))
    
    Tfinal = st.tau + 35

    if par_old.Dc == 10.0
        println("previous Dc = 10, nvm!")
        return nothing
    end

    output_dir = "/scratch03.local/gtucci/micro/julia/window_t_L_12.56_Nsites_3000.00_Tfinal_150.00_35_continued_Dc_10_2026-05-19_155535"
    # output_dir = par_old.output_dir
    
    par = Params(
        par_old.Dn1, par_old.Dn2, par_old.Dc * 0 + 10.0,
        par_old.gamma1, par_old.gamma2, par_old.kappa,
        par_old.μ, par_old.lambda1, par_old.lambda2,
        Tfinal,
        par_old.save_rate,
        par_old.save,
        output_dir
    )

    rhoc = sum(st.occc) / length(st.occc)
    Deff = par.Dn1 - par.lambda1 * (T1(rhoc, par.μ) - T0(rhoc, par.μ))
    @printf("Deff: %.6f\n", Deff)
    @printf("par.lambda1: %.6f\n", par.lambda1)
    @printf("par.rhoc: %.6f\n", rhoc)
    return nothing

    ν = par.lambda1 * T1(rhoc, par.μ) / (Deff + par.lambda1 * T1(rhoc, par.μ))
    @printf("nu: %.6f\n", ν)

    run_sim!(st, par)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()