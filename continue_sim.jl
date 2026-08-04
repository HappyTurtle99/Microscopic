include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/scratch.local/gtucci/micro/julia/homogenous_7000.00_sites_L_35.0Dc3.00kappa0.12_lambda5.00_2026-07-30_222801"

    par_old = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))
    
    Tfinal = 0.001
    
    output_dir = "/scratch.local/gtucci/micro/julia/homogenous_7000.00_sites_L_35.0Dc3.00kappa0.12_lambda5.00_2026-07-30_222801_extra_time_test"
    # output_dir = par_old.output_dir

    par = Params(
        par_old.Dn1, par_old.Dn2, par_old.Dc,
        par_old.gamma1, par_old.gamma2, par_old.kappa,
        par_old.μ, par_old.lambda1, par_old.lambda2,
        Tfinal,
        par_old.save_rate,
        par_old.save,
        output_dir
    )

    println("Copy 1: running for extra 100 time: /scratch.local/gtucci/micro/julia/homogenous_7000.00_sites_L_35.0Dc3.00kappa0.12_lambda5.00_2026-07-30_222801")

    run_sim!(st, par)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()
