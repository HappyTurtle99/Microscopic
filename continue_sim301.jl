include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/scratch03.local/gtucci/micro/julia/L_12.56_Dc_6.00_kappa_1.00_2026-06-03_162317"
    
    par_old = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))

    output_dir = "/scratch03.local/gtucci/micro/julia/L_12.56_Dc_6.00_kappa_1.00_contto_Tfinal115_2026-06-03_162317"
    Tfinal = st.tau + 80.0
    
    par = Params(
        par_old.Dn1, par_old.Dn2, par_old.Dc,
        par_old.gamma1, par_old.gamma2, par_old.kappa,
        par_old.μ, par_old.lambda1, par_old.lambda2,
        Tfinal,
        par_old.save_rate,
        par_old.save,
        output_dir
    )
    print(par.Dc)
    
    Tfinal = st.tau + 135.0
    
    println("running now", par.Dc / par.Dn1, par.kappa)

    run_sim!(st, par)
    dir_local = @sprintf("cluster_data/05jun/%.2f/%.2f", par.Dc / par.Dn1, par.kappa)
    save_sim_dir(st, par, dir_local)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()
