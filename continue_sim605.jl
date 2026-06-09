include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/scratch03.local/gtucci/micro/julia/L_12.56_Dc_6.00_kappa_0.50_2026-06-03_162317"
    
    par = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))
    
    Tfinal = st.tau + 80.0
    
    output_dir = "/scratch03.local/gtucci/micro/julia/L_12.56_Dc_6.00_kappa_0.50_contto_Tfinal115_2026-06-03_162317"
    
    println("running now")
    

    run_sim!(st, par)
    dir_local = @sprintf("cluster_data/05jun/%.2f/%.2f", par.Dc / par.Dn1, par.kappa)
    save_sim_dir(st, par, dir_local)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()
