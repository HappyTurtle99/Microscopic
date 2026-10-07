include("MicroSimFast.jl")
using .MicroSimFast
using Printf
using Dates
using Serialization

function main()
    path = "/scratch03.local/gtucci/micro/julia/homogenous_5950.00_sites_L_35.0Dc3.00kappa0.12_lambda5.00_2026-09-24_145947"
    
    par = deserialize(joinpath(path, "Params.bin"))
    st = deserialize(joinpath(path, "SimState.bin"))

    println(st.tau)
    Tfinal = st.tau + 0.15

    output_dir = "$(path)_Tfinal$(Tfinal)"

    par_new = Params(
        par.Dn1, par.Dn2, par.Dc,
        par.gamma1, par.gamma2, par.kappa,
        par.μ, par.lambda1, par.lambda2,
        Tfinal, # 
        par.save_rate,
        par.save,
        output_dir # (changed)
    )

    println(output_dir)

    run_sim!(st, par_new)

    println("Job finished! Number of sites: ", length(st.occc), " Time: ", Tfinal)
end

main()
