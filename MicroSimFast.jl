module MicroSimFast

using Random
using SpecialFunctions
using Statistics
using Printf
using Serialization
using NPZ

export Params, RandomDict, SimState,
      initialize_state, run_sim!, save_sim,
     T0, T1, meso_avg, history_array
     
# ----------------------------
# Utilities
# ----------------------------

@inline periodic_index(i::Int, n::Int) = mod(i - 1, n) + 1
@inline f(c::Real, μ::Real=1.0) = max(0.0, tanh(μ * c))

function besseli_scaled_asymptotic(ν, z; terms=3)
    t = 1 / (8z)
    μ = 4ν^2

    # build series iteratively
    term = 1.0
    sum = 1.0

    for k in 1:terms
        term *= -(μ - (2k - 1)^2) * t / k
        sum += term
    end

    return sum / sqrt(2π*z)
end

function p(γ, m)
    z = 2.0 * γ
    if z < 80
        return exp(-z) * besseli(m, z)
    else
        return besseli_scaled_asymptotic(m, z; terms=4)
    end
end

function T0(γ::Real, μ::Real=1.0; m::Int=70)
    s = 0.0
    for i in -m:m
        s += f(i, μ) * p(γ, i)
    end
    return s
end

function T1(γ::Real, μ::Real=1.0; m::Int=70)
    s = 0.0
    for i in -m:m
        s += i * f(i, μ) * p(γ, i)
    end
    return s
end

function report_progress(tau::Float64, T::Float64, next_progress::Float64, start_time::Float64)
    progress = tau / T
    elapsed = time() - start_time
    while progress >= next_progress && next_progress <= 1.0
        @printf("%d%% completed (tau = %.4f / %.4f) | elapsed: %.2f s\n",
                round(Int, 100 * next_progress), tau, T, elapsed)
        next_progress += 0.05
    end
    return next_progress
end

# Fenwick Tree / Binary Indexed Tree
#
# Supports:
#   - point updates
#   - prefix sums
#   - range sums
#   - weighted random sampling

mutable struct FenwickTree
    tree::Vector{Float64}
    n::Int
end

# Create an empty Fenwick tree of size n
function FenwickTree(n::Int)
    return FenwickTree(zeros(Float64, n), n)
end

# Build Fenwick tree from an array
function FenwickTree(arr::AbstractVector{<:Real})
    n = length(arr)

    ft = FenwickTree(n)

    for i in 1:n
        update!(ft, i, Float64(arr[i]))
    end

    return ft
end

# Change element at i to new_val
function update!(ft::FenwickTree, i::Int, delta::Float64)
    while i <= ft.n
        ft.tree[i] += delta
        i += i & (-i)
    end

    return nothing
end

# Compute prefix sum from 1 to i
function prefix_sum(ft::FenwickTree, i::Int)

    result = 0.0

    while i > 0
        result += ft.tree[i]
        i -= i & (-i)
    end

    return result
end

# Compute sum over range [left, right]
function range_sum(ft::FenwickTree, left::Int, right::Int)

    if left > right
        return 0.0
    end

    if left == 1
        return prefix_sum(ft, right)
    end

    return prefix_sum(ft, right) - prefix_sum(ft, left - 1)
end

# Total sum of all elements
function total_sum(ft::FenwickTree)
    return prefix_sum(ft, ft.n)
end

# Find smallest index such that prefix_sum(index) > target
#
# Assumes all stored values are nonnegative
function find_prefix(ft::FenwickTree, target::Float64)

    idx = 0
    bit = 1

    # Find largest power of 2 <= n
    while bit <= ft.n
        bit <<= 1
    end

    bit >>= 1

    current_sum = 0.0
    while bit > 0 

        next_idx = idx + bit

        if next_idx <= ft.n && current_sum + ft.tree[next_idx] <= target # bit: 1, will be true if epsilon>0, false otherwise
            current_sum += ft.tree[next_idx] 
            idx = next_idx
        end

        bit >>= 1
    end

    return idx + 1
end

# Weighted random sample
#
# Returns an index with probability proportional
# to its stored weight
function sample(ft::FenwickTree)

    total = total_sum(ft)

    if total <= 0
        error("Cannot sample from empty/nonpositive weights")
    end

    r = rand() * total

    return find_prefix(ft, r)
end

# ----------------------------
# RandomDict
# ----------------------------

mutable struct RandomDict
    d::Dict{Int,Int}                # particle_id => site
    keys::Vector{Int}               # active particle ids
    index::Dict{Int,Int}            # particle_id => index in keys
    val_keys::Dict{Int,Vector{Int}} # site => particle ids at site
    max_key::Int
end

RandomDict() = RandomDict(Dict{Int,Int}(),
                          Int[],
                          Dict{Int,Int}(),
                          Dict{Int,Vector{Int}}(),
                          0)

function Base.copy(rd::RandomDict)
    return RandomDict(
        copy(rd.d),
        copy(rd.keys),
        copy(rd.index),
        Dict(k => copy(v) for (k, v) in rd.val_keys),
        rd.max_key
    )
end

function insert!(rd::RandomDict, key::Int, value::Int)
    if !haskey(rd.d, key)
        push!(rd.keys, key)
        rd.index[key] = length(rd.keys)
    end

    old_val = get(rd.d, key, nothing)
    rd.d[key] = value

    if old_val !== nothing
        vec = rd.val_keys[old_val]
        pos = findfirst(==(key), vec)
        pos !== nothing && deleteat!(vec, pos)
        isempty(vec) && delete!(rd.val_keys, old_val)
    end

    push!(get!(rd.val_keys, value, Int[]), key)
    rd.max_key = max(rd.max_key, key)
    return rd
end

function insert_value!(rd::RandomDict, value::Int)
    key = rd.max_key + 1
    insert!(rd, key, value)
    return key
end

function Base.delete!(rd::RandomDict, key::Int)
    haskey(rd.d, key) || return rd # shouldn't be needed.

    value = rd.d[key]
    vec = rd.val_keys[value]
    pos = findfirst(==(key), vec)
    pos !== nothing && deleteat!(vec, pos)
    isempty(vec) && Base.delete!(rd.val_keys, value)

    idx = rd.index[key]
    last_key = rd.keys[end]
    rd.keys[idx] = last_key
    rd.index[last_key] = idx
    pop!(rd.keys)

    Base.delete!(rd.index, key)
    Base.delete!(rd.d, key)
    return rd
end

@inline function get_random(rd::RandomDict)
    key = rand(rd.keys)
    return key, rd.d[key]
end

@inline function random_key_from_val(rd::RandomDict, value::Int)
    return rand(rd.val_keys[value])
end

# ----------------------------
# Parameters and state
# ----------------------------

struct Params
    Dn1::Float64
    Dn2::Float64
    Dc::Float64
    gamma1::Float64
    gamma2::Float64
    kappa::Float64
    μ::Float64
    lambda1::Float64
    lambda2::Float64
    Tfinal::Float64
    save_rate::Int
    save::Bool
    output_dir::String
end

mutable struct SimState
    occ1::Vector{Int}
    occ2::Vector{Int}
    occc::Vector{Int}

    dict1::RandomDict
    dict2::RandomDict
    dictc::RandomDict

    chemo_rates1::Vector{Float64}
    chemo_rates2::Vector{Float64}
    chemo_rate_tree1::FenwickTree
    chemo_rate_tree2::FenwickTree

    N1::Int
    N2::Int
    Nc::Int

    tau::Float64
    counter::Int
    cache::Vector{ComplexF64}

    occ_history::Vector{Array{Int,2}}
    tau_history::Vector{Float64}
end

# ----------------------------
# Initialization
# ----------------------------

function initialize_field(N::Int; drho::Real=0.0, rho::Real=1.0, hom::Bool=false, k::Int=1)
    occ = Vector{Int}(undef, N)

    if hom
        fill!(occ, round(Int, rho))
    else
        kk =  k * (2π) / N
        for i in 1:N
            if rho != 0
                occ[i] = floor(Int, rho * (1 + (drho / (rho)) * real(exp(im * (i - 1) * kk))))
            else
                occ[i] = 0
            end
        end
    end

    d = RandomDict()
    p_index = 0
    for i in eachindex(occ)
        for _ in 1:occ[i]
            insert!(d, p_index, i)
            p_index += 1
        end
    end

    return occ, d
end

function initialize_chemo_rates(occupancies_n::Vector{Int},
                                occupancies_c::Vector{Int},
                                μ::Float64,
                                attractive::Bool)
    n = length(occupancies_c)
    rates = zeros(Float64, 2n)
    rates_tree = FenwickTree(2n)
    for i in 1:n
        get_new_chemo_rates!(rates, rates_tree, occupancies_n, occupancies_c, i, μ, attractive)
    end
    return rates, rates_tree
end

function initialize_state(Nsites::Int, rhon1::Int, rhon2::Int, rhoc::Int, μ::Float64; drho1::Real=0.0, drho2::Real=0.0, drhoc::Real=0.0)
    occ1, dict1 = initialize_field(Nsites; drho=drho1, rho=rhon1)
    occ2, dict2 = initialize_field(Nsites; drho=drho2, rho=rhon2)
    occc, dictc = initialize_field(Nsites; drho=drhoc, rho=rhoc)

    rates1, tree1 = initialize_chemo_rates(occ1, occc, μ, true)
    rates2, tree2 = initialize_chemo_rates(occ2, occc, μ, false)

    init_snapshot = vcat(permutedims(occ1), permutedims(occ2), permutedims(occc))

    return SimState(
        occ1, occ2, occc,
        dict1, dict2, dictc,
        rates1, rates2, tree1, tree2,
        sum(occ1), sum(occ2), sum(occc),
        0.0, 0, ComplexF64[0 + 0im, 0 + 0im],
        [init_snapshot], [0.0]
    )
end


# ----------------------------
# Chemotactic rates
# ----------------------------

function get_new_chemo_rates!(rates::Vector{Float64},
                              rates_tree::FenwickTree,
                              occupancies_n::Vector{Int},
                              occupancies_c::Vector{Int},
                              index::Int,
                              μ::Float64,
                              attractive::Bool)
    n = length(occupancies_c)

    l = occupancies_c[periodic_index(index - 1, n)] - occupancies_c[index]
    r = occupancies_c[periodic_index(index + 1, n)] - occupancies_c[index]

    if !attractive
        l = -l
        r = -r
    end

    left_pos  = 2(index - 1) + 1
    right_pos = 2(index - 1) + 2

    if occupancies_n[index] == 0
        update!(rates_tree, left_pos, 0.0 - rates[left_pos])
        rates[left_pos] = 0.0
        update!(rates_tree, right_pos, 0.0 - rates[right_pos])
        rates[right_pos] = 0.0
    else
        new_rate = f(l, μ) * occupancies_n[index]
        update!(rates_tree, left_pos, new_rate - rates[left_pos])
        rates[left_pos] = new_rate

        new_rate = f(r, μ) * occupancies_n[index]
        update!(rates_tree, right_pos, new_rate - rates[right_pos])
        rates[right_pos] = new_rate
    end

    left_site = periodic_index(index - 1, n)
    left_site_right_pos = 2(left_site - 1) + 2
    left_site_left_pos = 2(left_site - 1) + 1
    if occupancies_n[left_site] == 0
        update!(rates_tree, left_site_right_pos, 0.0 - rates[left_site_right_pos])
        rates[left_site_right_pos] = 0.0
        update!(rates_tree, left_site_left_pos, 0.0 - rates[left_site_left_pos])
        rates[left_site_left_pos] = 0.0
    else
        new_rate = f(-l, μ) * occupancies_n[left_site]
        update!(rates_tree, left_site_right_pos, new_rate - rates[left_site_right_pos])
        rates[left_site_right_pos] = new_rate
    end

    right_site = periodic_index(index + 1, n)
    right_site_left_pos = 2(right_site - 1) + 1
    right_site_right_pos = 2(right_site - 1) + 2

    if occupancies_n[right_site] == 0
        update!(rates_tree, right_site_left_pos, 0.0 - rates[right_site_left_pos])
        rates[right_site_left_pos] = 0.0
        update!(rates_tree, right_site_right_pos, 0.0 - rates[right_site_right_pos])
        rates[right_site_right_pos] = 0.0
    else
        new_rate = f(-r, μ) * occupancies_n[right_site]
        update!(rates_tree, right_site_left_pos, new_rate - rates[right_site_left_pos])
        rates[right_site_left_pos] = new_rate
    end

    return nothing
end

function get_new_chemo_rates_diff!(rates::Vector{Float64},
                                   rates_tree::FenwickTree,
                                   occupancies_n::Vector{Int},
                                   occupancies_c::Vector{Int},
                                   indexi::Int,
                                   indexf::Int,
                                   μ::Float64,
                                   attractive::Bool)
    get_new_chemo_rates!(rates, rates_tree, occupancies_n, occupancies_c, indexi, μ, attractive)
    get_new_chemo_rates!(rates, rates_tree, occupancies_n, occupancies_c, indexf, μ, attractive)
    return nothing
end

# ----------------------------
# Elementary events
# ----------------------------

function evap!(occupancies::Vector{Int}, d::RandomDict, p_index::Int)
    index = d.d[p_index]
    occupancies[index] -= 1
    Base.delete!(d, p_index)
    return ComplexF64[-1 + 0im, 0 + 0im], index, index
end

function create!(occupancies_c::Vector{Int}, dict_c::RandomDict, dict_n::RandomDict, p_index::Int)
    index = dict_n.d[p_index]
    occupancies_c[index] += 1
    insert_value!(dict_c, index)
    return ComplexF64[1 + 0im, 0 + 0im], index, index
end

function diff!(occupancies::Vector{Int}, d::RandomDict, p_index::Int)
    index = d.d[p_index]
    occupancies[index] > 0 || error("Improbable diff at empty site $index")

    if rand() >= 0.5
        target = periodic_index(index + 1, length(occupancies))
        dcache = ComplexF64[0 + 0im, 1 + 0im]
    else
        target = periodic_index(index - 1, length(occupancies))
        dcache = ComplexF64[0 + 0im, -1 + 0im]
    end

    occupancies[index] -= 1
    occupancies[target] += 1
    Base.delete!(d, p_index)
    insert!(d, p_index, target) # CHECK how this mechanism works, check class to see if this line does the right thing. (checked)

    return dcache, index, target
end

function get_random_chemo_index(d::RandomDict, chemo_rate_tree::FenwickTree)
    idx = sample(chemo_rate_tree)
    site = div(idx - 1, 2) + 1
    direction = mod(idx - 1, 2) # 0 left, 1 right

    haskey(d.val_keys, site) || error("chemotactic jump from empty site $site")
    return site, direction
end

function chemo!(occupancies::Vector{Int}, d::RandomDict, chemo_rate_tree::FenwickTree)
    index, direction = get_random_chemo_index(d, chemo_rate_tree)

    if direction == 0
        target = periodic_index(index - 1, length(occupancies))
        dcache = ComplexF64[0 + 0im, -1im]
    else
        target = periodic_index(index + 1, length(occupancies))
        dcache = ComplexF64[0 + 0im, 1im]
    end

    occupancies[index] -= 1
    occupancies[target] += 1

    p_index = random_key_from_val(d, index)
    Base.delete!(d, p_index)
    insert!(d, p_index, target) 

    return dcache, index, target
end

# ----------------------------
# One Gillespie step
# ----------------------------

function timestep!(st::SimState, par::Params)
    total_chemo_rate1 = total_sum(st.chemo_rate_tree1)
    total_chemo_rate2 = total_sum(st.chemo_rate_tree2)

    total_rate = par.Dn1 * st.N1 +
                 par.Dn2 * st.N2 +
                 par.Dc  * st.Nc +
                 par.gamma1 * st.N1 +
                 par.gamma2 * st.N2 +
                 par.kappa  * st.Nc +
                 par.lambda1 * total_chemo_rate1 +
                 par.lambda2 * total_chemo_rate2

    r1 = rand()
    r2 = rand()
    dtau = log(1 / r1) / total_rate

    if r2 <= par.Dn1 * st.N1 / total_rate
        p_index, index = get_random(st.dict1)
        _, i0, i1 = diff!(st.occ1, st.dict1, p_index)
        get_new_chemo_rates_diff!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, i1, par.μ, true)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2) / total_rate
        p_index, index = get_random(st.dict2)
        _, i0, i1 = diff!(st.occ2, st.dict2, p_index)
        get_new_chemo_rates_diff!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, i1, par.μ, false)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2 + par.Dc * st.Nc) / total_rate
        p_index, index = get_random(st.dictc)
        _, i0, i1 = diff!(st.occc, st.dictc, p_index)
        get_new_chemo_rates_diff!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, i1, par.μ, true)
        get_new_chemo_rates_diff!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, i1, par.μ, false)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2 + par.Dc * st.Nc + par.gamma1 * st.N1) / total_rate
        p_index, index = get_random(st.dict1)
        dcache, i0, _ = create!(st.occc, st.dictc, st.dict1, p_index)
        st.cache .+= dcache
        st.Nc += 1
        get_new_chemo_rates!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, par.μ, true)
        get_new_chemo_rates!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, par.μ, false)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2 + par.Dc * st.Nc + par.gamma1 * st.N1 + par.gamma2 * st.N2) / total_rate
        p_index, index = get_random(st.dict2)
        dcache, i0, _ = create!(st.occc, st.dictc, st.dict2, p_index)
        st.cache .+= dcache
        st.Nc += 1
        get_new_chemo_rates!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, par.μ, true)
        get_new_chemo_rates!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, par.μ, false)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2 + par.Dc * st.Nc + par.gamma1 * st.N1 + par.gamma2 * st.N2 + par.kappa * st.Nc) / total_rate
        p_index, index = get_random(st.dictc)
        dcache, i0, _ = evap!(st.occc, st.dictc, p_index)
        st.cache .+= dcache
        st.Nc -= 1
        get_new_chemo_rates!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, par.μ, true)
        get_new_chemo_rates!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, par.μ, false)

    elseif r2 <= (par.Dn1 * st.N1 + par.Dn2 * st.N2 + par.Dc * st.Nc + par.gamma1 * st.N1 + par.gamma2 * st.N2 + par.kappa * st.Nc + par.lambda1 * total_chemo_rate1) / total_rate
        _, i0, i1 = chemo!(st.occ1, st.dict1, st.chemo_rate_tree1)
        get_new_chemo_rates_diff!(st.chemo_rates1, st.chemo_rate_tree1, st.occ1, st.occc, i0, i1, par.μ, true)

    else
        _, i0, i1 = chemo!(st.occ2, st.dict2, st.chemo_rate_tree2)
        get_new_chemo_rates_diff!(st.chemo_rates2, st.chemo_rate_tree2, st.occ2, st.occc, i0, i1, par.μ, false)
    end

    st.tau += dtau
    st.counter += 1
    return dtau
end

# ----------------------------
# Run
# ----------------------------

function snapshot(st::SimState)
    return vcat(permutedims(st.occ1), permutedims(st.occ2), permutedims(st.occc))
end

function run_sim!(st::SimState, par::Params)
    start_time = time()
    next_progress = 0.05

    while st.tau <= par.Tfinal
        timestep!(st, par)
        next_progress = report_progress(st.tau, par.Tfinal, next_progress, start_time)

        if st.counter % par.save_rate == 0
            push!(st.occ_history, snapshot(st))
            push!(st.tau_history, st.tau)
        end

        # rebuild the tree to take away floating point error (very rare) and prevent freak sampling event
        if st.counter % 3_000_000 == 0
            st.chemo_rate_tree1 = FenwickTree(st.chemo_rates1)
            st.chemo_rate_tree2 = FenwickTree(st.chemo_rates2)
        end
    end

    push!(st.occ_history, snapshot(st))
    push!(st.tau_history, st.tau)

    println("Number of timesteps completed: $(st.counter), Nc: $(st.Nc)")

    if par.save
        save_sim(st, par)
    end#

    return st
end

# ----------------------------
# Save / post-processing
# ----------------------------

function history_array(st::SimState)
    nt = length(st.occ_history)
    nfields, nsites = size(st.occ_history[1])
    out = Array{Int}(undef, nt, nfields, nsites)
    for t in 1:nt
        out[t, :, :] .= st.occ_history[t]
    end
    return out
end 

function save_sim(st::SimState, par::Params)
    outdir = par.output_dir
    if isdir(outdir)
        outdir *= "_new"
    end
    mkpath(outdir)

    serialize(joinpath(outdir, "Params.bin"), par)
    serialize(joinpath(outdir, "SimState.bin"), st)
    npzwrite(joinpath(outdir, "occupancies_t.npy"), history_array(st))
    npzwrite(joinpath(outdir, "tau_t.npy"), collect(st.tau_history))
    npzwrite(joinpath(outdir, "cache.npy"), st.cache)
    npzwrite(joinpath(outdir, "chemo_rates1.npy"), st.chemo_rates1)
    npzwrite(joinpath(outdir, "chemo_rates2.npy"), st.chemo_rates2)

    println("Saved to: $outdir")
end

function meso_avg(occupancies::AbstractVector, w::Int)
    n = length(occupancies)
    out = zeros(Float64, n)
    for i in 1:n
        s = 0.0
        for j in -w:w
            s += occupancies[periodic_index(i + j, n)]
        end
        out[i] = s / (2w + 1)
    end
    return out
end

end