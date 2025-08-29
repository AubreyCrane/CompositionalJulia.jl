module JSSPLayeredPartitioning

using Graphs
using GraphPlot
using JuMP
using GLPK
using Random
using Colors
using Cairo
using Fontconfig
using Compose
using CairoMakie

# Local partition for subsets
struct LocalPartition
    k::Int
    ℓ_S::Vector{Int}
end

# Add jobs list to graph-like structure
struct JSSPGraph
    V::Vector{Int}
    E_c::Vector{Tuple{Int,Int}}
    E_d::Vector{Tuple{Int,Int}}
    machine::Vector{Int}
    jobs::Vector{Vector{Int}}
end

function compute_component_graph(local_part::LocalPartition, edges_S::Vector{Tuple{Int,Int}})
    C = DiGraph(local_part.k)
    for (a, b) in edges_S
        i, j = local_part.ℓ_S[a], local_part.ℓ_S[b]
        if i != j
            add_edge!(C, i, j)
        end
    end
    return C
end

function find_backward_edge(C::DiGraph)
    for e in edges(C)
        i, j = src(e), dst(e)
        if i > j
            return (i, j)
        end
    end
    return nothing
end

function refine_partition!(local_part::LocalPartition, conflict_local::Tuple{Int,Int}, G, S)
    a, b = conflict_local
    old_layer = local_part.ℓ_S[b]
    new_layer = local_part.k + 1
    machine_b = G.machine[S[b]]
    T = [i for i in 1:length(S) if G.machine[S[i]] == machine_b]
    for i in T
        local_part.ℓ_S[i] = new_layer
    end
    setfield!(local_part, :k, new_layer)
    println("Refined layer $old_layer -> $new_layer for machine $machine_b in subset $S")
end

function solve_layer(layer::Int, ℓ::Vector{Int}, G)
    V_layer = [v for v in 1:length(ℓ) if ℓ[v] == layer]
    E_d_layer = [(u, v) for (u, v) in G.E_d if ℓ[u] == layer && ℓ[v] == layer]
    E_c_layer = [(u, v) for (u, v) in G.E_c if ℓ[u] == layer && ℓ[v] == layer]
    if isempty(E_d_layer)
        return Dict{Tuple{Int,Int},Int}()
    end
    model = Model(GLPK.Optimizer)
    edge_indices = Dict((u, v) => i for (i, (u, v)) in enumerate(E_d_layer))
    @variable(model, x[1:length(E_d_layer)], Bin)
    @variable(model, rank[V_layer] >= 0)
    for (u, v) in E_c_layer
        @constraint(model, rank[u] + 1 <= rank[v])
    end
    for (idx, (u, v)) in enumerate(E_d_layer)
        @constraint(model, rank[u] + 1 <= rank[v] + 1000 * (1 - x[idx]))
        @constraint(model, rank[v] + 1 <= rank[u] + 1000 * x[idx])
    end
    optimize!(model)
    orientation = Dict{Tuple{Int,Int},Int}()
    for (uv, idx) in edge_indices
        dir = value(x[idx]) > 0.5 ? 1 : -1
        orientation[uv] = dir
    end
    return orientation
end

function compose_solutions(local_sols::Vector{Dict{Tuple{Int,Int},Int}}, ℓ::Vector{Int}, G)
    Gσ = DiGraph(length(ℓ))
    for (u, v) in G.E_c
        add_edge!(Gσ, u, v)
    end
    for layer_sols in local_sols
        for ((u, v), dir) in layer_sols
            if dir == 1
                add_edge!(Gσ, u, v)
            else
                add_edge!(Gσ, v, u)
            end
        end
    end
    return Gσ
end

function compute_levels(S, E_c, idx)
    n_S = length(S)
    in_degree = zeros(Int, n_S)
    for (u, v) in E_c
        if u in S && v in S
            in_degree[idx[v]] += 1
        end
    end
    queue = [i for i in 1:n_S if in_degree[i] == 0]
    levels = fill(0, n_S)
    current_level = 1
    while !isempty(queue)
        next_queue = Int[]
        for i in queue
            levels[i] = current_level
            for j in 1:n_S
                if (S[i], S[j]) in E_c
                    in_degree[j] -= 1
                    if in_degree[j] == 0
                        push!(next_queue, j)
                    end
                end
            end
        end
        queue = next_queue
        current_level += 1
    end
    return levels
end

function partition_subset(S::Vector{Int}, current_layer::Int, ℓ::Vector{Int}, G, max_ops::Int)
    if length(S) <= max_ops
        for v in S
            ℓ[v] = current_layer
        end
        println("Assigned subset $S to layer $current_layer")
        return current_layer + 1
    else
        n_S = length(S)
        idx = Dict(S[i] => i for i in 1:n_S)
        inv_idx = Dict(i => S[i] for i in 1:n_S)
        levels = compute_levels(S, G.E_c, idx)
        max_level = maximum(levels)
        local_layers = [[inv_idx[i] for i in 1:n_S if levels[i] == j] for j in 1:max_level]
        for layer_ops in local_layers
            if length(layer_ops) <= max_ops
                for v in layer_ops
                    ℓ[v] = current_layer
                end
                println("Assigned subset $layer_ops to layer $current_layer")
                current_layer += 1
            else
                num_groups = ceil(Int, length(layer_ops) / max_ops)
                for g in 1:num_groups
                    group = layer_ops[(g-1)*max_ops + 1 : min(g*max_ops, length(layer_ops))]
                    for v in group
                        ℓ[v] = current_layer
                    end
                    println("Assigned subset $group to layer $current_layer")
                    current_layer += 1
                end
            end
        end
        return current_layer
    end
end

function layered_partitioning(G; max_ops=Inf, verbose=true)
    n = length(G.machine)
    ℓ = zeros(Int, n)
    V = collect(1:n)
    if verbose
        println("Starting partitioning with max_ops=$max_ops")
    end
    current_layer = partition_subset(V, 1, ℓ, G, max_ops)
    k = maximum(ℓ)
    if verbose
        println("Final partition: k=$k, layers=$ℓ")
    end
    local_sols = [solve_layer(i, ℓ, G) for i in 1:k]
    Gσ = compose_solutions(local_sols, ℓ, G)
    if verbose
        println("Acyclicity: $(!is_cyclic(Gσ) ? "Acyclic" : "Cyclic")")
    end
    return ℓ, k, Gσ
end

"""
Existing graph visualizer (GraphPlot)
"""
function visualize_graph(G, ℓ=nothing; filename="jssp_graph.svg")
    g = DiGraph(length(G.machine))
    for (u, v) in G.E_c
        add_edge!(g, u, v)
    end
    for (u, v) in G.E_d
        add_edge!(g, u, v)
        add_edge!(g, v, u)
    end
    node_labels = ["Op$u (M$(G.machine[u]))" for u in 1:length(G.machine)]
    nodefillc = ℓ === nothing ? fill(RGB(0.7,0.7,0.9), length(G.machine)) :
        [rand(RGB) for _ in 1:length(G.machine)]
    if ℓ !== nothing
        unique_colors = [rand(RGB) for _ in 1:maximum(ℓ)]
        nodefillc = [unique_colors[ℓ[i]] for i in 1:length(ℓ)]
    end
    p = gplot(g, nodelabel=node_labels, nodefillc=nodefillc, layout=spring_layout)
    if endswith(lowercase(filename), ".png")
        draw(PNG(filename, 800px, 800px), p)
    else
        draw(SVG(filename, 800px, 800px), p)
    end
    println("Graph visualization saved to $filename")
end

"""
New disjunctive graph visualizer with proper layers and colors
"""
function visualize_disjunctive_graph(G::JSSPGraph, ℓ::Vector{Int}; filename="jssp_partitioned_graph.png")
    n = length(G.V)
    layers = maximum(ℓ)
    pos = Dict{Int, Tuple{Float32,Float32}}()
    for (j_idx, job) in enumerate(G.jobs)
        for (t_idx, v) in enumerate(job)
            pos[v] = (j_idx + 0.1*(t_idx-1), -ℓ[v])
        end
    end
    fig = Figure(size = (900, 600), backgroundcolor = :white)  # changed from resolution to size
    ax  = Axis(fig[1,1]; title="JSSP Disjunctive Graph", xlabel="Jobs", ylabel="Layers")
    cmap = cgrad(:viridis, layers)

    # Conjunctive arcs
    for (u,v) in G.E_c
        lines!(ax, [pos[u][1], pos[v][1]], [pos[u][2], pos[v][2]], color=:black, linewidth=2)
        # Use arrows! for arrowheads
        arrows!(ax, [Point(pos[u]...)], [Point(pos[v]...) .- Point(pos[u]...)], arrowsize=15, linewidth=2, color=:black)
    end

    # Disjunctive arcs
    for (u,v) in G.E_d
        if ℓ[u] == ℓ[v]
            lines!(ax, [pos[u][1], pos[v][1]], [pos[u][2], pos[v][2]], color=:red, linestyle=:dash)
        else
            from,to = ℓ[u] < ℓ[v] ? (u,v) : (v,u)
            lines!(ax, [pos[from][1], pos[to][1]], [pos[from][2], pos[to][2]], color=RGBAf(1,0,0,0.6), linestyle=:dash)
            arrows!(ax, [Point(pos[from]...)], [Point(pos[to]...) .- Point(pos[from]...)], arrowsize=15, color=RGBAf(1,0,0,0.6))
        end
    end

    # Nodes
    for v in G.V
        scatter!(ax, [pos[v][1]], [pos[v][2]], color=cmap[ℓ[v]], markersize=18)
        text!(ax, string("Op", v, "\nM", G.machine[v]), position=pos[v], align=(:center,:center), fontsize=12)
    end

    hidespines!(ax); hidedecorations!(ax)
    CairoMakie.save(filename, fig)
    println("✔ Disjunctive graph saved to $filename")
end


function run_test()
    jobs = [[1,2,3],[4,5,6]]
    G = JSSPGraph(
        1:6,
        [(1, 2), (2, 3), (4, 5), (5, 6)],
        [(1, 5), (2, 4), (3, 6)],
        [1, 2, 3, 2, 1, 3],
        jobs
    )
    println("Test case: 2 jobs × 3 ops, max_ops=2")
    ℓ, k, Gσ = layered_partitioning(G, max_ops=2)
    visualize_disjunctive_graph(G, ℓ, filename="jssp_partitioned_graph.png")
    return ℓ, k, Gσ
end

run_test()

end # module
