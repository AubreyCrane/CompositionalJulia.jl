### layered_partition.jl
# Skeleton implementation of recursive layered partitioning for disjunctive JSSP
# using Catlab/Graphs or LightGraphs (replace with your chosen graph package).

# === Data Structures ===
# ℓ: Vector{Int} mapping each vertex index to a layer index
# Alternatively, Dict{Vertex,Int} for arbitrary Vertex types

struct Partition
    ℓ::Vector{Int}          # layer index per vertex id
    k::Int                  # number of layers
end

# Disjunctive edges as pairs of vertex indices
const DisEdge = Tuple{Int,Int}

# Component graph representation (could use LightGraphs.DiGraph)
using LightGraphs

# === Key Functions ===

"""
compute_component_graph(ℓ, edges) -> DiGraph
Constructs the layer-level DAG C(ℓ) on 1:k from disjunctive edges.
"""
function compute_component_graph(part::Partition, edges::Vector{DisEdge})
    k = part.k
    C = DiGraph(k)
    for (u,v) in edges
        i = part.ℓ[u]; j = part.ℓ[v]
        if i != j && i < j
            add_edge!(C, i, j)
        elseif i != j && j < i
            add_edge!(C, j, i)
        end
    end
    return C
end

"""
find_backward_edge(C) -> Edge or nothing
Scans C for any edge i->j with i >= j.
"""
function find_backward_edge(C::DiGraph)
    for e in edges(C)
        u,v = src(e), dst(e)
        if u >= v
            return (u,v)
        end
    end
    return nothing
end

"""
refine_partition!(part, conflict::DisEdge)
Splits the layer containing one endpoint of the conflict into a new layer.
"""
function refine_partition!(part::Partition, conflict::DisEdge)
    u,v = conflict
    # choose endpoint w to move (e.g., u)
    w = u
    # assign new layer index
    new_layer = part.k + 1
    push!(part.ℓ, 0) # extend ℓ if needed
    part.ℓ[w] = new_layer
    part.k += 1
    return part
end

"""
solve_layer(layer::Int, part, g)
Extracts subgraph for 'layer' and solves its local orientation.
Returns a Dict{DisEdge,Direction} mapping.
"""
function solve_layer(layer::Int, part::Partition, g)
    # TODO: extract nodes where ℓ[node]==layer
    # TODO: build local subgraph and call CP/MIP solver
    return Dict{DisEdge,Symbol}()  # e.g., :u_to_v or :v_to_u
end

"""
compose_solutions(local_sols::Vector, part::Partition)
Combines local orientations with forward cross-layer rules.
"""
function compose_solutions(local_sols, part::Partition)
    global_orient = Dict{DisEdge,Symbol}()
    # fill from local_sols for intra-layer
    # enforce u->v when ℓ[u]<ℓ[v]
    return global_orient
end

# === API Sketch for Catlab/Graphs ===
# - Tag ACSet vertices with 'layer' attribute
# - Use 'parts(g, :DisArc)' to iterate disjunctive arcs
# - filter by layer and rebuild sub-ACSet

# === Example Usage ===
# 1. Initialize Partition
# part = Partition(fill(1, nv), initial_k)
# 2. loop:
#    C = compute_component_graph(part, disedges)
#    conf = find_backward_edge(C)
#    isnothing(conf) && break
#    refine_partition!(part, conf)
# 3. for each layer i in 1:part.k: local_sols[i] = solve_layer(i, part, g)
# 4. global_orient = compose_solutions(local_sols, part)
