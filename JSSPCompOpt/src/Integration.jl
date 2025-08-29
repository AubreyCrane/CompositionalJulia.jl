module Integration

using Graphs
using ..Partitioning: JSSPGraph
using Catlab.ACSetInterface: nparts, parts

# Compose solutions from all layers
function compose_solutions(g::JSSPGraph, local_sols::Vector{Dict{Tuple{Int,Int},Tuple{Int,Int}}})
    n = nparts(g, :Task)
    Gσ = DiGraph(n)
    for arc in parts(g, :ConArc)
        u, v = g[arc, :src_con], g[arc, :tgt_con]
        add_edge!(Gσ, u, v)
    end
    for sol in local_sols
        for (edge, (from, to)) in sol
            add_edge!(Gσ, from, to)
        end
    end
    return Gσ
end

# Test function
function test_integration()
    g = JSSPGraph{Float64,String,Int}()
    tasks = add_parts!(g, :Task, 4, layer=[1,1,2,2])
    local_sols = [
        Dict((1,2) => (1,2)),  # Layer 1: 1→2
        Dict((3,4) => (3,4))   # Layer 2: 3→4
    ]
    Gσ = compose_solutions(g, local_sols)
    println("Test Integration - Global orientation edges: $(edges(Gσ))")
    println("Is acyclic: $(!is_cyclic(Gσ))")
end

export compose_solutions, test_integration

end # module