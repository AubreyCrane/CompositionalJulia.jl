module Partitioning

using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphs
using Graphs
import Graphs
using JuMP, GLPK
import JuMP
using Random

# Define the JSSP schema with layer attribute
@present JSSPSchema(FreeSchema) begin
    Float64::AttrType
    String::AttrType
    Int::AttrType

    Task::Ob
    Machine::Ob
    ConArc::Ob
    DisArc::Ob

    src_con::Hom(ConArc, Task)
    tgt_con::Hom(ConArc, Task)
    src_dis::Hom(DisArc, Task)
    tgt_dis::Hom(DisArc, Task)
    task_machine::Hom(Task, Machine)

    proc_time::Attr(Task, Float64)
    machine_name::Attr(Machine, String)
    task_name::Attr(Task, String)
    layer::Attr(Task, Int)
end

@acset_type JSSPGraph(JSSPSchema, index=[:src_con, :tgt_con, :src_dis, :tgt_dis, :task_machine])

"""
    create_scaled_jssp(; num_jobs, tasks_per_job, num_machines, seed) -> JSSPGraph

Create a scaled JSSP instance with specified parameters.
"""
function create_scaled_jssp(; num_jobs=6, tasks_per_job=3, num_machines=4, seed=123)
    Random.seed!(seed)
    g = JSSPGraph{Float64,String,Int}()

    num_tasks = num_jobs * tasks_per_job
    task_names = ["J$(j)T$(t)" for j in 1:num_jobs for t in 1:tasks_per_job]
    proc_times = rand(1.0:0.1:5.0, num_tasks)
    task_ids = add_parts!(g, :Task, num_tasks, 
        proc_time=proc_times, 
        task_name=task_names,
        layer=fill(1, num_tasks))

    machine_ids = add_parts!(g, :Machine, num_machines, machine_name=["M$i" for i in 1:num_machines])

    for task in task_ids
        set_subpart!(g, task, :task_machine, machine_ids[rand(1:num_machines)])
    end

    for j in 1:num_jobs
        for t in 1:(tasks_per_job-1)
            task_idx = (j-1) * tasks_per_job + t
            add_part!(g, :ConArc, src_con=task_ids[task_idx], tgt_con=task_ids[task_idx+1])
        end
    end

    for t1 in parts(g, :Task), t2 in parts(g, :Task)
        if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
            add_part!(g, :DisArc, src_dis=t1, tgt_dis=t2)
            add_part!(g, :DisArc, src_dis=t2, tgt_dis=t1)
        end
    end

    return g
end

"""
    compute_component_graph(g::JSSPGraph) -> DiGraph

Construct the layer-DAG C(ℓ) with edges i→j where ∃(u,v)∈E_c with ℓ(u)=i < ℓ(v)=j.
"""
function compute_component_graph(g::JSSPGraph)
    k = maximum(g[:layer])
    C = DiGraph(k)
    for arc in parts(g, :ConArc)
        u = g[arc, :src_con]
        v = g[arc, :tgt_con]
        i, j = g[u, :layer], g[v, :layer]
        if i < j
            Graphs.add_edge!(C, i, j)
        end
    end
    return C
end

"""
    find_backward_edge(C::DiGraph) -> Union{Tuple{Int,Int},Nothing}

Return the first edge (i→j) in C with i ≥ j, or nothing if none exist.
"""
function find_backward_edge(C::DiGraph)
    for e in Graphs.edges(C)
        i, j = Graphs.src(e), Graphs.dst(e)
        if i >= j
            return (i, j)
        end
    end
    return nothing
end

"""
    refine_partition!(g::JSSPGraph, conflict_arc::Int)

Split the layer of the source task in a conjunctive conflict, updating layers.
"""
function refine_partition!(g::JSSPGraph, conflict_arc::Int)
    u = g[conflict_arc, :src_con]
    v = g[conflict_arc, :tgt_con]
    old_layer = g[u, :layer]
    new_layer = maximum(g[:layer]) + 1
    machine = g[u, :task_machine]
    for t in parts(g, :Task)
        if g[t, :task_machine] == machine
            set_subpart!(g, t, :layer, new_layer)
        end
    end
    println("Refined layer $old_layer -> $new_layer for machine $(g[machine, :machine_name])")
end

"""
    solve_layer(layer::Int, g::JSSPGraph) -> Dict{Tuple{Int,Int},Tuple{Int,Int}}

Solve the induced subgraph G[ℓ⁻¹(layer)] using MILP, returning oriented disjunctive edges.
"""
function solve_layer(layer::Int, g::JSSPGraph)
    tasks = [t for t in parts(g, :Task) if g[t, :layer] == layer]
    con_arcs = [a for a in parts(g, :ConArc) if g[g[a, :src_con], :layer] == layer && g[g[a, :tgt_con], :layer] == layer]
    dis_arcs = [a for a in parts(g, :DisArc) if g[g[a, :src_dis], :layer] == layer && g[g[a, :tgt_dis], :layer] == layer]

    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, start_time[tasks] >= 0)
    x = Dict(a => @variable(model, binary=true) for a in dis_arcs)

    for arc in con_arcs
        src, tgt = g[arc, :src_con], g[arc, :tgt_con]
        @constraint(model, start_time[tgt] >= start_time[src] + g[src, :proc_time])
    end

    M = 1000.0
    for arc in dis_arcs
        u, v = g[arc, :src_dis], g[arc, :tgt_dis]
        @constraint(model, start_time[u] + g[u, :proc_time] <= start_time[v] + M * (1 - x[arc]))
        @constraint(model, start_time[v] + g[v, :proc_time] <= start_time[u] + M * x[arc])
    end

    @variable(model, C_max >= 0)
    for t in tasks
        @constraint(model, C_max >= start_time[t] + g[t, :proc_time])
    end
    @objective(model, Min, C_max)

    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return Dict((g[a, :src_dis], g[a, :tgt_dis]) => value(x[a]) > 0.5 ? (g[a, :src_dis], g[a, :tgt_dis]) : (g[a, :tgt_dis], g[a, :src_dis]) for a in dis_arcs)
    else
        error("No optimal solution for layer $layer")
    end
end

"""
    compose_solutions(g::JSSPGraph, local_sols::Vector{Dict}) -> DiGraph

Combine intra-layer orientations with conjunctive edges, orienting cross-layer edges forward.
"""
function compose_solutions(g::JSSPGraph, local_sols::Vector{Dict{Tuple{Int,Int},Tuple{Int,Int}}})
    n = nparts(g, :Task)
    Gσ = DiGraph(n)
    for arc in parts(g, :ConArc)
        u, v = g[arc, :src_con], g[arc, :tgt_con]
        Graphs.add_edge!(Gσ, u, v)
    end
    for sol in local_sols
        for (edge, (from, to)) in sol
            Graphs.add_edge!(Gσ, from, to)
        end
    end
    return Gσ
end

"""
    layered_partitioning(g::JSSPGraph) -> Tuple{Vector{Int},Int,DiGraph}

Run the recursive layered-partitioning algorithm, returning the partition, number of layers, and global orientation.
"""
function layered_partitioning(g::JSSPGraph)
    # Initialize all tasks in layer 1
    for t in parts(g, :Task)
        set_subpart!(g, t, :layer, 1)
    end

    # Driver loop
    while true
        C = compute_component_graph(g)
        backward = find_backward_edge(C)
        if isnothing(backward)
            break
        end
        # Find a conjunctive arc causing the backward edge
        for arc in parts(g, :ConArc)
            u, v = g[arc, :src_con], g[arc, :tgt_con]
            if g[u, :layer] > g[v, :layer]
                refine_partition!(g, arc)
                break
            end
        end
    end

    k = maximum(g[:layer])
    local_sols = [solve_layer(i, g) for i in 1:k]
    Gσ = compose_solutions(g, local_sols)
    return g[:layer], k, Gσ
end

"""
    run_test() -> Nothing

Test the algorithm on a 2-jobs × 3-ops example.
"""
function run_test()
    g = JSSPGraph{Float64,String,Int}()
    tasks = add_parts!(g, :Task, 6, 
        task_name=["J1T1", "J1T2", "J1T3", "J2T1", "J2T2", "J2T3"],
        proc_time=[1.0, 2.0, 1.5, 1.2, 1.8, 1.0],
        layer=fill(1, 6))
    machines = add_parts!(g, :Machine, 3, machine_name=["M1", "M2", "M3"])
    set_subpart!(g, tasks, :task_machine, [1, 2, 3, 2, 1, 3])

    add_part!(g, :ConArc, src_con=tasks[1], tgt_con=tasks[2])
    add_part!(g, :ConArc, src_con=tasks[2], tgt_con=tasks[3])
    add_part!(g, :ConArc, src_con=tasks[4], tgt_con=tasks[5])
    add_part!(g, :ConArc, src_con=tasks[5], tgt_con=tasks[6])

    for t1 in tasks, t2 in tasks
        if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
            add_part!(g, :DisArc, src_dis=t1, tgt_dis=t2)
            add_part!(g, :DisArc, src_dis=t2, tgt_dis=t1)
        end
    end

    ℓ, k, Gσ = layered_partitioning(g)
    println("Layers: $ℓ, Number of layers: $k")
    println("Global orientation is acyclic: $(!is_cyclic(Gσ))")
end

# Export functions for integration
export JSSPGraph, create_scaled_jssp, layered_partitioning, run_test

end # module

# Run the test
Partitioning.run_test()