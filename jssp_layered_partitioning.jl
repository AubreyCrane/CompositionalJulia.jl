module JSSPLayeredPartitioning

using Catlab
using Catlab.CategoricalAlgebra
using Graphs
using JuMP
using GLPK
using Random

# === Schema Definition ===
# Define the JSSP schema with task names and layer attribute
Catlab.@present JSSPSchema(FreeSchema) begin
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

Catlab.@acset_type JSSPGraph(JSSPSchema, index=[:src_con, :tgt_con, :src_dis, :tgt_dis, :task_machine])

# === Graph Creation ===
"""
    create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed) -> JSSPGraph

Create a JSSP instance with the specified number of jobs, tasks per job, and machines.
- `num_jobs`: Number of jobs.
- `tasks_per_job`: Number of tasks per job.
- `num_machines`: Number of machines.
- `seed`: Random seed for reproducibility.

Returns a JSSPGraph instance with conjunctive and disjunctive edges.
"""
function create_jssp_instance(num_jobs::Int, tasks_per_job::Int, num_machines::Int, seed::Int=123)
    Random.seed!(seed)
    g = JSSPGraph{Float64, String, Int}()

    # Creating tasks with processing times and names
    num_tasks = num_jobs * tasks_per_job
    task_names = ["J$(j)T$(t)" for j in 1:num_jobs for t in 1:tasks_per_job]
    proc_times = rand(1.0:0.1:5.0, num_tasks)
    task_ids = Catlab.add_parts!(g, :Task, num_tasks,
        proc_time=proc_times,
        task_name=task_names,
        layer=fill(1, num_tasks))

    # Creating machines
    machine_ids = Catlab.add_parts!(g, :Machine, num_machines, machine_name=["M$i" for i in 1:num_machines])
``
    # Assigning tasks to machines randomly
    for task in task_ids
        Catlab.set_subpart!(g, task, :task_machine, machine_ids[rand(1:num_machines)])
    end

    # Adding conjunctive arcs (job precedence)
    for j in 1:num_jobs
        for t in 1:(tasks_per_job-1)
            task_idx = (j-1) * tasks_per_job + t
            Catlab.add_part!(g, :ConArc, src_con=task_ids[task_idx], tgt_con=task_ids[task_idx+1])
        end
    end

    # Adding disjunctive arcs (machine conflicts)
    for t1 in Catlab.parts(g, :Task), t2 in Catlab.parts(g, :Task)
        if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
            Catlab.add_part!(g, :DisArc, src_dis=t1, tgt_dis=t2)
            Catlab.add_part!(g, :DisArc, src_dis=t2, tgt_dis=t1)
        end
    end

    return g
end

# === Partitioning ===
"""
    compute_component_graph(g::JSSPGraph) -> Graphs.DiGraph

Construct the layer-DAG C(l) with edges i->j where there exists a conjunctive arc (u,v) with l(u)=i < l(v)=j,
or a crossing disjunctive conflict that requires layer separation.
"""
function compute_component_graph(g::JSSPGraph)
    k = maximum(g[:layer])
    C = Graphs.DiGraph(k)
    # Adding edges from conjunctive arcs
    for arc in Catlab.parts(g, :ConArc)
        u = g[arc, :src_con]
        v = g[arc, :tgt_con]
        i, j = g[u, :layer], g[v, :layer]
        if i < j
            Graphs.add_edge!(C, i, j)
        end
    end
    # Detecting crossing disjunctive arcs or cycles
    for arc1 in Catlab.parts(g, :DisArc)
        u1, v1 = g[arc1, :src_dis], g[arc1, :tgt_dis]
        i1, j1 = g[u1, :layer], g[v1, :layer]
        if i1 == j1
            for arc2 in Catlab.parts(g, :DisArc)
                if arc1 != arc2
                    u2, v2 = g[arc2, :src_dis], g[arc2, :tgt_dis]
                    i2, j2 = g[u2, :layer], g[v2, :layer]
                    if i2 == j2 && has_crossing_conflict(g, u1, v1, u2, v2)
                        Graphs.add_edge!(C, i1, j1)
                    end
                end
            end
        end
    end
    return C
end

"""
    has_crossing_conflict(g::JSSPGraph, u1::Int, v1::Int, u2::Int, v2::Int) -> Bool

Check if two disjunctive arc pairs (u1,v1) and (u2,v2) form a crossing or cycle.
"""
function has_crossing_conflict(g::JSSPGraph, u1::Int, v1::Int, u2::Int, v2::Int)
    if g[u1, :task_machine] == g[u2, :task_machine] && g[u1, :task_machine] == g[v1, :task_machine]
        dis_arcs = [(u1, v1), (v1, u1), (u2, v2), (v2, u2)]
        visited = Set{Int}()
        rec_stack = Set{Int}()
        for (s, t) in dis_arcs
            if !in(s, visited) && dfs_cycle_detect_subgraph(g, s, t, visited, rec_stack, dis_arcs)
                return true
            end
        end
    end
    return false
end

"""
    dfs_cycle_detect_subgraph(g::JSSPGraph, v::Int, target::Int, visited::Set{Int}, rec_stack::Set{Int}, edges::Vector{Tuple{Int,Int}}) -> Bool

Helper function to detect cycles in a subgraph defined by disjunctive edges.
"""
function dfs_cycle_detect_subgraph(g::JSSPGraph, v::Int, target::Int, visited::Set{Int}, rec_stack::Set{Int}, edges::Vector{Tuple{Int,Int}})
    push!(visited, v)
    push!(rec_stack, v)
    for (s, t) in edges
        if s == v && !in(t, visited)
            if dfs_cycle_detect_subgraph(g, t, target, visited, rec_stack, edges)
                return true
            end
        elseif s == v && t == target && in(target, rec_stack)
            return true
        end
    end
    pop!(rec_stack, v, 0)
    return false
end

"""
    has_conflict(g::JSSPGraph, u::Int, v::Int) -> Bool

Check if two tasks u and v conflict (e.g., are assigned to the same machine).
"""
function has_conflict(g::JSSPGraph, u::Int, v::Int)
    return g[u, :task_machine] == g[v, :task_machine]
end

"""
    refine_partition!(g::JSSPGraph, conflict_arc::Int, granularity_factor::Float64=1.0)

Refine the layer assignment for tasks involved in a conflict, with adjustable granularity.
"""
function refine_partition!(g::JSSPGraph, conflict_arc::Int, granularity_factor::Float64=1.0)
    u = g[conflict_arc, :src_con]
    v = g[conflict_arc, :tgt_con]
    machine = g[u, :task_machine]
    current_layer = g[u, :layer]
    
    if has_conflict(g, u, v) && g[u, :layer] == g[v, :layer]
        base_new_layer = maximum(g[:layer])
        conflicting_tasks = Set{Int}()
        for t in Catlab.parts(g, :Task)
            if g[t, :task_machine] == machine && has_conflict(g, u, t)
                push!(conflicting_tasks, t)
            end
        end
        num_conflicts = length(conflicting_tasks)
        new_layer = base_new_layer + floor(Int, granularity_factor * num_conflicts)
        
        for t in conflicting_tasks
            if g[t, :layer] == current_layer
                Catlab.set_subpart!(g, t, :layer, new_layer)
            end
        end
        println("Refined layer $current_layer -> $new_layer for machine $(g[machine, :machine_name]) due to conflict involving tasks $u and $v")
    end
end

"""
    layered_partitioning(g::JSSPGraph, max_ops_per_layer::Int=Inf, granularity_factor::Float64=1.0) -> Tuple{Vector{Int},Int}

Perform layered partitioning on the JSSP graph with flexible granularity.
"""
function layered_partitioning(g::JSSPGraph, max_ops_per_layer::Int=Inf, granularity_factor::Float64=1.0)
    for t in Catlab.parts(g, :Task)
        Catlab.set_subpart!(g, t, :layer, 1)
    end

    max_iterations = 10
    resolved_conflicts = Set{Tuple{Int,Int}}()
    last_layer_count = 0
    for iter in 1:max_iterations
        C = compute_component_graph(g)
        backward = find_backward_edge(C)
        if isnothing(backward)
            break
        end
        found_conflict = false
        for arc in Catlab.parts(g, :ConArc)
            u, v = g[arc, :src_con], g[arc, :tgt_con]
            if g[u, :layer] > g[v, :layer]
                refine_partition!(g, arc, granularity_factor)
                found_conflict = true
                break
            end
        end
        if !found_conflict
            for arc1 in Catlab.parts(g, :DisArc)
                u, v = g[arc1, :src_dis], g[arc1, :tgt_dis]
                conflict_pair = (min(u, v), max(u, v))
                if g[u, :layer] == g[v, :layer] && !(conflict_pair in resolved_conflicts)
                    for arc2 in Catlab.parts(g, :DisArc)
                        if arc1 != arc2
                            u2, v2 = g[arc2, :src_dis], g[arc2, :tgt_dis]
                            if has_crossing_conflict(g, u, v, u2, v2)
                                refine_partition!(g, arc1, granularity_factor)
                                push!(resolved_conflicts, conflict_pair)
                                found_conflict = true
                                break
                            end
                        end
                    end
                    if found_conflict
                        break
                    end
                end
            end
        end
        current_layer_count = maximum(g[:layer])
        if current_layer_count == last_layer_count
            break
        end
        last_layer_count = current_layer_count
        if iter == max_iterations
            println("Warning: Reached maximum iterations without resolving all backward edges.")
        end
    end

    k = maximum(g[:layer])
    return g[:layer], k
end

"""
    find_backward_edge(C::Graphs.DiGraph) -> Union{Tuple{Int,Int}, Nothing}

Find a backward edge in a directed graph (returns a tuple (u, v) or nothing if acyclic).
"""
function find_backward_edge(C::Graphs.DiGraph)
    for v in Graphs.vertices(C)
        for w in Graphs.outneighbors(C, v)
            visited = Set{Int}()
            stack = [w]
            while !isempty(stack)
                u = pop!(stack)
                if u == v
                    return (v, w)
                end
                if !(u in visited)
                    push!(visited, u)
                    append!(stack, Graphs.outneighbors(C, u))
                end
            end
        end
    end
    return nothing
end

# === Optimization ===
"""
    solve_layer(layer::Int, g::JSSPGraph) -> Dict{Tuple{Int,Int},Tuple{Int,Int}}

Solve the induced subgraph G[l^-1(layer)] using MILP, returning oriented disjunctive edges.
"""
function solve_layer(layer::Int, g::JSSPGraph)
    tasks = [t for t in Catlab.parts(g, :Task) if g[t, :layer] == layer]
    con_arcs = [a for a in Catlab.parts(g, :ConArc) if g[g[a, :src_con], :layer] == layer && g[g[a, :tgt_con], :layer] == layer]
    dis_arcs = [a for a in Catlab.parts(g, :DisArc) if g[g[a, :src_dis], :layer] == layer && g[g[a, :tgt_dis], :layer] == layer]

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
    if JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        return Dict((g[a, :src_dis], g[a, :tgt_dis]) => JuMP.value(x[a]) > 0.5 ? (g[a, :src_dis], g[a, :tgt_dis]) : (g[a, :tgt_dis], g[a, :src_dis]) for a in dis_arcs)
    else
        error("No optimal solution for layer $layer")
    end
end

# === Integration ===
"""
    compose_solutions(g::JSSPGraph, local_sols::Vector{Dict}) -> Graphs.DiGraph

Combine intra-layer orientations with conjunctive edges, orienting cross-layer edges forward.
"""
function compose_solutions(g::JSSPGraph, local_sols::Vector{Dict{Tuple{Int,Int},Tuple{Int,Int}}})
    n = Catlab.nparts(g, :Task)
    Gσ = Graphs.DiGraph(n)
    for arc in Catlab.parts(g, :ConArc)
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

# === Utilities ===
"""
    is_acyclic(G::Graphs.DiGraph) -> Bool

Check if the directed graph G is acyclic, ensuring a valid schedule.
"""
function is_acyclic(G::Graphs.DiGraph)
    return !Graphs.is_cyclic(G)
end

# === Testing ===
"""
    test_create_jssp()

Test the creation of a JSSP instance with a small example.
"""
function test_create_jssp(num_jobs=2, tasks_per_job=3, num_machines=3, seed=123)
    g = create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed)
    println("Test Create JSSP - Number of tasks: $(Catlab.nparts(g, :Task)), Number of machines: $(Catlab.nparts(g, :Machine))")
end

"""
    test_layered_partitioning()

Test the layered partitioning function with a small example.
"""
function test_layered_partitioning(max_ops_per_layer=typemax(Int), num_jobs=2, tasks_per_job=3, num_machines=3, seed=123)
    g = create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed)
    ℓ, k = layered_partitioning(g, max_ops_per_layer)
    println("Test Layered Partitioning - Layers: $ℓ, Number of layers: $k")
end

"""
    test_solve_layer()

Test solving a single layer with a small example.
"""
function test_solve_layer(max_ops_per_layer=typemax(Int), num_jobs=2, tasks_per_job=3, num_machines=3, seed=123)
    g = create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed)
    layered_partitioning(g, max_ops_per_layer)
    sol = solve_layer(1, g)
    println("Test Solve Layer - Solution for layer 1: $sol")
end

"""
    test_compose_solutions()

Test composing solutions from layers with a small example.
"""
function test_compose_solutions(max_ops_per_layer=typemax(Int), num_jobs=2, tasks_per_job=3, num_machines=3, seed=123)
    g = create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed)
    layered_partitioning(g, max_ops_per_layer)
    local_sols = [solve_layer(i, g) for i in 1:maximum(g[:layer])]
    Gσ = compose_solutions(g, local_sols)
    println("Test Compose Solutions - Global orientation is acyclic: $(is_acyclic(Gσ))")
end

# Exporting functions for use in run_experiment.jl
export create_jssp_instance, layered_partitioning, solve_layer, compose_solutions, is_acyclic, test_create_jssp, test_layered_partitioning, test_solve_layer, test_compose_solutions

end # module