using JuMP, GLPK

# Define a struct for JSSP instances (unchanged)
struct JSSPInstance
    tasks::Vector{String}               # List of task IDs
    conjunctive_arcs::Vector{Tuple{String, String}}  # Precedence constraints (u, v)
    machine_tasks::Dict{String, Vector{String}}      # Machine to tasks requiring it
    processing_times::Dict{String, Float64}          # Task to processing time
end

# Function to compose two JSSP instances along shared machines (unchanged)
function compose_jssp(inst1::JSSPInstance, inst2::JSSPInstance, shared_machines::Vector{String})
    tasks = unique(vcat(inst1.tasks, inst2.tasks))
    conjunctive_arcs = unique(vcat(inst1.conjunctive_arcs, inst2.conjunctive_arcs))
    machine_tasks = Dict{String, Vector{String}}()
    all_machines = union(keys(inst1.machine_tasks), keys(inst2.machine_tasks))
    for machine in all_machines
        tasks1 = get(inst1.machine_tasks, machine, String[])
        tasks2 = get(inst2.machine_tasks, machine, String[])
        machine_tasks[machine] = unique(vcat(tasks1, tasks2))
    end
    processing_times = merge(inst1.processing_times, inst2.processing_times)
    return JSSPInstance(tasks, conjunctive_arcs, machine_tasks, processing_times)
end

# Corrected function to solve a JSSP instance using MILP
function solve_jssp(inst::JSSPInstance)
    model = Model(GLPK.Optimizer)
    
    # Variables: start times for each task
    @variable(model, t[inst.tasks] >= 0)
    
    # Binary variables for disjunctive constraints, stored in a dictionary
    # CHANGE: Use anonymous variables to avoid naming conflicts
    x = Dict{Tuple{String, String}, VariableRef}()
    for (machine, tasks) in inst.machine_tasks
        if length(tasks) > 1
            for i in 1:length(tasks)-1
                for j in i+1:length(tasks)
                    u, v = tasks[i], tasks[j]
                    x[(u, v)] = @variable(model, binary=true)  # Anonymous binary variable
                end
            end
        end
    end
    
    # Constraints: precedence
    for (u, v) in inst.conjunctive_arcs
        @constraint(model, t[v] >= t[u] + inst.processing_times[u])
    end
    
    # Constraints: disjunctive (non-overlapping tasks on same machine)
    M = 1000.0  # Big-M value; adjust if needed
    for (machine, tasks) in inst.machine_tasks
        if length(tasks) > 1
            for i in 1:length(tasks)-1
                for j in i+1:length(tasks)
                    u, v = tasks[i], tasks[j]
                    # CHANGE: Use x[(u, v)] instead of x[u, v]
                    @constraint(model, t[u] + inst.processing_times[u] <= t[v] + M * (1 - x[(u, v)]))
                    @constraint(model, t[v] + inst.processing_times[v] <= t[u] + M * x[(u, v)])
                end
            end
        end
    end
    
    # Objective: minimize makespan
    @variable(model, C_max >= 0)
    for task in inst.tasks
        @constraint(model, C_max >= t[task] + inst.processing_times[task])
    end
    @objective(model, Min, C_max)
    
    # Solve
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        return value(C_max), Dict(task => value(t[task]) for task in inst.tasks)
    else
        error("No optimal solution found.")
    end
end

# Test with small instances (unchanged)
inst1 = JSSPInstance(
    ["T1", "T2"],
    [("T1", "T2")],
    Dict("M1" => ["T1", "T2"]),
    Dict("T1" => 3.0, "T2" => 2.0)
)

inst2 = JSSPInstance(
    ["T3", "T4"],
    [("T3", "T4")],
    Dict("M2" => ["T3", "T4"]),
    Dict("T3" => 4.0, "T4" => 1.0)
)

# Compose instances (no shared machines for simplicity)
composed = compose_jssp(inst1, inst2, String[])
makespan, start_times = solve_jssp(composed)

println("Makespan: ", makespan)
println("Start times: ", start_times)



###### Disjunctive Graph Construction Visual ######


using Graphs

function build_disjunctive_graph(instance::JSSPInstance)
    task_ids = instance.tasks
    task_to_id = Dict(task => i for (i, task) in enumerate(task_ids))
    G = SimpleGraph(length(task_ids))

    for tasks in values(instance.machine_tasks)
        for i in 1:length(tasks)-1
            for j in i+1:length(tasks)
                add_edge!(G, task_to_id[tasks[i]], task_to_id[tasks[j]])
            end
        end
    end

    return G, task_ids
end

function print_graph_edges(G::SimpleGraph, task_ids::Vector{String}, machine_tasks::Dict{String, Vector{String}})
    println("\n========== Disjunctive Graph (Edge List) ==========")
    for e in edges(G)
        u = task_ids[src(e)]
        v = task_ids[dst(e)]
        # FIX: collect keys into a vector to make them iterable
        machine = findfirst(m -> (u in machine_tasks[m] && v in machine_tasks[m]), collect(keys(machine_tasks)))
        println("  $u ⋄ $v  (machine: $machine)")
    end
    println("===================================================\n")
end

G, task_ids = build_disjunctive_graph(composed)
print_graph_edges(G, task_ids, composed.machine_tasks)
