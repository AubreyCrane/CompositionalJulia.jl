module Optimizer

using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphs
using JuMP: Model, @variable, @constraint, @objective, optimize!, value, termination_status
using GLPK

using Random

# === Schema ===
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
    task_name::Attr(Task, String)
    layer::Attr(Task, Int)  # Used later by Partitioning.jl

    machine_name::Attr(Machine, String)
end

@acset_type JSSPGraph(JSSPSchema, index=[:src_con, :tgt_con, :src_dis, :tgt_dis, :task_machine])

# === Graph Construction ===
function create_scaled_jssp(; num_jobs=6, tasks_per_job=3, num_machines=4, seed=123)
    Random.seed!(seed)
    g = JSSPGraph{Float64,String,Int}()

    # Tasks
    num_tasks = num_jobs * tasks_per_job
    task_names = ["J$(j)T$(t)" for j in 1:num_jobs for t in 1:tasks_per_job]
    proc_times = rand(1.0:0.1:5.0, num_tasks)
    task_ids = add_parts!(g, :Task, num_tasks;
                          proc_time=proc_times,
                          task_name=task_names,
                          layer=fill(1, num_tasks))

    # Machines
    machine_ids = add_parts!(g, :Machine, num_machines;
                             machine_name=["M$(i)" for i in 1:num_machines])

    # Assign tasks to machines randomly
    for t in task_ids
        set_subpart!(g, t, :task_machine, machine_ids[rand(1:num_machines)])
    end

    # Conjunctive arcs (job order)
    for j in 1:num_jobs
        for t in 1:(tasks_per_job-1)
            task_idx = (j-1)*tasks_per_job + t
            add_part!(g, :ConArc, src_con=task_ids[task_idx], tgt_con=task_ids[task_idx+1])
        end
    end

    # Disjunctive arcs (machine conflicts)
    for t1 in task_ids, t2 in task_ids
        if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
            add_part!(g, :DisArc, src_dis=t1, tgt_dis=t2)
            add_part!(g, :DisArc, src_dis=t2, tgt_dis=t1)
        end
    end

    return g
end

# === Solver ===
function solve_jssp(g::JSSPGraph)
    tasks = parts(g, :Task)
    model = Model(GLPK.Optimizer)

    @variable(model, start_time[tasks] >= 0)

    # Binary variables for machine ordering
    x = Dict{Tuple{Int,Int}, VariableRef}()
    for t1 in tasks, t2 in tasks
        if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
            x[(t1, t2)] = @variable(model, binary=true)
        end
    end

    # Precedence constraints
    for arc in parts(g, :ConArc)
        s, t = g[arc, :src_con], g[arc, :tgt_con]
        @constraint(model, start_time[t] >= start_time[s] + g[s, :proc_time])
    end

    # Disjunctive constraints
    M = 1000.0
    for (t1, t2) in keys(x)
        @constraint(model, start_time[t1] + g[t1, :proc_time] <= start_time[t2] + M * (1 - x[(t1, t2)]))
        @constraint(model, start_time[t2] + g[t2, :proc_time] <= start_time[t1] + M * x[(t1, t2)])
    end

    @variable(model, C_max >= 0)
    for t in tasks
        @constraint(model, C_max >= start_time[t] + g[t, :proc_time])
    end
    @objective(model, Min, C_max)

    optimize!(model)
    return value(C_max), Dict(t => value(start_time[t]) for t in tasks)
end

# === Visualization (Text-based) ===
function visualize_jssp(g::JSSPGraph, solution::Tuple{Float64, Dict{Int, Float64}})
    println("=== JSSP Instance ===")
    for t in parts(g, :Task)
        println("Task $(g[t,:task_name]): Machine=$(g[g[t,:task_machine], :machine_name]), Layer=$(g[t,:layer])")
    end
    println("\nMakespan = ", round(solution[1], digits=2))
end

export JSSPGraph, create_scaled_jssp, solve_jssp, visualize_jssp

end # module
