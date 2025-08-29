using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphs
using JuMP, GLPK
using Random

# Define the JSSP schema with task names
@present JSSPSchema(FreeSchema) begin
    Float64::AttrType
    String::AttrType

    # Define the objects for tuple G(V, Ec, Ed, M) with dummy src and tgt listed in the morphisms section
    Task::Ob
    Machine::Ob
    ConArc::Ob
    DisArc::Ob

    # Define the morphisms for conjunctive and disjunctive arcs
    src_con::Hom(ConArc, Task)
    tgt_con::Hom(ConArc, Task)
    src_dis::Hom(DisArc, Task)
    tgt_dis::Hom(DisArc, Task)
    task_machine::Hom(Task, Machine)
    
    # attributres for processing times and names the A in ACSet :P which adds data to the objects (where decorated cospans will come into play)
    proc_time::Attr(Task, Float64)
    machine_name::Attr(Machine, String)
    task_name::Attr(Task, String)
end

# Generate the ACSet type
@acset_type JSSPGraph(JSSPSchema, index=[:src_con, :tgt_con, :src_dis, :tgt_dis, :task_machine])

# Create a scaled JSSP instance (6 jobs, 3 tasks per job, 4 machines)
function create_scaled_jssp(; num_jobs=6, tasks_per_job=3, num_machines=4, seed=123)
    Random.seed!(seed) # Ensure reproducibility
    g = JSSPGraph{Float64,String}()
    
    # Adding tasks for each job
    num_tasks = num_jobs * tasks_per_job
    task_names = ["J$(j)T$(t)" for j in 1:num_jobs for t in 1:tasks_per_job]
    proc_times = rand(1.0:0.1:5.0, num_tasks) # Random processing times between 1.0 and 5.0
    task_ids = add_parts!(g, :Task, num_tasks, 
        proc_time=proc_times, 
        task_name=task_names)
    
    # Adding machines
    machine_names = ["M$(i)" for i in 1:num_machines]
    machine_ids = add_parts!(g, :Machine, num_machines, machine_name=machine_names)
    
    # Assigning tasks to machines randomly
    for (i, task) in enumerate(task_ids)
        machine_id = machine_ids[rand(1:num_machines)]
        set_subpart!(g, task, :task_machine, machine_id)
    end
    
    # Adding conjunctive arcs (task precedence within jobs)
    for j in 1:num_jobs
        for t in 1:(tasks_per_job-1)
            task_idx = (j-1) * tasks_per_job + t
            next_task_idx = task_idx + 1
            add_part!(g, :ConArc, src_con=task_ids[task_idx], tgt_con=task_ids[next_task_idx])
        end
    end
    
    # Adding disjunctive arcs (machine conflicts)
    for t1 in parts(g, :Task)
        for t2 in parts(g, :Task)
            if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
                add_part!(g, :DisArc, src_dis=t1, tgt_dis=t2)
                add_part!(g, :DisArc, src_dis=t2, tgt_dis=t1)
            end
        end
    end
    
    return g
end

# Solve the JSSP instance using MILP
function solve_jssp(g::JSSPGraph)
    model = JuMP.Model(GLPK.Optimizer)
    
    tasks = parts(g, :Task)
    
    @variable(model, start_time[tasks] >= 0)
    
    x = Dict{Tuple{Int,Int}, VariableRef}()
    for t1 in tasks
        for t2 in tasks
            if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
                x[(t1, t2)] = @variable(model, binary=true)
            end
        end
    end
    
    for arc in parts(g, :ConArc)
        src, tgt = g[arc, :src_con], g[arc, :tgt_con]
        @constraint(model, start_time[tgt] >= start_time[src] + g[src, :proc_time])
    end
    
    M = 1000.0
    for t1 in tasks
        for t2 in tasks
            if t1 < t2 && haskey(x, (t1, t2))
                @constraint(model, start_time[t1] + g[t1, :proc_time] <= start_time[t2] + M * (1 - x[(t1, t2)]))
                @constraint(model, start_time[t2] + g[t2, :proc_time] <= start_time[t1] + M * x[(t1, t2)])
            end
        end
    end
    
    @variable(model, C_max >= 0)
    for task in tasks
        @constraint(model, C_max >= start_time[task] + g[task, :proc_time])
    end
    @objective(model, Min, C_max)
    
    optimize!(model)
    
    if termination_status(model) == MOI.OPTIMAL
        return value(C_max), Dict(task => value(start_time[task]) for task in tasks)
    else
        error("No optimal solution found.")
    end
end

# Visualize the JSSP graph and solution in the terminal
function visualize_jssp(g::JSSPGraph, solution::Tuple{Float64, Dict{Int, Float64}}=(-1.0, Dict()))
    println("=== JSSP Disjunctive Graph ===")
    
    println("\nTasks:")
    for t in parts(g, :Task)
        name = g[t, :task_name]
        proc_time = g[t, :proc_time]
        machine_id = g[t, :task_machine]
        machine_name = g[machine_id, :machine_name]
        start_time = haskey(solution[2], t) ? round(solution[2][t], digits=2) : "N/A"
        println("  $name: Processing Time = $proc_time, Machine = $machine_name, Start Time = $start_time")
    end
    
    println("\nConjunctive Edges (Precedence):")
    con_arcs = parts(g, :ConArc)
    if isempty(con_arcs)
        println("  None")
    else
        for arc in con_arcs
            src = g[arc, :src_con]
            tgt = g[arc, :tgt_con]
            src_name = g[src, :task_name]
            tgt_name = g[tgt, :task_name]
            println("  $src_name → $tgt_name")
        end
    end
    
    println("\nDisjunctive Edges (Machine Conflicts):")
    machine_conflicts = Dict{String, Vector{Tuple{String,String}}}()
    for t1 in parts(g, :Task)
        for t2 in parts(g, :Task)
            if t1 < t2 && g[t1, :task_machine] == g[t2, :task_machine]
                m = g[g[t1, :task_machine], :machine_name]
                push!(get!(machine_conflicts, m, []), (g[t1, :task_name], g[t2, :task_name]))
            end
        end
    end
    if isempty(machine_conflicts)
        println("  None")
    else
        for (machine, pairs) in machine_conflicts
            println("  On $machine:")
            for (t1, t2) in pairs
                println("    $t1 -- $t2")
            end
        end
    end
    
    if solution[1] >= 0
        println("\nGantt Chart (Text-Based):")
        max_time = ceil(Int, solution[1])
        for t in parts(g, :Task)
            name = g[t, :task_name]
            start = haskey(solution[2], t) ? solution[2][t] : 0.0
            duration = g[t, :proc_time]
            timeline = repeat(" ", floor(Int, start)) * repeat("#", floor(Int, duration))
            println("  $name: [$timeline]")
        end
    end
    
    if solution[1] >= 0
        println("\nMakespan: ", round(solution[1], digits=2))
    end
    println("============================")
end

# Test the scaled instance, solve, and visualize
g = create_scaled_jssp()
solution = solve_jssp(g)
visualize_jssp(g, solution)