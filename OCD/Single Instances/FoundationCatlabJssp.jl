using Catlab, Catlab.CategoricalAlgebra, Catlab.Graphs

# Define the JSSP schema with task names
@present JSSPSchema(FreeSchema) begin
    Float64::AttrType
    String::AttrType
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
end

# Generate the ACSet type
@acset_type JSSPGraph(JSSPSchema, index=[:src_con, :tgt_con, :src_dis, :tgt_dis, :task_machine])

# Create a minimal JSSP instance
function create_minimal_jssp()
    g = JSSPGraph{Float64,String}()
    
    # Add tasks with processing times and names
    task_ids = add_parts!(g, :Task, 4, proc_time=[3.0, 2.0, 4.0, 1.0], task_name=["J1T1", "J1T2", "J2T1", "J2T2"])
    
    # Add machines with names
    machine_ids = add_parts!(g, :Machine, 2, machine_name=["M1", "M2"])
    
    # Assign tasks to machines
    set_subpart!(g, task_ids[1], :task_machine, machine_ids[1])  # J1T1 -> M1
    set_subpart!(g, task_ids[2], :task_machine, machine_ids[2])  # J1T2 -> M2
    set_subpart!(g, task_ids[3], :task_machine, machine_ids[1])  # J2T1 -> M1
    set_subpart!(g, task_ids[4], :task_machine, machine_ids[2])  # J2T2 -> M2
    
    # Add conjunctive arcs (job precedence)
    add_part!(g, :ConArc, src_con=task_ids[1], tgt_con=task_ids[2])  # J1T1 -> J1T2
    add_part!(g, :ConArc, src_con=task_ids[3], tgt_con=task_ids[4])  # J2T1 -> J2T2
    
    # Add disjunctive arcs (machine conflicts)
    add_part!(g, :DisArc, src_dis=task_ids[1], tgt_dis=task_ids[3])  # J1T1 <-> J2T1 on M1
    add_part!(g, :DisArc, src_dis=task_ids[3], tgt_dis=task_ids[1])  # J2T1 <-> J1T1 on M1
    add_part!(g, :DisArc, src_dis=task_ids[2], tgt_dis=task_ids[4])  # J1T2 <-> J2T2 on M2
    add_part!(g, :DisArc, src_dis=task_ids[4], tgt_dis=task_ids[2])  # J2T2 <-> J1T2 on M2
    
    return g
end

# Visualize the JSSP graph in the terminal
function visualize_jssp(g::JSSPGraph)
    println("=== JSSP Disjunctive Graph ===")
    
    # Tasks with processing times and machine assignments
    println("\nTasks:")
    for t in parts(g, :Task)
        name = g[t, :task_name]
        proc_time = g[t, :proc_time]
        machine_id = g[t, :task_machine]
        machine_name = g[machine_id, :machine_name]
        println("  $name: Processing Time = $proc_time, Machine = $machine_name")
    end
    
    # Conjunctive edges (directed)
    println("\nConjunctive Edges (Precedence):")
    # CHANGE: Iterate over ConArc indices to avoid BoundsError
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
    
    # Disjunctive edges (undirected, grouped by machine)
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
    
    println("============================")
end

# Test the instance and visualization
g = create_minimal_jssp()
visualize_jssp(g)