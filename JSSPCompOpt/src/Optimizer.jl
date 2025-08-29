module Optimizer

using JuMP, GLPK
using ..Partitioning: JSSPGraph
using Catlab.ACSetInterface: parts

# Solve a single layer using MILP
function solve_layer(layer::Int, g::JSSPGraph)
    tasks = [t for t in parts(g, :Task) if g[t, :layer] == layer]
    con_arcs = [a for a in parts(g, :ConArc) if g[g[a, :src_con], :layer] == layer && g[g[a, :tgt_con], :layer] == layer]
    dis_arcs = [a for a in parts(g, :DisArc) if g[g[a, :src_dis], :layer] == layer && g[g[a, :tgt_dis], :layer] == layer]

    model = JuMP.Model(GLPK.Optimizer)
    @variable(model, start_time[t in tasks] >= 0)
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

# Test function
function test_optimizer()
    g = JSSPGraph{Float64,String,Int}()
    tasks = add_parts!(g, :Task, 2, task_name=["T1", "T2"], proc_time=[1.0, 1.0], layer=[1, 1])
    machines = add_parts!(g, :Machine, 1, machine_name=["M1"])
    set_subpart!(g, tasks, :task_machine, [1, 1])
    add_part!(g, :DisArc, src_dis=tasks[1], tgt_dis=tasks[2])
    add_part!(g, :DisArc, src_dis=tasks[2], tgt_dis=tasks[1])
    sol = solve_layer(1, g)
    println("Test Optimizer - Solution for layer 1: $sol")
end

export solve_layer, test_optimizer

end # module