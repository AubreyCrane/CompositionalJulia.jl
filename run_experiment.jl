include("jssp_layered_partitioning.jl")
using .JSSPLayeredPartitioning
using Graphs
using GraphPlot
using Compose
using Cairo

# Importing compute_component_graph if defined in JSSPLayeredPartitioning
import .JSSPLayeredPartitioning: compute_component_graph

# Importing JSSPGraph type for type annotations in visualization functions
import .JSSPLayeredPartitioning: JSSPGraph

# === Experiment Parameters ===
num_jobs = 2
tasks_per_job = 3
num_machines = 2
seed = 123
max_ops_per_layer = typemax(Int)
visualize = true

# === Visualization Functions ===
"""
    visualize_jssp_graph(g::JSSPGraph, filename::String)

Visualize the initial JSSP graph with tasks, conjunctive arcs (solid directed edges), and disjunctive arcs (dashed undirected edges).
"""
function visualize_jssp_graph(g::JSSPGraph, filename::String)
    plot_graph = Graphs.SimpleDiGraph(Catlab.nparts(g, :Task))
    for arc in Catlab.parts(g, :ConArc)
        Graphs.add_edge!(plot_graph, g[arc, :src_con], g[arc, :tgt_con])
    end
    for arc in Catlab.parts(g, :DisArc)
        u, v = g[arc, :src_dis], g[arc, :tgt_dis]
        Graphs.add_edge!(plot_graph, u, v)
        Graphs.add_edge!(plot_graph, v, u)
    end

    node_labels = [g[t, :task_name] for t in Catlab.parts(g, :Task)]
    node_colors = fill("lightblue", Catlab.nparts(g, :Task))

    p = GraphPlot.gplot(plot_graph, nodelabel=node_labels, nodefillc=node_colors, layout=GraphPlot.spring_layout)
    Compose.draw(Compose.PNG(filename, 800px, 800px), p)
    println("Initial JSSP graph saved to $filename")
end

"""
    visualize_partitioned_graph(g::JSSPGraph, filename::String)

Visualize the partitioned JSSP graph with tasks colored by their layer.
"""
function visualize_partitioned_graph(g::JSSPGraph, filename::String)
    plot_graph = Graphs.SimpleDiGraph(Catlab.nparts(g, :Task))
    for arc in Catlab.parts(g, :ConArc)
        Graphs.add_edge!(plot_graph, g[arc, :src_con], g[arc, :tgt_con])
    end
    for arc in Catlab.parts(g, :DisArc)
        u, v = g[arc, :src_dis], g[arc, :tgt_dis]
        Graphs.add_edge!(plot_graph, u, v)
        Graphs.add_edge!(plot_graph, v, u)
    end

    node_labels = [g[t, :task_name] for t in Catlab.parts(g, :Task)]
    layers = g[:layer]
    unique_layers = unique(layers)
    palette = ["lightblue", "lightgreen", "lightpink", "lightyellow", "orange", "violet", "lightgray", "gold", "cyan", "salmon"]
    color_map = Dict(l => palette[mod1(i, length(palette))] for (i, l) in enumerate(unique_layers))
    node_colors = [color_map[l] for l in layers]

    p = GraphPlot.gplot(plot_graph, nodelabel=node_labels, nodefillc=node_colors, layout=GraphPlot.spring_layout)
    Compose.draw(Compose.PNG(filename, 800px, 800px), p)
    println("Partitioned JSSP graph saved to $filename")
end

"""
    visualize_component_graph(C::Graphs.DiGraph, filename::String)

Visualize the component graph showing dependencies between layers.
"""
function visualize_component_graph(C::Graphs.DiGraph, filename::String)
    node_labels = ["Layer $i" for i in 1:Graphs.nv(C)]
    p = GraphPlot.gplot(C, nodelabel=node_labels, layout=GraphPlot.spring_layout)
    Compose.draw(Compose.PNG(filename, 400px, 400px), p)
    println("Component graph saved to $filename")
end

"""
    visualize_schedule(Gσ::Graphs.DiGraph, g::JSSPGraph, filename::String)

Visualize the final schedule with all edges oriented.
"""
function visualize_schedule(Gσ::Graphs.DiGraph, g::JSSPGraph, filename::String)
    node_labels = [g[t, :task_name] for t in 1:Graphs.nv(Gσ)]
    p = GraphPlot.gplot(Gσ, nodelabel=node_labels, layout=GraphPlot.spring_layout)
    Compose.draw(Compose.PNG(filename, 800px, 800px), p)
    println("Final schedule saved to $filename")
end

# === Run Full Experiment ===
"""
    run_full_experiment(num_jobs, tasks_per_job, num_machines, seed, max_ops_per_layer, visualize)

Run a complete JSSP experiment with the specified parameters and optional visualizations.
"""
function run_full_experiment(num_jobs, tasks_per_job, num_machines, seed, max_ops_per_layer, visualize::Bool=false)
    g = create_jssp_instance(num_jobs, tasks_per_job, num_machines, seed)
    if visualize
        visualize_jssp_graph(g, "initial_jssp_graph.png")
    end
    ℓ, k = layered_partitioning(g, max_ops_per_layer)
    if visualize
        visualize_partitioned_graph(g, "partitioned_jssp_graph.png")
        C = compute_component_graph(g)
        visualize_component_graph(C, "component_graph.png")
    end
    local_sols = [solve_layer(i, g) for i in 1:k]
    Gσ = compose_solutions(g, local_sols)
    if visualize
        visualize_schedule(Gσ, g, "final_schedule.png")
    end
    println("Experiment with $num_jobs jobs, $tasks_per_job tasks/job, $num_machines machines")
    println("Number of layers: $k")
    println("Global orientation is acyclic: $(is_acyclic(Gσ))")
end

# Executing the experiment
run_full_experiment(num_jobs, tasks_per_job, num_machines, seed, max_ops_per_layer, visualize)

# === Run Tests ===
println("\nRunning Component Tests:")
test_create_jssp(num_jobs, tasks_per_job, num_machines, seed)
test_layered_partitioning(max_ops_per_layer, num_jobs, tasks_per_job, num_machines, seed)
test_solve_layer(max_ops_per_layer, num_jobs, tasks_per_job, num_machines, seed)
test_compose_solutions(max_ops_per_layer, num_jobs, tasks_per_job, num_machines, seed)