include("../src/Partitioning.jl")
include("../src/Optimizer.jl")
include("../src/Integration.jl")

using .Partitioning
using .Optimizer
using .Integration
using Graphs

# Run a full experiment
function run_full_experiment(num_jobs, tasks_per_job, num_machines)
    g = Partitioning.create_scaled_jssp(num_jobs=num_jobs, tasks_per_job=tasks_per_job, num_machines=num_machines)
    ℓ, k, _ = layered_partitioning(g)
    local_sols = [solve_layer(i, g) for i in 1:k]
    Gσ = compose_solutions(g, local_sols)
    println("Experiment with $num_jobs jobs, $tasks_per_job tasks/job, $num_machines machines")
    println("Number of layers: $k")
    println("Global orientation is acyclic: $(!is_cyclic(Gσ))")
end

# Run a small experiment
run_full_experiment(2, 3, 3) # jobs, tasks per job, machines

# Test individual components
println("\nTesting Partitioning:")
Partitioning.test_partitioning()

println("\nTesting Optimizer:")
Optimizer.test_optimizer()

println("\nTesting Integration:")
Integration.test_integration()