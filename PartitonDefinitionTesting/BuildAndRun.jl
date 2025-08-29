
# BuildAndRun.jl: include(raw"C:\Users\willc\.julia\dev\CompositionalJulia\PartitonDefinitionTesting\OpenJSSPCore.jl"); include(raw"C:\Users\willc\.julia\dev\CompositionalJulia\PartitonDefinitionTesting\BuildAndRun.jl")
# User-facing file: build a disjunctive graph, call the math/core, and visualize.
isfile(joinpath(@__DIR__, "OpenJSSPCore.jl")) || error("OpenJSSPCore.jl not next to BuildAndRun.jl")

include(joinpath(@__DIR__, "OpenJSSPCore.jl"))
using .OpenJSSPCore

# ---------- Example builder API ----------

function build_tiny_example()
  G = JSSPGraph()

  # Machines
  mNone = add_machine!(G; label="None")
  mA = add_machine!(G; label="A")
  mB = add_machine!(G; label="B")

  # Ops
  s  = add_op!(G; label="src", mach=mNone, is_src=true)
  a1 = add_op!(G; label="A1",  mach=mA)
  a2 = add_op!(G; label="A2",  mach=mA)
  b1 = add_op!(G; label="B1",  mach=mB)
  b2 = add_op!(G; label="B2",  mach=mB)
  t  = add_op!(G; label="tgt", mach=mNone, is_tgt=true)

  # Precedence (DAG)
  add_prec!(G, s, a1)
  add_prec!(G, a1, a2)
  add_prec!(G, s, b1)
  add_prec!(G, b1, b2)
  add_prec!(G, a2, t)
  add_prec!(G, b2, t)

  # Conflicts (machine disjunctives)
  add_conf!(G, a1, b1)
  add_conf!(G, a2, b2)

  return G
end

# ---------- Run: decompose once and visualize ----------

function main(; outdir::AbstractString="out", λ::Real=1.0, μ::Real=0.0)
  isdir(outdir) || mkpath(outdir)

  # 1) Build or load your JSSP graph
  G = build_tiny_example()

  # 2) Try a single admissible cut (BFS-based proposal)
  cut = best_admissible_cut(G; λ=λ, μ=μ)
  partition = Dict{Int,Int}()
  interface = Int[]
  if cut === nothing
    println("No admissible or beneficial cut found. Visualizing full graph.")
  else
    (G1, C, G2, U) = cut
    println("Cut utility U = ", U)
    interface = C
    # For visualization: color A/B by pulling back the split
    # We can approximate A as the non-interface vertices of G1.
    Aset = Set{Int}()
    # Recover which original vertices ended up in G1 (we don't return maps in this simple version)
    # As a proxy, color interface as 0, and randomly split remaining by BFS parity for the picture.
    # If you want exact coloring, extend best_admissible_cut to return the A-side set.
    # Here we just color interface as 0 and leave others white unless you prefer a heuristic.
    # We'll instead skip partition colors unless you adapt core to return A.
  end

  # 3) Visualize
  dot_all = joinpath(outdir, "graph.dot")
  write_dot(G, dot_all; partition=nothing, interface=interface)
  println("Wrote DOT to: ", dot_all)
  println("Render with: dot -Tsvg $(dot_all) -o $(joinpath(outdir, "graph.svg"))")

  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
