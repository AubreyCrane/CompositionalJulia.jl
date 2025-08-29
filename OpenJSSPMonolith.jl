##############################
# OpenJSSP_Monolith.jl (no-attr_types mapping)
# Single-file, self-contained pipeline for disjunctive graphs (JSSP) with
# - Catlab ACSet schema
# - Friendly DSL to edit machines/jobs/ops
# - Admissible cut proposal + utility
# - DOT visualization
# Run: julia OpenJSSP_Monolith.jl
##############################

module OpenJSSPMono

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.ACSetInterface
using Catlab.Theories

# =========================
#  Schema & ACSet type
# =========================

@present SchJSSP(FreeSchema) begin
  (Op, Mach, Prec, Conf)::Ob
  src_prec::Hom(Prec, Op)
  tgt_prec::Hom(Prec, Op)
  src_conf::Hom(Conf, Op)
  tgt_conf::Hom(Conf, Op)
  mach::Hom(Op, Mach)

  # Attribute types declared here; mapped implicitly by Catlab (no attr_types kw)
  Str::AttrType
  BoolT::AttrType

  # Attributes
  olabel::Attr(Op, Str)    # human-readable op name
  is_src::Attr(Op, BoolT)
  is_tgt::Attr(Op, BoolT)

  mlabel::Attr(Mach, Str)  # human-readable machine name
end

# NOTE: Older Catlab/ACSets versions do not accept the attr_types/attrtype keyword in @acset_type.
# We therefore omit it here; attributes will be stored without explicit Julia type mapping.
@acset_type JSSPGraph(SchJSSP,
  index = [:src_prec, :tgt_prec, :src_conf, :tgt_conf, :mach]
)

# =========================
#  Model wrapper & DSL
# =========================

mutable struct JSSPModel
  G::JSSPGraph
  mid::Dict{Symbol,Int}      # machine name -> id
  oid::Dict{Symbol,Int}      # op name -> id
  rid::Dict{Int,Symbol}      # reverse: id -> op name
  jid::Dict{Symbol,Vector{Symbol}}  # job name -> op symbols
end

jssp() = JSSPModel(JSSPGraph(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Dict{Int,Symbol}(), Dict{Symbol,Vector{Symbol}}())

# ---- Helpers ----
nv(G::JSSPGraph) = nparts(G, :Op)
nprec(G::JSSPGraph) = nparts(G, :Prec)
nconf(G::JSSPGraph) = nparts(G, :Conf)
nmach(G::JSSPGraph) = nparts(G, :Mach)

function _ensure_machine!(M::JSSPModel, name::Symbol)
  if haskey(M.mid, name); return M.mid[name]; end
  m = add_part!(M.G, :Mach)
  set_attr!(M.G, :mlabel, m, String(name))
  M.mid[name] = m
  return m
end

function _ensure_op_name_free!(M::JSSPModel, name::Symbol)
  haskey(M.oid, name) && error("Operation name $name already exists.")
end

# ---- DSL: machines ----
function mach!(M::JSSPModel, name::Symbol)
  _ensure_machine!(M, name)
  return M
end

# ---- DSL: operations ----
function op!(M::JSSPModel, name::Symbol; machine::Symbol, job::Union{Symbol,Nothing}=nothing, src::Bool=false, tgt::Bool=false)
  _ensure_op_name_free!(M, name)
  mid = _ensure_machine!(M, machine)
  v = add_part!(M.G, :Op)
  set_attr!(M.G, :olabel, v, String(name))
  set_attr!(M.G, :is_src, v, src)
  set_attr!(M.G, :is_tgt, v, tgt)
  set_subpart!(M.G, :mach, v, mid)
  M.oid[name] = v
  M.rid[v] = name
  if job !== nothing
    vec = get!(M.jid, job, Vector{Symbol}())
    push!(vec, name)
  end
  return M
end

# ---- DSL: precedence & conflicts ----
function prec!(M::JSSPModel, u::Symbol, v::Symbol)
  haskey(M.oid,u) || error("Unknown op: $u")
  haskey(M.oid,v) || error("Unknown op: $v")
  e = add_part!(M.G, :Prec)
  set_subpart!(M.G, :src_prec, e, M.oid[u])
  set_subpart!(M.G, :tgt_prec, e, M.oid[v])
  return M
end

function conf!(M::JSSPModel, u::Symbol, v::Symbol)
  haskey(M.oid,u) || error("Unknown op: $u")
  haskey(M.oid,v) || error("Unknown op: $v")
  u2, v2 = (M.oid[u] <= M.oid[v]) ? (M.oid[u], M.oid[v]) : (M.oid[v], M.oid[u])
  e = add_part!(M.G, :Conf)
  set_subpart!(M.G, :src_conf, e, u2)
  set_subpart!(M.G, :tgt_conf, e, v2)
  return M
end

# Convenience: define a job as a chain of ops on a machine
function job!(M::JSSPModel, jobname::Symbol, ops::Vector{Symbol}; machine::Symbol)
  _ensure_machine!(M, machine)
  for (i,opname) in pairs(ops)
    op!(M, opname; machine=machine, job=jobname)
    if i > 1
      prec!(M, ops[i-1], ops[i])
    end
  end
  return M
end

# Chain precedence across any provided sequence (possibly spanning machines)
function chain!(M::JSSPModel, seq::Vector{Symbol})
  for i in 2:length(seq)
    prec!(M, seq[i-1], seq[i])
  end
  return M
end

# Mark special ops
src!(M::JSSPModel, opname::Symbol) = (set_attr!(M.G, :is_src, M.oid[opname], true); M)
tgt!(M::JSSPModel, opname::Symbol) = (set_attr!(M.G, :is_tgt, M.oid[opname], true); M)

# ---- Low-level getters ----
function prec_edges(G::JSSPGraph)
  [(subpart(G,:src_prec,e), subpart(G,:tgt_prec,e)) for e in 1:nprec(G)]
end
function conf_edges(G::JSSPGraph)
  [(subpart(G,:src_conf,e), subpart(G,:tgt_conf,e)) for e in 1:nconf(G)]
end
machid(G::JSSPGraph, v::Int) = subpart(G, :mach, v)
olabel(G::JSSPGraph, v::Int) = get_attr(G, :olabel, v)

# =========================
#  Interface, split, admissibility
# =========================

# Build induced subgraph on a subset of vertices
function _induce(G::JSSPGraph, subset::Set{Int})
  H = JSSPGraph()
  # clone machines
  mach_map = Dict{Int,Int}()
  for m in 1:nmach(G)
    m2 = add_part!(H, :Mach)
    set_attr!(H, :mlabel, m2, get_attr(G,:mlabel,m))
    mach_map[m] = m2
  end
  # add vertices
  vmap = Dict{Int,Int}()
  for v in subset
    v2 = add_part!(H, :Op)
    set_attr!(H,:olabel,v2, get_attr(G,:olabel,v))
    set_attr!(H,:is_src,v2, get_attr(G,:is_src,v))
    set_attr!(H,:is_tgt,v2, get_attr(G,:is_tgt,v))
    set_subpart!(H,:mach,v2, mach_map[machid(G,v)])
    vmap[v] = v2
  end
  # internal precedence
  for (u,v) in prec_edges(G)
    if (u in subset) && (v in subset)
      e = add_part!(H,:Prec)
      set_subpart!(H,:src_prec,e, vmap[u])
      set_subpart!(H,:tgt_prec,e, vmap[v])
    end
  end
  # internal conflicts
  for (u,v) in conf_edges(G)
    if (u in subset) && (v in subset)
      e = add_part!(H,:Conf)
      u2, v2 = (vmap[u] <= vmap[v]) ? (vmap[u], vmap[v]) : (vmap[v], vmap[u])
      set_subpart!(H,:src_conf,e, u2)
      set_subpart!(H,:tgt_conf,e, v2)
    end
  end
  return H, vmap
end

# Split by interface C and A-side assignment
function split_by_interface(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int})
  Cset = Set(C); Aset = Set(A_side)
  allV = Set(1:nv(G))
  Bset = setdiff(allV, union(Cset, Aset))
  G1, map1 = _induce(G, union(Cset, Aset))
  G2, map2 = _induce(G, union(Cset, Bset))
  return G1, map1, G2, map2
end

# Admissibility checks (A1: precedence; A2: machine conflicts)
function _check_A1_precedence(G::JSSPGraph, C::Set{Int}, A::Set{Int}, B::Set{Int})
  for (u,v) in prec_edges(G)
    if (u in A && !(u in C)) && (v in B && !(v in C))
      return false
    end
  end
  return true
end

function _check_A2_machine(G::JSSPGraph, C::Set{Int}, A::Set{Int}, B::Set{Int})
  for (u,v) in conf_edges(G)
    if (u in A && !(u in C)) && (v in B && !(v in C))
      return false
    end
  end
  return true
end

function admissible_cut(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int})
  Cset = Set(C); Aset = Set(A_side)
  allV = Set(1:nv(G))
  @assert isempty(setdiff(Aset, setdiff(allV, Cset))) "A_side ⊆ V \\ C required"
  Bset = setdiff(allV, union(Cset, Aset))
  _check_A1_precedence(G, Cset, Aset, Bset) && _check_A2_machine(G, Cset, Aset, Bset)
end

# =========================
#  Diagnostics & utility
# =========================

function effort(G::JSSPGraph; α=1.6, κV=1.0, κP=0.2, κC=0.4)
  κV * (nv(G)^α) + κP * nprec(G) + κC * nconf(G)
end

function coupling_on_interface(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int}; ρ=1.0)
  Cset = Set(C); Aset = Set(A_side)
  allV = Set(1:nv(G))
  Bset = setdiff(allV, union(Cset, Aset))
  cross_prec = 0
  for (u,v) in prec_edges(G)
    if (u in Aset && v in Bset) || (u in Bset && v in Aset)
      cross_prec += 1
    end
  end
  cross_conf = 0
  for (u,v) in conf_edges(G)
    if (u in Aset && v in Bset) || (u in Bset && v in Aset)
      cross_conf += 1
    end
  end
  return cross_prec + ρ*cross_conf
end

overhead(C::Vector{Int}; τ0=0.0, τ1=0.5) = τ0 + τ1*length(C)
bound_slack_estimate(::JSSPGraph, ::Vector{Int}, ::Vector{Int}) = 0.0

function cut_utility(G::JSSPGraph, G1::JSSPGraph, C::Vector{Int}, G2::JSSPGraph; λ::Real=1.0, μ::Real=0.0)
  E = effort(G) - (effort(G1) + effort(G2) + overhead(C))
  K = coupling_on_interface(G, C, Int[])
  B = 0.0
  return E - λ*K - μ*B
end

# =========================
#  Interface proposal (simple, reliable)
# =========================

# BFS from a marked src (if any), else vertex 1
function _bfs_levels(G::JSSPGraph)
  roots = Int[ v for v in 1:nv(G) if get_attr(G,:is_src,v) ]
  root = isempty(roots) ? 1 : first(roots)
  succs = Dict{Int,Vector{Int}}(v=>Int[] for v in 1:nv(G))
  for (u,v) in prec_edges(G); push!(succs[u], v); end
  level = Dict{Int,Int}(v=>typemax(Int) for v in 1:nv(G))
  level[root]=0
  queue = [root]
  while !isempty(queue)
    u = popfirst!(queue)
    for v in succs[u]
      if level[v] == typemax(Int)
        level[v] = level[u]+1
        push!(queue, v)
      end
    end
  end
  return level
end

# Grow interface C until (A1)-(A2) are satisfied
function _grow_interface!(G::JSSPGraph, C::Vector{Int}, A::Vector{Int})
  Cset = Set(C); Aset = Set(A)
  allV = Set(1:nv(G))
  changed = true
  while changed
    changed = false
    Bset = setdiff(allV, union(Cset, Aset))
    for (u,v) in prec_edges(G)
      if (u in Aset && v in Bset) || (u in Bset && v in Aset)
        if !(u in Cset); push!(C,u); Cset = Set(C); changed=true; end
        if !(v in Cset); push!(C,v); Cset = Set(C); changed=true; end
      end
    end
    for (u,v) in conf_edges(G)
      if (u in Aset && v in Bset) || (u in Bset && v in Aset)
        if !(u in Cset); push!(C,u); Cset = Set(C); changed=true; end
        if !(v in Cset); push!(C,v); Cset = Set(C); changed=true; end
      end
    end
  end
  # dedup
  unique!(C)
  return C, A
end

# One-shot best admissible cut (single proposal)
function best_admissible_cut(G::JSSPGraph; λ::Real=1.0, μ::Real=0.0)
  level = _bfs_levels(G)
  A = Int[]; B = Int[]
  for v in 1:nv(G)
    ℓ = get(level, v, 0)
    if isodd(ℓ); push!(A,v) else push!(B,v) end
  end
  C = Int[]
  _grow_interface!(G, C, A)
  if !admissible_cut(G, C, A); return nothing; end
  G1, _, G2, _ = split_by_interface(G, C, A)
  U = cut_utility(G, G1, C, G2; λ=λ, μ=μ)
  return (G1, C, G2, U, A)  # return A for coloring
end

# =========================
#  DOT Visualization
# =========================

"""
write_dot(G, path; partition, interface)
- precedence: solid arrows
- conflicts: dashed, dir=none
- interface: doublecircle
- src/tgt colored
- partition coloring: A=1 (blue), B=2 (gold)
"""
function write_dot(G::JSSPGraph, path::AbstractString;
                   partition::Union{Nothing,Dict{Int,Int}}=nothing,
                   interface::Vector{Int}=Int[])
  io = open(path, "w")
  try
    println(io, "digraph G {")
    println(io, "  rankdir=LR; node [shape=circle, style=filled, fillcolor=white];")
    I = Set(interface)
    for v in 1:nv(G)
      label = get_attr(G, :olabel, v)
      isS = get_attr(G, :is_src, v)
      isT = get_attr(G, :is_tgt, v)
      shape = (v in I) ? "doublecircle" : "circle"
      fill = "white"
      if partition !== nothing && haskey(partition, v)
        p = partition[v]
        fill = (p == 1) ? "lightblue" : (p == 2 ? "lightgoldenrod" : "white")
      end
      if isS; fill = "palegreen"; end
      if isT; fill = "lightcoral"; end
      println(io, "  v$v [label=\"$(label)\", shape=$shape, fillcolor=$fill];")
    end
    for (u,v) in prec_edges(G)
      println(io, "  v$u -> v$v [color=black, style=solid];")
    end
    for (u,v) in conf_edges(G)
      println(io, "  v$u -> v$v [color=gray40, style=dashed, dir=none];")
    end
    println(io, "}")
  finally
    close(io)
  end
  return path
end

# Convenience: partition dictionary from A/C sets
function _partition_from_A_C(n::Int, A::Vector{Int}, C::Vector{Int})
  P = Dict{Int,Int}()
  Aset, Cset = Set(A), Set(C)
  for v in 1:n
    if v in Cset
      P[v] = 0
    elseif v in Aset
      P[v] = 1
    else
      P[v] = 2
    end
  end
  return P
end

# =========================
#  Demo & CLI
# =========================

function demo_model()
  M = jssp()
  mach!(M, :A); mach!(M, :B); mach!(M, :None)  # None for src/tgt if desired

  # Jobs on machines
  job!(M, :J1, [:A1,:A2]; machine=:A)
  job!(M, :J2, [:B1,:B2]; machine=:B)

  # Cross job precedence
  op!(M, :src; machine=:None, src=true)
  op!(M, :tgt; machine=:None, tgt=true)
  chain!(M, [:src,:A1,:A2,:tgt])
  chain!(M, [:src,:B1,:B2,:tgt])

  # Disjunctive conflicts (same machines across jobs)
  conf!(M, :A1, :B1)
  conf!(M, :A2, :B2)

  return M
end

function main(; outdir::AbstractString="out", λ::Real=1.0, μ::Real=0.0)
  isdir(outdir) || mkpath(outdir)
  # Build or edit your model here.
  M = demo_model()
  G = M.G

  # Try a single admissible cut
  cut = best_admissible_cut(G; λ=λ, μ=μ)
  if cut === nothing
    println("No admissible or beneficial cut found. Visualizing full graph only.")
    write_dot(G, joinpath(outdir, "G.dot"))
  else
    (G1, C, G2, U, A) = cut
    println("Cut utility U = ", U)
    P = _partition_from_A_C(nv(G), A, C)
    write_dot(G,  joinpath(outdir, "G.dot");  partition=P, interface=C)
    write_dot(G1, joinpath(outdir, "G1.dot"))
    write_dot(G2, joinpath(outdir, "G2.dot"))
    println("Wrote DOT files G.dot, G1.dot, G2.dot in: ", outdir)
  end
  println("Render with Graphviz, e.g.: dot -Tsvg out/G.dot -o out/G.svg")
end

# Execute if called as a script
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end

end # module
