
module OpenJSSPCore
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.ACSetInterface
# (Optional, but helps in some Catlab versions:)
using Catlab.Theories

# =========================
#  A. Schema & Basic Ops
# =========================

@present SchJSSP(FreeSchema) begin
  (Op, Mach, Prec, Conf)::Ob
  src_prec::Hom(Prec, Op)
  tgt_prec::Hom(Prec, Op)
  src_conf::Hom(Conf, Op)
  tgt_conf::Hom(Conf, Op)
  mach::Hom(Op, Mach)

  # --- Declare attribute TYPES inside the schema ---
  Str::AttrType
  BoolT::AttrType

  # --- Attributes use those AttrTypes, not Julia types ---
  olabel::Attr(Op, Str)
  is_src::Attr(Op, BoolT)
  is_tgt::Attr(Op, BoolT)
  mlabel::Attr(Mach, Str)
end

# Map AttrTypes -> Julia types when creating the ACSet type.
# Catlab keyword name has changed across versions; try attrtype first, attr_types if needed.

try
  @acset_type JSSPGraph(SchJSSP,
    index = [:src_prec, :tgt_prec, :src_conf, :tgt_conf, :mach],
    attrtype = Dict(:Str => String, :BoolT => Bool)
  )
catch err
  # Fallback for alternate keyword in some Catlab versions
  @acset_type JSSPGraph(SchJSSP,
    index = [:src_prec, :tgt_prec, :src_conf, :tgt_conf, :mach],
    attr_types = Dict(:Str => String, :BoolT => Bool)
  )
end

# Helpers: counts
nv(G::JSSPGraph) = nparts(G, :Op)
nprec(G::JSSPGraph) = nparts(G, :Prec)
nconf(G::JSSPGraph) = nparts(G, :Conf)
nmach(G::JSSPGraph) = nparts(G, :Mach)

# Add parts
function add_machine!(G::JSSPGraph; label::AbstractString="M")
  m = add_part!(G, :Mach)
  set_attr!(G, :mlabel, m, String(label))
  return m
end

function add_op!(G::JSSPGraph; label::AbstractString="op", mach::Int, is_src::Bool=false, is_tgt::Bool=false)
  v = add_part!(G, :Op)
  set_attr!(G, :olabel, v, String(label))
  set_attr!(G, :is_src, v, is_src)
  set_attr!(G, :is_tgt, v, is_tgt)
  set_subpart!(G, :mach, v, mach)
  return v
end

function add_prec!(G::JSSPGraph, u::Int, v::Int)
  e = add_part!(G, :Prec)
  set_subpart!(G, :src_prec, e, u)
  set_subpart!(G, :tgt_prec, e, v)
  return e
end

function add_conf!(G::JSSPGraph, u::Int, v::Int)
  u2, v2 = (u <= v) ? (u, v) : (v, u)  # canonical orientation for undirected edge
  e = add_part!(G, :Conf)
  set_subpart!(G, :src_conf, e, u2)
  set_subpart!(G, :tgt_conf, e, v2)
  return e
end

# Accessors
function get_prec_edges(G::JSSPGraph)
  [(subpart(G,:src_prec,e), subpart(G,:tgt_prec,e)) for e in 1:nprec(G)]
end
function get_conf_edges(G::JSSPGraph)
  [(subpart(G,:src_conf,e), subpart(G,:tgt_conf,e)) for e in 1:nconf(G)]
end
function get_mach(G::JSSPGraph, v::Int)
  subpart(G, :mach, v)
end

# =========================
#  B. Interfaces & Splits
# =========================

"""
    make_interface(G, C::Vector{Int})

Return the list of interface vertices `C` (as a Set), and two convenience maps:
- `is_interface(v)`: Bool
- `label(v)`: String
We do not create a separate ACSet for L(C) because in this simplified pipeline, we only need
the vertex set; the interface is edgeless by construction.
"""
function make_interface(G::JSSPGraph, C::Vector{Int})
  Cset = Set(C)
  is_interface = v->(v in Cset)
  label = v->get_attr(G, :olabel, v)
  return Cset, is_interface, label
end

"""
    split_by_interface(G, C::Vector{Int}, A_side::Vector{Int})

Given an interface C and a side assignment A for non-interface vertices,
build induced subgraphs G1 and G2 that both include C, with:
- internal precedence edges kept,
- internal conflict edges kept,
- cross-side edges dropped (must be admissible to recover by pushout later).

Returns (G1, map1, G2, map2) where map* are Dicts from original op id -> new id (or 0 if absent).
"""
function split_by_interface(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int})
  Cset = Set(C)
  Aset = Set(A_side)
  allV = collect(1:nv(G))
  # Sanity: A should be subset of V\C
  @assert isempty(setdiff(Aset, setdiff(Set(allV), Cset))) "A_side must be subset of V \\ C"
  Bset = Set(setdiff(allV, union(Cset, Aset))) # remaining vertices go to B-side

  function induce(subset::Set{Int})
    H = JSSPGraph()
    # Clone machines
    mach_map = Dict{Int,Int}()
    for m in 1:nmach(G)
      m2 = add_machine!(H; label=get_attr(G,:mlabel,m))
      mach_map[m] = m2
    end
    # Add vertices in `subset`
    vmap = Dict{Int,Int}()
    for v in subset
      v2 = add_op!(H; label=get_attr(G,:olabel,v),
                      mach=mach_map[get_mach(G,v)],
                      is_src=get_attr(G,:is_src,v),
                      is_tgt=get_attr(G,:is_tgt,v))
      vmap[v] = v2
    end
    # Add internal precedence edges
    for (u,v) in get_prec_edges(G)
      if (u in subset) && (v in subset)
        add_prec!(H, vmap[u], vmap[v])
      end
    end
    # Add internal conflict edges
    for (u,v) in get_conf_edges(G)
      if (u in subset) && (v in subset)
        add_conf!(H, vmap[u], vmap[v])
      end
    end
    return H, vmap
  end

  G1, map1 = induce(union(Cset, Aset))
  G2, map2 = induce(union(Cset, Bset))
  return G1, map1, G2, map2
end

# =========================
#  C. Admissibility Checks
# =========================

"""
(A1) Precedence-respecting interface:
No precedence edge may cross from A\\C to B\\C.
"""
function check_A1_precedence(G::JSSPGraph, C::Set{Int}, A::Set{Int}, B::Set{Int})
  for (u,v) in get_prec_edges(G)
    if (u in A && !(u in C)) && (v in B && !(v in C))
      return false
    end
  end
  return true
end

"""
(A2) Machine-consistency:
No conflict (disjunctive) edge may cross from A\\C to B\\C.
"""
function check_A2_machine(G::JSSPGraph, C::Set{Int}, A::Set{Int}, B::Set{Int})
  for (u,v) in get_conf_edges(G)
    if (u in A && !(u in C)) && (v in B && !(v in C))
      return false
    end
  end
  return true
end

"""
Admissible cut predicate for a proposed (C, A-side) partition.
"""
function admissible_cut(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int})
  Cset = Set(C)
  Aset = Set(A_side)
  allV = Set(1:nv(G))
  @assert isempty(setdiff(Aset, setdiff(allV, Cset))) "A_side ⊆ V \\ C required"
  Bset = setdiff(allV, union(Cset, Aset))
  check_A1_precedence(G, Cset, Aset, Bset) && check_A2_machine(G, Cset, Aset, Bset)
end

# =========================
#  D. Diagnostics & Utility
# =========================

"""
Effort model: simple superlinear in |V| with linear edge terms.
"""
function effort(G::JSSPGraph; α=1.6, κV=1.0, κP=0.2, κC=0.4)
  κV * (nv(G)^α) + κP * nprec(G) + κC * nconf(G)
end

"""
Coupling on interface: count of cross edges implied by the partition
(i.e., edges that forced inclusion into C). We compute it directly from (C,A,B).
"""
function coupling_on_interface(G::JSSPGraph, C::Vector{Int}, A_side::Vector{Int}; ρ=1.0)
  Cset, Aset = Set(C), Set(A_side)
  allV = Set(1:nv(G))
  Bset = setdiff(allV, union(Cset, Aset))
  cross_prec = 0
  for (u,v) in get_prec_edges(G)
    if (u in Aset && v in Bset) || (u in Bset && v in Aset)
      cross_prec += 1
    end
  end
  cross_conf = 0
  for (u,v) in get_conf_edges(G)
    if (u in Aset && v in Bset) || (u in Bset && v in Aset)
      cross_conf += 1
    end
  end
  return cross_prec + ρ*cross_conf
end

"""
Overhead model: affine in |C|.
"""
overhead(C::Vector{Int}; τ0=0.0, τ1=0.5) = τ0 + τ1*length(C)

"""
Bound slack estimate: placeholder 0.0 (hook for later calibration).
"""
bound_slack_estimate(::JSSPGraph, ::Vector{Int}, ::Vector{Int}) = 0.0

"""
Utility of cut:
U = [effort saved] - λ * K∂ - μ * B∂
"""
function cut_utility(G::JSSPGraph, G1::JSSPGraph, C::Vector{Int}, G2::JSSPGraph;
                     λ::Real=1.0, μ::Real=0.0)
  E = effort(G) - (effort(G1) + effort(G2) + overhead(C))
  K = coupling_on_interface(G, C, Int[])
  B = bound_slack_estimate(G, C, Int[])
  return E - λ*K - μ*B
end

# =========================
#  E. Simple Interface Proposals
# =========================

"""
propose_interface_by_bfs(G; from_src=true)

Compute a crude A/B split by precedence BFS levels, then set C to be endpoints
of crossing edges and iterate until admissible. Returns (C, A_side).
"""
function propose_interface_by_bfs(G::JSSPGraph; from_src::Bool=true)
  # Pick source or target if marked; otherwise fall back to 1
  roots = [v for v in 1:nv(G) if get_attr(G,:is_src,v)]
  roots = (from_src && !isempty(roots)) ? roots : ([v for v in 1:nv(G) if get_attr(G,:is_tgt,v)])
  root = isempty(roots) ? 1 : first(roots)

  # Build precedence adjacency (directed)
  succs = Dict{Int,Vector{Int}}()
  for v in 1:nv(G); succs[v] = Int[]; end
  for (u,v) in get_prec_edges(G); push!(succs[u], v); end

  # BFS levels
  level = Dict{Int,Int}(); for v in 1:nv(G); level[v] = typemax(Int); end
  queue = [root]; level[root]=0
  while !isempty(queue)
    u = popfirst!(queue)
    for v in succs[u]
      if level[v] == typemax(Int)
        level[v] = level[u] + 1
        push!(queue, v)
      end
    end
  end
  # Partition by parity of level (fallback level=0 if unreachable)
  A = Int[]; B = Int[]
  for v in 1:nv(G)
    ℓ = (level[v] == typemax(Int)) ? 0 : level[v]
    if isodd(ℓ); push!(A, v) else push!(B, v) end
  end
  C = Int[]
  # Grow interface until admissible by moving endpoints of crossing edges into C
  function grow_until_admissible!(C::Vector{Int}, A::Vector{Int})
    Cset = Set(C)
    Aset = Set(A)
    allV = Set(1:nv(G))
    changed = true
    while changed
      changed = false
      Bset = setdiff(allV, union(Cset, Aset))
      # Precedence crossings
      for (u,v) in get_prec_edges(G)
        if (u in Aset && v in Bset) || (u in Bset && v in Aset)
          if !(u in Cset); push!(C, u); push!(C, v); Cset = Set(C); changed = true; end
        end
      end
      # Conflict crossings
      for (u,v) in get_conf_edges(G)
        if (u in Aset && v in Bset) || (u in Bset && v in Aset)
          if !(u in Cset); push!(C, u); push!(C, v); Cset = Set(C); changed = true; end
        end
      end
    end
  end
  grow_until_admissible!(C, A)
  # Ensure uniqueness
  C = collect(Set(C))
  return C, A
end

"""
best_admissible_cut(G; λ=1.0, μ=0.0)

Currently: single-shot proposal via BFS-based interface growth.
Returns (G1, C, G2, U) or nothing if inadmissible or non-positive utility.
"""
function best_admissible_cut(G::JSSPGraph; λ::Real=1.0, μ::Real=0.0)
  C, A = propose_interface_by_bfs(G)
  if !admissible_cut(G, C, A); return nothing; end
  G1, _, G2, _ = split_by_interface(G, C, A)
  U = (effort(G) - (effort(G1)+effort(G2)+overhead(C))) - λ*coupling_on_interface(G, C, A) - μ*0.0
  return (G1, C, G2, U)
end

# =========================
#  F. DOT Visualization
# =========================

"""
write_dot(G, path; partition=nothing, interface=Int[])

Render a digraph in DOT:
- precedence edges: solid, arrowed
- conflict edges: dashed, dir=none
- interface vertices: doublecircle
- src/tgt colored distinctly
- optional partition coloring for A/B (1/2)
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
    # precedence edges
    for (u,v) in get_prec_edges(G)
      println(io, "  v$u -> v$v [color=black, style=solid];")
    end
    # conflict edges (undirected look)
    for (u,v) in get_conf_edges(G)
      println(io, "  v$u -> v$v [color=gray40, style=dashed, dir=none];")
    end
    println(io, "}")
  finally
    close(io)
  end
  return path
end

end # module
