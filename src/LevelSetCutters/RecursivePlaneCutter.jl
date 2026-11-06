## RecursivePlaneCutter: the level-set cutter built on the recursive cell cutter
#
# Background cells cut by some level set are split into simplices using `_initial_sub_triangulation`. 
# Every simplex is then cut by all the level sets at once by `CellCutter`.
#
# Zeros: a level-set value with |φ| ≤ atol is zero. Level set k is shifted symbolically by
# -εᵏ (ε → 0⁺), so zeros are IN. The shift is global, so neighbouring cells agree, 
# and the powers are distinct, so no D+1 planes meet at a point.

"""
    struct RecursivePlaneCutter <: Cutter
    RecursivePlaneCutter(; atol=1e-12)

Cutter for `DiscreteGeometry` and `AnalyticalGeometry`. 
Every cut cell is split into simplices, then every simplex is cut 
by all the level sets at once by the recursive `CellCutter`.
The parameter `atol` gives an absolute zero tolerance for level-set values.

## Usage

    cut(RecursivePlaneCutter(),background::DiscreteModel,geom)
    cut_facets(RecursivePlaneCutter(),background::DiscreteModel,geom)
    compute_bgcell_to_inoutcut(RecursivePlaneCutter(),background::DiscreteModel,geom)
    compute_bgfacet_to_inoutcut(RecursivePlaneCutter(),background::DiscreteModel,geom)

It also cuts a simplicial (immersed) model by the cells of a background model:

    cut(RecursivePlaneCutter(),bgmodel::DiscreteModel,imodel::DiscreteModel) -> ImmersedDiscreteModel
"""
struct RecursivePlaneCutter <: Cutter
  atol::Float64
end

RecursivePlaneCutter(; atol::Real=1e-12) = RecursivePlaneCutter(Float64(atol))

function cut(cutter::RecursivePlaneCutter,background::DiscreteModel,geom)
  data = _cut(cutter,get_grid(background),geom,Val(true))
  EmbeddedDiscretization(background, data..., geom)
end

function cut_facets(cutter::RecursivePlaneCutter,background::DiscreteModel,geom)
  D = num_cell_dims(background)
  grid = Grid(ReferenceFE{D-1},background)
  data = _cut(cutter,grid,geom,Val(false))
  EmbeddedFacetDiscretization(background, data..., geom)
end

function compute_bgcell_to_inoutcut(cutter::RecursivePlaneCutter,model::DiscreteModel,geom)
  _compute_bgcell_to_inoutcut(cutter,get_grid(model),geom)
end

function compute_bgfacet_to_inoutcut(cutter::RecursivePlaneCutter,model::DiscreteModel,geom)
  D = num_cell_dims(model)
  _compute_bgcell_to_inoutcut(cutter,Grid(ReferenceFE{D-1},model),geom)
end

## Level-set values: one vector per unique leaf of the CSG tree, snapped (|φ| ≤ atol → 0)

function _levelset_values(cutter::RecursivePlaneCutter,grid::Grid,geom::AnalyticalGeometry)
  _, oid_to_ls = _find_unique_leaves(get_tree(geom))
  ls_to_point_to_value, _ = _levelset_values(cutter,grid,discretize(geom,grid))
  ls_to_point_to_value, oid_to_ls
end

function _levelset_values(cutter::RecursivePlaneCutter,grid::Grid,geom::DiscreteGeometry)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(get_tree(geom))
  snap(v) = abs(v) <= cutter.atol ? zero(v) : v
  map(values -> map(snap,values),ls_to_point_to_value), oid_to_ls
end

function _compute_bgcell_to_inoutcut(cutter::RecursivePlaneCutter,grid::Grid,geom)
  ugrid = UnstructuredGrid(grid)
  ls_to_point_to_value, oid_to_ls = _levelset_values(cutter,ugrid,geom)
  cell_to_points = get_cell_node_ids(ugrid)
  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(cell_to_points,point_to_value) for point_to_value in ls_to_point_to_value]
  function conversion(data)
    f,name,meta = data
    ls_to_cell_to_inoutcut[oid_to_ls[objectid(f)]], name, meta
  end
  compute_inoutcut(replace_data(identity,conversion,get_tree(geom)))
end

## The cut

# `with_boundary`: also the subfacets (the pieces of the level sets), for `cut`
function _cut(cutter::RecursivePlaneCutter,grid::Grid,geom,with_boundary::Val)
  ugrid = UnstructuredGrid(grid)
  ls_to_point_to_value, oid_to_ls = _levelset_values(cutter,ugrid,geom)
  subsimplices, ls_to_subpoint_to_value, ls_to_bgcell_to_inoutcut =
    _initial_sub_triangulation(ugrid,ls_to_point_to_value)
  N = length(ls_to_point_to_value)
  out = _cut_simplices(cutter,subsimplices,ls_to_subpoint_to_value,Val(N),with_boundary)
  (ls_to_bgcell_to_inoutcut, out..., oid_to_ls)
end

# Cuts every simplex of `st` (its points carry the level-set values) by the N level sets.
# Returns the subcells and their IN/OUT per level set, and, with boundary, the subfacets
# (normals towards φ > 0) and their IN/OUT/INTERFACE per level set.
function _cut_simplices(
  cutter::RecursivePlaneCutter,
  st::SubCellData{Dc,Dp,T},
  ls_to_point_to_value,
  ::Val{N},
  ::Val{with_boundary}) where {Dc,Dp,T,N,with_boundary}

  cellcutter = CellCutter(Simplex(Val(Dc)),N;atol=cutter.atol)
  F = num_facets(cellcutter.refcell)
  dims = with_boundary ? (Dc-1,Dc) : (Dc,)
  extractor = CellFaceExtractor(dims=dims,planes=F+1:F+N) # (Dc-1)-faces on level sets only
  # Level set k is shifted by -εᵏ
  shifts = [SVector{N,Float64}(ntuple(j -> j == k ? -1.0 : 0.0,Val(N))) for k in 1:N]

  φs = ntuple(k -> zeros(T,Dc+1),Val(N))
  ccache = return_cache(cellcutter,shifts,φs...)
  fcache = return_cache(extractor,ccache)

  point_to_coords = Point{Dp,T}[]
  point_to_rcoords = Vector{eltype(st.point_to_rcoords)}()
  cells = (; points=Int32[], bgcell=Int32[], ls_to_inout=[Int8[] for k in 1:N])
  facets = (; points=Int32[], bgcell=Int32[], normal=VectorValue{Dp,T}[], ls_to_inout=[Int8[] for k in 1:N])

  for s in 1:length(st.cell_to_points)
    spoints = st.cell_to_points[s]
    for k in 1:N, (i,p) in enumerate(spoints)
      φs[k][i] = ls_to_point_to_value[k][p]
    end
    c = evaluate!(ccache,cellcutter,shifts,φs...)
    cf = evaluate!(fcache,extractor,c)

    # Points: the vertices of the cut simplex, mapped affinely from the reference simplex
    X = st.point_to_coords[spoints]
    R = st.point_to_rcoords[spoints]
    offset = Int32(length(point_to_coords))
    for v in 1:num_vertices(c)
      ξ = get_vertex_coordinates(c,v)
      push!(point_to_coords,X[1] + sum(j -> ξ[j]*(X[j+1]-X[1]),1:Dc))
      push!(point_to_rcoords,R[1] + sum(j -> ξ[j]*(R[j+1]-R[1]),1:Dc))
    end
    bgcell = st.cell_to_bgcell[s]

    # Subcells: the leaves, positively oriented in the reference simplex (as the simplex)
    for f in 1:num_faces(cf,Dc)
      phase = get_face_phase(cf,Dc,f)
      for sim in get_face_simplices(cf,Dc,f)
        verts = get_simplex_vertices(cf,Dc,sim)
        _push_oriented!(cells.points,c,verts,offset)
        push!(cells.bgcell,bgcell)
        for k in 1:N
          push!(cells.ls_to_inout[k],phase[k])
        end
      end
    end

    # Subfacets: the (Dc-1)-faces on a level set, normal (physical) towards φ > 0
    if with_boundary
      J = TensorValue{Dp,Dc}(ntuple(i -> (X[(i-1)÷Dp+2]-X[1])[(i-1)%Dp+1],Val(Dp*Dc)))
      invJt = transpose(inv(J))
      for f in 1:num_faces(cf,Dc-1)
        phase = get_face_phase(cf,Dc-1,f)
        n = push_normal(invJt,get_face_normal(cf,f))
        for sim in get_face_simplices(cf,Dc-1,f)
          for v in get_simplex_vertices(cf,Dc-1,sim)
            push!(facets.points,offset+v)
          end
          push!(facets.bgcell,bgcell)
          push!(facets.normal,n)
          for k in 1:N
            push!(facets.ls_to_inout[k],phase[k])
          end
        end
      end
    end
  end

  _table(points,n) = Table(points,collect(Int32,1:n:length(points)+1))
  subcells = SubCellData(
    _table(cells.points,Dc+1),cells.bgcell,point_to_coords,point_to_rcoords
  )
  if !with_boundary 
    return (subcells, cells.ls_to_inout)
  end
  subfacets = SubFacetData(
    _table(facets.points,Dc),facets.normal,facets.bgcell,point_to_coords,point_to_rcoords
  )
  return (subcells, cells.ls_to_inout, subfacets, facets.ls_to_inout)
end

# The vertices of a Dc-simplex of the cut cell, swapping the last two if it is negatively
# oriented in reference coordinates
function _push_oriented!(points,c::CutCell{D},verts,offset) where D
  x1 = get_vertex_coordinates(c,verts[1])
  J = TensorValue{D,D}(ntuple(i -> (get_vertex_coordinates(c,verts[(i-1)÷D+2])-x1)[(i-1)%D+1],Val(D*D)))
  flip = det(J) < 0
  for (i,v) in enumerate(verts)
    j = flip && i >= D ? (i == D ? D+1 : D) : i
    push!(points,offset+verts[j])
  end
end

## Cutting an immersed model by a background model

"""
    cut(cutter::RecursivePlaneCutter, bgmodel::DiscreteModel, imodel::DiscreteModel) -> ImmersedDiscreteModel

Cut the simplicial model `imodel` by the cells of `bgmodel`. 
Results in an `ImmersedDiscreteModel` with the cutfaces of every dimension 0:D.
The immersed model must lie completely inside the background model.
The background mode MUST have affine cell maps, to guarantee convexity.
"""
function cut(cutter::RecursivePlaneCutter,bgmodel::DiscreteModel{D,D},imodel::DiscreteModel{D,D}) where D
  igrid, bggrid = get_grid(imodel), get_grid(bgmodel)
  p, bg_p = only(get_polytopes(imodel)), only(get_polytopes(bgmodel))
  @notimplementedif !is_simplex(p) "the immersed model must be simplicial"
  @notimplementedif !(eltype(get_cell_map(bggrid)) <: AffineField) """
    The background cell maps must be `AffineField`s (e.g. Cartesian or simplicial grids): every \
    background cell is cut by its facet planes and the cutfaces are mapped back by its inverse \
    cell map. Got cell maps of type $(eltype(get_cell_map(bggrid)))."""
  bg_X, i_X = get_cell_coordinates(bggrid), get_cell_coordinates(igrid)
  bg_maps = get_cell_map(bggrid)
  icell_to_bgcells = _icell_to_bgcells(i_X,bg_X)

  topo = get_grid_topology(imodel)
  c_to_f = [get_faces(topo,D,d) for d in 0:D]
  f_to_v = [get_faces(topo,d,0) for d in 0:D]
  owner = [map(first,get_faces(topo,d,D)) for d in 0:D]
  # local d-face of T on the facets of T a cutface lies on (its planes 1:D+1)
  lface = [Dict(sort(collect(Int32,fs)) => lf for (lf,fs) in enumerate(get_faces(p,d,D-1))) for d in 0:D-1]
  push!(lface,Dict(Int32[] => 1))

  N = num_facets(bg_p)
  F = num_facets(p)
  cellcutter = CellCutter(p,N;atol=cutter.atol,phases=ntuple(_ -> IN,N))
  φs = ntuple(_ -> zeros(D+1),N)
  shifts = zeros(SVector{D,Float64},N)
  c = return_cache(cellcutter,shifts,φs...)
  extractor = CellFaceExtractor(;planes=1:F) # faces on facets of T only
  cf = return_cache(extractor,c)

  P = Point{D,Float64}
  Group = Tuple{Int32,Vector{Vector{Int}},Vector{P}}
  groups = [[Group[] for _ in 1:num_faces(topo,d)] for d in 0:D]
  for icell in 1:num_cells(imodel)
    X = i_X[icell]
    x(ξ) = X[1] + sum(j -> ξ[j]*(X[j+1]-X[1]),1:D)
    for bg in icell_to_bgcells[icell]
      for (k,h) in enumerate(_cell_halfspaces(bg_p,bg_maps[bg],bg_X[bg]))
        shifts[k] = -SVector(Tuple(h.normal))
        for v in 1:D+1
          φs[k][v] = h(X[v])
        end
      end
      evaluate!(cf,extractor,evaluate!(c,cellcutter,shifts,φs...))
      for d in 0:D, f in 1:num_faces(cf,d)
        face = c_to_f[d+1][icell][lface[d+1][get_face_planes(cf,d,f)]]
        owner[d+1][face] == icell || continue
        verts = collect(get_face_vertices(cf,d,f))
        simplices = [[findfirst(==(v),verts) for v in get_simplex_vertices(cf,d,s)] for s in get_face_simplices(cf,d,f)]
        push!(groups[d+1][face],(Int32(bg),simplices,[x(get_vertex_coordinates(cf,v)) for v in verts]))
      end
    end
  end

  bg_invmaps = lazy_map(inverse_map,get_cell_map(bggrid))
  i_invmaps = lazy_map(inverse_map,get_cell_map(igrid))
  c_to_v = get_faces(topo,D,0)
  iface_to_cutfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,D+1)
  cutfaces = ntuple(Val(D+1)) do k
    data, f_to_cf = _immersed_face_data(Val(k-1),Val(D),groups[k],f_to_v[k],c_to_v,owner[k],i_invmaps,bg_invmaps)
    iface_to_cutfaces[k] = f_to_cf
    data
  end
  metadata = (; atol=cutter.atol, candidates=:bounding_boxes, degeneracies=:translated_background)
  return ImmersedDiscreteModel(imodel,bgmodel,cutfaces,iface_to_cutfaces,metadata)
end

"""
    cut(bgmodel::DiscreteModel, imodel::DiscreteModel) -> ImmersedDiscreteModel

`cut(RecursivePlaneCutter(), bgmodel, imodel)`.
"""
cut(bgmodel::DiscreteModel,imodel::DiscreteModel) = cut(RecursivePlaneCutter(),bgmodel,imodel)

# For each immersed cell, return the background cell candidates whose 
# bounding box intersects the immersed cell's bounding box.
# Similar to what STLCutters does for unstructured grids.
function _icell_to_bgcells(i_X::AbstractArray,bg_X::AbstractArray)
  δ = 1e-3
  function bbox(X)
    pmin, pmax = reduce((a,b) -> min.(a,b),X), reduce((a,b) -> max.(a,b),X)
    Δ = (pmax-pmin)*δ
    pmin-Δ, pmax+Δ
  end
  bg_boxes = map(bbox,bg_X)
  i_boxes = map(bbox,i_X)
  pmin = reduce((a,b) -> min.(a,b),first.(bg_boxes))
  pmax = reduce((a,b) -> max.(a,b),last.(bg_boxes))
  h = maximum(((a,b),) -> maximum(Tuple(b-a)),bg_boxes)
  D = length(pmin)
  n = ntuple(k -> max(1,ceil(Int,(pmax[k]-pmin[k])/h)),D)
  bin(x) = CartesianIndex(ntuple(k -> clamp(floor(Int,(x[k]-pmin[k])/h)+1,1,n[k]),D))
  bins(a,b) = bin(a):bin(b)
  lids = LinearIndices(n)
  bin_to_bgcells = [Int32[] for _ in 1:prod(n)]
  for (bg,(a,b)) in enumerate(bg_boxes), I in bins(a,b)
    push!(bin_to_bgcells[lids[I]],bg)
  end
  overlap((a,b),(c,d)) = all(k -> a[k] ≤ d[k] && c[k] ≤ b[k],1:D)
  data, ptrs, cands = Int32[], Int32[1], Int32[]
  for ibox in i_boxes
    empty!(cands)
    for I in bins(ibox...), bg in bin_to_bgcells[lids[I]]
      overlap(ibox,bg_boxes[bg]) && push!(cands,bg)
    end
    append!(data,unique!(sort!(cands)))
    push!(ptrs,length(data)+1)
  end
  Table(data,ptrs)
end

# Half-spaces whose intersection is the background cell with polytope `p`, affine cell map
# `cmap` and vertices `X`, one per facet, in Gridap's local facet order: `Plane`s with the
# outward unit normal (Gridap's reference facet normals pushed forward, as for boundary
# triangulations), so that the value is the signed distance.
function _cell_halfspaces(p::Polytope{D},cmap::AffineField,X::AbstractVector{<:Point{D}}) where D
  invJt = pinvJt(cmap.gradient)
  map(get_facet_normal(p),get_faces(p,D-1,0)) do n̂, fv
    n = push_normal(invJt,n̂)
    Plane(n,-(n⋅X[first(fv)]))
  end
end

# ImmersedFaceData from the pieces of every iface. `fgroups[f]`: the pieces of d-iface f,
# (bgcell, simplices as local point ids, points).
function _immersed_face_data(::Val{d},::Val{D},fgroups,f_to_v,c_to_v,owner,i_invmaps,bg_invmaps) where {d,D}
  P = Point{D,Float64}
  npoints = sum(g -> sum(((_,_,X),) -> length(X),g;init=0),fgroups;init=0)
  ncells = sum(g -> sum(((_,simplices,_),) -> length(simplices),g;init=0),fgroups;init=0)
  coords, rcoords = Vector{P}(undef,npoints), Vector{P}(undef,npoints)
  ircoords = Vector{Point{d,Float64}}(undef,npoints)
  data = Vector{Int32}(undef,ncells*(d+1))
  cell_to_bgcell = Vector{Int32}(undef,ncells)
  ptrs = zeros(Int32,length(fgroups)+1)
  o, c = 0, 0
  for (f,groups) in enumerate(fgroups)
    icell = owner[f]
    # local vertex (in the owner icell) of every vertex of the iface, in topology order
    lv = [findfirst(==(v),c_to_v[icell]) for v in f_to_v[f]]
    for (bg,simplices,X) in groups
      for (k,x) in enumerate(X)
        coords[o+k] = x
        rcoords[o+k] = evaluate(bg_invmaps[bg],x)
        λ = _barycentric(evaluate(i_invmaps[icell],x))
        ircoords[o+k] = Point{d,Float64}(ntuple(k -> λ[lv[k+1]],Val(d)))
      end
      for t in simplices
        data[c*(d+1)+1:(c+1)*(d+1)] .= t .+ o
        c += 1
        cell_to_bgcell[c] = bg
      end
      o += length(X)
      ptrs[f+1] += length(simplices)
    end
  end
  length_to_ptrs!(ptrs)
  cell_to_points = Table(data,collect(Int32,1:d+1:ncells*(d+1)+1))
  subfaces = SubFaceData{d}(cell_to_points,cell_to_bgcell,coords,rcoords)
  ImmersedFaceData(subfaces,ircoords), Table(Int32.(1:ncells),ptrs)
end

# Barycentric coordinates of the reference point ξ in Gridap's reference simplex
_barycentric(ξ::Point{D}) where D = (1-sum(ξ),ntuple(k -> ξ[k],Val(D))...)
