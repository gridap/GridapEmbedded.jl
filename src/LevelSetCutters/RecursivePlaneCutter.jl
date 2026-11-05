## RecursivePlaneCutter: the level-set cutter built on the recursive cell cutter
#
# The background cells cut by some level set are split into simplices exactly as for
# `LevelSetCutter` (`_initial_sub_triangulation`: Gridap's `simplexify` of the reference cell,
# level-set values at the vertices). Every simplex is then cut by all the level sets at once by
# `CellCutter` (P1 level sets are affine on a simplex, so their zero sets are planes): the
# leaves of its tree are the subcells, the faces of the leaves lying on a level set are the
# subfacets. The result is the same `EmbeddedDiscretization` / `EmbeddedFacetDiscretization`
# as `LevelSetCutter`'s, so the CSG machinery works unchanged.
#
# Zeros: a level-set value with |φ| ≤ atol is zero. Level set k is shifted symbolically by
# -εᵏ (ε → 0⁺), so zeros are IN, as `LevelSetCutter`'s φ ≤ 0. The shift is global, so
# neighbouring cells agree, and the powers are distinct, so no D+1 planes meet at a point
# (no `DegenerateCutError`); zero-measure pieces can appear. The nodal values are snapped
# (|φ| ≤ atol → 0) before anything else, so that the IN/OUT/CUT labels of the background
# cells use the same rule.

"""
    struct RecursivePlaneCutter <: Cutter

Cutter for `DiscreteGeometry` and `AnalyticalGeometry` (as `LevelSetCutter`): every cut cell
is split into simplices and every simplex is cut by all the level sets at once by the
recursive `CellCutter`.

    RecursivePlaneCutter(; atol=1e-12)

`atol`: absolute zero tolerance of the level-set values (level sets are meant to be distance
functions).

## Usage

    cut(RecursivePlaneCutter(),background::DiscreteModel,geom)
    cut_facets(RecursivePlaneCutter(),background::DiscreteModel,geom)
    compute_bgcell_to_inoutcut(RecursivePlaneCutter(),background::DiscreteModel,geom)
    compute_bgfacet_to_inoutcut(RecursivePlaneCutter(),background::DiscreteModel,geom)
"""
struct RecursivePlaneCutter <: Cutter
  atol::Float64
end

RecursivePlaneCutter(; atol::Real=1e-12) = RecursivePlaneCutter(Float64(atol))

function cut(cutter::RecursivePlaneCutter,background::DiscreteModel,geom)
  data = _rp_cut(cutter,get_grid(background),geom,Val(true))
  EmbeddedDiscretization(background, data..., geom)
end

function cut_facets(cutter::RecursivePlaneCutter,background::DiscreteModel,geom)
  D = num_cell_dims(background)
  grid = Grid(ReferenceFE{D-1},background)
  data = _rp_cut(cutter,grid,geom,Val(false))
  EmbeddedFacetDiscretization(background, data..., geom)
end

function compute_bgcell_to_inoutcut(cutter::RecursivePlaneCutter,model::DiscreteModel,geom)
  _rp_compute_bgcell_to_inoutcut(cutter,get_grid(model),geom)
end

function compute_bgfacet_to_inoutcut(cutter::RecursivePlaneCutter,model::DiscreteModel,geom)
  D = num_cell_dims(model)
  _rp_compute_bgcell_to_inoutcut(cutter,Grid(ReferenceFE{D-1},model),geom)
end

## Level-set values: one vector per unique leaf of the CSG tree, snapped (|φ| ≤ atol → 0)

function _rp_levelset_values(cutter::RecursivePlaneCutter,grid::Grid,geom::AnalyticalGeometry)
  _, oid_to_ls = _find_unique_leaves(get_tree(geom))
  ls_to_point_to_value, _ = _rp_levelset_values(cutter,grid,discretize(geom,grid))
  ls_to_point_to_value, oid_to_ls
end

function _rp_levelset_values(cutter::RecursivePlaneCutter,grid::Grid,geom::DiscreteGeometry)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(get_tree(geom))
  snap(v) = abs(v) <= cutter.atol ? zero(v) : v
  map(values -> map(snap,values),ls_to_point_to_value), oid_to_ls
end

function _rp_compute_bgcell_to_inoutcut(cutter::RecursivePlaneCutter,grid::Grid,geom)
  ugrid = UnstructuredGrid(grid)
  ls_to_point_to_value, oid_to_ls = _rp_levelset_values(cutter,ugrid,geom)
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
function _rp_cut(cutter::RecursivePlaneCutter,grid::Grid,geom,with_boundary::Val)
  ugrid = UnstructuredGrid(grid)
  ls_to_point_to_value, oid_to_ls = _rp_levelset_values(cutter,ugrid,geom)
  subsimplices, ls_to_subpoint_to_value, ls_to_bgcell_to_inoutcut =
    _initial_sub_triangulation(ugrid,ls_to_point_to_value)
  N = length(ls_to_point_to_value)
  out = _rp_cut_simplices(cutter,subsimplices,ls_to_subpoint_to_value,Val(N),with_boundary)
  (ls_to_bgcell_to_inoutcut, out..., oid_to_ls)
end

# Cuts every simplex of `st` (its points carry the level-set values) by the N level sets.
# Returns the subcells and their IN/OUT per level set, and, with boundary, the subfacets
# (normals towards φ > 0) and their IN/OUT/INTERFACE per level set.
function _rp_cut_simplices(
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
        _rp_push_oriented!(cells.points,c,verts,offset)
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
        n = invJt⋅get_face_normal(cf,f)
        n = n/norm(n)
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
  subcells = SubCellData(_table(cells.points,Dc+1),cells.bgcell,point_to_coords,point_to_rcoords)
  with_boundary || return (subcells, cells.ls_to_inout)
  subfacets = SubFacetData(_table(facets.points,Dc),facets.normal,facets.bgcell,point_to_coords,point_to_rcoords)
  (subcells, cells.ls_to_inout, subfacets, facets.ls_to_inout)
end

# The vertices of a Dc-simplex of the cut cell, swapping the last two if it is negatively
# oriented in reference coordinates
function _rp_push_oriented!(points,c::CutCell{D},verts,offset) where D
  x1 = get_vertex_coordinates(c,verts[1])
  J = TensorValue{D,D}(ntuple(i -> (get_vertex_coordinates(c,verts[(i-1)÷D+2])-x1)[(i-1)%D+1],Val(D*D)))
  flip = det(J) < 0
  for (i,v) in enumerate(verts)
    j = flip && i >= D ? (i == D ? D+1 : D) : i
    push!(points,offset+verts[j])
  end
end
