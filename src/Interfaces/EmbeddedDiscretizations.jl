
struct EmbeddedDiscretization{Dp,T} <: GridapType
  bgmodel::DiscreteModel
  ls_to_bgcell_to_inoutcut::Vector{Vector{Int8}}
  subcells::SubTriangulation{Dp,T}
  ls_to_subcell_to_inout::Vector{Vector{Int8}}
  subfacets::FacetSubTriangulation{Dp,T}
  ls_to_subfacet_to_inout::Vector{Vector{Int8}}
  oid_to_ls::Dict{UInt,Int}
  geo::CSG.Geometry
end

function DiscreteModel(cut::EmbeddedDiscretization)
  DiscreteModel(cut,cut.geo)
end

function DiscreteModel(cut::EmbeddedDiscretization,geo)
  DiscreteModel(cut,geo,(IN,CUT))
end

function DiscreteModel(cut::EmbeddedDiscretization,name::String,in_or_out)
  geo = get_geometry(cut.geo,name)
  DiscreteModel(cut,geo,in_or_out)
end

function DiscreteModel(cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out)
  pred = i-> i in in_or_out
  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  cell_list = findall(pred, bgcell_to_inoutcut)
  DiscreteModel(cut.bgmodel,cell_list)
end

function compute_bgcell_to_inoutcut(cut::EmbeddedDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_bgcell_to_inoutcut[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

function compute_inoutcut(a::Leaf)
  first(a.data)
end

function compute_inoutcut(a::Node)
  cell_to_inoutcut_1 = compute_inoutcut(a.leftchild)
  cell_to_inoutcut_2 = compute_inoutcut(a.rightchild)
  op = first(a.data)
  if op  == :∪
    return _compute_inoutcut_union.(cell_to_inoutcut_1,cell_to_inoutcut_2)
  elseif op == :∩
    return _compute_inoutcut_intersection.(cell_to_inoutcut_1,cell_to_inoutcut_2)
  elseif op == :-
    return _compute_inoutcut_setdiff.(cell_to_inoutcut_1,cell_to_inoutcut_2)
  else
    @error "operation $op not implemented"
  end
end

function compute_inoutcut(a::UnaryNode)
  cell_to_inoutcut_1 = compute_inoutcut(a.leftchild)
  op = first(a.data)
  if op  == :!
    return _compute_inoutcut_complementary.(cell_to_inoutcut_1)
  else
    @error "operation $op not implemented"
  end
end

function _compute_inoutcut_union(inout_1,inout_2)
  inout_12 = (inout_1,inout_2)
  if (inout_1==OUT) && (inout_2==OUT)
    Int8(OUT)
  elseif (inout_12 == (CUT,CUT)) || (inout_12 == (CUT,OUT)) || (inout_12 == (OUT,CUT))
    Int8(CUT)
  else
    Int8(IN)
  end
end

function _compute_inoutcut_intersection(inout_1,inout_2)
  inout_12 = (inout_1,inout_2)
  if (inout_1==IN) && (inout_2==IN)
    Int8(IN)
  elseif (inout_12 == (CUT,CUT)) || (inout_12 == (CUT,IN)) || (inout_12 == (IN,CUT))
    Int8(CUT)
  else
    Int8(OUT)
  end
end

function _compute_inoutcut_setdiff(inout_1,inout_2)
  inout_12 = (inout_1,inout_2)
  if (inout_1==IN) && (inout_2==OUT)
    Int8(IN)
  elseif (inout_12 == (CUT,CUT)) || (inout_12 == (CUT,OUT)) || (inout_12 == (IN,CUT))
    Int8(CUT)
  else
    Int8(OUT)
  end
end

function _compute_inoutcut_complementary(inout_1)
  if (inout_1==OUT)
    Int8(IN)
  elseif (inout_1==IN)
    Int8(OUT)
  else
    Int8(CUT)
  end
end

function compute_inoutboundary(a::Leaf)
  d = first(a.data)
  d, ones(Int8,length(d))
end

function compute_inoutboundary(a::Node)
  cell_to_inoutcut_1, orientation1 = compute_inoutboundary(a.leftchild)
  cell_to_inoutcut_2, orientation2 = compute_inoutboundary(a.rightchild)
  op = first(a.data)
  if op  == :∪
    r1 = _compute_inoutcut_union.(cell_to_inoutcut_1,cell_to_inoutcut_2)
    r2 = _compute_orientation_union.(cell_to_inoutcut_1,cell_to_inoutcut_2,orientation1,orientation2)
    return r1, r2
  elseif op == :∩
    r1 = _compute_inoutcut_intersection.(cell_to_inoutcut_1,cell_to_inoutcut_2)
    r2 = _compute_orientation_intersect.(cell_to_inoutcut_1,cell_to_inoutcut_2,orientation1,orientation2)
    return r1, r2
  elseif op == :-
    r1 = _compute_inoutcut_setdiff.(cell_to_inoutcut_1,cell_to_inoutcut_2)
    r2 = _compute_orientation_setdiff.(cell_to_inoutcut_1,cell_to_inoutcut_2,orientation1,orientation2)
    return r1, r2
  else
    @error "operation $op not implemented"
  end
end

function compute_inoutboundary(a::UnaryNode)
  cell_to_inoutboundary_1, orientation1 = compute_inoutboundary(a.leftchild)
  op = first(a.data)
  if op  == :!
    r1 = _compute_inoutcut_complementary.(cell_to_inoutboundary_1)
    r2 = _compute_orientation_complementary.(cell_to_inoutboundary_1,orientation1)
    return r1, r2
  else
    @error "operation $op not implemented"
  end
end

function _compute_orientation_union(iob1,iob2,d1,d2)
  if iob1 == INTERFACE
    d1
  elseif iob2 == INTERFACE
    d2
  else
    d1
  end
end

const _compute_orientation_intersect = _compute_orientation_union

function _compute_orientation_setdiff(iob1,iob2,d1,d2)
  if iob1 == INTERFACE
    d1
  elseif iob2 == INTERFACE
    -d2
  else
    d1
  end
end

function _compute_orientation_complementary(iob1,d1)
  if iob1 == INTERFACE
    -d1
  else
    d1
  end
end

function Triangulation(cut::EmbeddedDiscretization)
  Triangulation(cut,cut.geo)
end

function Triangulation(cut::EmbeddedDiscretization,geo)
  Triangulation(cut,geo,(CUTIN,IN))
end

function Triangulation(cut::EmbeddedDiscretization,name::String,in_or_out)
  geo = get_geometry(cut.geo,name)
  Triangulation(cut,geo,in_or_out)
end

function Triangulation(cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out::Tuple)
  trian1 = Triangulation(cut,geo,in_or_out[1])
  trian2 = Triangulation(cut,geo,in_or_out[2])
  lazy_append(trian1,trian2)
end

function Triangulation(cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out::CutInOrOut)
  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  subcell_to_inoutcut = reindex(bgcell_to_inoutcut,cut.subcells.cell_to_bgcell)
  subcell_to_inout = compute_subcell_to_inout(cut,geo)
  mask = apply( (a,b) -> a==CUT && b==in_or_out.in_or_out, subcell_to_inoutcut, subcell_to_inout   )
  newsubcells = findall(mask)
  st = SubTriangulation(cut.subcells,newsubcells)
  SubTriangulationWrapper(st)
end

function Triangulation(cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out::Integer)
  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  cell_to_mask = collect(Bool,bgcell_to_inoutcut .== in_or_out)
  trian = Triangulation(cut.bgmodel)
  RestrictedTriangulation(trian,cell_to_mask)
end

function compute_subcell_to_inout(cut::EmbeddedDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    subcell_to_inout = cut.ls_to_subcell_to_inout[ls]
    subcell_to_inout, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

function compute_inout(a::Leaf)
  first(a.data)
end

function compute_inout(a::Node)
  cell_to_inout_1 = compute_inout(a.leftchild)
  cell_to_inout_2 = compute_inout(a.rightchild)
  op = first(a.data)
  if op  == :∪
    return _compute_inout_union.(cell_to_inout_1,cell_to_inout_2)
  elseif op == :∩
    return _compute_inout_intersection.(cell_to_inout_1,cell_to_inout_2)
  elseif op == :-
    return _compute_inout_setdiff.(cell_to_inout_1,cell_to_inout_2)
  else
    @error "operation $op not implemented"
  end
end

function compute_inout(a::UnaryNode)
  cell_to_inout_1 = compute_inout(a.leftchild)
  op = first(a.data)
  if op  == :!
    return _compute_inout_complementary.(cell_to_inout_1)
  else
    @error "operation $op not implemented"
  end
end

function _compute_inout_union(inout_1,inout_2)
  if (inout_1==IN) || (inout_2==IN)
    Int8(IN)
  else
    Int8(OUT)
  end
end

function _compute_inout_intersection(inout_1,inout_2)
  if (inout_1==IN) && (inout_2==IN)
    Int8(IN)
  else
    Int8(OUT)
  end
end

function _compute_inout_setdiff(inout_1,inout_2)
  if (inout_1==IN) && (inout_2==OUT)
    Int8(IN)
  else
    Int8(OUT)
  end
end

function _compute_inout_complementary(inout_1)
  if (inout_1==IN)
    Int8(OUT)
  else
    Int8(IN)
  end
end

function EmbeddedBoundary(cut::EmbeddedDiscretization)
  EmbeddedBoundary(cut,cut.geo)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,name::String)
  geo2 = get_geometry(cut.geo,name)
  EmbeddedBoundary(cut,geo2)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,geo::CSG.Geometry)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_subfacet_to_inout[ls]
    cell_to_inoutcut, name, meta
  end

  tree = get_tree(geo)
  newtree = replace_data(identity,conversion,tree)
  subfacet_to_inoutcut, orientation = compute_inoutboundary(newtree)
  newsubfacets = findall(subfacet_to_inoutcut .== INTERFACE)
  neworientation = orientation[newsubfacets]
  fst = FacetSubTriangulation(cut.subfacets,newsubfacets,neworientation)
  FacetSubTriangulationWrapper(fst)

end

function EmbeddedBoundary(cut::EmbeddedDiscretization,name1::String,name2::String)
  geo1 = get_geometry(cut.geo,name1)
  geo2 = get_geometry(cut.geo,name2)
  EmbeddedBoundary(cut,geo1,geo2)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,geo1::CSG.Geometry,geo2::CSG.Geometry)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_subfacet_to_inout[ls]
    cell_to_inoutcut, name, meta
  end

  tree1 = get_tree(geo1)
  tree2 = get_tree(geo2)
  newtree1 = replace_data(identity,conversion,tree1)
  newtree2 = replace_data(identity,conversion,tree2)
  subfacet_to_inoutcut1, orientation = compute_inoutboundary(newtree1)
  subfacet_to_inoutcut2 = compute_inoutcut(newtree2)
  mask = apply( (i,j)->(i==INTERFACE) && (j==INTERFACE), subfacet_to_inoutcut1, subfacet_to_inoutcut2 )
  newsubfacets = findall( mask )
  neworientation = orientation[newsubfacets]
  fst = FacetSubTriangulation(cut.subfacets,newsubfacets,neworientation)
  FacetSubTriangulationWrapper(fst)

end

function GhostSkeleton(cut::EmbeddedDiscretization)
  GhostSkeleton(cut,cut.geo)
end

function GhostSkeleton(cut::EmbeddedDiscretization,geo)
  GhostSkeleton(cut,geo,IN)
end

function GhostSkeleton(cut::EmbeddedDiscretization,name::String,in_or_out)
  geo = get_geometry(cut.geo,name)
  GhostSkeleton(cut,geo,in_or_out)
end

function GhostSkeleton(cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out)

  @assert in_or_out in (IN,OUT)
  cell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  model = cut.bgmodel
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  facet_to_cells = Table(get_faces(topo,D-1,D))
  facet_to_mask = fill(false,length(facet_to_cells))
  _fill_ghost_skeleton_mask!(facet_to_mask,facet_to_cells,cell_to_inoutcut,in_or_out)

  SkeletonTriangulation(model,facet_to_mask)
end

function _fill_ghost_skeleton_mask!(facet_to_mask,facet_to_cells::Table,cell_to_inoutcut,in_or_out)

  nfacets = length(facet_to_cells)
  for facet in 1:nfacets
    a = facet_to_cells.ptrs[facet]
    b = facet_to_cells.ptrs[facet+1]
    ncells_around = b-a
    ncells_around_cut = 0
    ncells_around_active = 0
    for cell_around in 1:ncells_around
      cell = facet_to_cells.data[a-1+cell_around]
      inoutcut = cell_to_inoutcut[cell]
      if (inoutcut == CUT)
        ncells_around_cut += 1
      end
      if (inoutcut == CUT) || (inoutcut == in_or_out)
        ncells_around_active += 1
      end
    end
    if (ncells_around_cut >0) && (ncells_around_active == 2)
      facet_to_mask[facet] = true
    end
  end

end

