
struct EmbeddedDiscretization{Dp,T} <: GridapType
  bgmodel::DiscreteModel
  ls_to_bgcell_to_inoutcut::Vector{Vector{Int8}}
  subcells::SubTriangulation{Dp,T}
  ls_to_cell_to_inout::Vector{Vector{Int8}}
  ls_to_subfacets::Vector{FacetSubTriangulation{Dp,T}}
  ls_to_ls_to_facet_to_inout::Vector{Vector{Vector{Int8}}}
  oid_to_ls::Dict{UInt,Int}
  geo::CSG.Geometry
end

function DiscreteModel(cut::EmbeddedDiscretization)
  DiscreteModel(cut,cut.geo)
end

function DiscreteModel(cut::EmbeddedDiscretization,geo)
  DiscreteModel(cut,geo,(IN,CUT))
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
  else
    @error "operation $op not implemented"
  end
end

function _compute_inoutcut_union(inout_1,inout_2)
  inout_12 = (inout_1,inout_2)
  if (inout_1==IN) || (inout_2==IN)
    Int8(IN)
  elseif (inout_12 == (CUT,CUT)) || (inout_12 == (CUT,OUT)) || (inout_12 == (IN,OUT))
    Int8(CUT)
  else
    Int8(OUT)
  end
end

function _compute_inoutcut_intersection(inout_1,inout_2)
  inout_12 = (inout_1,inout_2)
  if (inout_1==OUT) || (inout_2==OUT)
    Int8(OUT)
  elseif (inout_12 == (CUT,CUT)) || (inout_12 == (CUT,IN)) || (inout_12 == (IN,CUT))
    Int8(CUT)
  else
    Int8(IN)
  end
end






