
function aggregate(strategy,cut::EmbeddedDiscretization)
  aggregate(strategy,cut,cut.geo)
end

function aggregate(strategy,cut::EmbeddedDiscretization,geo)
  aggregate(strategy,cut,geo,IN)
end

function aggregate(strategy,cut::EmbeddedDiscretization,name::String,in_or_out)
  geo = get_geometry(cut.geo,name)
  aggregate(strategy,cut,geo,in_or_out)
end

function aggregate(strategy,cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out)
  @abstractmethod
end

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  cut_facets::EmbeddedFacetDiscretization,
  name::String,
  in_or_out)

  geo = get_geometry(cut.geo,name)
  aggregate(strategy,cut,cut_facets,geo,in_or_out)
end

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  cut_facets::EmbeddedFacetDiscretization)

  aggregate(strategy,cut,cut_facets,cut.geo,IN)
end

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  cut_facets::EmbeddedFacetDiscretization,
  geo::CSG.Geometry)

  aggregate(strategy,cut,cut_facets,geo,IN)
end

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  facet_to_inoutcut::AbstractVector)

  aggregate(strategy,cut,cut.geo,IN,facet_to_inoutcut)
end

function aggregate(
  strategy,
  cut::EmbeddedDiscretization,
  cut_facets::EmbeddedFacetDiscretization,
  geo::CSG.Geometry,
  in_or_out)

  @abstractmethod
end

struct AggregateCutCellsByThreshold
  threshold::Float64
  AggregateCutCellsByThreshold(x) =
    ( x < 0.0 ) || ( x > 1.0 ) ? error("Agg. threshold must be in [0,1]") : new(x)
end

AggregateAllCutCells() = AggregateCutCellsByThreshold(1.0)

function aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::EmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out)

  facet_to_inoutcut = compute_bgfacet_to_inoutcut(cut.bgmodel,geo)
  _aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

function aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::EmbeddedDiscretization,
  cut_facets::EmbeddedFacetDiscretization,
  geo::CSG.Geometry,
  in_or_out)

  facet_to_inoutcut = compute_bgfacet_to_inoutcut(cut_facets,geo)
  _aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

function aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::EmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out,
  facet_to_inoutcut::AbstractVector)

  _aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

function _aggregate_by_threshold(threshold,cut,geo,loc,facet_to_inoutcut)
  @assert loc in (IN,OUT)

  cutinorout = loc == IN ? (CUT_IN,IN) : (CUT_OUT,OUT)
  trian = Triangulation(cut,cutinorout,geo)
  model = cut.bgmodel
  bgtrian = get_triangulation(model)
  cell_to_cut_meas = get_cell_measure(trian,bgtrian)
  cell_to_meas = get_cell_measure(bgtrian)
  cell_to_unit_cut_meas = lazy_map(/,cell_to_cut_meas,cell_to_meas)

  cell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)

  cell_to_coords = get_cell_coordinates(bgtrian)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  cell_to_faces = get_faces(topo,D,D-1)
  face_to_cells = get_faces(topo,D-1,D)

  _aggregate_by_threshold_barrier(
    threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
    loc,cell_to_coords,cell_to_faces,face_to_cells)
end

function _aggregate_by_threshold_barrier(
  threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
  loc,cell_to_coords,cell_to_faces,face_to_cells)

  n_cells = length(cell_to_unit_cut_meas)
  cell_to_cellin = zeros(Int32,n_cells)
  cell_to_touched = fill(false,n_cells)

  for cell in 1:n_cells
    if cell_to_unit_cut_meas[cell] ≥ threshold
      cell_to_cellin[cell] = cell
      cell_to_touched[cell] = true
    end
  end

  c1 = array_cache(cell_to_faces)
  c2 = array_cache(face_to_cells)
  c3 = array_cache(cell_to_coords)
  c4 = array_cache(cell_to_coords)

  max_iters = 20

  all_aggregated = false
  for iter in 1:max_iters
    all_aggregated = true
    for cell in 1:n_cells
      if ! cell_to_touched[cell] && cell_to_inoutcut[cell] == CUT
        neigh_cell = _find_best_neighbor(
          c1,c2,c3,c4,cell,
          cell_to_faces,
          face_to_cells,
          cell_to_coords,
          cell_to_touched,
          cell_to_cellin,
          facet_to_inoutcut,
          loc)
        if neigh_cell > 0
          cellin = cell_to_cellin[neigh_cell]
          cell_to_cellin[cell] = cellin
        else
          all_aggregated = false
        end
      end
    end
    if all_aggregated
      break
    end
    _touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
  end

  @assert all_aggregated

  cell_to_cellin
end

function _find_best_neighbor(
  c1,c2,c3,c4,cell,
  cell_to_faces,
  face_to_cells,
  cell_to_coords,
  cell_to_touched,
  cell_to_cellin,
  facet_to_inoutcut,
  loc)

  faces = getindex!(c1,cell_to_faces,cell)
  dmin = Inf
  T = eltype(eltype(face_to_cells))
  best_neigh_cell = zero(T)
  i_to_coords = getindex!(c3,cell_to_coords,cell)
  for face in faces
    inoutcut = facet_to_inoutcut[face]
    if  inoutcut != CUT && inoutcut != loc
      continue
    end
    neigh_cells = getindex!(c2,face_to_cells,face)
    for neigh_cell in neigh_cells
      if neigh_cell != cell && cell_to_touched[neigh_cell]
        cellin = cell_to_cellin[neigh_cell]
        j_to_coords = getindex!(c4,cell_to_coords,cellin)
        d = 0.0
        for p in i_to_coords
          for q in j_to_coords
            d = max(d,Float64(norm(p-q)))
          end
        end
        if (1.0+1.0e-9)*d < dmin
          dmin = d
          best_neigh_cell = neigh_cell
        end
      end
    end
  end
  best_neigh_cell
end

function _touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
  for (cell,cellin) in enumerate(cell_to_cellin)
    if cellin > 0
      cell_to_touched[cell] = true
    end
  end
end

function color_aggregates(cell_to_cellin,model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  cell_to_faces = get_faces(topo,D,D-1)
  face_to_cells = get_faces(topo,D-1,D)
  _color_aggregates_barrier(cell_to_cellin,cell_to_faces,face_to_cells)
end

function _color_aggregates_barrier(cell_to_cellin,cell_to_faces,face_to_cells)

  n_cells = length(cell_to_cellin)
  cell_to_isin = fill(false,n_cells)
  for cellin in cell_to_cellin
    if cellin>0
      cell_to_isin[cellin] = true
    end
  end

  n_incell = count(cell_to_isin)
  g = SimpleGraph(n_incell)

  incell_to_cell = findall(cell_to_isin)
  cell_to_incell = zeros(Int,n_cells)
  cell_to_incell[cell_to_isin] .= 1:n_incell

  c1 = array_cache(cell_to_faces)
  c2 = array_cache(face_to_cells)

  for cell in 1:n_cells
    cellin = cell_to_cellin[cell]
    if cellin > 0
      incell = cell_to_incell[cellin]
      faces = getindex!(c1,cell_to_faces,cell)
      for face in faces
        neigh_cells = getindex!(c2,face_to_cells,face)
        for neigh_cell in neigh_cells
          neigh_cellin = cell_to_cellin[neigh_cell]
          if neigh_cellin > 0
            neigh_incell = cell_to_incell[neigh_cellin]
            if neigh_incell != incell
              add_edge!(g,incell,neigh_incell)
            end
          end
        end
      end
    end
  end

  c = greedy_color(g)

  incell_to_nslaves = zeros(Int,n_incell)
  for cell in 1:n_cells
    cellin = cell_to_cellin[cell]
    if cellin > 0
      incell = cell_to_incell[cellin]
      incell_to_nslaves[incell] += 1
    end
  end

  incell_to_color = c.colors
  cell_to_color = zeros(Int,n_cells)
  for cell in 1:n_cells
    cellin = cell_to_cellin[cell]
    if cellin > 0
      incell = cell_to_incell[cellin]
      if incell_to_nslaves[incell] > 1
        color = incell_to_color[incell]
        cell_to_color[cell] = color
      end
    end
  end

  cell_to_color
end
