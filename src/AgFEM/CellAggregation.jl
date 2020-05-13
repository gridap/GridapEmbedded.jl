
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

struct AggregateAllCutCells end

function aggregate(strategy::AggregateAllCutCells,cut::EmbeddedDiscretization,geo::CSG.Geometry,in_or_out)
  cell_to_inoutcut = compute_bgcell_to_inoutcut(cut,geo)
  _compute_aggregates(cell_to_inoutcut,cut.bgmodel,in_or_out)
end

function _compute_aggregates(cell_to_inoutcut,model,loc)
  @assert loc in (IN,OUT)
  trian = get_triangulation(model)
  cell_to_coords = get_cell_coordinates(trian)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  cell_to_faces = get_faces(topo,D,D-1)
  face_to_cells = get_faces(topo,D-1,D)
  _compute_aggregates_barrier(
    cell_to_inoutcut,loc,cell_to_coords,cell_to_faces,face_to_cells)
end

function _compute_aggregates_barrier(
  cell_to_inoutcut,loc,cell_to_coords,cell_to_faces,face_to_cells)

  n_cells = length(cell_to_inoutcut)
  cell_to_cellin = zeros(Int32,n_cells)
  cell_to_touched = fill(false,n_cells)

  for (cell, inoutcut) in enumerate(cell_to_inoutcut)
    if inoutcut == loc
      cell_to_cellin[cell] = cell
      cell_to_touched[cell] = true
    end
  end

  c1 = array_cache(cell_to_faces)
  c2 = array_cache(face_to_cells)

  max_iters = 20

  all_aggregated = false
  for iter in 1:max_iters
    all_aggregated = true
    for (cell, inoutcut) in enumerate(cell_to_inoutcut)
      if inoutcut == CUT && ! cell_to_touched[cell]
        neigh_cell = _find_best_neighbor(c1,c2,cell,cell_to_faces,face_to_cells,cell_to_touched)
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

function _find_best_neighbor(c1,c2,cell,cell_to_faces,face_to_cells,cell_to_touched)
  #TODO skip faces that are out
  #TODO take into account cell coordinates
  faces = getindex!(c1,cell_to_faces,cell)
  for face in faces
    neigh_cells = getindex!(c2,face_to_cells,face)
    for neigh_cell in neigh_cells
      if neigh_cell != cell && cell_to_touched[neigh_cell]
        return neigh_cell
      end
    end
  end
  T = eltype(eltype(face_to_cells))
  zero(T)
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

