module TriangulationsTests

using GridapEmbedded
using GridapEmbedded: IN, OUT, CUT
using Gridap.Arrays

@inline function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector,
  cell::Integer)

  case = compute_case(cell_to_points,point_to_value,cell)
  table.case_to_inoutcut[case]
end

@inline function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  ls_to_point_to_value::Vector{<:AbstractVector},
  cell::Integer)

  for point_to_value in ls_to_point_to_value
    if OUT == compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
      return OUT
    end
  end
  for point_to_value in ls_to_point_to_value
    if CUT == compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
      return CUT
    end
  end
  return IN
end

function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector)

  ncells = length(cell_to_points)
  cell_to_inoutcut = zeros(Int8,ncells)
  for cell in 1:ncells
    cell_to_inoutcut[cell] = compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
  end
  cell_to_inoutcut
end


using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded

table = LookupTable(QUAD)

domain = (0,1,0,1)
n = 50
partition =(n,n)
grid = UnstructuredGrid(CartesianGrid(domain,partition))

point_to_coords = get_node_coordinates(grid)
cell_to_points = get_cell_nodes(grid)

const R1 = 0.7
ls1(x) = x[1]^2 + x[2]^2 - R1

const R2 = 0.7
ls2(x) = x[1]^2 + (x[2]-1)^2 - R1

point_to_ls1 = ls1.(point_to_coords)
point_to_ls2 = ls2.(point_to_coords)

ls_to_point_to_value = [point_to_ls1, point_to_ls2]

cell_to_inoutcut = compute_in_out_or_cut(table,cell_to_points,ls_to_point_to_value)

using Gridap.Visualization
write_vtk_file(grid,"grid",celldata=["inoutcut"=>cell_to_inoutcut])

cutcell_to_cell = findall(cell_to_inoutcut .== CUT)
cutgrid = simplexify(GridPortion(grid,cutcell_to_cell))
write_vtk_file(cutgrid,"cutgrid")




end #module
