using Gridap
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Test

# Background Mesh
domain = (-0.01,1.01,-0.01,1.01)
n_x = 2
n_t = 2
partition = (n_x,n_t)
bgmodel = CartesianDiscreteModel(domain,partition)

# domain
geo = quadrilateral(x0=Point(0,0),d1=VectorValue(1,0),d2=VectorValue(0,1))
cutgeo = cut(bgmodel,geo)
model = DiscreteModel(cutgeo)

# AgFEM
#stategy = AggregateSpaceCutCells()
#aggregates = aggregatespace(strategy,cutgeo)


trian = get_triangulation(model)
cell_to_coords = get_cell_coordinates(trian)
topo = get_grid_topology(model)
D = num_cell_dims(model)

n_cell = num_cells(topo)
filter = [false,false,true,true]
filters = fill(filter,n_cell)

cell_to_faces = get_faces(topo,D,D-1)
face_to_cells = get_faces(topo,D-1,D)

# Filter Cell_to_Faces
ptrs = collect(1:n_cell)
a = CompressedArray(cell_to_faces,ptrs)
fa = FilteredCellArray(a,filters)
@test fa == [[3,4],[4,7],[9,10],[10,12]]

# Extracting only the faces to be aggregated
face = get_faces(topo,D-1,D-1)
n_faces = length(face)
filter_f = [false]
filters_f = fill(filter_f,n_faces)

for i in 1:n_faces
    for j in 1:length(fa)
        for k in 1:2
            if face[i][1]==fa[j][k]
                filters_f[i] = [true]
            end
        end
    end
end

ptrs_f = collect(1:n_faces)
b = CompressedArray(face,ptrs_f)
fb = FilteredCellArray(b,filters_f)
@test fb == [[],[],[3],[4],[],[],[7],[],[9],[10],[],[12]]

# Filter Faces_to_cells
ptrs_fc = collect(1:n_faces)
c = CompressedArray(face_to_cells,ptrs_fc)

filters_fc = deepcopy(filters_f)

for i in 1:n_faces
    l = length(face_to_cells[i])
    println(l)
    for j in 2:l
        append!(filters_fc[i],filters_f[i][1])
        println(filters_fc[i])
    end
end

@test filters_fc == [[false],[false,false],[true],[true,true],[false],[false,false],[true],[false],[true],[true,true],[false],[true]]

fc = FilteredCellArray(c,filters_fc)
