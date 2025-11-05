"""
    struct NonConformingGridTopology{Dc,Dp,T,O} <: GridTopology{Dc,Dp}
      # private fields
    end
"""
struct NonConformingGridTopology{Dc,Dp,T,O,A,B,C,D,E,F,G} <: GridTopology{Dc,Dp}
  ugt::UnstructuredGridTopology{Dc,Dp,T,O}
  ncg::NonConformingGlue{Dc,A,B,C,D,E,F,G}
  coarse_cell_to_hanging_faces::Table{Int32,Int32}
  hanging_faces_to_coarse_cell::Table{Int32,Int32}
end

# Implementation of abstract API

OrientationStyle(
  ::Type{<:NonConformingGridTopology{Dc,Dp,T,O,A,B,C,D,E,F,G}}) where 
  {Dc,Dp,T,O,A,B,C,D,E,F,G} = O()

RegularityStyle(
  ::Type{<:NonConformingGridTopology{Dc,Dp,T,O,A,B,C,D,E,F,G}}) where 
  {Dc,Dp,T,O,A,B,C,D,E,F,G} = Irregular()

get_vertex_coordinates(ncgt::NonConformingGridTopology) = ncgt.ugt.vertex_coordinates

get_cell_type(ncgt::NonConformingGridTopology) = ncgt.ugt.cell_type

get_polytopes(ncgt::NonConformingGridTopology) = collect1d(ncgt.ugt.polytopes)

# Can I do this?

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O,A,B,C,D,E,F,G}, 
                   ::Val{Dc}, ::Val{Dc-1}) where {Dc,Dp,T,O,A,B,C,D,E,F,G}
  
end

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O,A,B,C,D,E,F,G}, 
                   ::Val{Dc-1}, ::Val{Dc}) where {Dc,Dp,T,O,A,B,C,D,E,F,G}
  
end

# What other procs do I need? Review CA (e.g., facet_to_inoutcut)