module GridapEmbedded

using Libdl


#int qh_new_qhull(int dim, int numpoints, coordT *points, boolT ismalloc,
#                 char *qhull_cmd, FILE *outfile, FILE *errfile) 
export qh_new_qhull
export qh_triangulate
export qh_memfreeshort
export qh_freeqhull
export qh_pointid

include("load.jl")
include("bindings.jl")

end # module
