
#int qh_new_qhull(int dim, int numpoints, coordT *points, boolT ismalloc,
#                 char *qhull_cmd, FILE *outfile, FILE *errfile) 
function qh_new_qhull(  dim::Int32,
                        numpoints::Int32,
                        points::Vector{Float64},ismalloc::Bool,
                        qhull_cmd::String)
  @check_if_loaded
  outfile = ccall(:fopen, Ptr{Cvoid}, (Ptr{UInt8}, Ptr{UInt8}), "out", "rw+")
  errfile = ccall(:fopen, Ptr{Cvoid}, (Ptr{UInt8}, Ptr{UInt8}), "err", "rw+")
  ccall(qh_new_qhull_ptr[], 
        Cint, (
            Cint, 
            Cint,
            Ptr{Cdouble},
            Cuint,
            Ptr{UInt8},
            Ptr{Cvoid},
            Ptr{Cvoid}
            ), 
        dim, 
        numpoints, 
        points, 
        ismalloc, 
        qhull_cmd, 
        outfile, 
        errfile)
end


#void qh_triangulate(void /*qh facet_list*/) 
function qh_triangulate()
    @check_if_loaded
    ccall(qh_triangulate_ptr[], Cvoid, ())
end

#void qh_memfreeshort(int *curlong, int *totlong) 
function qh_memfreeshort(curlong::Int32, totlong::Int32)
    @check_if_loaded
    ccall(qh_memfreeshort_ptr[], 
        Cvoid, (
            Ptr{Cint}, 
            Ptr{Cint}
            ), 
        Ref(Cint(curlong)), 
        Ref(Cint(totlong)))
end

#void qh_freeqhull(boolT allmem) 
#define boolT unsigned int
#define qh_ALL      True
#define False 0
#define True 1
function qh_freeqhull(allmem::Bool)
  @check_if_loaded
  ccall(qh_freeqhull_ptr[], Cvoid, (Ptr{Cuint},), Ref(Cuint(allmem)))
end

#int qh_pointid(pointT *point) 
#define pointT   coordT
#define coordT   realT
#define realT   double
function qh_pointid(point::Vector{Float64})
  @check_if_loaded
  ccall(qh_pointid_ptr[], Cint, (Ptr{Cdouble},), point)
end


