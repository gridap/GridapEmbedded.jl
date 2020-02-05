deps_jl = joinpath(@__DIR__, "..", "deps", "deps.jl")

if !isfile(deps_jl)
  error("Package GridapEmbedded not installed properly.")
end

include(deps_jl)

const qh_pointid_ptr       = Ref{Ptr}()
const qh_freeqhull_ptr     = Ref{Ptr}()
const qh_memfreeshort_ptr  = Ref{Ptr}()
const qh_new_qhull_ptr     = Ref{Ptr}()
const qh_triangulate_ptr   = Ref{Ptr}()
const QHULL_LOADED         = Ref(false)

function __init__()
    if QHULL_FOUND
        flags = Libdl.RTLD_LAZY | Libdl.RTLD_DEEPBIND | Libdl.RTLD_GLOBAL
        QHULL = Libdl.dlopen(QHULL_LIB_PATH, flags)

        # Initialization / Finalization
        GridapEmbedded.qh_pointid_ptr[]      = Libdl.dlsym(QHULL,:qh_pointid)
        GridapEmbedded.qh_freeqhull_ptr[]    = Libdl.dlsym(QHULL,:qh_freeqhull)
        GridapEmbedded.qh_memfreeshort_ptr[] = Libdl.dlsym(QHULL,:qh_memfreeshort)
        GridapEmbedded.qh_new_qhull_ptr[]    = Libdl.dlsym(QHULL,:qh_new_qhull)
        GridapEmbedded.qh_triangulate_ptr[]  = Libdl.dlsym(QHULL,:qh_triangulate)

        QHULL_LOADED[] = true
    end
end


macro check_if_loaded()
  quote
    if ! QHULL_LOADED[]
      error("QHULL is not properly loaded")
    end
  end
end


