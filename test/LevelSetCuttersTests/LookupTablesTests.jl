module LookupTablesTests

using Test
using Gridap

using GridapEmbedded.LevelSetCutters: compute_case, LookupTable

@test compute_case([0,0,0]) == 1
@test compute_case([1,0,0]) == 2
@test compute_case([1,0,1]) == 6
@test compute_case([1,1,1]) == 8

table = LookupTable(SEGMENT)
table = LookupTable(TRI)
table = LookupTable(QUAD)
table = LookupTable(HEX)
table = LookupTable(TET)

end # module
