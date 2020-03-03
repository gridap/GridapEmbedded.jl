module LookupTablesTests

using Test
using GridapEmbedded
using Gridap

@test num_cases(3) == 8
@test compute_case([0,0,0]) == 1
@test compute_case([1,0,0]) == 2
@test compute_case([1,0,1]) == 6
@test compute_case([1,1,1]) == 8

table = LookupTable(TRI)
table = LookupTable(QUAD)
table = LookupTable(HEX)
table = LookupTable(TET)

end # module
