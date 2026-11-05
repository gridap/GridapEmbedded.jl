module LookupTablesTests

using Test
using Gridap

using GridapEmbedded.Interfaces: Simplex, CUT
using GridapEmbedded.LevelSetCutters: compute_case, LookupTable, get_lookup_table

@test compute_case([0,0,0]) == 1
@test compute_case([1,0,0]) == 2
@test compute_case([1,0,1]) == 6
@test compute_case([1,1,1]) == 8

table = LookupTable(SEGMENT)
table = LookupTable(TRI)
table = LookupTable(QUAD)
table = LookupTable(HEX)
table = LookupTable(TET)

# Cached tables: one per polytope, also for polytopes built anew (Polytope hashes by identity)
@test get_lookup_table(TRI) === get_lookup_table(Simplex(Val{2}()))
@test get_lookup_table(Simplex(Val{1}())) === get_lookup_table(Simplex(Val{1}()))
@test get_lookup_table(TET) !== get_lookup_table(TRI)
tables = fetch.([Threads.@spawn get_lookup_table(QUAD) for _ in 1:4])
@test all(t -> t === get_lookup_table(QUAD), tables)

# Points of a case: the vertices, then one per cut edge; the subcells use all of them
for p in (SEGMENT,TRI,TET)
  table = get_lookup_table(p)
  for case in 1:2^num_vertices(p)
    @test table.case_to_inoutcut[case] == CUT || isempty(table.case_to_cut_edges[case])
    npoints = num_vertices(p) + length(table.case_to_cut_edges[case])
    @test sort(unique(reduce(vcat,table.case_to_subcell_to_points[case]))) == 1:npoints
  end
end

end # module
