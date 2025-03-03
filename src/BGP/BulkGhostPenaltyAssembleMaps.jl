# TO-DO: Better name?
struct BulkGhostPenaltyAssembleLhsMap{A} <: Gridap.Fields.Map
    agg_cells_lhs_contribs::A
end

function _get_rank(::Type{Array{T,N}}) where {T,N}
  N
end

function Gridap.Fields.return_cache(m::BulkGhostPenaltyAssembleLhsMap,cells)
    cache_unassembled_lhs=array_cache(m.agg_cells_lhs_contribs)
    T=eltype(m.agg_cells_lhs_contribs)
    evaluate_result=Gridap.Arrays.CachedArray(eltype(T),_get_rank(T))
    cache_unassembled_lhs,evaluate_result
end

function Gridap.Fields.evaluate!(cache,m::BulkGhostPenaltyAssembleLhsMap,cells)
    cache_unassembled_lhs,result=cache
    contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,1)

    Gridap.Arrays.setsize!(result,size(contrib))
    result.array .= 0.0
    for (i,cell) in enumerate(cells)
        contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,cell)
        result.array .+= contrib
    end
    result.array
end
struct BulkGhostPenaltyAssembleRhsMap{A,B} <: Gridap.Fields.Map
    agg_cells_local_dof_ids::A
    agg_cells_rhs_contribs::B
end

function Gridap.Fields.return_cache(m::BulkGhostPenaltyAssembleRhsMap,aggregate_local_cells)
  cache_agg_cells_local_dof_ids=array_cache(m.agg_cells_local_dof_ids)
  cache_unassembled_rhs=array_cache(m.agg_cells_rhs_contribs)
  evaluate_result=Gridap.Arrays.CachedArray(eltype(eltype(m.agg_cells_rhs_contribs)),2)
  cache_agg_cells_local_dof_ids,cache_unassembled_rhs,evaluate_result
end

function Gridap.Fields.evaluate!(cache,m::BulkGhostPenaltyAssembleRhsMap,aggregate_local_cells)
  cache_agg_cells_local_dof_ids,cache_unassembled_rhs,result=cache
  contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,1)

  max_local_dof_id=-1
  for (i,cell) in enumerate(aggregate_local_cells)
    current_cell_local_dof_ids = getindex!(cache_agg_cells_local_dof_ids,m.agg_cells_local_dof_ids,cell)
    for local_dof in current_cell_local_dof_ids
        max_local_dof_id=max(max_local_dof_id,local_dof)
    end
  end

  Gridap.Arrays.setsize!(result,(size(contrib,1),max_local_dof_id))

  result.array .= 0.0
  for (i,cell) in enumerate(aggregate_local_cells)
        current_cell_local_dof_ids = getindex!(cache_agg_cells_local_dof_ids,m.agg_cells_local_dof_ids,cell)
        contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,cell)
        for (j,local_dof) in enumerate(current_cell_local_dof_ids)
          result.array[:,local_dof] += contrib[:,j]
        end
  end
  result.array
end

struct BulkGhostPenaltyAssembleRhsFEFunctionMap{A} <: Gridap.Fields.Map
    agg_cells_rhs_contribs::A
end

function Gridap.Fields.return_cache(m::BulkGhostPenaltyAssembleRhsFEFunctionMap,aggregate_local_cells)
  cache_unassembled_rhs=array_cache(m.agg_cells_rhs_contribs)
  evaluate_result=Gridap.Arrays.CachedArray(eltype(eltype(m.agg_cells_rhs_contribs)),1)
  cache_unassembled_rhs,evaluate_result
end

function Gridap.Fields.evaluate!(cache,m::BulkGhostPenaltyAssembleRhsFEFunctionMap,aggregate_local_cells)
  cache_unassembled_rhs,result=cache
  contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,1)
  Gridap.Arrays.setsize!(result,(size(contrib,1),))
  result.array .= 0.0
  for (i,cell) in enumerate(aggregate_local_cells)
        contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,cell)
        result.array .+= contrib
  end
  result.array
end