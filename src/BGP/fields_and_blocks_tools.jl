function _restrict_to_block(cell_dof_ids::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}}, blockid)
    map=cell_dof_ids.maps.value 
    @assert length(map.size)==1
    @assert blockid >= 1
    @assert blockid <= map.size[1]
    cell_dof_ids.args[blockid]
end 

function _get_single_field_fe_basis(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.single_field
end
function _get_single_field_fe_basis(a)
    a
end
function _is_multifield_fe_basis_component(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    true
end
function _is_multifield_fe_basis_component(a)
    false
end
function _nfields(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.nfields
end
function _fieldid(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.fieldid
end
