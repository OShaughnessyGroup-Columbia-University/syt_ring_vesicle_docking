
function a = get_vertex_area(id,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri)

	a = 0;
	for i = 1 : v_nbr_tri_cnt(id)
		j = v_nbr_tri_seq_list(id,i);
		a = a + atri(j);
	end
	
	a = a/3.0;
end



