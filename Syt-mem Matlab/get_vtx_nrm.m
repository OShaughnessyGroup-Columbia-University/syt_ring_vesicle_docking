
function nrm = get_vtx_nrm(id,v_nbr_tri_cnt,v_nbr_tri_seq_list,tri_nrm,atri)
% returns normal direction of i-th vertex

	nrm = [0, 0, 0];
	for i = 1 : v_nbr_tri_cnt(id)	% for each connected face
		j = v_nbr_tri_seq_list(id,i);
		nrm = nrm + tri_nrm(j,:)*atri(j);
	end
	
	nrm = nrm/norm(nrm);
end
