
function nrm = get_vtx_nrm_slow(id,v_nbr_face,tri_nrm,atri)
% returns normal direction of i-th vertex

	nrm = [0, 0, 0];
	for i = find(v_nbr_face(id,:))	% for each connected face
		nrm = nrm + tri_nrm(i,:)*atri(i);	% area-weighted normal
	end
	nrm = nrm/norm(nrm);
end