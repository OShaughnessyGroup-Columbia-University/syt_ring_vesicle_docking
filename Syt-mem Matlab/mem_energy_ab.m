function [tri_nrm_new,atri_new,vtx_nrm_new,area_vertex_new,mean_curv_vtx_new,Ebend_vtx_new,de_t,de_b] = ...
	mem_energy_ab(x,y,z,znew,n,ntri,tri_nrm_sign,tri_nrm,atri,vtx_nrm,area_vertex,triver,...
	tens,isetn,v_nbr_tri_cnt,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list,kappadb)
% calculate change in membrane tension and bending energy due to move (x,y,znew) for all membrane nodes

	% ====== udpate membrane area ======
	d_area = 0;
	for i = 1 : ntri	% going through all faces
		[tri_nrm_new(i,:), atri_new(i)] = get_tri_nrm_area(i,triver,x,y,znew,tri_nrm_sign);
		da_tri_i = atri_new(i) - atri(i);
		d_area = d_area + da_tri_i;
	end
	
	de_t = d_area * tens;	% change in tension energy
    	
	% ====== membrane bending ======
	eb_old = 0;
	eb_new = 0;

	for i = 1 : n
		if isetn(i)
			continue;
		end
	
		%vtx_nrm_new(i,:) = get_vtx_nrm(i,v_nbr_face,tri_nrm_new,atri_new);
		%area_vertex_new(i) = get_vertex_area(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri_new);
		vtx_nrm_new(i,:) = get_vtx_nrm(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,tri_nrm_new,atri_new);
		area_vertex_new(i) = get_vertex_area(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri_new);
		
		mean_curv_vtx(i) = get_curvature_vtx(i,x,y,z,vtx_nrm,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list);
		mean_curv_vtx_new(i) = get_curvature_vtx(i,x,y,znew,vtx_nrm_new,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list);
	  
		Ebend_vtx(i) = kappadb * mean_curv_vtx(i)^2 * area_vertex(i);
		Ebend_vtx_new(i) = kappadb * mean_curv_vtx_new(i)^2 * area_vertex_new(i);
	    		
		eb_old = eb_old + Ebend_vtx(i);
		eb_new = eb_new + Ebend_vtx_new(i);
	end
	
	de_b = eb_new - eb_old;
	
end