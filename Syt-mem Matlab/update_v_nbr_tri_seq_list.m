function [v_nbr_tri_seq_list, v_nbr_tri_cnt] = update_vnbr_tri_seq_list(id,triver,ntri,v_nbr_vtx_seq_list,v_nbr_cnt);
% update the list of faces around vertex-i in a ccw manner

	v_nbr_tri_seq_list = zeros(1,10);
	v_nbr_tri_cnt = 0;
	
	for i = 1 : v_nbr_cnt	% for each 3 vertecies: (i,j,j+1)
		j = i + 1;
		if j > v_nbr_cnt
			j = 1;
		end
		
		vtri = [id, v_nbr_vtx_seq_list(id,i), v_nbr_vtx_seq_list(id,j)];
		vtri = sort(vtri);
		
		k = 1;
		flg = 0;
		while flg==0 && k<=ntri
			%if vtri == triver(k,:)
			if vtri(1) == triver(k,1) && vtri(2) == triver(k,2) && vtri(3) == triver(k,3)
				flg = 1;
			end
			k = k + 1;
		end
		
		if flg==1	% has a matching triangle
			v_nbr_tri_cnt = v_nbr_tri_cnt + 1;
			v_nbr_tri_seq_list(1,v_nbr_tri_cnt) = k-1;
		end
	end
	
end