
function [area_i, v_nbr_face_i] = update_nbr_face_of_mem_vertex(imem,ntri,bond,triver,atri)

	v_nbr_face_i = zeros(1,ntri);	% reset imem-th row
	
	m1 = ismember(triver,imem);	% find imem in triver list, returns a nx3 matrix
	m2 = sum(m1,2);	% reduce m1 to a nx1 column
	tris = find(m2);	% triangle indecies of non-zero elements
	
	atmp = 0;
	for i = 1 : length(tris)
		j = tris(i);	% id of triangle
		v_nbr_face_i(j) = 1;
		atmp = atmp + atri(j);
	end
	
	area_i = atmp/3.0;
end



