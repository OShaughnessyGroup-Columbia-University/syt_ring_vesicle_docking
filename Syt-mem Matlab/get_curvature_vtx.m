
function H = get_curvature_vtx(id,x,y,z,vtx_nrm,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list)
% calculate the mean curvature at vertex(id) using Watanabe method

	dr1 = zeros(3, 1);
	dr2 = zeros(3, 1);
	ki = zeros(v_nbr_vtx_cnt(id), 1);
	phi = zeros(v_nbr_vtx_cnt(id), 1);
	
	for i = 1 : v_nbr_vtx_cnt(id)
		i1 = i + 1;
		if i1 > v_nbr_vtx_cnt(id)
			i1 = 1;
		end
		
		j = v_nbr_vtx_seq_list(id, i);	% j is the i-th neighbor of id
		k = v_nbr_vtx_seq_list(id, i1);	% k i sthe (i+1)-th neighbor of id
		
		dr1(1) = x(j) - x(id);
		dr1(2) = y(j) - y(id);
		dr1(3) = z(j) - z(id);
		dr2(1) = x(k) - x(id);
		dr2(2) = y(k) - y(id);
		dr2(3) = z(k) - z(id);
		
		d2 = dr1(1)^2 + dr1(2)^2 + dr1(3)^2;
		if d2 == 0
			d2 = 1e-9;
		end
		
		%ki(i) = -2*dot(vtx_nrm(id,:), dr1)/d2;	% negative sign s.t. convex is positive
		%ki(i) = -2*sum(vtx_nrm(id,:) .* dr1')/d2;	% faster than dot
		ki(i) = -2*(vtx_nrm(id,1)*dr1(1)+vtx_nrm(id,2)*dr1(2)+vtx_nrm(id,3)*dr1(3))/d2;
		phi(i) = getangle(dr1, dr2);	% 0<=phi<=pi is between id_(i) and id_(i+1)
	end
	
	H = 0;
	
	for i = 1 : v_nbr_vtx_cnt(id)
		i1 = i - 1;
		if i1 < 1
			i1 = v_nbr_vtx_cnt(id);
		end
		
		H = H + ki(i) * (phi(i) + phi(i1));
	end
	
	H = H/4/pi;
end




function a = getangle(r1, r2)
	
	%nrm1 = r1/norm(r1);
	%nrm2 = r2/norm(r2);
	%dotp = dot(nrm1, nrm2);
	%dotp = sum(nrm1 .* nrm2);	% faster than dot()
	
	nrm1 = r1/sqrt(r1(1)^2+r1(2)^2+r1(3)^2);
	nrm2 = r2/sqrt(r2(1)^2+r2(2)^2+r2(3)^2);
	dotp = nrm1(1)*nrm2(1)+nrm1(2)*nrm2(2)+nrm1(3)*nrm2(3);
	
	if dotp > 1
		dotp = 1;
	elseif dotp < -1
		dotp = -1;
	end
	
	a = acos(dotp);
end

