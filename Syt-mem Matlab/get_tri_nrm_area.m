
function [nrm, area] = get_tri_nrm_area(i,triver,x,y,z,tri_nrm_sign)
% returns normal direction of i-th triangle

	v1 = triver(i,1);
	v2 = triver(i,2);
	v3 = triver(i,3);
	
	dr12 = [x(v2)-x(v1), y(v2)-y(v1), z(v2)-z(v1)];	% vector v1->v2
	dr13 = [x(v3)-x(v1), y(v3)-y(v1), z(v3)-z(v1)];	% vector v1->v3
	
	%nrm = cross(dr12, dr13);	% normal is always defined as dr12 x dr13
	nrm = [dr12(2)*dr13(3)-dr12(3)*dr13(2), dr12(3)*dr13(1)-dr12(1)*dr13(3), dr12(1)*dr13(2)-dr12(2)*dr13(1)];
	lnrm = norm(nrm);
	nrm = tri_nrm_sign(i) * nrm / lnrm;
	area = 0.5*lnrm;
end

