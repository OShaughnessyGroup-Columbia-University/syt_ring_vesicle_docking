function [dn,nrm,mid] = getSMdistPlus(id,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,triver,ntri)
% smallest distance between KKKK of syt-id to membrane, also direction and vertices ID
	
	mid = zeros(3,1);
	nrm = zeros(3,1);
	
	dr = zeros(3,1);
	d2min = 1e6;
	dn = 1e6;
	
	% for each membrane vertex
	for i = 1 : n
		dr(1) = r4k(id,1) - x(i);
		dr(2) = r4k(id,2) - y(i);
		dr(3) = r4k(id,3) - z(i);
		d2 = dr(1)^2 + dr(2)^2 + dr(3)^2;
		
		if d2<d2min
			dotp = dr(1)*vtx_nrm(i,1) + dr(2)*vtx_nrm(i,2) + dr(3)*vtx_nrm(i,3);
			if dotp>0
				d2min = d2;
				dn = dotp;
				for j = 1 : 3
					nrm(j) = dr(j)/sqrt(d2);
					mid(j) = i;
				end
			end
		end
	end
	
	for i = 1 : ntri
		for j = 1 : 3
			dr(j) = r4k(id,j) - tricnt(i,j);
		end
		d2 = dr(1)^2 + dr(2)^2 + dr(3)^2;
		
		if d2<d2min
			dotp = dr(1)*tri_nrm(i,1) + dr(2)*tri_nrm(i,2) + dr(3)*tri_nrm(i,3);
			if dotp>0
				d2min = d2;
				dn = dotp;
				for j = 1 : 3
					nrm(j) = dr(j)/sqrt(d2);
					mid(j) = triver(i,j);
				end
			end
		end
	end
	
end