
function e = sytmem_att_2(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn)
% syt-mem attraction energy from area-node interaction

	e = 0;
	nrm = [0, 0, 0];
	
	for i = 1 : n
		if isetn(i)
			continue;
		end
		for j = 1 : Nsyt
			dr1 = x(i) - r4k(j,1);
			dr2 = y(i) - r4k(j,2);
			dr3 = z(i) - r4k(j,3);
			
			if abs(dr1)<dESmax && abs(dr2)<dESmax && abs(dr3)<dESmax
				d = sqrt(dr1*dr1+dr2*dr2+dr3*dr3);
				if d < 1e-9;
					d = 1e-9;
				end
				nrm(1) = dr1/d;
				nrm(2) = dr2/d;
				nrm(3) = dr3/d;
				%cs = dot(nrm, vtx_nrm(i,:));
				%cs = sum(nrm .* vtx_nrm(i,:));	% faster than dot()
				cs = nrm(1)*vtx_nrm(i,1) + nrm(2)*vtx_nrm(i,2) + nrm(3)*vtx_nrm(i,3);
				cs = abs(cs);
				eij = getEsmPF(d,cs,area_vertex(i),LDebye);	% eij>0
				e = e + eij;
			end
			
		end
	end
	
	e = Eescoeff * e;	% Eescoeff<0
end


