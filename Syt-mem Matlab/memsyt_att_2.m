
function ems = memsyt_att_2(xv,yv,zv,av,nv,Nsyt,r4k,dESmax,Eescoeff,LDebye)
% mem-syt attraction energy from area-node interaction

	e = 0;
	nrm = [0, 0, 0];
	LDebye0 = 0.7;
	
	for i = 1 : Nsyt
		dr1 = r4k(i,1) - xv;
		dr2 = r4k(i,2) - yv;
		dr3 = r4k(i,3) - zv;
		
		if abs(dr1)<dESmax && abs(dr2)<dESmax && abs(dr3)<dESmax
			d = sqrt(dr1*dr1+dr2*dr2+dr3*dr3);
			if d < 1e-9;
				d = 1e-9;
			end
			nrm(1) = dr1/d;
			nrm(2) = dr2/d;
			nrm(3) = dr3/d;
			%cs = dot(nrm, nv);
			%cs = sum(nrm .* nv);	% faster than dot()
			cs = nrm(1)*nv(1) + nrm(2)*nv(2) + nrm(3)*nv(3);
			cs = abs(cs);
			eij = getEsmPF(d,cs,av,LDebye);	% eij>0
			e = e + eij;
		end
	end
	e = Eescoeff * e;	% Eescoeff<0
end


