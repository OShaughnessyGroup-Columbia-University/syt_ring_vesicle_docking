
function e = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff)
% syt-mem attraction energy from node-node interaction

	e = 0;
	
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
				rcutoff = interpolate_func(r_cutoff,area_vertex(i),amin_cutoff,da_cutoff,n_cutoff);
				d = max(d, rcutoff);
				eij = area_vertex(i) / d * exp(-d/LDebye);	% Eescoeff<0
				e = e + eij;
			end
			
		end	% end of j-loop
	end	% end of i-loop
	
	e = Eescoeff * e;
end




% function e = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,r_4k_mem_min)

	% e = 0;
	% for i = 1 : Nsyt
		% for j = 1 : n
			% dr1 = r4k(i,1) - x(j);
			% dr2 = r4k(i,2) - y(j);
			% dr3 = r4k(i,3) - z(j);
			
			% if abs(dr1)<dESmax && abs(dr2)<dESmax && abs(dr3)<dESmax
				% d = sqrt(dr1*dr1+dr2*dr2+dr3*dr3);
				% %d = max(d, 0.1);
				% d = max(d, r_4k_mem_min);
				% eij = Eescoeff * area_vertex(j) / d * exp(-d/LDebye);	% Eescoeff<0
				% e = e + eij;
			% end
			
		% end
	% end
% end

