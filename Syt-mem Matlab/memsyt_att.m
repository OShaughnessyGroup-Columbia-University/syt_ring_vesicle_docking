
function ems = memsyt_att(xv,yv,zv,av,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff)
% mem-syt attraction energy from node-node interaction

	ems = 0;
	LDebye0 = 0.7;
	for i = 1 : Nsyt
		dr1 = r4k(i,1) - xv;
		dr2 = r4k(i,2) - yv;
		dr3 = r4k(i,3) - zv;
		
		if abs(dr1)<dESmax && abs(dr2)<dESmax && abs(dr3)<dESmax
			d = sqrt(dr1*dr1+dr2*dr2+dr3*dr3);
			rcutoff = interpolate_func(r_cutoff,av,amin_cutoff,da_cutoff,n_cutoff);
			d = max(d, rcutoff);
			eij = av / d * exp(-d/LDebye);	% Eescoeff<0
			ems = ems + eij;
		end
	end
	ems = Eescoeff * ems;
end



% function ems = memsyt_att(xv,yv,zv,av,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_4k_mem_min)

	% ems = 0;
	% for i = 1 : Nsyt
		% dr1 = r4k(i,1) - xv;
		% dr2 = r4k(i,2) - yv;
		% dr3 = r4k(i,3) - zv;
		
		% if abs(dr1)<dESmax && abs(dr2)<dESmax && abs(dr3)<dESmax
			% d = sqrt(dr1*dr1+dr2*dr2+dr3*dr3);
			% %d = max(d, 0.1);
			% d = max(d, r_4k_mem_min);
			% eij = Eescoeff * av / d * exp(-d/LDebye);	% Eescoeff<0
			% ems = ems + eij;
		% end
	% end
% end


