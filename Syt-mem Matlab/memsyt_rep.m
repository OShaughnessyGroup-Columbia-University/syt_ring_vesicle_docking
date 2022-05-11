
function e = memsyt_rep(xi,yi,zi,Nsyt,rsyt,rsyta,rsytb,Rc2b,Rc2b2,Krepulsion)
% syt-mem steric repulsion for membrane imem
	
	e = 0;
	
	for i = 1 : Nsyt
		% -------- center sphere --------
		dx = abs(xi-rsyt(i,1));
		dy = abs(yi-rsyt(i,2));
		dz = abs(zi-rsyt(i,3));
		
		if dx<Rc2b && dy<Rc2b && dz<Rc2b
			dr2 = dx*dx+dy*dy+dz*dz;
			if dr2<Rc2b2	% collision in new position
				e = e + (Rc2b-sqrt(dr2))^2;
			end
		end
		
		% -------- side sphere-a --------
		dx = abs(xi-rsyta(i,1));
		dy = abs(yi-rsyta(i,2));
		dz = abs(zi-rsyta(i,3));
		
		if dx<Rc2b && dy<Rc2b && dz<Rc2b
			dr2 = dx*dx+dy*dy+dz*dz;
			if dr2<Rc2b2	% collision in new position
				e = e + (Rc2b-sqrt(dr2))^2;
			end
		end
		
		% -------- side sphere-b --------
		dx = abs(xi-rsytb(i,1));
		dy = abs(yi-rsytb(i,2));
		dz = abs(zi-rsytb(i,3));
		
		if dx<Rc2b && dy<Rc2b && dz<Rc2b
			dr2 = dx*dx+dy*dy+dz*dz;
			if dr2<Rc2b2	% collision in new position
				e = e + (Rc2b-sqrt(dr2))^2;
			end
		end
	end	% end of i-loop
	
	e = 0.5 * Krepulsion * e;
end

