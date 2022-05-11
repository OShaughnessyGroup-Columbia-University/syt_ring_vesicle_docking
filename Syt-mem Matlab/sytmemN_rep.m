
function e = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion)
% syt-mem steric repulsion for all syts
	
	e = 0;
	
	for i = 1 : Nsyt	% for all syts
		for j = 1 : n	% for all membrane nodes
			% -------- center sphere --------
			dx = abs(rsyt(i,1)-x(j));
			dy = abs(rsyt(i,2)-y(j));
			dz = abs(rsyt(i,3)-z(j));
			
			if dx<Rc2b && dy<Rc2b && dz<Rc2b
				dr2 = dx*dx+dy*dy+dz*dz;
				if dr2<Rc2b2	% collision in new position
					e = e + (Rc2b-sqrt(dr2))^2;
				end
			end
			
			% -------- side sphere-a --------
			dx = abs(rsyta(i,1)-x(j));
			dy = abs(rsyta(i,2)-y(j));
			dz = abs(rsyta(i,3)-z(j));
			
			if dx<Rc2b && dy<Rc2b && dz<Rc2b
				dr2 = dx*dx+dy*dy+dz*dz;
				if dr2<Rc2b2	% collision in new position
					e = e + (Rc2b-sqrt(dr2))^2;
				end
			end
			
			% -------- side sphere-b --------
			dx = abs(rsytb(i,1)-x(j));
			dy = abs(rsytb(i,2)-y(j));
			dz = abs(rsytb(i,3)-z(j));
			
			if dx<Rc2b && dy<Rc2b && dz<Rc2b
				dr2 = dx*dx+dy*dy+dz*dz;
				if dr2<Rc2b2	% collision in new position
					e = e + (Rc2b-sqrt(dr2))^2;
				end
			end
		end	% end of j-loop
	end	% end of i-loop
	
	e = 0.5 * Krepulsion * e;
end

