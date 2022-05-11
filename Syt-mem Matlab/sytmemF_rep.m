
function e = sytmemF_rep(Nsyt,rsyt,rsyta,rsytb,tricnt,ntri,Rc2b,Rc2b2,Krepulsion)
% syt-mem steric repulsion for all syts
	
	e = 0;
	
	for i = 1 : Nsyt	% for all syts
		for j = 1 : ntri	% for all membrane triangles
			% -------- center sphere --------
			dx = abs(rsyt(i,1)-tricnt(j,1));
			dy = abs(rsyt(i,2)-tricnt(j,2));
			dz = abs(rsyt(i,3)-tricnt(j,3));
			
			if dx<Rc2b && dy<Rc2b && dz<Rc2b
				dr2 = dx*dx+dy*dy+dz*dz;
				if dr2<Rc2b2	% collision in new position
					e = e + (Rc2b-sqrt(dr2))^2;
				end
			end
			
			% -------- side sphere-a --------
			dx = abs(rsyta(i,1)-tricnt(j,1));
			dy = abs(rsyta(i,2)-tricnt(j,2));
			dz = abs(rsyta(i,3)-tricnt(j,3));
			
			if dx<Rc2b && dy<Rc2b && dz<Rc2b
				dr2 = dx*dx+dy*dy+dz*dz;
				if dr2<Rc2b2	% collision in new position
					e = e + (Rc2b-sqrt(dr2))^2;
				end
			end
			
			% -------- side sphere-b --------
			dx = abs(rsytb(i,1)-tricnt(j,1));
			dy = abs(rsytb(i,2)-tricnt(j,2));
			dz = abs(rsytb(i,3)-tricnt(j,3));
			
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

