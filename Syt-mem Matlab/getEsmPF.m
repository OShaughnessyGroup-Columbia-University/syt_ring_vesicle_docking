
function e = getEsmPF(r,cs,a,LDebye)
% full approximation for syt-mem interaction from area-node

	radius = sqrt(a/pi);
	r0 = 1.0*radius;
	
	dr = 0.1*r0;
	r1 = r0 - dr/2;
	r2 = r0 + dr/2;
	
	if r<=r1
		e = getEsmNear(r,cs,radius,LDebye);
	elseif r>=r2
		e = getEsmFar(r,cs,radius,LDebye);
	else
		e1 = getEsmNear(r,cs,radius,LDebye)*(1-(r-r1)/dr);	% linear transition
		e2 = getEsmFar(r,cs,radius,LDebye)*(r-r1)/dr;
		if e1 < 0
			e1 = 0;
		end
		if e2 < 0
			e2 = 0;
		end
		e = e1 + e2;
	end
	
	% common prefactor for near & far is q*alpha=ESM0/LDebye0/2/pi
	% e>0 and decreases monotonically with r, 
	
end


