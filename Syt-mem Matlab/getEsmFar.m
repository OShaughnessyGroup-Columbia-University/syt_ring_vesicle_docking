
function e = getEsmFar(r,cs,radius,LDebye)
% far-field approximation for syt-mem interaction from area-node
	
	radius2 = radius*radius;	% radius = effective radius for membrane patch
	radius4 = radius2*radius2;
	radius6 = radius4*radius2;
	
	LDebye2=LDebye*LDebye;
	LDebye3=LDebye2*LDebye;
	LDebye4=LDebye3*LDebye;
	
	exp_r_L = exp(-r/LDebye);
		
	r2 = r*r;
	r3 = r2*r;
	r4 = r3*r;
	r5 = r4*r;
	r6 = r5*r;
	
	cs2 = cs*cs;
	cs4 = cs2*cs2;
	sn2 = 1-cs2;
	sn4 = sn2*sn2;
	
	Z = pi*radius2;	% without the DebyeMemCoeff (alpha) term
	Q = -pi/4*radius4;
	F = pi/8*radius6;
	
	c1 = Z+Q/2/LDebye2*(cs2-1)+F/24/LDebye4*sn4;
	c2 = Q/2/LDebye*(3*cs2-1)-F/3/LDebye3*(sn2-5/32.0*sn4);
	c3 = Q/2*(3*cs2-1)+F/3/LDebye2*(1-6*sn2+45/8.0*sn4);
	c4 = F/8/LDebye*(35*cs4-30*cs2+3);
	c5 = F/8*(35*cs4-30*cs2+3);
	
	e = c1/r + c2/r2 + c3/r3 + c4/r4 + c5/r5;
	e = e*exp_r_L;	% without the DebyeMemCoeff or DebyeVesCoeff factor
	
end


