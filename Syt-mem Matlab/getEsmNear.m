
function e = getEsmNear(r,cs,radius,LDebye)
% near-field approximation for syt-mem interaction from area-node

	radius2 = radius*radius;	% radius = effective radius for membrane patch
	radius3 = radius2*radius;
	radius5 = radius3*radius2;
	radius7 = radius5*radius2;
	
	LDebye2=LDebye*LDebye;
	LDebye3=LDebye2*LDebye;
	LDebye4=LDebye3*LDebye;
	LDebye5=LDebye4*LDebye;
	LDebye6=LDebye5*LDebye;
	LDebye7=LDebye6*LDebye;
	LDebye8=LDebye7*LDebye;
	
	exp_radius_L = exp(-radius/LDebye);
		
	r2 = r*r;
	r3 = r2*r;
	r4 = r3*r;
	r5 = r4*r;
	r6 = r5*r;
	r7 = r6*r;
	r8 = r7*r;
	%r9 = r8*r;
	%r10 = r9*r;
	
	cs2 = cs*cs;
	cs3 = cs2*cs;
	cs4 = cs3*cs;
	cs5 = cs4*cs;
	cs6 = cs5*cs;
	cs7 = cs6*cs;
	cs8 = cs7*cs;
	%cs9 = cs8*cs;
	%cs10 = cs9*cs;
	%cs11 = cs10*cs;
	
	% coefficients of series expansion of potential: exp[-z/D]-exp[-sqrt(z^2+R^2)/D] = c0+c1*r+c2*r^2+...
	c0 = 1 - exp_radius_L;
	c1 = -1.0/LDebye;
	c2 = 0.5*(radius+exp_radius_L*LDebye)/radius/LDebye2;
	c3 = -1/6.0/LDebye3;
	c4 = (radius3 - 3*exp_radius_L*LDebye2*(radius+LDebye))/24.0/radius3/LDebye4;
	c5 = -1/120.0/LDebye5;
	c6 = (radius5+15*exp_radius_L*LDebye3*(radius2+3*radius*LDebye+3*LDebye2))/720.0/radius5/LDebye6;
	c7 = -1/5040.0/LDebye7;
	c8 = (radius7-105*exp_radius_L*LDebye4*(radius3+6*radius2*LDebye+15*radius*LDebye2+15*LDebye3))/40320.0/radius7/LDebye8;
	%c9 = -1/362880.0/LDebye9;
	
	% Legendre polynomials
	p0 = 1;
	p1 = cs;
	p2 = (3*cs2-1)/2.0;
	p3 = (5*cs3-3*cs)/2.0;
	p4 = (35*cs4-30*cs2+3)/8.0;
	p5 = (63*cs5-70*cs3+15*cs)/8.0;
	p6 = (231*cs6-315*cs4+105*cs2-5)/16.0;
	p7 = (429*cs7-693*cs5+315*cs3-35*cs)/16.0;
	p8 = (6435*cs8-12012*cs6+6930*cs4-1260*cs2+35)/128.0;
	%p9 = (12155*cs9-25740*cs7+18018*cs5-4620*cs3+315*cs)/128.0;
	%p10 = (46189*cs10-109395*cs8+90090*cs6-30030*cs4+3465*cs2-63)/256.0;
	%p11 = (88179*cs11-230945*cs9+218790*cs7-90090*cs5+15015*cs3-693*cs)/256.0;
	
	%pot = c0*p0 + c1*p1*r + c2*p2*r2 + c3*p3*r3 + ...
	%pot *= Pidb*DebyeMemCoeff*LDebye;
	
	e = c0*p0 + c1*p1*r + c2*p2*r2 + c3*p3*r3 + c4*p4*r4 + c5*p5*r5 + c6*p6*r6 + c7*p7*r7 + c8*p8*r8;
	e = e*2*pi*LDebye;	% without the DebyeMemCoeff (alpha) or DebyeVesCoeff factor
	
end


