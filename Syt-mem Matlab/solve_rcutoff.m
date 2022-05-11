
function rcut = solve_rcutoff(a,d,LDebye,LDebye0,LDebyepidb)
% calculate the min distance between mem vertex to KKKK

	% 1 node
	%fun = @(x) exp(-x/LDebye)/x-LDebyepidb/a;

	% 1-ring neighbor
	%fun = @(x) exp(-x/LDebye)/x+6*exp(-sqrt(x^2+d^2)/LDebye)/sqrt(x^2+d^2)-LDebyepidb/a;
	
	% 3-ring neighbor
	%fun = @(x) exp(-x/LDebye)/x+6*exp(-sqrt(x^2+d^2)/LDebye)/sqrt(x^2+d^2)+6*exp(-sqrt(x^2+3*d^2)/LDebye)/sqrt(x^2+3*d^2)+6*exp(-sqrt(x^2+4*d^2)/LDebye)/sqrt(x^2+4*d^2)-LDebyepidb/a;
	
	%in-plane off-center 1-ring neighbor
	%fun = @(x) exp(-x/LDebye)/x + 2*exp(-sqrt(x^2+d^2-sqrt(3)*x*d)/LDebye)/sqrt(x^2+d^2-sqrt(3)*x*d) + 2*exp(-sqrt(x^2+d^2)/LDebye)/sqrt(x^2+d^2) + 2*exp(-sqrt(x^2+d^2+sqrt(3)*x*d)/LDebye)/sqrt(x^2+d^2+sqrt(3)*x*d) - LDebyepidb/a;
	
	%in-plane off-center 3-ring neighbor
	%fun = @(x) exp(-x/LDebye)/x + 2*exp(-sqrt(x^2+d^2-sqrt(3)*x*d)/LDebye)/sqrt(x^2+d^2-sqrt(3)*x*d) + 2*exp(-sqrt(x^2+d^2)/LDebye)/sqrt(x^2+d^2) + 2*exp(-sqrt(x^2+d^2+sqrt(3)*x*d)/LDebye)/sqrt(x^2+d^2+sqrt(3)*x*d) +6*exp(-sqrt(x^2+3*d^2)/LDebye)/sqrt(x^2+3*d^2)+6*exp(-sqrt(x^2+4*d^2)/LDebye)/sqrt(x^2+4*d^2)-LDebyepidb/a;
	
	fun = @(x) func(x,a,d,LDebye,LDebye0,LDebyepidb);
	
	rcut = fzero(fun,0.1);	% find the root of function fun
	
	if isnan(rcut) || rcut<1e-3
		rcut = fzero(fun,0.2);	% try different initial guesses
	end
	if isnan(rcut) || rcut<1e-3
		rcut = fzero(fun,0.5);
	end
	if isnan(rcut) || rcut<1e-3
		rcut = fzero(fun,1);
	end
	if isnan(rcut) || rcut<1e-3
		rcut = fzero(fun,2);
	end
	if isnan(rcut) || rcut<1e-3
		rcut = fzero(fun,5);
	end
	
end




function y = func(x,a,d,LDebye,LDebye0,LDebyepidb)
	d2 = d*d;
	x2 = x*x;
	
	% use discrete calc. to find potential of nodes inside the 2 rings
	y1 = exp(-x/LDebye)/x;	% center node
	
	y2 = exp(-sqrt(x2+d2-sqrt(3)*x*d)/LDebye)/sqrt(x2+d2-sqrt(3)*x*d);	% 1st ring neighbors in the same plane
	y3 = exp(-sqrt(x2+d2)/LDebye)/sqrt(x2+d2);
	y4 = exp(-sqrt(x2+d2+sqrt(3)*x*d)/LDebye)/sqrt(x2+d2+sqrt(3)*x*d);
	
	y5 = exp(-sqrt(x2+3*d2)/LDebye)/sqrt(x2+3*d2);	% 2nd ring neighbors
	y6 = exp(-sqrt(x2+4*d2)/LDebye)/sqrt(x2+4*d2);
	
	% use continuum approx. to estimate potential of nodes outside the 2 rings
	rin = sqrt(19*a/pi);	% radius of the 2 inner rings
	
	y = y1 + 2*y2 + 2*y3 + 2*y4 + 6*y5 + 6*y6 - LDebyepidb/a*(1-exp(-rin/LDebye));
end