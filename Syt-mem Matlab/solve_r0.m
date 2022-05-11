
function r0 = solve_r0(ESM,LDebye,a)
% calculate the min distance between mem vertex to KKKK
	fun = @(x) x*exp(x/LDebye)-a*ESM/2/pi/LDebye;
	r0 = fzero(fun,0.5);
end


