function y = interpolate_func(func,x,xmin,dx,n)
	i = floor((x-xmin)/dx);
	if i<=1
		y = func(1);
	elseif i >= n
		y = func(n);
	else
		x1 = xmin+i*dx;
		frac = (x-x1)/dx;
		y = func(i) + (func(i+1)-func(i))*frac;
	end
end

