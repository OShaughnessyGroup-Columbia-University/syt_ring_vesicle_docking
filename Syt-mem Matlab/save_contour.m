
%----------------- controls -----------------
n_bin = 51;
r_bin = zeros(n_bin,1);
z_bin = zeros(n_bin,1);
rc = [0 0 0];

%----------------- binning -----------------

rmax = sidehalf;
dr = rmax/(n_bin-1);

for i = 1 : Nsyt
	rc(1,:) = rc(1,:) + rsyt(i,:);
end

rc = rc/Nsyt;	% ring center

for i = 1 : n	% for each vertex
	ri = sqrt((x(i)-rc(1,1))^2 + (y(i)-rc(1,2))^2);
	j = floor(ri/dr+1.5);
	if j <n_bin
		r_bin(j) = r_bin(j) + 1;
		z_bin(j) = z_bin(j) + z(i);
	end
end

for i = 1 : n_bin
	if r_bin(i)>0
		z_bin(i) = z_bin(i)/r_bin(i);	% average z
	end
	r_bin(i) = (i-1)*dr;
end


r_bin_chop = r_bin(z_bin>0);	% remove empty elements
z_bin_chop = z_bin(z_bin>0);


%----------------- save membrane contour -----------------

rs = sqrt(rsyt(1,1)^2 + rsyt(1,2)^2);

fn = sprintf('contour_%d.dat', loop);
fid = fopen(fn, 'w');
fprintf(fid, 'r\tz\n');
for i = 1 : length(r_bin_chop)
	if i == 1
		fprintf(fid, '%f\t%f\t%f\t%f\n', r_bin_chop(i), z_bin_chop(i), rs, Zsyt);
	else
		fprintf(fid, '%f\t%f\n', r_bin_chop(i), z_bin_chop(i));
	end
end
fclose(fid);


