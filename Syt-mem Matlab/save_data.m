
% ====== syt-mem attraction ======
ees = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
if Extra_PIP2==1
	d = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% ave. distance
	ees = ees + sytpip2_att(d,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
end


% ====== membrane bending ======
ebend = 0;

for i = 1 : n
	if isetn(i)
		continue;
	end
	%vtx_nrm(i,:) = get_vtx_nrm(i,v_nbr_face,tri_nrm,atri);
	vtx_nrm(i,:) = get_vtx_nrm(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,tri_nrm,atri);
	mean_curv_vtx(i) = get_curvature_vtx(i,x,y,z,vtx_nrm,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list);
	Ebend_vtx(i) = kappadb * mean_curv_vtx(i)^2 * area_vertex(i);
	ebend = ebend + Ebend_vtx(i);
end

% ====== membrane tension ======
area = 0;
for i = 1:size(triver,1)
    dx12 = x(triver(i,1)) - x(triver(i,2));
    dy12 = y(triver(i,1)) - y(triver(i,2));
    dz12 = z(triver(i,1)) - z(triver(i,2));
    dr12 = [dx12; dy12; dz12];
    dx13 = x(triver(i,1)) - x(triver(i,3));
    dy13 = y(triver(i,1)) - y(triver(i,3));
    dz13 = z(triver(i,1)) - z(triver(i,3));
    dr13 = [dx13; dy13; dz13];
    %area = area + .5 * norm(cross(dr12,dr13));
    crs = [dr12(2)*dr13(3)-dr12(3)*dr13(2), dr12(3)*dr13(1)-dr12(1)*dr13(3), dr12(1)*dr13(2)-dr12(2)*dr13(1)];
    area = area + .5 * norm(crs);
    
end
etens = tens * (area - area_mem_0);

% ====== syt-mem repulsion ======
erep = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion);

% ====== total energy ======
energy = ees + ebend + etens + erep;

% ====== height ======
zm = max(z);
zrel = zm - Zsyt + Rc2b;

% display
fprintf('%d: T = %.3g, zsyt = %.3g, zmax = %.3g, zrel = %.3g, Ees = %.3g, Ebend = %.3g, Etens = %.3g, Etot = %.3g\n', ...
	loop, temperature, Zsyt, zm, zrel, ees, ebend, etens, energy);

fprintf('\tAccept = (%.3g, %.3g, %.3g, %.3g, %.3g, %.3g)\n', sucrate1, sucrate2, sucrate3, sucrate4, sucrate5, sucrate6);


% save file
fid = fopen('summary.dat', 'a');
%fprintf(fid, '%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n', loop, temperature, ees, ebend, etens, energy, zm, zrel, area);
fprintf(fid, '%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n', ...
	loop, temperature, zm, zrel, ees, ebend, etens, energy, ees/Nsyt, ebend/Nsyt, (ees+ebend)/Nsyt);
fclose(fid);


fid = fopen('acceptance.dat', 'a');
fprintf(fid, '%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n', loop, sucrate1, sucrate2, sucrate3, sucrate4, sucrate5, sucrate6);
fclose(fid);


