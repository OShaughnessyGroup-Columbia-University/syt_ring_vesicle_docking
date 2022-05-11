
%----------------- controls -----------------
plot_mem = 3;	% 1=2D plot, 2=3D scatter, 3=3D surf
mark_mem = 0;	%	mark nodes
show_nrm = 0;	% show normals

plot_syt = 2;	% 1=only show KKKK, 2=show both KKKK & C2B
mark_syt = 0;

%----------------- plot membrane -----------------

if plot_mem == 1	% 2D plot
	[ix, jx] = find(bond);
	for ii = 1:length(ix)
	    iplot = ix(ii);
	    jplot = jx(ii);
	    plot([x(iplot),x(jplot)],[y(iplot), y(jplot)],'Color',[0 1 0]);
	    hold on;
	end
	
elseif plot_mem == 2	% 3D line-by-line mesh
	intx = x;
	inty = y;
	intz = z;
	intx(isetn) = [];
	inty(isetn) = [];
	intz(isetn) = [];
	scatter3(intx, inty, intz, 20, 'blue', 'fill');
	hold on;
	
	% plot the external vertices 
	etnx = x(isetn);
	etny = y(isetn);
	etnz = z(isetn);
	scatter3(etnx, etny, etnz, 20, 'red', 'fill');
	hold on;
	
%	% plot the rim of nanodisc, if applicable
%	if exist('indfixed', 'var')
%	    etnx = x(indfixed);
%	    etny = y(indfixed);
%	    etnz = z(indfixed);
%	    scatter3(etnx, etny, etnz, 50, 'black', 'fill');
%	end
	
	% plot crosslinkers as green lines
	% pairs of indices of crosslinkers
	[ix, jx] = find(bond);
	for ii = 1:length(ix)
	    iplot = ix(ii);
	    jplot = jx(ii);
	    plot3([x(iplot),x(jplot)],[y(iplot), y(jplot)],[z(iplot),z(jplot)],'Color',[0 1 0]);
	end
	

elseif plot_mem == 3	% 3D triangle mesh
	m = trimesh(triver,x,y,z);
	set(m, 'LineWidth',1.5);
	hold on;
end

%----------------- mark membrane -----------------
if mark_mem
	for i = 1 : n	% for each vertex
		str = sprintf('%d', i);
		%str = sprintf('%.3g', mean_curv_vtx(i));
		%str = sprintf('%.2f', area_vertex(i));
		%str = sprintf('%.2f, %.2f', x(i), y(i));
		text(x(i), y(i), str, 'HorizontalAlignment', 'center', 'color', 'b', 'Fontsize', 10);
		hold on;
	end
	
	for i = 1 : ntri	% for each face
		str = sprintf('%d', i);
		text(tricnt(i,1), tricnt(i,2), str, 'HorizontalAlignment', 'center', 'color', 'g', 'Fontsize', 10);
		hold on;
	end
end

%----------------- show membrane normal -----------------

if show_nrm
	r1 = tricnt;
	r2 = nv;
	for i = 1 : ntri
		quiver3(r1(i,1),r1(i,2),r1(i,3),r2(i,1),r2(i,2),r2(i,3));
	end
	hold on;
end

%----------------- plot syt ring -----------------
rad_4k = 0.25*Rc2b;
color_sf = [0 1 0];	% syt face
color_se = [0 0.75 0];	% syt edge
color_kf = [0.3 0.3 1];	% 4K face
color_ke = [0.3 0.3 1];	% 4k edge

if plot_syt == 1
	[xk,yk,zk] = sphere;
	for i = 1:Nsyt
	    etnx = r4k(i,1);
	    etny = r4k(i,2);
	    etnz = r4k(i,3);
			sk = surf(rad_4k*xk+etnx,rad_4k*yk+etny,rad_4k*zk+etnz);
			set(sk,'FaceColor',color_kf,'FaceAlpha',1,'EdgeColor',color_ke,'CDataMapping','direct');
	end
elseif plot_syt == 2
	[xs,ys,zs] = sphere;
	[xsa,ysa,zsa] = sphere;
	[xsb,ysb,zsb] = sphere;
	[xk,yk,zk] = sphere;
	for i = 1:Nsyt
	    etnx = rsyt(i,1);
	    etny = rsyt(i,2);
	    etnz = rsyt(i,3);
	    
	    etnxa = rsyta(i,1);
	    etnya = rsyta(i,2);
	    etnza = rsyta(i,3);
	    
	    etnxb = rsytb(i,1);
	    etnyb = rsytb(i,2);
	    etnzb = rsytb(i,3);
	    
	    etnxk = r4k(i,1);
	    etnyk = r4k(i,2);
	    etnzk = r4k(i,3);
	    
	    s = surf(xs*Rc2b+etnx,ys*Rc2b+etny,zs*Rc2b+etnz);
	    sa = surf(xsa*Rc2b+etnxa,ysa*Rc2b+etnya,zsa*Rc2b+etnza);
	    sb = surf(xsb*Rc2b+etnxb,ysb*Rc2b+etnyb,zsb*Rc2b+etnzb);
			set(s, 'FaceColor',color_sf,'FaceAlpha',1,'EdgeColor',color_se,'CDataMapping','direct','LineWidth',1.5);
			set(sa,'FaceColor',color_sf,'FaceAlpha',1,'EdgeColor',color_se,'CDataMapping','direct','LineWidth',1.5);
			set(sb,'FaceColor',color_sf,'FaceAlpha',1,'EdgeColor',color_se,'CDataMapping','direct','LineWidth',1.5);
	    
	    sk = surf(rad_4k*xk+etnxk,rad_4k*yk+etnyk,rad_4k*zk+etnzk);
			set(sk,'FaceColor',color_kf,'FaceAlpha',1,'EdgeColor',color_ke,'CDataMapping','direct');

	end
end

hold on;


%----------------- end -----------------

axis equal
hold off