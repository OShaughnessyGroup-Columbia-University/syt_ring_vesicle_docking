%clear
fprintf('Initializing ...\n');

% random number generator
rng('shuffle');	%% random number generator
S = rng;

%------------------ switches -----------------
Test_Mem = 0;	% test a cylindrical membrane
SytMemAtt = 1;	% syt-mem attraction mode: 1=point-point, 0=point-face

Mem_Flat = 0;	% initial membrane shape if Test_Mem=0: 1=flat, 0=step
Hmax_ini = 3.0;	% initial max height in nm, if Mem_Flat=0

ShiftSyt = 1;	% shift syt with membrane when moving membrane
Extra_PIP2 = 1;	% have extra PIP2 (1-on-1 syt-pip2) or not

%------------------ parameters -----------------
side = 30;	% size of the membrane in nm
dmin = 1;	% hard ball diameter of vertices

kT = 4.1;	% kT in pN*nm

kappa = 25;	% bending modulus in kT
tens = 0.01/kT;	% ring tension in terms of kT / nm^2

Nsyt = 10;	% number of syts in a ring
Asm0 = 0;	% location of KKKK on syt ring in degrees (0 is at 3 O'clock, 90 is at 6 O'clock)

ESM0 = 8.2;	% Syt-mem binding energy in kT for PS=25%, PIP2=0, [salt]=140 uM
ESP0 = 11.0;	% 1-to-1 syt-PIP2 binding energy in kT (phys. salt, see Bogaart...Jahn, JBC 2012), used if Extra_PIP2=1

PS = 25;	% PS%
PIP2 = 2;	% PIP2%, assuming 1PIP2=4PS

Salt = 140.0;	% current salt concentration in uM
LDebye0 = 0.7;	% Debye length in nm when [Salt]=140 uM

Ess = 10.0;	% syt-syt binding energy in kT
Rsyt0 = 14;	% spontaneous radius of syt oligomer (through center) in nm
LP = 220;	% persistence length of syt oligomer in nm

Rc2b = 1.5;	% radius of C2B in nm
Lc2b = 5;	% length of C2B in nm

Krepulsion = 1e4;	% a large number for collision penalty

%------------------ initial membrane shape -----------------
Rdome = side/2;	% dome radius

% for Test_Mem=1:
MemAng = pi;	% MemAng=pi for half a cylinder
Rcylinder = side/2/sin(MemAng/2);
Rcylinder2 = Rcylinder*Rcylinder;
Hcylinder = Rcylinder*cos(MemAng/2);

%------------------ general variables -----------------
r2d = 180/pi;
d2r = pi/180;
Asm = Asm0*d2r;

dmax = sqrt(3) * dmin;	% maximal tether length
dmin2 = dmin * dmin;
dmax2 = dmax * dmax;
maxd = 0.06 * dmin;	% max step size (0.06)	% test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lte = 1.05*dmin;	% tentetive length of tether, nm

atrimin = sqrt(3)/4 * dmin^2;
atrimax = sqrt(3)/4 * dmax^2;
atriave = 0.5*(atrimin+atrimax);
avtxmin = 6*atrimin/3.0;
avtxmax = 6*atrimax/3.0;

sidehalf = 0.5*side;
diaghalf = sqrt(2)*sidehalf;

Rc2b2 = Rc2b*Rc2b;
Rring = Nsyt*Lc2b/2/pi;
Rring2 = Rring*Rring;
Rring4k = Rring - Rc2b;
Rring4k2 = Rring4k*Rring4k;

LDebye = LDebye0*sqrt(140.0/Salt);	% actual Debye length
LDebyepidb = LDebye*2*pi;
LDebye2Pidb = LDebye*LDebye*2*pi;

ESM0eff = ESM0*(PS+4*PIP2)/25.0;	% contribution from both PS and PIP2
ESM = ESM0eff*LDebye/LDebye0;	% effective ESM

dESmax = 4.0*LDebye + dmax;
dESmax2 = dESmax*dESmax;

%Eescoeff = -ESM0eff/LDebye0/2/pi;	% prefactor of E_pp = ESM0eff*dA/(2pi)/LDebye0*exp(-r/LDebye)/r, netagive!
Eescoeff = -ESM0eff/2/pi/LDebye0;	% prefactor for syt-mem attraction energy, Ei(r)=Eescoeff*Area/r*exp(-r/LDebye), negative!

kappadb = 2.0*kappa;

sucrate = 0;
sucrate1 = 0;
sucrate2 = 0;
sucrate3 = 0;
sucrate4 = 0;
sucrate5 = 0;
sucrate6 = 0;


%------------------ membrane -----------------

% estimate an upper bound of # vertices
fprintf(' * Membrane nodes...\n');

% membrane area
maxarea = side * side;
% lower bound of area of each triangle
atr = 0.7 * lte * lte * sqrt(3) / 4;
% upper bound of # vertices
maxnver = round(maxarea / atr);

% preallocate
x = zeros(maxnver,1);
y = x;
z = x;
bond = sparse(maxnver,maxnver);
%bond = zeros(maxnver,maxnver);
isetn = false(maxnver,1);

% separate membrane into lines in the x direction
ysep = sqrt(3) / 2 * lte;	% separation of lines in the y direction
nyline = round(side / ysep);	% number of lines in the y direction
ysep = side / (nyline - 1);	% recalibrate separation of lines in the y direction
nxmax = round(side / lte);	% number of points in the x direction
xsep = side / (nxmax - 1);	% recalibrate separation in the x direction


% lay down points on each line
ind = 1;	% index of vertex
thisy = 0;	% go through each line

% for each y-line
for iline = 1:nyline
	% determine number of vertices on this line
	if mod(iline, 2) == 0
		nx = nxmax;
	else
		nx = nxmax - 1;
	end
	
	% lay down vertices on this line
	% x of first vertex
	if mod(iline,2) == 0	%; extra semicolon. -Jie
		thisx = 0;
	else
		thisx = xsep / 2;
	end
        
	% for each x-point
	for i = 1:nx
		if i == 1 || i == nx || iline == 1 || iline == nyline		% external points
			isetn(ind) = true;                
		end
		
		x(ind) = thisx;
		y(ind) = thisy;
		
		%====================================
		
		if Test_Mem == 1
			dx = thisx - side/2;
			z2 = Rcylinder2-dx*dx;
			if z2 < 0
				z2 = 0;
			end
			z(ind) = sqrt(z2)-Hcylinder;
		else
			if Mem_Flat==1
				z(ind) = 0; % initially a flat membrane
				%z(ind) = 1.0*rand;
			else
				dx = thisx - side/2;
				dy = thisy - side/2;
				dr2 = dx*dx+dy*dy;
				dr = sqrt(dr2);
%				if dr2 > Rring4k2
%					z(ind) = 0;
%				else
%					z(ind) = 0.9*Rc2b;	% step
%					%z(ind) = sqrt(Rring4k2-dr2);	% sphere
%				end
				
				if dr>Rdome
					z(ind) = 0;
				else
					zcos = 0.5*(1+cos(pi*sqrt(dr2)/Rdome));
					z(ind) = Hmax_ini*zcos^1;	% test: starting with a smooth shape !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				end
			end
		end
		%====================================
		
		% bond with vertices on the same line
		if i ~= 1
			bond(ind,ind-1) = 1;	% connecting ind to its left neighbor
			bond(ind-1,ind) = 1;
		end
		
		% bond with vertices on the previous line
		if iline ~= 1
			if mod(iline, 2) == 1	% if this vertex is on odd lines
				bond(ind,ind-nxmax) = 1;	% connecting ind to its lowerleft neighbor
				bond(ind-nxmax,ind) = 1;
				bond(ind,ind-nxmax+1) = 1;	% connecting ind to its lowerright neighbor
				bond(ind-nxmax+1,ind) = 1;
			else % if this vertex is on even lines
				if i ~= 1
					bond(ind,ind-nxmax) = 1;	% lowerleft neighbor
					bond(ind-nxmax,ind) = 1;                            
				end
				if i ~= nxmax
					bond(ind,ind-nxmax+1) = 1;	% lowerright neighbor
					bond(ind-nxmax+1,ind) = 1;
				end
			end
		end
		
   	ind = ind + 1;
		thisx = thisx + xsep;
		
	end	% end of i-loop
	
	thisy = thisy + ysep;
end	% end of iline

% clear zeros at the end of th, r and z, and give bond the correct size
x(ind:end) = [];
y(ind:end) = [];
z(ind:end) = [];

%% recenter x and y at the origin
x = x - mean(x);
y = y - mean(y);


bond = bond(1:ind-1, 1:ind-1);
isetn(ind:end) = [];

n = ind-1;	% number of vertices
netn = sum(isetn);	% number of external vertices
max_min_lbond;	% check maximum and minimum tether length

if ~Test_Mem
	if lmax > dmax || lmin < dmin
		error ('max or min length wrong')    
	end
end


%% identify triangles and calculater their unit normal vectors
fprintf(' * Membrane faces...\n');

% preallocate
triver = zeros(maxnver, 3);	% id of 3 vertices of triangles
tricnt = zeros(maxnver, 3);	% center position of triangles
nv = zeros(maxnver, 3);	% normal of vertices
atri = zeros(maxnver, 1);	% area of triangles

area_mem_0 = 0;

ind = 1;
% go through each vertex
for i = 1 : n
    % go through each neighbor of the first vertex
    nbr = find(bond(i,:));
    nbr(nbr < i) = [];	% to avoid double counting, nbr always has a higher index than i
    for j = nbr
        % find common neighbors of these two vertices
        nbr2 = find(and(bond(i,:), bond(j,:)));
        % there should be at most 2 common neighbors
        if length(nbr2) > 2
            error ('3 or more common neighbors')
        end
        % avoid double counting, nbr2 always has a higher index than j
        nbr2(nbr2 < j) = [];
        % go through each nbr2
        for k = nbr2
            % store vertex indices
            triver(ind,:) = [i, j, k];	% i<j<k
            % find center
            xc = (x(i) + x(j) + x(k))/3.0;
            yc = (y(i) + y(j) + y(k))/3.0;
            zc = (z(i) + z(j) + z(k))/3.0;
            tricnt(ind,:) = [xc, yc, zc];
            % find normal vector
            r1 = [x(i) - x(j), y(i) - y(j), z(i) - z(j)];
            r2 = [x(i) - x(k), y(i) - y(k), z(i) - z(k)];
            thisnv = cross(r1,r2);
            %thisnv = [r1(2)*r2(3)-r1(3)*r2(2), r1(3)*r2(1)-r1(1)*r2(3), r1(1)*r2(2)-r1(2)*r2(1)];
            % make its length 1
            thisatri = .5 * norm(thisnv);
            thisnv = thisnv / norm(thisnv);
            % store unit normal vector
            atri(ind) = thisatri;
            nv(ind,:) = thisnv;
            
            % calculate area of flat membrane
            r1(3)=0;
            r2(3)=0;
            area_mem_0 = area_mem_0 + 0.5 * norm(cross(r1,r2));
            
            % increase ind by 1
            ind = ind + 1;
        end
    end
end

% cut unnecessary lines
triver(ind:end,:) = [];
tricnt(ind:end,:) = [];
nv(ind:end,:) = [];
atri(ind:end,:) = [];
ntri = ind - 1;	% number of triangles

atriave = mean(atri);
dave = sqrt(4*atriave/sqrt(3));

%% calculate the neighboring relation between triangles
trinbr = zeros(ntri);
for i = 1:ntri-1
    for j = i+1:ntri
        if triver(i,1) == triver(j,1)
            if triver(i,2) == triver(j,2) || triver(i,3) == triver(j,2) || triver(i,3) == triver(j,3)
                trinbr(i,j) = 1;
                trinbr(j,i) = 1;              
            end
        elseif triver(i,2) == triver(j,1)
            if triver(i,3) == triver(j,2) || triver(i,3) == triver(j,3)
                trinbr(i,j) = 1;
                trinbr(j,i) = 1; 
            end
        elseif triver(i,2) == triver(j,2) && triver(i,3) == triver(j,3)
            trinbr(i,j) = 1;
            trinbr(j,i) = 1; 
        end
    end
end

trinbr = sparse(trinbr);


% get neighboring faces of vertex i
v_nbr_face = zeros(n, ntri);	% indicating which faces are around vertex v (row=vertex#, column=face#)
area_vertex = zeros(n, 1);
for i = 1 : n	% for each vertex
	[area_vertex(i), v_nbr_face(i,:)] = update_nbr_face_of_mem_vertex(i,ntri,bond,triver,atri);
end
v_nbr_face = sparse(v_nbr_face);	% un-sorted list of neighboring faces
area_vertex_max = max(area_vertex);
es_vertex = zeros(n, 1);


%====== Watanabe method ======
fprintf(' * Membrane normal...\n');

% get normal of each face
tri_nrm = zeros(ntri,3);	% normal of each triangle
tri_nrm_sign = ones(ntri,1);	% stores the signs to adjust normal direction

for i = 1 : ntri	% for each triangle
	tri_nrm(i,3) = 1;	% default value
	v1 = triver(i,1);
	v2 = triver(i,2);
	v3 = triver(i,3);
	
	dr12 = [x(v2)-x(v1), y(v2)-y(v1), z(v2)-z(v1)];	% vector v1->v2
	dr13 = [x(v3)-x(v1), y(v3)-y(v1), z(v3)-z(v1)];	% vector v1->v3
	
	nrm = cross(dr12,dr13);	% normal is always defined as dr12 x dr13
	nrm = nrm/norm(nrm);
	if nrm(3)<0
		tri_nrm_sign(i) = -1;
	end
	
	tri_nrm(i,:) = tri_nrm_sign(i) * nrm;	% normal of triangle i
end


% get normal of each vertex
vtx_nrm = zeros(n, 3);	% normal of each vertex
for i = 1 : n	% for each vertex
	if isetn(i)	% excluding boundary nodes
		vtx_nrm(i,3) = 1;	% up
	else
		vtx_nrm(i,:) = get_vtx_nrm_slow(i,v_nbr_face,tri_nrm,atri);
		%vtx_nrm(i,:) = get_vtx_nrm(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,tri_nrm,atri);
	end
end


fprintf(' * Membrane list...\n');
v_nbr_vtx_cnt = zeros(n, 1);	% number of vertices around each vertex
v_nbr_tri_cnt = zeros(n, 1);	% number of triangles around each vertex
v_nbr_vtx_seq_list = zeros(n, 10);	% max number of vertices around a vertex is 10, record vertex # sequentially in ccw direction
v_nbr_tri_seq_list = zeros(n, 10);	% max number of faces around a vertex is 10, record face # sequentially in ccw direction

for i = 1 : n
	v_nbr_vtx_cnt(i) = length(find(bond(i,:)));	% typical value is 6
	v_nbr_tri_cnt(i) = length(find(v_nbr_face(i,:)));	% typical value is 6
end

for i = 1 : n
	vnrm = vtx_nrm(i,:);
	v_nbr_vtx_seq_list(i,:) = update_v_nbr_vtx_seq_list(i,x,y,z,vnrm,bond,v_nbr_vtx_cnt(i));
	[v_nbr_tri_seq_list(i,:), v_nbr_tri_cnt(i)] = update_v_nbr_tri_seq_list(i,triver,ntri,v_nbr_vtx_seq_list,v_nbr_vtx_cnt(i));
end


% mean curvature & bending energy at each vertex
fprintf(' * Membrane curvature...\n');
mean_curv_vtx = zeros(n, 1);
Ebend_vtx = zeros(n, 1);
for i = 1 : n
	if isetn(i)
		mean_curv_vtx(i) = 0;
		Ebend_vtx(i) = 0;
	else
		mean_curv_vtx(i) = get_curvature_vtx(i,x,y,z,vtx_nrm,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list);
		Ebend_vtx(i) = kappadb * mean_curv_vtx(i)^2 * area_vertex(i);
	end
end

%====== end of Watanabe method ======



%% calculate mean r squared in x and y direction for all external vertices
r2bar = mean(x(isetn) .* x(isetn) + y(isetn) .* y(isetn));


%------------------ Syt -----------------

fprintf(' * Syt...\n');
Lsytring = Nsyt*Lc2b;	% perimeter of syt ring through center in nm
Rsyt = Lsytring/(2*pi);	% radius of syt ring through center in nm
R4k = Rsyt - Rc2b;	% radius of polylysine patch in ring in nm

Zsyt = Rc2b;	% initial z value of syt ring (no tilting motion afterwards)
%Zsyt = Rc2b + 3;	% test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if Mem_Flat == 0
	zcos = 0.5*(1+cos(pi*Rring/Rdome));
	Zsyt = Zsyt + Hmax_ini*zcos^1;	% syt on top of membrane !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

th0 = 2*pi*rand;
dth = 2*pi/Nsyt;
e1_syt = zeros(Nsyt,3);	% local e1 vector of each syt (radial direction)
e2_syt = zeros(Nsyt,3); % local e2 vector of each syt (tangential direction)
e3_syt = [0 0 1];

rsyt = zeros(Nsyt,3);	% position of C2B center
rsyta = zeros(Nsyt,3);	% end caps
rsytb = zeros(Nsyt,3);
r4k0 = zeros(3,1);
r4k = zeros(Nsyt,3);	% position of polylysine patch
dz4k = zeros(Nsyt,1);

dr = Rc2b*(1-cos(Asm));
dh =-Rc2b*sin(Asm);
	
for i = 1 : Nsyt
	th = th0 + (i-1)*dth;
	e1_syt(i,:) = [ cos(th), sin(th), 0];	% radial direction, e1 x e2 = e3 = nz
	e2_syt(i,:) = [-sin(th), cos(th), 0];	% tangential direction ccw
	
	rsyt(i,1) = Rsyt*e1_syt(i,1);	% position of each C2B center
	rsyt(i,2) = Rsyt*e1_syt(i,2);
	rsyt(i,3) = Zsyt;
	
	rsyta(i,:) = rsyt(i,:) - (Lc2b/2-Rc2b)*e2_syt(i,:);
	rsytb(i,:) = rsyt(i,:) + (Lc2b/2-Rc2b)*e2_syt(i,:);
	
	r4k0(1) = R4k*e1_syt(i,1);	% position of 3 O'clock on each syt
	r4k0(2) = R4k*e1_syt(i,2);
	r4k0(3) = Zsyt;
	
	r4k(i,1) = r4k0(1) + dr*e1_syt(i,1);
	r4k(i,2) = r4k0(2) + dr*e1_syt(i,2);
	r4k(i,3) = r4k0(3) + dh*e3_syt(3);
	
	dz4k(i) = r4k(i,3) - rsyt(i,3);
end


%------------------ extra syt-pip2 interaction -----------------

Q_KKKK=4;	% charge of KKKK
Q_PIP2=4;	% charge of PIP2
e_charge=1.6e-19;	% e in Coulomb

Inv4Pie0e = 1.0/(4*pi*8.85e-12*80);
ESP=ESP0*LDebye/LDebye0;	% syt-pip2 binding in current salt (in kT)

QsytQpip2_4pie0eSI = Q_KKKK*Q_PIP2*(e_charge^2)*Inv4Pie0e;	% in SI unit: J*m
QsytQpip2_4pie0e = QsytQpip2_4pie0eSI * 1e30 / kT;	% in kT*nm, >0
r0_PIP2 = LDebye * lambertw(QsytQpip2_4pie0e/ESP/LDebye);	% cut-off distance for KKKK-PIP2 interaction in nm (E=ESP at this distance)


%------------------ interpolation function -----------------

fprintf(' * Interpolation...\n');
n_cutoff = 201;
amin_cutoff = 0.5*avtxmin;
amax_cutoff = 2*avtxmax;	% up to 8 neighbors
da_cutoff = (amax_cutoff - amin_cutoff)/(n_cutoff-1);
r_cutoff = zeros(n_cutoff,1);	% product log
for i = 1 : n_cutoff
	ai = amin_cutoff+(i-1)*da_cutoff;
	r_cutoff(i) = solve_rcutoff(ai,dave,LDebye,LDebye0,LDebyepidb);
end


%------------------ Display -----------------

rcutoffx = interpolate_func(r_cutoff,area_vertex_max,amin_cutoff,da_cutoff,n_cutoff);

fprintf('ESM = %.3g,  LDebye = %.3g,  dave = %.3g, rcutoff = %.3g, rcutoff_SP = %.3g\n',...
	ESM, LDebye, dave, rcutoffx, r0_PIP2);


%% save file
%fid = fopen('e.dat', 'w');
%
%m = 20;
%dzz = Rc2b/m;
%
%for i = 1 : m
%	zz = i*dzz;
%	e = memsyt_att(r4k(1,1),r4k(1,2),zz,atriave,Nsyt,r4k,dESmax,Eescoeff,LDebye,rcutoff);
%	fprintf(fid, '%f\t%f\n', Rc2b-zz, e);
%	fprintf('%f\t%f\n', Rc2b-zz, e);
%end
%
%fclose(fid);



