tic
% monte carlo
% initialize random number generator
rng(S);   % set rng to different initial state each time (by default it goes to same state every time matlab is started)

big = 1e8;	% the factor to enlarge (z-h) when a point goes from external to interal

% relative probability of each kind of move
p1 = 0; % change the position of one vertex
p2 = 0; % flip one bond
p3 = 1; % move syt ring
p4 = 1; % scale the height of all membrane nodes according to height
p5 = 0;	% change the height of all membrane nodes according to power law
p6 = 1;	% change the height of all membrane nodes based on Fourier modes

% normalize and get the cumulative probability
p_sum = p1 + p2 + p3 + p4 + p5 + p6;
p1 = p1 / p_sum;
p2 = p2 / p_sum + p1;
p3 = p3 / p_sum + p2;
p4 = p4 / p_sum + p3;
p5 = p5 / p_sum + p4;
p6 = p6 / p_sum + p5;

% successful steps, used to calculate success rate
nsuc = 0;
ntry1 = 0;
ntry2 = 0;
ntry3 = 0;
ntry4 = 0;
ntry5 = 0;
ntry6 = 0;

nsuc1 = 0;
nsuc2 = 0;
nsuc3 = 0;
nsuc4 = 0;
nsuc5 = 0;
nsuc6 = 0;

% initialize total step count
if ~exist('totstep','var')
    totstep = 0;
end
% initialize total energy
if ~exist('tote','var')
    tote = 0;
end

% initialize count of saved configurations
if ~exist('savecount','var')
    savecount = 0;
end
% initialize indices of fixed vertices
if ~exist('indfixed','var')
    indfixed = [];
end
% current temperature, unit: room temperature (300K)
if ~exist('temperature', 'var')
    temperature = 1;
else
    temperature = temperature * annealing_factor;
end

% minimal total energy in the last cycle
min_tote = tote;


for ismtep = 1 : nstep
    pick = rand;
    %fprintf('(%d)\n', ismtep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if pick < p1    %% move the position of one vertex
        %imove = randi(n-length(indfixed));     % pick one vertex at random
        imove = randi(n);     % pick one vertex at random
        
        while (isetn(imove))	% if at boundary
        	%imove = randi(n-length(indfixed));
        	imove = randi(n);
        end
        
%=====================

%fprintf('%d\t', imove);
				
        % pick a random new position
        xnew = x(imove) + (2*rand - 1) * maxd;
        ynew = y(imove) + (2*rand - 1) * maxd;
        znew = z(imove) + (2*rand - 1) * maxd;
        
        
        % check that the new position is at least dmin away from ALL beads
        xynbr = and(abs(x-xnew) < dmin, abs(y-ynew) < dmin);
        xyznbr = and(xynbr, abs(z-znew) < dmin);
        xyznbr(imove) = 0;
        if any(xyznbr)
            dx = x(xyznbr) - xnew;
            dy = y(xyznbr) - ynew;
            dz = z(xyznbr) - znew;
            dl2 = dx.*dx + dy.*dy + dz.*dz;
            if any(dl2 < dmin2)
                continue;
            end
        end
        
        % check that the new position is at most dmax away from TETHERED
        % NEIGHBORS
        nbr = find(bond(imove,:));
        dx = x(nbr) - xnew;
        dy = y(nbr) - ynew;
        dz = z(nbr) - znew;
        dl2 = dx.*dx + dy.*dy + dz.*dz;
        if any(dl2 > dmax2)
            continue;
        end
        
%=====================
				xtemp = x;
				ytemp = y;
				ztemp = z;
				xtemp(imove) = xnew;	% copy x,y,z to xtemp,ytemp,ztemp, except for point at imove
				ytemp(imove) = ynew;
				ztemp(imove) = znew;
				
        % energy change determines if this proposed step is accepted
        
        %------ area, tension ------
        d_area = 0;
        da_vertex_i = zeros(n,1);	% change in area of each vertex
        atri_new = atri;
        area_vertex_new = area_vertex;
        tri_nrm_new = tri_nrm;
        
        for i = 1 : v_nbr_vtx_cnt(imove)	% for the i-th neighboring vertex, imove-i-(i+1) makes the i-th neighboring face
        	i1 = i + 1;
        	if i1 > v_nbr_vtx_cnt(imove)
        		i1 = 1;
        	end
        	
        	j = v_nbr_vtx_seq_list(imove,i);		% j & k are the IDs of two neighboring vertices, in ccw direction about imove
        	k = v_nbr_vtx_seq_list(imove,i1);
        	itri = v_nbr_tri_seq_list(imove,i);	% ID of i-th face consisting imove-i-(i+1)
        	
        	[tri_nrm_new(itri,:), atri_new(itri)] = get_tri_nrm_area(itri,triver,xtemp,ytemp,ztemp,tri_nrm_sign);
        	
        	d_area_tri_i = atri_new(itri) - atri(itri);
        	da_vtx = d_area_tri_i / 3.0;
        	
        	da_vertex_i(imove) = da_vertex_i(imove) + da_vtx;
        	da_vertex_i(j) = da_vertex_i(j) + da_vtx;
        	da_vertex_i(k) = da_vertex_i(k) + da_vtx;
        	
        	d_area = d_area + d_area_tri_i;	% change in area of all triangles
        end
        
				de_t = d_area * tens;	% change in tension energy
        
        %------- bending -------
        vtx_nrm_new = vtx_nrm;
				mean_curv_vtx_new = mean_curv_vtx;	% make a copy
				Ebend_vtx_new = Ebend_vtx;
				
				de_b = 0;
				for i = 0 : v_nbr_vtx_cnt(imove)
					if i == 0
						j = imove;
					else
						j = v_nbr_vtx_seq_list(imove,i);
					end
					%vtx_nrm_new(j,:) = get_vtx_nrm(j,v_nbr_face,tri_nrm_new,atri_new);	% new normal of vertex
					vtx_nrm_new(j,:) = get_vtx_nrm(j,v_nbr_tri_cnt,v_nbr_tri_seq_list,tri_nrm_new,atri_new);
					
					%area_vertetx_new(j) = get_vertex_area(j,v_nbr_face,atri_new);	% new area of vertex
					area_vertex_new(j) = area_vertex(j) + da_vertex_i(j);
					mean_curv_vtx_new(j) = get_curvature_vtx(j,xnew,ynew,znew,vtx_nrm_new,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list);
					Ebend_vtx_new(j) = kappadb * mean_curv_vtx_new(j)^2 * area_vertex_new(j);
					de_b = de_b + Ebend_vtx_new(j) - Ebend_vtx(j);
				end
        
        %------- electrostatic -------
				es_old = 0;
				es_new = 0;
				for i = 0 : v_nbr_vtx_cnt(imove)
					if i == 0
						j = imove;
					else
						j = v_nbr_vtx_seq_list(imove,i);
					end
					
					if SytMemAtt==1
						es_old_i = memsyt_att(x(j),y(j),z(j),area_vertex(j),Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
						es_new_i = memsyt_att(xtemp(j),ytemp(j),ztemp(j),area_vertex_new(j),Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
					else
						es_old_i = memsyt_att_2(x(j),y(j),z(j),area_vertex(j),vtx_nrm(j),Nsyt,r4k,dESmax,Eescoeff,LDebye);
						es_new_i = memsyt_att_2(xtemp(j),ytemp(j),ztemp(j),area_vertex_new(j),vtx_nrm_new(j),Nsyt,r4k,dESmax,Eescoeff,LDebye);
					end
					
					es_old = es_old + es_old_i;	% sum
					es_new = es_new + es_new_i;
				end
				
				de_s = es_new - es_old;	% change in syt-mem ES energy
        
%				% find closest triangles, indct contains the list of adjacent triangle id
%				indct = find(sum(ismember(triver, imove),2));	% sum(A,dim) sum along the dim-th dimension of A (1=along columns, 2=along rows)
%				
%				% recalculate areas for these triangles
%				newarea = zeros(10,1);
%				d_area = 0;
%				
%				da_vertex_i = zeros(n,1);	% change in area of each vertex
%			
%				for itri = 1:length(indct)	% go through each triangle that contains imove
%					ii = indct(itri);
%					ijk = triver(ii,:);	% find vertex indices
%					ijk(ijk == imove) = [];		% remove the one that equals imove
%					
%					i = ijk(1);		% i and j are the two indices that are not equal to imove
%					j = ijk(2);                    
%				
%					r1 = [x(i) - x(j), y(i) - y(j), z(i) - z(j)];	% find new normal vector
%					r2 = [x(i) - xnew, y(i) - ynew, z(i) - znew];
%					thisnv = cross(r1,r2);
%					this_area = .5 * norm(thisnv);
%					newarea(itri) = this_area;	% store new nv and new area
%					
%					d_area_tri_i = this_area - atri(ii);	% change in area of triangle-itri
%					da_vtx = d_area_tri_i/3.0;	% change in area for each vertex in triangle-itri
%					
%					da_vertex_i(i) = da_vertex_i(i) + da_vtx;
%					da_vertex_i(j) = da_vertex_i(j) + da_vtx;
%					da_vertex_i(imove) = da_vertex_i(imove) + da_vtx;
%                    
%					d_area = d_area + d_area_tri_i;	% change in area of all triangles
%				end
%				
%				de_t = d_area * tens;	% change in tension energy
%        
%        
%				xtemp = x;
%				ytemp = y;
%				ztemp = z;
%				xtemp(imove) = xnew;	% copy x,y,z to xtemp,ytemp,ztemp, except for point at imove
%				ytemp(imove) = ynew;
%				ztemp(imove) = znew;
%				sumdnn = 0;	% calculate the total change in n1n2 due to this move
%				for itri = 1:length(indct)	% go through each triangle
%					ii = indct(itri);
%					for jj = find(trinbr(ii,:))	% go through its neighbors
%						nn_old = get_nn(ii,jj,triver,x,y,z);
%						nn_new = get_nn(ii,jj,triver,xtemp,ytemp,ztemp);
%						dnn = nn_new - nn_old;	% change in n1n2
%						jtri = find(indct == jj);	% check if jj is also a closest neighbor
%						if ~isempty(jtri)	% if jj is also a closest neighbor
%							dnn = .5 * dnn; % half to avoid double counting
%						end
%						sumdnn = sumdnn + dnn;	% accumulate
%					end
%				end
%				
%				de_b = - kappa * sumdnn;	% change in bending energy is - kappa * change in n1n2
%				
%				
%				% syt-mem electrostatic energy
%				% for vertex imove
%				es_old_sum = 0;
%				es_new_sum = 0;
%				
%				area_new = area_vertex(imove) + da_vertex_i(imove);
%				es_old_i = memsyt_att(x(imove),y(imove),z(imove),area_vertex(imove),Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
%				es_new_i = memsyt_att(xnew,ynew,znew,area_new,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
%				
%				es_old_sum = es_old_sum + es_old_i;	% sum
%				es_new_sum = es_new_sum + es_new_i;
%				
%				nbr = find(bond(imove,:));	% all neighboring vertices
%				for i = nbr
%					area_new = area_vertex(i) + da_vertex_i(i);
%					es_old_i = memsyt_att(x(i),y(i),z(i),area_vertex(i),Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
%					es_new_i = memsyt_att(x(i),y(i),z(i),area_new,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
%					es_old_sum = es_old_sum + es_old_i;	% sum
%					es_new_sum = es_new_sum + es_new_i;
%				end
%				
%				de_s = es_new_sum - es_old_sum;	% change in syt-mem ES energy
				
				
				% syt-mem steric repulsion
				e_rep_old = memsyt_rep(x(imove),y(imove),z(imove),Nsyt,rsyt,rsyta,rsytb,Rc2b,Rc2b2,Krepulsion);
				e_rep_new = memsyt_rep(xnew,ynew,znew,Nsyt,rsyt,rsyta,rsytb,Rc2b,Rc2b2,Krepulsion);
				
				de_r = e_rep_new - e_rep_old;	% change in syt-mem repulsion energy
				

				de = de_b + de_t + de_s + de_r;	% change in total energy
				
				ntry1 = ntry1 + 1;
				
				if de > 0 	% decide if this move should be accepted
					if rand > exp(-de/temperature)
						continue;
					end
				end
        

        % update position and successful steps
        x(imove) = xnew;
        y(imove) = ynew;
        z(imove) = znew;
        
	      %atri(indct) = newarea(1:length(indct));
%        for i = 1 : length(indct)	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%        	atri(indct(i)) = newarea(i);
%        end

				% update triangles
				for i = 1 : v_nbr_tri_cnt(imove)	% for each neighboring triangle
					j = v_nbr_tri_seq_list(imove,i);
					atri(j) = atri_new(j);
					tri_nrm(j) = tri_nrm_new(j);
					
        	i1 = triver(j,1);
        	i2 = triver(j,2);
        	i3 = triver(j,3);
        	xc = (x(i1) + x(i2) + x(i3))/3.0;
        	yc = (y(i1) + y(i2) + y(i3))/3.0;
        	zc = (z(i1) + z(i2) + z(i3))/3.0;
        	tricnt(fii,:) = [xc, yc, zc];
				end


        % update vertices
				for i = 0 : v_nbr_vtx_cnt(imove)	% for each neighboring vertex
					if i == 0
						j = imove;
					else
						j = v_nbr_vtx_seq_list(imove,i);
					end
					area_vertex(j) = get_vertex_area(j,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri);
					vtx_nrm(j) = vtx_nrm_new(j);
					mean_curv_vtx(j) = mean_curv_vtx_new(j);
					Ebend_vtx(j) = Ebend_vtx_new(j);
				end
        
        
%        % update ES energy and area per vertex
%        area_vertex(imove) = get_vertex_area(imove,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri);
%        %area_vertex(imove) = area_vertex(imove) + da_vertex_i(imove);	% this generates cumulative error
%        
%        nbr = find(bond(imove,:));
%        for i = nbr
%        	area_vertex(i) = get_vertex_area(i,v_nbr_tri_cnt,v_nbr_tri_seq_list,atri);
%        	%area_vertex(i) = area_vertex(i) + da_vertex_i(i);	% this generates cumulative error
%        end
%        
%        % update face centers
%        fnbr = find(v_nbr_face(imove,:));
%        for fii = fnbr
%        	i = triver(fii,1);
%        	j = triver(fii,2);
%        	k = triver(fii,3);
%        	xc = (x(i) + x(j) + x(k))/3.0;
%        	yc = (y(i) + y(j) + y(k))/3.0;
%        	zc = (z(i) + z(j) + z(k))/3.0;
%        	tricnt(fii,:) = [xc, yc, zc];
%        end
%        
        
        nsuc = nsuc + 1;
        nsuc1 = nsuc1 + 1;
        
%         if isetn(imove)
%             % change in r2bar
%             r2bar = r2bar * (netn-1) + x(imove)*x(imove) + y(imove)*y(imove);
%             r2bar = r2bar / netn;
%         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif pick < p2    %% flip one tether
        % pick a pair of triangles
        [rtri, ctri] = find(trinbr);
        ind = randi(length(rtri));
        i = rtri(ind);
        j = ctri(ind);
        % far-end vertices (not shared between triangles i and j)
        iver = triver(i,:);
        jver = triver(j,:);
        fei = iver(~ismember(iver,jver));
        fej = jver(~ismember(jver,iver));
        % check if the distance between fei and fej is within tether length
        dx = x(fei) - x(fej);
        dy = y(fei) - y(fej);
        dz = z(fei) - z(fej);
        distsq = dx*dx + dy*dy + dz*dz;
        
        if distsq > dmax2 || distsq < dmin2
           continue;
        end
        
        % this operation should not link two external vertices together
        if isetn(fei) && isetn(fej)
            continue;
        end
        
        % check if both shared vertices are connected to at most 8 
        % vertices (they gain one neighbor due to flipping, and must have
        % < 10 neighbors after flipping)
        if sum(bond(fei,:)) > 8 || sum(bond(fej,:)) > 8
        	continue;
        end
        
        % shared vertices
        sv = iver;
        sv(sv == fei) = [];
        sv1 = sv(1);
        sv2 = sv(2);
        
        %if ismember(sv1, indfixed) && ismember(sv2, indfixed)	% not && but || !!!!!!!!!!!!!!!!!!!!!!!
        if (isetn(sv1) || isetn(sv2)) && (isetn(fei) || isetn(fej))	% excluding triangles with 2 boundary vertices
					continue;
        end
            
        % check if both shared vertices are connected to at least 4
        % vertices (they lose one neighbor due to flipping, and must have 3
        % neighbors after flipping)
        if sum(bond(sv1,:)) < 4 || sum(bond(sv2,:)) < 4
            continue;
        end
        
%         % compute new normal vectors
%         r1 = [x(fei) - x(sv1), y(fei) - y(sv1), z(fei) - z(sv1)];
%         r2 = [x(fej) - x(sv1), y(fej) - y(sv1), z(fej) - z(sv1)];
%         newnv1 = cross(r1,r2);
%         newnv1 = newnv1 / norm(newnv1);
%         r1 = [x(fei) - x(sv2), y(fei) - y(sv2), z(fei) - z(sv2)];
%         r2 = [x(fej) - x(sv2), y(fej) - y(sv2), z(fej) - z(sv2)];
%         newnv2 = cross(r1,r2);
%         newnv2 = newnv2 / norm(newnv2);

        % old neighbors of i
        inbr = find(trinbr(i,:));
        inbr(inbr == j) = [];
        %temp = triver(inbr(1),:);
        if any(triver(inbr(1),:) == sv2)
            inbr = inbr([2,1]);
        end
        jnbr = find(trinbr(j,:));
        jnbr(jnbr == i) = [];
        
        if any(triver(jnbr(1),:) == sv2)
            jnbr = jnbr([2,1]);
        end
        % normal vectors of inbr and jnbr. If # nbr < 2, fill with zeros
        nvinbr = zeros(2,3);
        nvjnbr = nvinbr;
        nvinbr(1:length(inbr),:) = nv(inbr,:);
        nvjnbr(1:length(jnbr),:) = nv(jnbr,:);
        % old n1n2, before flipping the bond
        oldnn = get_nn(i,j,triver,x,y,z);
        oldnn = oldnn + get_nn(i,inbr(1),triver,x,y,z);
        oldnn = oldnn + get_nn(i,inbr(2),triver,x,y,z);
        oldnn = oldnn + get_nn(j,jnbr(1),triver,x,y,z);
        oldnn = oldnn + get_nn(j,jnbr(2),triver,x,y,z);
        % new neighbors of the new triangle containing vertex sv1
        temp_triver = triver;
        temp_triver(i,:) = [fei, fej, sv1];	% new triangle i' is fei-fej-sv1
        temp_triver(j,:) = [fei, fej, sv2];	% new triangle j' is fei-fej-sv2
        newnn = get_nn(i,j,temp_triver,x,y,z);
        newnn = newnn + get_nn(i,inbr(1),temp_triver,x,y,z);
        newnn = newnn + get_nn(i,jnbr(1),temp_triver,x,y,z);
        newnn = newnn + get_nn(j,inbr(2),temp_triver,x,y,z);
        newnn = newnn + get_nn(j,jnbr(2),temp_triver,x,y,z);
        
        % change in bending energy is - kappa * change in n1n2
        de_b = kappa * (oldnn - newnn);
        
        new_atri_i = get_atri(i,temp_triver,x,y,z);	% area of new triangle i'
        new_atri_j = get_atri(j,temp_triver,x,y,z);	% area of new triangle j'
        
        % change in area
        d_area = new_atri_i + new_atri_j - atri(i) - atri(j); % change in area
        de_t = tens * d_area;
        
        % change in mem-syt attraction
        da_fei = (-atri(i) + new_atri_i + new_atri_j)/3.0;	% change in the area of each vertex
        da_fej = (-atri(j) + new_atri_i + new_atri_j)/3.0;
        da_sv1 = (-atri(i) - atri(j) + new_atri_i)/3.0;
        da_sv2 = (-atri(i) - atri(j) + new_atri_j)/3.0;
        
        % change in mem-syt energy is proportional to change in area of each vertex
        de_s_fei = memsyt_att(x(fei),y(fei),z(fei),da_fei,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
        de_s_fej = memsyt_att(x(fej),y(fej),z(fei),da_fej,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
        de_s_sv1 = memsyt_att(x(sv1),y(sv1),z(sv1),da_sv1,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
        de_s_sv2 = memsyt_att(x(sv2),y(sv2),z(sv2),da_sv2,Nsyt,r4k,dESmax,Eescoeff,LDebye,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
        
%        de_s_fei = es_vertex(fei)*da_fei/area_vertex(fei);	% ES energy is proportional to area_vertex
%        de_s_fej = es_vertex(fej)*da_fej/area_vertex(fej);	% problem: es_vertex() may not be the latest value
%        de_s_sv1 = es_vertex(sv1)*da_sv1/area_vertex(sv1);
%        de_s_sv2 = es_vertex(sv2)*da_sv2/area_vertex(sv2);
        de_s = de_s_fei + de_s_fej + de_s_sv1 + de_s_sv2;
        
        % total change
        de = de_b + de_t + de_s;
        
				ntry2 = ntry2 + 1;
				
        % decide if this move should be accepted
        if de > 0 
        	if rand > exp(-de/temperature)	% reject
            continue;
          end
        end
        
        % update bond
        bond(sv1, sv2) = 0;
        bond(sv2, sv1) = 0;
        bond(fei, fej) = 1;
        bond(fej, fei) = 1;
        
        %fprintf('flipping %d-%d to %d-%d\n', sv1,sv2,fei,fej);
        
        % update triver 
        triver(i,:) = [fei, fej, sv1];
        triver(j,:) = [fei, fej, sv2];
        
        % update trinbr
        trinbr(i,inbr(2)) = false;
        trinbr(inbr(2),i) = false;
        trinbr(i,jnbr(1)) = true;
        trinbr(jnbr(1),i) = true;
        trinbr(j,jnbr(1)) = false;
        trinbr(jnbr(1),j) = false;
        trinbr(j,inbr(2)) = true;
        trinbr(inbr(2),j) = true;
        
        % update neighboring faces & area
        [area_vertex(fei), v_nbr_face(fei,:)] = update_nbr_face_of_mem_vertex(fei,ntri,bond,triver,atri);
        [area_vertex(fej), v_nbr_face(fej,:)] = update_nbr_face_of_mem_vertex(fej,ntri,bond,triver,atri);
        [area_vertex(sv1), v_nbr_face(sv1,:)] = update_nbr_face_of_mem_vertex(sv1,ntri,bond,triver,atri);
        [area_vertex(sv2), v_nbr_face(sv2,:)] = update_nbr_face_of_mem_vertex(sv2,ntri,bond,triver,atri);
        
        atri(i) = new_atri_i;
        atri(j) = new_atri_j;
        
        % update 2 face centers
        for nf = 1 : 2
        	if nf == 1
        		fij = i;
        	else
        		fij = j;
        	end
        	ii = triver(fij,1);
        	jj = triver(fij,2);
        	kk = triver(fij,3);
        	xc = (x(ii) + x(jj) + x(kk))/3.0;
        	yc = (y(ii) + y(jj) + y(kk))/3.0;
        	zc = (z(ii) + z(jj) + z(kk))/3.0;
        	tricnt(fij,:) = [xc, yc, zc];
        end
        
%         trinbr(i,:) = zeros(1,ntri);
%         trinbr(:,i) = zeros(ntri,1);
%         trinbr(j,:) = zeros(1,ntri);
%         trinbr(:,j) = zeros(ntri,1);
%         trinbr(i,[j, newnbr1]) = ones(1,length(newnbr1) + 1);
%         trinbr([j, newnbr1],i) = ones(length(newnbr1) + 1,1);
%         trinbr(j,[i, newnbr2]) = ones(1,length(newnbr2) + 1);
%         trinbr([i, newnbr2],j) = ones(length(newnbr2) + 1,1);
%         % update nv
%         nv(i,:) = newnv1;
%         nv(j,:) = newnv2;
        % update successful steps

        nsuc = nsuc + 1;
        nsuc2 = nsuc2 + 1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif pick < p3    %% move syt ring
    	
    	dz = (1 - 2*rand) * maxd;
    	Zsytnew = Zsyt + dz;
    	
%    	dr = zeros(3,1);
%    	dr(1) = (1 - 2*rand) * maxd;
%    	dr(2) = (1 - 2*rand) * maxd;
%    	dr(3) = (1 - 2*rand) * maxd;
%    	Zsytnew = Zsyt + dr(3);
    	
    	
    	for i = 1 : Nsyt
    		for j = 1 : 2
		    	rsytnew(i,j) = rsyt(i,j);
		    	rsytanew(i,j) = rsyta(i,j);
		    	rsytbnew(i,j) = rsytb(i,j);
		    	r4knew(i,j) = r4k(i,j);
		    end
		    
    		rsytnew(i,3) = Zsytnew;	% shift in dz
    		rsytanew(i,3) = Zsytnew;
    		rsytbnew(i,3) = Zsytnew;
    		r4knew(i,3) = Zsytnew + dz4k(i);
    		
%    		for j = 1 : 3
%    			rsytnew(i,j) = rsyt(i,j) + dr(j);	% shift in dz
%	    		rsytanew(i,j) = rsyta(i,j) + dr(j);
%	    		rsytbnew(i,j) = rsytb(i,j) + dr(j);
%	    		r4knew(i,j) = r4k(i,j) + dr(j);
%	    	end
    	end
    	if SytMemAtt==1
				es_old = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
				es_new = sytmem_att(Nsyt,n,r4knew,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
			else
				es_old = sytmem_att_2(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn);
				es_new = sytmem_att_2(Nsyt,n,r4knew,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn);
			end
			
			if Extra_PIP2==1	% extra 1-to-1 syt-PIP2 interaction
				d1 = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% old ave. distance
				d2 = getSMdistAve(Nsyt,r4knew,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% new ave. distance
				esp_old = sytpip2_att(d1,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				esp_new = sytpip2_att(d2,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				es_old = es_old + esp_old;
				es_new = es_new + esp_new;
			end
			
    	de_s = es_new - es_old;
			
			
			er_old = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion);
			er_old = er_old + sytmemF_rep(Nsyt,rsyt,rsyta,rsytb,tricnt,ntri,Rc2b,Rc2b2,Krepulsion);
			er_new = sytmemN_rep(Nsyt,rsytnew,rsytanew,rsytbnew,x,y,z,n,Rc2b,Rc2b2,Krepulsion);
			er_new = er_new + sytmemF_rep(Nsyt,rsytnew,rsytanew,rsytbnew,tricnt,ntri,Rc2b,Rc2b2,Krepulsion);
    	
    	de_r = er_new - er_old;
    	
    	de = de_s + de_r;
    	
			ntry3 = ntry3 + 1;
			
    	if de > 0
    		if rand > exp(-de/temperature)
    			continue;
    		end
    	end
    	
    	% keep trial value
    	Zsyt = Zsytnew;
		 	for i = 1 : Nsyt
		    rsyt(i,3) = rsytnew(i,3);	% keep change
		   	rsyta(i,3) = rsytanew(i,3);
		    rsytb(i,3) = rsytbnew(i,3);
    		r4k(i,3) = r4knew(i,3);

%				for j = 1 : 3
%		    	rsyt(i,j) = rsyt(i,j) - dr(j);	% shift in dz
%					rsyta(i,j) = rsyta(i,j) - dr(j);
%					rsytb(i,j) = rsytb(i,j) - dr(j);
%			  	r4k(i,j) = r4k(i,j) - dr(j);
%			 	end
			end
    	
      nsuc = nsuc + 1;
      nsuc3 = nsuc3 + 1;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif pick < p4    %% scale the height of all mem nodes and syts
    	
    	dlt_z = 0.05;
    	fct = 1 + (1-2*rand)*dlt_z;
    	
    	% ====== shifting all syts ======
    	rsytnew = rsyt;	% make a copy
    	rsytanew = rsyta;
    	rsytbnew = rsytb;
    	r4knew = r4k;
    	
    	if ShiftSyt==0	% not moving syt
	    	Zsytnew = Zsyt;
	    else	% moving syt with membrane
    		Zsytnew = fct*Zsyt;
	    	for i = 1 : Nsyt
	    		rsytnew(i,3) = Zsytnew;
	    		rsytanew(i,3) = Zsytnew;
	    		rsytbnew(i,3) = Zsytnew;
    			r4knew(i,3) = Zsytnew + dz4k(i);
	    	end
			end
    	
    	% ====== shifting membrane ======
    	znew = z;
    	for i = 1 : n
    		if ~isetn(i)
    			znew(i) = fct*z(i);	% temporary move
    		end
    	end
    	
    	tricntnew = tricnt;	% coord. of mem face center
    	for i = 1 : ntri
    		i1 = triver(i,1);
    		i2 = triver(i,2);
    		i3 = triver(i,3);
    		tricntnew(i, 3) = (znew(i1) + znew(i2) + znew(i3))/3.0;	% only update zc
    	end
    		
    	
			% ====== membrane tension and bending ======
    	tri_nrm_new = tri_nrm;
    	atri_new = atri;
    	vtx_nrm_new = vtx_nrm;
    	area_vertex_new = area_vertex;
    	mean_curv_vtx_new = zeros(n, 1);
    	Ebend_vtx_new = zeros(n, 1);
    	
    	[tri_nrm_new,atri_new,vtx_nrm_new,area_vertex_new,mean_curv_vtx_new,Ebend_vtx_new,de_t,de_b] = ...
				mem_energy_ab(x,y,z,znew,n,ntri,tri_nrm_sign,tri_nrm,atri,vtx_nrm,area_vertex,triver,...
				tens,isetn,v_nbr_tri_cnt,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list,kappadb);
			
			% ====== syt-mem electrostatic energy ======
			if SytMemAtt==1
				es_old = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
				es_new = sytmem_att(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
			else
				es_old = sytmem_att_2(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn);
				es_new = sytmem_att_2(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,vtx_nrm_new,LDebye,isetn);
			end
			
			if Extra_PIP2==1	% extra 1-to-1 syt-PIP2 interaction
				d1 = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% old ave. distance
				d2 = getSMdistAve(Nsyt,r4knew,x,y,znew,vtx_nrm_new,n,tricnt,tri_nrm_new,ntri);	% new ave. distance
				esp_old = sytpip2_att(d1,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				esp_new = sytpip2_att(d2,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				es_old = es_old + esp_old;
				es_new = es_new + esp_new;
			end
			
			de_s = es_new - es_old;	% change in syt-mem ES energy
			
			% ====== syt-mem repulsion ======
			e_rep_old = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_old = e_rep_old + sytmemF_rep(Nsyt,rsyt,rsyta,rsytb,tricnt,ntri,Rc2b,Rc2b2,Krepulsion);
			
			e_rep_new = sytmemN_rep(Nsyt,rsytnew,rsytanew,rsytbnew,x,y,znew,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_new = e_rep_new + sytmemF_rep(Nsyt,rsytnew,rsytanew,rsytbnew,tricntnew,ntri,Rc2b,Rc2b2,Krepulsion);
			
			de_r = e_rep_new - e_rep_old;	% change in syt-mem repulsion energy
			
			de = de_t + de_b + de_s + de_r;	% change in total energy
			
	    %fprintf('%.3g\t%.3g\t%.3g\t%.3g\n', de_b, de_t, de_s, de_r);	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			ntry4 = ntry4 + 1;
			
			if de > 0 	% decide if this move should be accepted
				if rand > exp(-de/temperature)
					continue;
				end
			end
	    	
	    % update position and successful steps
	   	Zsyt = Zsytnew;
			for i = 1 : Nsyt
			  rsyt(i,3) = rsytnew(i,3);	% keep change
		   	rsyta(i,3) = rsytanew(i,3);
			  rsytb(i,3) = rsytbnew(i,3);
	   		r4k(i,3) = r4knew(i,3);
	   	end
	    
	    z = znew;
	    tri_nrm = tri_nrm_new;
	    atri = atri_new;
	    vtx_nrm = vtx_nrm_new;
	    area_vertex = area_vertex_new;
	    mean_curv_vtx = mean_curv_vtx_new;
	    Ebend_vtx = Ebend_vtx_new;
	    
	    % update face center
			for i = 1 : ntri
				i1 = triver(i,1);
				i2 = triver(i,2);
				i3 = triver(i,3);
				%xc = (x(i1) + x(i2) + x(i3))/3.0;	% no change
				%yc = (y(i1) + y(i2) + y(i3))/3.0;
				zc = (z(i1) + z(i2) + z(i3))/3.0;
				%tricnt(i,:) = [xc, yc, zc];
				tricnt(i,3) = zc;
			end
    	
      nsuc = nsuc + 1;
      nsuc4 = nsuc4 + 1;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif pick < p5    %% change the height of all membrane nodes according to power law
    	
    	rm = sidehalf;
    	%rm = diaghalf;
    	
    	%dlt_z0 = maxd * (1 - 2*rand);	% z-shift at center
    	dlt_z0 = 0.1 * (1 - 2*rand);	% z-shift at center
    	
    	pwr = 0.5*(4)^rand;	% a random number between 0.5 and 2, on log scale
    	
    	% ====== shifting all syts ======
    	rsytnew = rsyt;	% make a copy
    	rsytanew = rsyta;
    	rsytbnew = rsytb;
    	r4knew = r4k;
    	
    	if ShiftSyt==0	% not moving syt
    		Zsytnew = Zsyt;
			else	% moving syt with membrane
	    	r = sqrt(rsyt(1,1)^2 + rsyt(1,2)^2);
	    	if r >= rm
	    		dlt_zs = 0;
	    	else
	    		dlt_zs = dlt_z0*(1-r/rm)^pwr;
	    	end
	    	
	    	Zsytnew = Zsyt + dlt_zs;	% temporary move
	    	for i = 1 : Nsyt
	    		rsytnew(i,3) = Zsytnew;
	    		rsytanew(i,3) = Zsytnew;
	    		rsytbnew(i,3) = Zsytnew;
    			r4knew(i,3) = Zsytnew + dz4k(i);
	    	end
	    end
    	
    	% ====== shifting membrane ======
    	znew = z;
    	atri_new = atri;
    	area_vertex_new = area_vertex;
    	
    	for i = 1 : n
    		if ~isetn(i)
		    	r = sqrt(x(i)^2 + y(i)^2);
		    	if r >= rm
		    		dlt_z = 0;
		    	else
		    		dlt_z = dlt_z0*(1-r/rm)^pwr;
		    	end
		    	
    			znew(i) = z(i) + dlt_z;	% temporary move
    		end
    	end
    	
    	tricntnew = tricnt;	% coord. of mem face center
    	for i = 1 : ntri
    		i1 = triver(i,1);
    		i2 = triver(i,2);
    		i3 = triver(i,3);
    		tricntnew(i, 3) = (znew(i1) + znew(i2) + znew(i3))/3.0;	% only update zc
    	end
    	
			% ====== membrane tension and bending ======
    	tri_nrm_new = tri_nrm;
    	atri_new = atri;
    	vtx_nrm_new = vtx_nrm;
    	area_vertex_new = area_vertex;
    	mean_curv_vtx_new = zeros(n, 1);
    	Ebend_vtx_new = zeros(n, 1);
    	
    	[tri_nrm_new,atri_new,vtx_nrm_new,area_vertex_new,mean_curv_vtx_new,Ebend_vtx_new,de_t,de_b] = ...
				mem_energy_ab(x,y,z,znew,n,ntri,tri_nrm_sign,tri_nrm,atri,vtx_nrm,area_vertex,triver,...
				tens,isetn,v_nbr_tri_cnt,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list,kappadb);
			
			% ====== syt-mem electrostatic energy ======
			if SytMemAtt==1
				es_old = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
				es_new = sytmem_att(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
			else
				es_old = sytmem_att_2(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn);
				es_new = sytmem_att_2(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,vtx_nrm_new,LDebye,isetn);
			end
			
			if Extra_PIP2==1	% extra 1-to-1 syt-PIP2 interaction
				d1 = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% old ave. distance
				d2 = getSMdistAve(Nsyt,r4knew,x,y,znew,vtx_nrm_new,n,tricnt,tri_nrm_new,ntri);	% new ave. distance
				esp_old = sytpip2_att(d1,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				esp_new = sytpip2_att(d2,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				es_old = es_old + esp_old;
				es_new = es_new + esp_new;
			end
			
			de_s = es_new - es_old;	% change in syt-mem ES energy
			
			% ====== syt-mem repulsion ======
			e_rep_old = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_old = e_rep_old + sytmemF_rep(Nsyt,rsyt,rsyta,rsytb,tricnt,ntri,Rc2b,Rc2b2,Krepulsion);
			
			e_rep_new = sytmemN_rep(Nsyt,rsytnew,rsytanew,rsytbnew,x,y,znew,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_new = e_rep_new + sytmemF_rep(Nsyt,rsytnew,rsytanew,rsytbnew,tricntnew,ntri,Rc2b,Rc2b2,Krepulsion);
			
			de_r = e_rep_new - e_rep_old;	% change in syt-mem repulsion energy
			
			de = de_t + de_b + de_s + de_r;	% change in total energy
		
    %fprintf('%.3g\t%.3g\t%.3g\t%.3g\n', de_b, de_t, de_s, de_r);	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %fprintf('%.3g\t%.3g\n', e_rep_old, e_rep_new);
		
		ntry5 = ntry5 + 1;
		
		if de > 0 	% decide if this move should be accepted
			if rand > exp(-de/temperature)
				continue;
			end
		end
    	
    % update position and successful steps
   	Zsyt = Zsytnew;
		for i = 1 : Nsyt
		  rsyt(i,3) = rsytnew(i,3);	% keep change
	   	rsyta(i,3) = rsytanew(i,3);
		  rsytb(i,3) = rsytbnew(i,3);
   		r4k(i,3) = r4knew(i,3);
   	end
    
	  z = znew;
	  tri_nrm = tri_nrm_new;
	  atri = atri_new;
	  vtx_nrm = vtx_nrm_new;
	  area_vertex = area_vertex_new;
	  mean_curv_vtx = mean_curv_vtx_new;
	  Ebend_vtx = Ebend_vtx_new;
	    
    % update face center
		for i = 1 : ntri
			i1 = triver(i,1);
			i2 = triver(i,2);
			i3 = triver(i,3);
			%xc = (x(i1) + x(i2) + x(i3))/3.0;	% no change
			%yc = (y(i1) + y(i2) + y(i3))/3.0;
			zc = (z(i1) + z(i2) + z(i3))/3.0;
			%tricnt(i,:) = [xc, yc, zc];
			tricnt(i,3) = zc;
		end
    	
		nsuc = nsuc + 1;
		nsuc5 = nsuc5 + 1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   elseif pick < p6    %% change the height of all membrane nodes based on Fourier modes
   	
   	rm = sidehalf;
   	%rm = diaghalf;
   	
   	%----------------
%   	nmode = 3;	% number of Fourier modes
%   	ci = zeros(nmode,1);	% coefficient of Fourier modes
%   	for i = 1 : nmode
%   		ci(i) = 0.5 * maxd * (1 - 2*rand);	% 0.5 comes from 0.5*(1+cos(i*pi*r/rm))
%   	end
   	%----------------
		nmodemax = floor(0.1*rm/dmin);	% max # of modes (0.5*rm/dmin)
		nmode = floor(nmodemax*rand+0.99);	% pick a random mode
		%ci = 0.5 * maxd * (1 - 2*rand);	% 0.5 comes from 0.5*(1+cos(i*pi*r/rm))
		ci = 0.5 * 0.1 * (1 - 2*rand);	% 0.5 comes from 0.5*(1+cos(i*pi*r/rm))
   	%----------------
   	
   	% ====== shifting all syts ======
   	rsytnew = rsyt;	% make a copy
   	rsytanew = rsyta;
   	rsytbnew = rsytb;
   	r4knew = r4k;
   	
    if ShiftSyt==0	% not moving syt
	   	Zsytnew = Zsyt;
	  else	% moving syt with membrane
			r = sqrt(rsyt(1,1)^2 + rsyt(1,2)^2);
	   	dlt_zs = 0;
	   	if r < rm
	   	%----------------
	%	   	for i = 1 : nmode
	%				dlt_zs = dlt_zs + ci(i)*(1+cos(i*pi*r/rm));	% Fourier modes
	%			end
	   	%----------------
	   		dlt_zs = ci*(1+cos(nmode*pi*r/rm));
	   	%----------------
			end
	   	
	   	Zsytnew = Zsyt + dlt_zs;	% temporary move
	   	for i = 1 : Nsyt
	   		rsytnew(i,3) = Zsytnew;
	   		rsytanew(i,3) = Zsytnew;
	   		rsytbnew(i,3) = Zsytnew;
    		r4knew(i,3) = Zsytnew + dz4k(i);
	   	end
		end

    
   	% ====== shifting membrane ======
   	znew = z;
   	atri_new = atri;
   	area_vertex_new = area_vertex;
   	
   	for i = 1 : n
   		if ~isetn(i)
	    	r = sqrt(x(i)^2 + y(i)^2);
   			dlt_z = 0;
	    	if r < rm
   	%----------------
%		    	for j = 1 : nmode
%			    	dlt_z = dlt_z + ci(j)*(1+cos(j*pi*r/rm));
%			    end
   	%----------------
   			dlt_z = ci*(1+cos(nmode*pi*r/rm));
   	%----------------
			  end
		    
   			znew(i) = z(i) + dlt_z;	% temporary move
   		end
   	end
   	
    	tricntnew = tricnt;	% coord. of mem face center
    	for i = 1 : ntri
    		i1 = triver(i,1);
    		i2 = triver(i,2);
    		i3 = triver(i,3);
    		tricntnew(i, 3) = (znew(i1) + znew(i2) + znew(i3))/3.0;	% only update zc
    	end
    	
			% ====== membrane tension and bending ======
    	tri_nrm_new = tri_nrm;
    	atri_new = atri;
    	vtx_nrm_new = vtx_nrm;
    	area_vertex_new = area_vertex;
    	mean_curv_vtx_new = zeros(n, 1);
    	Ebend_vtx_new = zeros(n, 1);
    	
    	[tri_nrm_new,atri_new,vtx_nrm_new,area_vertex_new,mean_curv_vtx_new,Ebend_vtx_new,de_t,de_b] = ...
				mem_energy_ab(x,y,z,znew,n,ntri,tri_nrm_sign,tri_nrm,atri,vtx_nrm,area_vertex,triver,...
				tens,isetn,v_nbr_tri_cnt,v_nbr_vtx_cnt,v_nbr_tri_seq_list,v_nbr_vtx_seq_list,kappadb);
			
			% ====== syt-mem electrostatic energy ======
			if SytMemAtt==1
				es_old = sytmem_att(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
				es_new = sytmem_att(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,LDebye,isetn,r_cutoff,amin_cutoff,da_cutoff,n_cutoff);
			else
				es_old = sytmem_att_2(Nsyt,n,r4k,x,y,z,dESmax,Eescoeff,area_vertex,vtx_nrm,LDebye,isetn);
				es_new = sytmem_att_2(Nsyt,n,r4knew,x,y,znew,dESmax,Eescoeff,area_vertex_new,vtx_nrm_new,LDebye,isetn);
			end
			
			if Extra_PIP2==1	% extra 1-to-1 syt-PIP2 interaction
				d1 = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);	% old ave. distance
				d2 = getSMdistAve(Nsyt,r4knew,x,y,znew,vtx_nrm_new,n,tricnt,tri_nrm_new,ntri);	% new ave. distance
				esp_old = sytpip2_att(d1,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				esp_new = sytpip2_att(d2,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt);
				es_old = es_old + esp_old;
				es_new = es_new + esp_new;
			end
			
			de_s = es_new - es_old;	% change in syt-mem ES energy
			
			% ====== syt-mem repulsion ======
			e_rep_old = sytmemN_rep(Nsyt,rsyt,rsyta,rsytb,x,y,z,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_old = e_rep_old + sytmemF_rep(Nsyt,rsyt,rsyta,rsytb,tricnt,ntri,Rc2b,Rc2b2,Krepulsion);
			
			e_rep_new = sytmemN_rep(Nsyt,rsytnew,rsytanew,rsytbnew,x,y,znew,n,Rc2b,Rc2b2,Krepulsion);
			e_rep_new = e_rep_new + sytmemF_rep(Nsyt,rsytnew,rsytanew,rsytbnew,tricntnew,ntri,Rc2b,Rc2b2,Krepulsion);
			
			de_r = e_rep_new - e_rep_old;	% change in syt-mem repulsion energy
			
			de = de_t + de_b + de_s + de_r;	% change in total energy
		
    %fprintf('%.3g\t%.3g\t%.3g\t%.3g\n', de_b, de_t, de_s, de_r);	% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %fprintf('%.3g\t%.3g\n', e_rep_old, e_rep_new);
		
		ntry6 = ntry6 + 1;
		
		if de > 0 	% decide if this move should be accepted
			if rand > exp(-de/temperature)
				continue;
			end
		end
    
    %fprintf('.');
    
    % update position and successful steps
   	Zsyt = Zsytnew;
		for i = 1 : Nsyt
			rsyt(i,3) = rsytnew(i,3);	% keep change
			rsyta(i,3) = rsytanew(i,3);
			rsytb(i,3) = rsytbnew(i,3);
			r4k(i,3) = r4knew(i,3);
   	end
    
	  z = znew;
	  tri_nrm = tri_nrm_new;
	  atri = atri_new;
	  vtx_nrm = vtx_nrm_new;
	  area_vertex = area_vertex_new;
	  mean_curv_vtx = mean_curv_vtx_new;
	  Ebend_vtx = Ebend_vtx_new;
    
    % update face center
		for i = 1 : ntri
			i1 = triver(i,1);
			i2 = triver(i,2);
			i3 = triver(i,3);
			%xc = (x(i1) + x(i2) + x(i3))/3.0;	% no change
			%yc = (y(i1) + y(i2) + y(i3))/3.0;
			zc = (z(i1) + z(i2) + z(i3))/3.0;
			%tricnt(i,:) = [xc, yc, zc];
			tricnt(i,3) = zc;
		end
    	
    	
		nsuc = nsuc + 1;
		nsuc6 = nsuc6 + 1;
    	
	end %--------------------- end of if P1...P6 ---------------------
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%	if de > 1
%		if pick < p1
%			mcmove=1;
%		elseif pick < p2
%			mcmove=2;
%		elseif pick < p3
%			mcmove=3;
%		elseif pick < p4
%			mcmove=4;
%		elseif pick < p5
%			mcmove=5;
%		elseif pick < p6
%			mcmove=6;
%		else
%			mcmove=0;
%		end
%		fprintf('%d: dE = %.3g\n', mcmove, de);
%%		fprintf('er1 = %.3g, er2 = %.3g\n', er_old, er_new);
%%		fprintf('der = %.3g, des = %.3g, de = %.3g\n', de_r, de_s, de);
%	end
	
	tote = tote + de;
	if tote < min_tote
		min_tote = tote;
		best_x = x;
		best_y = y;
		best_z = z;
		best_atri = atri;
		best_bond = bond;
		best_triver = triver;
		best_trinbr = trinbr;
	end
		
		
end	% end of for ismtep


% update total steps
totstep = totstep + nstep;
% success rate
sucrate = nsuc / nstep;
sucrate1 = nsuc1 / (ntry1 + 1e-8);
sucrate2 = nsuc2 / (ntry2 + 1e-8);
sucrate3 = nsuc3 / (ntry3 + 1e-8);
sucrate4 = nsuc4 / (ntry4 + 1e-8);
sucrate5 = nsuc5 / (ntry5 + 1e-8);
sucrate6 = nsuc6 / (ntry6 + 1e-8);


te = toc;

% restore the system to the lowest energy state in this cycle
if annealing_factor < 1 && exist('best_x','var')
    tote = min_tote;
    x = best_x;
    y = best_y;
    z = best_z;
    atri = best_atri;
    bond = best_bond;
    triver = best_triver;
    trinbr = best_trinbr;
end

