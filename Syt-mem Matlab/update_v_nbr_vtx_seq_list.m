function v_nbr_vtx_seq_list = update_vnbr_vtx_seq_list(id,x,y,z,vnrm,bond,nnbr)
% update the list of faces around vertex-id in a ccw manner

	v_nbr_vtx_seq_list = zeros(1,10);
	
	nbrlist = find(bond(id,:));	% list of neighboring vertices
	%nbrlist
	
	% project all neighboring bonds to a plane perpendicular to vnrm
	edge = zeros(10,3);
	edgeperp = zeros(10,3);
	edgepara = zeros(10,3);
	for i = 1 : nnbr
		j = nbrlist(i);	% i-th neighboring vertex is #j
		edge(i,1) = x(j) - x(id);
		edge(i,2) = y(j) - y(id);
		edge(i,3) = z(j) - z(id);
		%edgeperp(i,:) = dot(edge(i,:), vnrm) * vnrm;	% perpendicular component of edge-i
		%edgeperp(i,:) = sum(edge(i,:) .* vnrm) * vnrm;	% perpendicular component of edge-i
		edgeperp(i,:) = (edge(i,1)*vnrm(1)+edge(i,2)*vnrm(2)+edge(i,3)*vnrm(3)) * vnrm;	% perpendicular component of edge-i
		edgepara(i,:) = edge(i,:) - edgeperp(i,:);	% parallel component of edge-i
	end
	
	%edge
	
	angle_nbr_vtx = 100*ones(10,1);
	angle_nbr_vtx(1) = 0;
	for i = 2 : nnbr	% 1st element is i=1
		angle_nbr_vtx(i) = getangle(edgepara(1,:), edgepara(i,:), vnrm);	% angle = [0, 2pi)
	end
	
	%angle_nbr_vtx
	
	[temp, list_order] = sort(angle_nbr_vtx,1);	% sort index
	for i = 1 : nnbr
		v_nbr_vtx_seq_list(i) = nbrlist(list_order(i));
	end
	
	v_nbr_vtx_seq_list(nnbr+1:end) = 0;

end




function a = getangle(r1,r2,vnrm)
% get the angle [0, 2pi) between r1 and r2 around axis vnrm

	nr1 = r1/norm(r1);
	nr2 = r2/norm(r2);
	
	%snx = cross(nr1,nr2);
	snx = [nr1(2)*nr2(3)-nr1(3)*nr2(2), nr1(3)*nr2(1)-nr1(1)*nr2(3), nr1(1)*nr2(2)-nr1(2)*nr2(1)];
	
	%sn = dot(snx, vnrm);
	%sn = sum(snx .* vnrm);
	sn = snx(1)*vnrm(1) + snx(2)*vnrm(2) + snx(3)*vnrm(3);
	
	if sn>1
		sn = 1;
	elseif sn<-1
		sn = -1;
	end
	a = asin(sn);
	
	%if(dot(nr1,nr2)>=0)
	%if(sum(nr1 .* nr2)>=0)
	if nr1(1)*nr2(1)+nr1(2)*nr2(2)+nr1(3)*nr2(3) >= 0
		cs = 1;
	else
		cs = -1;
	end
	
	if cs < 0	% 2nd & 3rd quadrants
		a = pi - a;
	elseif cs > 0	&& sn < 0 % 4th quadrant
		a = 2*pi + a;
	end

end
