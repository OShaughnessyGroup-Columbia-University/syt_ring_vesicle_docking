function atri_i = get_atri(i,triver,x,y,z)
    iver = triver(i,:);
        
    dx_1 = x(iver(1)) - x(iver(2));
    dy_1 = y(iver(1)) - y(iver(2));
    dz_1 = z(iver(1)) - z(iver(2));
%     dr_i = [dx_i; dy_i; dz_i];
    
    dx_2 = x(iver(1)) - x(iver(3));
    dy_2 = y(iver(1)) - y(iver(3));
    dz_2 = z(iver(1)) - z(iver(3));
%     dr_j = [dx_j; dy_j; dz_j];
    
%     ni = cross(dr_shared, dr_i);
    ni = [dy_1*dz_2-dz_1*dy_2; dz_1*dx_2-dx_1*dz_2; dx_1*dy_2-dy_1*dx_2];
    atri_i = .5 * norm(ni);
end