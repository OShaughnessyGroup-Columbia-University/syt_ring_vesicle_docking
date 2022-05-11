lmin = 1e6;
lmax = 0;
% find lmin
for i = 1:n
    xnew = x(i);
    ynew = y(i);
    znew = z(i);
    dx = x - xnew;
    dy = y - ynew;
    dz = z - znew;
    dl = sqrt(dx.*dx + dy.*dy + dz.*dz);
    dl(i) = [];
    dl = min(dl);
    if dl < lmin
        lmin = dl;
    end
end
for i = 1:n
    nbr = find(bond(i,:));
    for j = nbr
        dx = x(i) - x(j);
        dy = y(i) - y(j);
        dz = z(i) - z(j);
        lbond = sqrt(dx*dx + dy*dy + dz*dz);
        if lbond > lmax
            lmax = lbond;
        end
    end
end