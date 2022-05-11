function dave = getSMdistAve(Nsyt,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri)
% ave distance between KKKK of syt-id to membrane
	
	dave = 0;
	for i = 1 : Nsyt
		dave =  dave + getSMdist(i,r4k,x,y,z,vtx_nrm,n,tricnt,tri_nrm,ntri);
	end
	dave = dave / Nsyt;
end