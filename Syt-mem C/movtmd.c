
void TMDmotion(CHAIN *CH)
{
	int i;
	double dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i,dr)
#endif

		for(i=0; i<(p->nseg); i++) {	// motion of individual subunits
			vecprod(p->ftmd[i], MuTMDdt, dr);
			limitvector(dr, drmemmin);
			vecadd(p->rtmd[i], dr, p->rtmd[i]);
		}
		
		p = p->next;
	}
}


void TMDfluctuation(CHAIN *CH)
{
	int i;
	double dr[3];
	CHAIN *p;
	
#if OPENMP
	int tid;
#endif
	
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i,dr,tid)
#endif

		for(i=0; i<(p->nseg); i++) {
		#if OPENMP
			tid = omp_get_thread_num();
			dr[0]=dTMDnoise*(1-2*RNDMP(tid));
			dr[1]=dTMDnoise*(1-2*RNDMP(tid));
			dr[2]=dTMDnoise*(1-2*RNDMP(tid));
		#else
			dr[0]=dTMDnoise*(1-2*RND());
			dr[1]=dTMDnoise*(1-2*RND());
			dr[2]=dTMDnoise*(1-2*RND());
		#endif
			vecadd(p->rtmd[i], dr, p->rtmd[i]);
		}
		
		p = p->next;
	}
}



void getTMDVesPos(double rtmd[3], double rtmdnew[3])
{	// interpolate membrane to grid, run after getNbrPositions() !
	int i, j, ni, nj, flg, nnbr;
	int id, ii, jj;
	double tmdrc[3], tmdrcnew[3], d2, d2min;
	double r01[3], r02[3];
	double rnbr[NnbrMax][3], nrm[NnbrMax][3];
	VNODE *p;
	
	vecsub(rtmd, RCves, tmdrc);	// coord. of TMD in vesicle's frame
	//printf("%.3g, %.3g, %.3g\n", RCves[0], RCves[1], RCves[2]);
	
	d2min = inf;
	id = -1;
	for(i=0; i<NVnode; i++) {	// find the nearest vnode on vesicle
		p = &vnode[i];
		if(dotprod(p->rc, tmdrc)<0) continue;
		d2 = distance2(p->rc, tmdrc);
		if(d2<d2min) {
			d2min = d2;
			id = i;	// vnode-id is the nearest vnode to tmd
		}
	}
	if(id<0) {
		printf("Error in getTMDVesPos().\n");
		printf("rtmd = (%.3g, %.3g, %.3g)\n", rtmd[0], rtmd[1], rtmd[2]);
		exit(0);
	}
	
	p = &vnode[id];	// p is the nearest vesicle vnode to tmd
	nnbr = p->nnbr;
	for(i=0; i<nnbr; i++) {	// i-th neighbor
		j = (i+1)%nnbr;	// j is the (i+1)-th neighbor
		ni = p->nbr[i];	// id of i-th neighbor
		nj = p->nbr[j];	// id of (i+1)-th neighbor
		veccopy(vnode[ni].rc, rnbr[i]);	// location of the i-th neighbor in RC's frame
		
		vecsub(vnode[ni].rc, p->rc, r01);
		vecsub(vnode[nj].rc, p->rc, r02);
		crossprod(r01, r02, nrm[i]);
		normalize(nrm[i], nrm[i]);	// nrm[i] is the inward/outward normal of triangle id-ni-nj
		if(dotprod(nrm[i], p->rc)<0) vecprod(nrm[i], -1, nrm[i]);	// outward normal
	}
	
	flg=1;
	for(i=0; i<nnbr && flg; i++) {	// check if rtmd falls in r0-ri-rj triangle
		j = (i+1)%nnbr;
		if(checkface3D(tmdrc, p->rc, rnbr[i], rnbr[j], nrm[i], tmdrcnew)) { ii=i;	jj=j;	flg=0; }
	}
	
	if(flg) {	// do something or nothing if interpolation failed
		//printf("Wanring: tmd nbr not found in getTMDVesPos().\n");
		//printf("rtmd = (%.3g, %.3g, %.3g)\n", rtmd[0], rtmd[1], rtmd[2]);
		veccopy(p->rc, tmdrcnew);
	}
	
	vecadd(tmdrcnew, RCves, rtmdnew);
}





void MapTMD(CHAIN *CH)
{
	int i;
	//double dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i)
#endif

		for(i=0; i<(p->nseg); i++) {
			getTMDVesPos(p->rtmd[i], p->rtmd[i]);
			/*
			vecsub(p->rtmd[i], RCves, dr);	// in vesicle's frame
			mapsphere(Rv, dr);
			vecadd(dr, RCves, p->rtmd[i]);
			*/
		}
		
		p = p->next;
	}
}


void MotionTMD(CHAIN *CH)
{
	TMDmotion(CH);
	if(noise) TMDfluctuation(CH);
	MapTMD(CH);
}

