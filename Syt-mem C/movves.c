
void UpdateVesFaceCoord(VNODE *VN, VFACE *VF)	// copy vnode.r to vface.r
{
	int i, j;
	VFACE *f;
	
#if OPENMP
#pragma omp parallel for private(i,j,f)
#endif

	for(i=0; i<NVface; i++) {
		f = &VF[i];
		j = f->inode[0];
		veccopy(VN[j].r, f->r1);
		j = f->inode[1];
		veccopy(VN[j].r, f->r2);
		j = f->inode[2];
		veccopy(VN[j].r, f->r3);
		vecave3(f->r1, f->r2, f->r3, f->rclab);	// rc is the coord. of vface center in lab frame
		vecsub(f->rclab, RCves, f->rcrc);
	}
	
}



void getRCves(VNODE *VN)
{
	int i;
	double dr[3];
	
	//CalcVesAreaVolume();	// get vnode.area
	
	veczero(RCves);
	for(i=0; i<NVnode; i++) {
		vecprod(VN[i].r, VN[i].area, dr);	// area-weighted center
		vecadd(RCves, dr, RCves);
	}
	vecdiv(RCves, AreaVes, RCves);	// area-weighted center
}



void ShiftSytWithVes(CHAIN *CH, double dr[3])	// shift SYT+TMD together with vesicle by dr
{
	int i, n, nm1;
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		nm1 = n-1;
		
#if OPENMP
#pragma omp parallel for private(i)
#endif
		
		for(i=0; i<(p->nseg); i++) {
			vecadd(p->r[i], dr, p->r[i]);	// SYT
			vecadd(p->rtmd[i], dr, p->rtmd[i]);	// TMD
		}
		
		veccopy(p->r[0], p->rhead);
		veccopy(p->r[nm1], p->rtail);
		vecadd(p->rcenter, dr, p->rcenter);
		
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		
		p = p->next;
	}
}



void ShiftTMDWithVes(CHAIN *CH, double dr[3])	// shift TMD together with vesicle by dr
{
	int i;
	CHAIN *p;
	
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i)
#endif

		for(i=0; i<(p->nseg); i++) vecadd(p->rtmd[i], dr, p->rtmd[i]);	// TMD
		
		p = p->next;
	}
}



void RotateSytWithVes(CHAIN *CH, double axis[3], double sn, double cs)	// rotate SYT+TMD together with vesicle
{
	int i, n, nm1, nc;
	double dr[3];
	CHAIN *p;
	
	p = CH;
	
	while(p) {
		n = p->nseg;
		nm1 = n-1;
		
#if OPENMP
#pragma omp parallel for private(i,dr)
#endif
		for(i=0; i<n; i++) {
			// SYT
			vecsub(p->r[i], RCves, dr);
			RotMatrix3(sn, cs, axis, dr, dr);
			vecadd(RCves, dr, p->r[i]);
			
			RotMatrix3(sn, cs, axis, p->nx[i], p->nx[i]);
			RotMatrix3(sn, cs, axis, p->ny[i], p->ny[i]);
			RotMatrix3(sn, cs, axis, p->nz[i], p->nz[i]);
			RotMatrix3(sn, cs, axis, p->nnxt[i], p->nnxt[i]);
			RotMatrix3(sn, cs, axis, p->nxa[i], p->nxa[i]);
			RotMatrix3(sn, cs, axis, p->nya[i], p->nya[i]);
			RotMatrix3(sn, cs, axis, p->nza[i], p->nza[i]);
			
			// TMD
			vecsub(p->rtmd[i], RCves, dr);	// dr is rtmd[i] in vesicle's frame
			RotMatrix3(sn, cs, axis, dr, dr);
			vecadd(dr, RCves, p->rtmd[i]);	// back to lab frame
		}
		
		veccopy(p->r[0], p->rhead);
		veccopy(p->r[nm1], p->rtail);
		
		if(n>Nsytperring0 && !CLOSED) nc=Nsytperring0;	// get rcenter from the first nc syts
		else nc=max2(n,1);
		
		veczero(p->rcenter);
		for(i=0; i<nc; i++) vecadd(p->rcenter, p->r[i], p->rcenter);
		vecdiv(p->rcenter, 1.0*nc, p->rcenter);
		
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		
		p = p->next;
	}
}


void RotateTMDWithVes(CHAIN *CH, double axis[3], double sn, double cs)	// rotate TMD together with vesicle
{
	int i;
	double dr[3];
	CHAIN *p;
	
	p = CH;
	
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i,dr)
#endif
		for(i=0; i<(p->nseg); i++) {	// TMD
			vecsub(p->rtmd[i], RCves, dr);	// dr is rtmd[i] in vesicle's frame
			RotMatrix3(sn, cs, axis, dr, dr);
			vecadd(dr, RCves, p->rtmd[i]);	// back to lab frame
		}
		
		p = p->next;
	}
}



void VesTranslation(VNODE *VN, VFACE *VF, CHAIN *CH)
{
	int i;
	double df, dfn[3], dfp[3], drn[3], drp[3], dr[3];
	
#if OPENMP
#pragma omp parallel for private(i,df,dfn,dfp,drn,drp,dr)
#endif

	for(i=0; i<NVnode; i++) {
		//vecprod(VN[i].f, MuVesTdt, dr);
		df=dotprod(VN[i].f, VN[i].e3);
		vecprod(VN[i].e3, df, dfn);	// normal force
		vecsub(VN[i].f, dfn, dfp);	// parallel force
		
		vecprod(dfn, MuvesPerpdt, drn);	// normal displacement
		vecprod(dfp, MuvesParadt, drp);	// parallel displacement
		vecadd(drn, drp, dr);	// total displacement
		
		limitvector(dr, dRveseps);
		vecadd(VN[i].r, dr, VN[i].r);
	}
	
	veccopy(RCves, RCvesprev);
	getRCves(VN);
	vecsub(RCves, RCvesprev, dRCves);	// net change
	//*/
	
	for(i=0; i<NVnode; i++) vecsub(VN[i].r, RCves, VN[i].rc);	// update each vnode's rc for remapping
	for(i=0; i<NVface; i++) vecsub(VF[i].rclab, RCves, VF[i].rcrc);
	
	if(TMD) ShiftTMDWithVes(CH, dr);	// OPENMP inside
}



void VesRotation(VNODE *VN, CHAIN *CH)
{
	int i;
	double omg[3], da, cs, sn, axis[3];
	VNODE *n;
	
	vecprod(Tauves, MuVesTotRdt, omg);
	da = norm(omg);
	da = min2(da, 0.1);	// limit angle
	cs = cos(da);
	sn = sin(da);
	vecdiv(omg, max2(da,eps), axis);
	
	// rotate vesicle
	
#if OPENMP
#pragma omp parallel for private(i,n)
#endif

	for(i=0; i<NVnode; i++) {
		n = &VN[i];
		RotMatrix3(sn, cs, axis, n->rc, n->rc);	// rotate rc in vesicle's frame
		vecadd(n->rc, RCves, n->r);	// back to lab frame
	}
	
	if(TMD) RotateTMDWithVes(CH, axis, sn, cs);	// OPENMP inside
}



void VesNoise(VNODE *VN, CHAIN *CH)
{
	int i;
	double dr[3], axis[3];
	double da, sn, cs;
	VNODE *n;
	
	dr[0]=dr[1]=0;
	dr[2]=dVRnoise*(1-2*RND());	// only diffuse in z
	vecadd(RCves, dr, RCves);	// shift ves center
	
	getRndDir(axis);	// don't put this in OPENMP!
	da=dVAnoise*sqrt(3)*(1-2*RND());	// total angle in 3 directions (sqrt(3))
	sn=sin(da);
	cs=cos(da);
	
#if OPENMP
#pragma omp parallel for private(i,n)
#endif
	
	for(i=0; i<NVnode; i++) {
		n = &VN[i];
		RotMatrix3(sn, cs, axis, n->rc, n->rc);	// rc's rotation (rc doesn't shift)
		vecadd(n->rc, RCves, n->r);	// total
	}
	
	/*
	ShiftSytWithVes(CH, dr);	// OPENMP inside
	RotateSytWithVes(CH, axis, sn, cs);	// OPENMP inside
	*/
	if(TMD) {
		ShiftTMDWithVes(CH, dr);	// OPENMP inside
		RotateTMDWithVes(CH, axis, sn, cs);	// OPENMP inside
	}
}




void VesNodeTether(VNODE *VN)	// tethering potential
{
	int i, j, nnbr;
	double dr[3], d2, d, dd;
	VNODE *p, *q;
	
#if OPENMP
#pragma omp parallel for private(i,j,nnbr,dr,d2,d,dd,p,q)
#endif

	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		nnbr = p->nnbr;
		for(j=0; j<nnbr; j++) {
			q = &VN[p->nbr[j]];
			if(i > q->id) {	// to avoid double-counting
				vecsub(q->r, p->r, dr);	// dr is from p to q
				d2 = norm2(dr);
				// repulsion
				///*
				if(d2 < LBvesmin2) {	// p & q are in contact
					d = max2(sqrt(d2), eps);
					dd = (LBvesmin - d)/2.0;	// dd>0
					vecprod(dr, dd/d, dr);	// displacement for p
					vecsub(p->r, dr, p->r);	// push back
#if OPENMP	// be careful, use omp critical !!!
#pragma omp critical(VesNodeTether)	// share the same critical name as below
{
#endif
					vecadd(q->r, dr, q->r);
#if OPENMP
}
#endif
				}
				//*/
				/*
				d = max2(sqrt(d2), drmemeps);
				dd = drmemeps*exp(-10*d/drmem);
				vecprod(dr, dd/d, dr);	// displacement for p
				vecsub(p->r, dr, p->r);	// push back
				vecadd(q->r, dr, q->r);
				*/
				
				// attraction
				if(d2>LBvesmax2) {
					d = sqrt(d2);
					dd = (d - LBvesmax)/2.0;	// dd>0
					vecprod(dr, dd/d, dr);	// displacement for p
					vecadd(p->r, dr, p->r);
#if OPENMP	// be careful, use omp critical !!!
#pragma omp critical(VesNodeTether)	// share the same critical name as below
{
#endif
					vecsub(q->r, dr, q->r);
#if OPENMP
}
#endif
				}
			}
		}
	}
}






void VesShift(VNODE *VN, MNODE **MN, CHAIN *CH, int state)
{
	int i;
	double zmax, dr[3];
	VNODE *n;
	
	if(Membrane) GetMemMaxMin(MN);	// get Zmax
	else Zmax = 0;
	
	zmax = Zmax + Rv + Hvmax;
	
	if(state==1) {	// center x & y, limit z within tethering distance
		dr[0]=-RCves[0];
		dr[1]=-RCves[1];
		if(RCves[2] < zmax) dr[2] = 0;
		else dr[2] = zmax-RCves[2];	// too far away
	}
	else if(state==2) {	// center x & y, shift z according to force
		dr[0]=-RCves[0];	// center x & y
		dr[1]=-RCves[1];
		if(RCves[2] < zmax) {
			dr[2] = MuVesTotTdt*Fves[2];	// upper bound
			dr[2] -= dRCves[2];	// subtract the existing displacement from VesTranslation()
			dr[2] = sign(dr[2])*min2(fabs(dr[2]), drmemmin);
		}
		else dr[2] = zmax-RCves[2];	// too far away
	}
	else if(state==3) {	// shift x, y and z according to force
		vecprod(Fves, MuVesTotTdt, dr);
		if(RCves[2] > zmax) dr[2] = zmax - RCves[2];	// too far away
	}
	else if(state==4) {	// fix x and y when Membrane=0
		dr[0]=-RCves[0];
		dr[1]=-RCves[1];
		dr[2]=Fves[2]*MuVesTotTdt;
	}
	else if(state==5) {	// fix x, y and z when Membrane=0
		dr[0]=-RCves[0];
		dr[1]=-RCves[1];
		dr[2]=Rv+Hv-RCves[2];
	}
	else veczero(dr);	// not doing anything
	
	vecadd(RCves, dr, RCves);
	
#if OPENMP
#pragma omp parallel for private(i,n)
#endif

	for(i=0; i<NVnode; i++) {
		n = &VN[i];
		vecadd(n->r, dr, n->r);	// n->rc is not changing!
	}
	
	//ShiftSytWithVes(dr);    // OPENMP inside
	if(TMD) ShiftTMDWithVes(CH, dr);
}




void VesSytRotate(VNODE *VN, CHAIN *CH)	// run this when Membrane=0
{
	int i;
	double tau[3], omg[3];
	double da, cs, sn, axis[3];
	CHAIN *p;
	VNODE *n;
	
	veccopy(Tauves, tau);	// torque on vesicle
	p = CH;
	while(p) {
		vecadd(tau, p->tautot, tau);	// plus torque on all syts
		p = p->next;
	}
	
	vecprod(tau, -MuVesTotRdt, omg);	// cancel the net torque
	da = norm(omg);
	//da = min2(da, 0.1);	// limit angle
	cs = cos(da);
	sn = sin(da);
	vecdiv(omg, max2(da,eps), axis);
	
	// rotate vesicle
	
#if OPENMP
#pragma omp parallel for private(i,n)
#endif

	for(i=0; i<NVnode; i++) {
		n = &VN[i];
		RotMatrix3(sn, cs, axis, n->rc, n->rc);	// rotate rc in vesicle's frame
		vecadd(n->rc, RCves, n->r);	// back to lab frame
	}
	
	// rotate both TMD & SYT with it
	
	if(TMD) RotateTMDWithVes(CH, axis, sn, cs);	// OPENMP inside
	RotateSytWithVes(CH, axis, sn, cs);	// OPENMP inside
}




void MotionVes(VNODE *VN, VFACE *VF, MNODE **MN, CHAIN *CH)
{
	VesTranslation(VN, VF, CH);
	if(Membrane) VesRotation(VN, CH);
	
	GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
	VesNodeTether(VN);
	
#if FLPBND
	if(Cntflip_ves++ >= Nflip_ves) {
		GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
		UpdateVesFaceCoord(VN, VF);	// update face coord. don't forget !!!!!!!
		if(noise) VesBondFlip(CH,VN,VF,1);
		else VesBondFlip(CH,VN,VF,0);
		Cntflip_ves = 0;
	}
#endif
	
	GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
	if(Membrane) {
		VesShift(VN, MN, CH, 2);	// fixing (x,y), shifting z according to force
		//VesShift(VN, MN, CH, 3);	// shifting (x,y,z) according to force
	}
	else {
		VesShift(VN, MN, CH, 5);	// fixing (x,y,z)
		//VesSytRotate(VN, CH);
	}
	
	UpdateVesFaceCoord(VN, VF);	// after moving nodes, update face coord.

	if(noise) {
		VesNoise(VN, CH);
		GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
	}
	
	UpdateVesFaceCoord(VN, VF);	// after moving nodes, update face coord.
}


