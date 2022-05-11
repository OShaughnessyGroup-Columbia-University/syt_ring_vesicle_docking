
double getFattrSytVesPP1(double r, double a0)	// point-point attraction, f>0 for attraction
{
	double r02, sq, f;
	
	r02=a0/Pi;
	sq=sqrt(r02+r*r);
	f=DebyeVesCoeff*Pidb*(exp(-r/LDebye)-r/sq*exp(-sq/LDebye));	// f>0 for attraction
	return f;
}




double getFattrSytVesPP(double r, double a)	// point-point attraction w/ cut-off distance, f>0 for attraction
{
	double rcut, invr, f;
	
	rcut = interpolate_func(r_cutoff_sv, a, amin_cutoff_sv, da_cutoff_sv, n_cutoff_smv);
	if(r <= rcut) f = 0;
	else {
		invr = 1.0/r;
		f = DebyeVesCoeff * invr * (invr + LDebyeinv) * exp(-r/LDebye);
	}
	
	return f;
}


double getFattrSytVesPF(double r, double cs, double a)	// point-face attraction
{
	double radius;
	double dr, r0, r1, r2;
	double elec, f;
	
	radius = sqrt(a/Pi);
	r0 = 1.0*radius;
	
	dr = 0.1*r0;
	r1 = r0 - dr/2;
	r2 = r0 + dr/2;
	
	if(r<=r1) elec = max2(0, getEfieldNear(r, cs, radius));
	else if(r>=r2) elec = max2(0, getEfieldFar(r, cs, radius));
	else elec = max2(0, getEfieldNear(r, cs, radius))*(1-(r-r1)/dr) + max2(0, getEfieldFar(r, cs, radius))*(r-r1)/dr;	// linear transition
	
	f = DebyeVesCoeff*elec;	// f>0 for attraction
	return f;
}



void SytVesAttract(VNODE *n, CHAIN *CH)	// for each face vnode
{
	int i;
	double dr[3], r2, r, nr[3], f, fv[3], taus[3], tauv[3];
	CHAIN *p;
	
#if !ATTRCT
	double cs;
#endif

	p = CH;	// be careful with p when using openmp !
	while(p) {	// for each chain
		for(i=0; i<(p->nseg); i++) {	// for each segment-i
			vecsub(n->r, p->rsm[i], dr);
			r2=norm2(dr);
			if(r2<RSMmax2) {
				r=sqrt(r2);
				vecdiv(dr, r, nr);	// nr pointing from syt site to membrane grid
				
#if ATTRCT	// point-point attraction
				f=getFattrSytVesPP(r, n->area);	// f>0 for attraction
#else	// point-face attraction
				cs=fabs(dotprod(n->e3, nr));
				f=getFattrSytVesPF(r, cs, n->area);	// f>0 for attraction
#endif

				vecprod(nr, f, fv);	// force on syt-site
				vecsub(n->f, fv, n->f);	// force on vesicle vnode
				vecsub(n->fext, fv, n->fext);	// external force on vesicle vnode
				
				vecsub(p->rsm[i], p->r[i], dr);	// arm for syt's torque
				crossprod(dr, fv, taus);	// torque on syt
				crossprod(fv, n->rc, tauv);	// torque on vesicle (force is -fv)
				
#if OPENMP
#pragma omp critical(SytVesAttract)	// be careful, use omp critical !!!
{
#endif
				vecadd(p->f[i], fv, p->f[i]);	// force on syt
				vecadd(p->fsm[i], fv, p->fsm[i]);
				vecadd(p->fatt[i], fv, p->fatt[i]);
				TrackFsyt += fv[2];
				TrackFves -= fv[2];
				//vecsub(Fves, fv, Fves);
				vecadd(p->taum[i], taus, p->taum[i]);	// torque on syt
				vecadd(Tauves, tauv, Tauves);	// torque on vesicle
#if OPENMP
}
#endif
			}
		}	// end of i
		p = p->next;
	}	// end of p
}



void SytVesRepulseN(VNODE *n, CHAIN *CH)	// for each vnode
{
	int i, j;
	double dr[3], r2, r, f, fv[3], tau[3];
	CHAIN *p;
	
	p = CH;	// be careful with p when using openmp !
	while(p) {	// for each chain
		for(i=0; i<(p->nseg); i++) {	// for each segment-i
			
			for(j=0; j<4; j++) {	// for each C2AB bead-j
				vecsub(p->rC2AB[i][j], n->r, dr);	// pointing from vesicle vnode to syt center
				
				if(dotprod(dr, n->e3)>0) {	// correct side
					r2 = norm2(dr);
					
					if(r2 < drSVrepulsion2) {
						r = max2(sqrt(r2), eps);
						f = KV*(drSVrepulsion-r);	// f>0
						vecprod(dr, f/r, fv);	// force on syt
						crossprod(fv, n->rc, tau);	// torque on vesicle
						
#if OPENMP 
#pragma omp critical(SytVesRepulseN)	// be careful, use omp critical !!!
{
#endif
						vecadd(p->f[i], fv, p->f[i]);	// force on syt
						vecadd(p->frep[i][j], fv, p->frep[i][j]);
						TrackFsyt += fv[2];
						TrackFves -= fv[2];
						//vecsub(Fves, fv, Fves);
						vecadd(Tauves, tau, Tauves);
#if OPENMP
}
#endif
						vecsub(n->f, fv, n->f);	// force on vesicle vnode
						vecsub(n->fext, fv, n->fext);
					
					}
				}
			}
		}
		p = p->next;
	}
}



void SytVesRepulseF(VFACE *pf, CHAIN *CH)	// for each vface
{
	int i, j, k, nid;
	double dr[3], r2, r, f, fv[3], fv_3[3], tau[3];
	CHAIN *p;
	VNODE *n;
	
	p = CH;	// be careful with p when using openmp !
	while(p) {	// for each chain
		for(i=0; i<(p->nseg); i++) {	// for each segment-i
			
			for(j=0; j<4; j++) {	// for each C2AB bead-j
				vecsub(p->rC2AB[i][j], pf->rclab, dr);	// pointing from vesicle's face center to syt center
				
				if(dotprod(dr, pf->nrm)>0) {	// correct side
					r2 = norm2(dr);
					
					if(r2 < drSVrepulsion2) {
						r = max2(sqrt(r2), eps);
						f = KV*(drSVrepulsion-r);	// f>0
						vecprod(dr, f/r, fv);	// fv=KV*dr is force on syt
					
#if OPENMP
#pragma omp critical(SytVesRepulseF)	// be careful, use omp critical !!!
{
#endif
						vecadd(p->f[i], fv, p->f[i]);	// force on syt
						vecadd(p->frep[i][j], fv, p->frep[i][j]);
						vecdiv(fv, 3.0, fv_3);	// force on each vnode of face-i
						for(k=0; k<3; k++) {
							nid = pf->inode[k];
							n = &vnode[nid];
							vecsub(n->f, fv_3, n->f);	// force on each vertex of vesicle face
							vecsub(n->fext, fv_3, n->fext);
							
							crossprod(fv_3, n->rc, tau);	// torque on vesicle
							vecadd(Tauves, tau, Tauves);
						}
						TrackFsyt += fv[2];
						TrackFves -= fv[2];
						//vecsub(Fves, fv, Fves);
#if OPENMP
}
#endif
					}
				}
			}	// end of j
		}	// end of i
		p = p->next;
	}
}



double getWLCforce(double d)	// stretching force on WLC
{
	double r, f;
	r=d/Llinker0;
	if(r<rlinkercutoff) f=Flinker0*(0.25/pow(1-r,2)+r-0.25);	// Marko-Siggia model (see Yongli's science paper)
	else f=Flinkermax;
	return f;
}



void TMDVesForce(VNODE *VN, double rtmd[3], double fv[3], int nbrid[3], double fvnew[3][3], double nrmtmd[3])
{	// interpolate membrane to grid, run after getNbrPositions() !
	int i, j, ni, nj, flg, nnbr;
	int id, ii, jj;
	double tmdrc[3], tmdrcnew[3], d2, d2min;
	double r01[3], r02[3];
	double rc[3], rnbr[NnbrMax][3], nrm[NnbrMax][3];
	double a1, a2, a3, atot, w1, w2, w3;
	double nr1[3], nr2[3], nr3[3];
	VNODE *p;
	
	vecsub(rtmd, RCves, tmdrc);	// coord. of TMD in vesicle's frame
	
	d2min = inf;
	id = -1;
	for(i=0; i<NVnode; i++) {	// find the nearest vnode on vesicle
		p = &VN[i];
		if(dotprod(p->rc, tmdrc)<0) continue;
		d2 = distance2(p->rc, tmdrc);
		if(d2<d2min) {
			d2min = d2;
			id = i;	// vnode-id is the nearest vnode to tmd
		}
	}
	if(id<0) {
		printf("Error in TMDVesForce()\n");
		printf("RCves = (%.3g, %.3g, %.3g)\n", RCves[0], RCves[1], RCves[2]);
		printf("rtmd = (%.3g, %.3g, %.3g)\n", rtmd[0], rtmd[1], rtmd[2]);
		exit(0);
	}
	
	p = &VN[id];	// p is the nearest vesicle vnode to tmd
	nnbr = p->nnbr;
	veccopy(p->rc, rc);	// center vnode in vesicle's frame
	
	for(i=0; i<nnbr; i++) {	// i-th neighbor
		j = (i+1)%nnbr;	// j is the (i+1)-th neighbor
		ni = p->nbr[i];	// id of i-th neighbor
		nj = p->nbr[j];	// id of (i+1)-th neighbor
		
		veccopy(VN[ni].rc, rnbr[i]);	// location of the i-th neighbor in RC's frame
		
		vecsub(VN[ni].rc, rc, r01);
		vecsub(VN[nj].rc, rc, r02);
		crossprod(r01, r02, nrm[i]);
		normalize(nrm[i], nrm[i]);	// nrm[i] is the inward/outward normal of triangle id-ni-nj
	}
	
	flg=1;
	ii=jj=0;
	for(i=0; i<nnbr && flg; i++) {	// check if rtmd falls inside r0-ri-rj triangle
		j = (i+1)%nnbr;
		if(checkface3D(tmdrc, rc, rnbr[i], rnbr[j], nrm[i], tmdrcnew)) {
			ii = p->nbr[i];
			jj = p->nbr[j];
			flg = 0;
		}
	}
	
	if(flg==0) {
		a1 = getarea3db(tmdrcnew, VN[ii].rc, VN[jj].rc);
		a2 = getarea3db(tmdrcnew, VN[jj].rc, p->rc);
		a3 = getarea3db(tmdrcnew, rc, VN[ii].rc);
		
		atot = 	a1+a2+a3;	// total area x2
		w1 = a1/atot;
		w2 = a2/atot;
		w3 = a3/atot;
		
		nbrid[0] = id;
		nbrid[1] = ii;
		nbrid[2] = jj;
		
		vecprod(fv, w1, fvnew[0]);
		vecprod(fv, w2, fvnew[1]);
		vecprod(fv, w3, fvnew[2]);
		
		vecprod(VN[nbrid[0]].e3, w1, nr1);	// weighted average
		vecprod(VN[nbrid[1]].e3, w2, nr2);
		vecprod(VN[nbrid[2]].e3, w3, nr3);
		vecadd3(nr1, nr2, nr3, nrmtmd);
		normalize(nrmtmd, nrmtmd);	// tmd normal
		//veccopy(VN[id].e3, nrmtmd);	// test !!!!!!!!!!!!!!!!
	}
	else {	// do something or nothing if interpolation failed
		for(i=0; i<3; i++) {
			nbrid[i] = id;
			vecdiv(fv, 3.0, fvnew[i]);	// all force goes to vnode-id (3x w/ 1/3 each)
		}
		veccopy(VN[id].e3, nrmtmd);
	}
}



void SytTmdVesAttract(VNODE *VN, CHAIN *p, int i)
{
	int j, nbrid[3];
	double dr[3], d, nr[3], f, fv[3], fvnew[3][3], tau[3];
	VNODE *q;
	
	vecsub(p->rtmd[i], p->rlinker[i], dr);
	d=max2(norm(dr),eps);
	vecdiv(dr, d, nr);	// vector from syt to tmd
	
	f=getWLCforce(d);
	vecprod(nr, f, fv);
	vecadd(p->f[i], fv, p->f[i]);	// tmd force on syt
	vecadd(p->fatt[i], fv, p->fatt[i]);
	vecsub(p->ftmd[i], fv, p->ftmd[i]);	// opposite force on tmd
	
	TMDVesForce(VN, p->rtmd[i], fv, nbrid, fvnew, p->nrmtmd[i]);	// distribute tmd force to nearby vesicle nodes
	
#if OPENMP
#pragma omp critical(SytTmdVesAttract)
{
#endif
	for(j=0; j<3; j++) {
		q = &VN[nbrid[j]];	// be careful, use omp critical !!!
		vecsub(q->f, fvnew[j], q->f);	// opposite force on vesicle
		vecsub(q->fext, fvnew[j], q->fext);
		//TrackFves -= fvnew[j][2];
		
		crossprod(fvnew[j], q->rc, tau);	// torque on vesicle
		vecadd(Tauves, tau, Tauves);
	}
	TrackFsyt += fv[2];
	TrackFves -= fv[2];
	//vecsub(Fves, fv, Fves);
#if OPENMP
}
#endif
	
	crossprod(p->ny[i], fv, tau);
	vecprod(tau, Lc2ahalf, tau);
	vecadd(p->tautmd[i], tau, p->tautmd[i]);	// torque
}



void FTSytVes(CHAIN *CH, VNODE *VN, VFACE *VF)
{
	int i;
	CHAIN *p;

#if OPENMP
#pragma omp parallel for private(i) schedule(dynamic,10)
#endif
	for(i=0; i<NVnode; i++) SytVesAttract(&VN[i], CH);


#if OPENMP
#pragma omp parallel for private(i) schedule(dynamic,10)
#endif
	for(i=0; i<NVnode; i++) SytVesRepulseN(&VN[i], CH);	// syt-vnode repulsion


#if OPENMP
#pragma omp parallel for private(i) schedule(dynamic,10)
#endif
	for(i=0; i<NVface; i++) SytVesRepulseF(&VF[i], CH);	// syt-face_center repulsion


	if(TMD) {
		p = CH;
		while(p) {
#if OPENMP
#pragma omp parallel for private(i)
#endif
			for(i=0; i<(p->nseg); i++) {
				SytTmdVesAttract(VN, p, i);	// syt-vesicle
			}
			p = p->next;
		}
	}
	
}

