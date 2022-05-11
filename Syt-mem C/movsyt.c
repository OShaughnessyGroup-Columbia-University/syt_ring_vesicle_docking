
void PrepSytMotion(CHAIN *CH)
{
	CHAIN *p;
	p = CH;
	while(p) {
		veccopy(p->rcenter, p->rcenterprev);
		p = p->next;
	}
}


void StretchMotion(CHAIN *CH)	// syt-syt stretching
{
	int i, n;
	double dr[3], dmax;
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n<=1) {
			p = p->next;
			continue;
		}
		for(i=0; i<n; i++) {
			vecadd(p->dr1[i], p->dr2[i], dr);	// dr is between i & i+1
			
			dmax = 0.5*norm(dr);
		#if FRICTN
			vecprod(dr, KSMudt, dr);	// converting deformation to displacement
		#else
			vecprod(dr, KSMu0dt, dr);
		#endif
			limitvector(dr, dmax);
			
			//vecprod(dr, 0.25, dr);
			
			vecadd(p->r[i], dr, p->r[i]);
		}
		p = p->next;
	}
}


void BendingMotion(CHAIN *CH)	// update syt's orientation (nx, ny, nz)
{
	int i, n;
	double tq, axis[3], da, sn, cs;
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		// don't ignore n=1 since there is still non-zero taum
		
		for(i=0; i<n; i++) {
			tq=norm(p->tauorient[i]);
			vecdiv(p->tauorient[i], max2(tq,eps), axis);
			/*
			if(vecdivx(p->tauorient[i], max2(tq,eps), axis)==0) {
				printf("Error in BendingMotion: r=(%.4g, %.4g, %.4g), x=%.4g\n", p->tauorient[i][0], p->tauorient[i][1], p->tauorient[i][2], max2(tq,eps));
				exit(0);
			}
			*/
		#if FRICTN
			da=Murotdt*tq;
		#else
			da=Murot0dt*tq;
		#endif
			da=min2(da, tq/KBsytmax);
			sn=sin(da);
			cs=cos(da);
			
			RotMatrix3(sn, cs, axis, p->nx[i], p->nx[i]);
			RotMatrix3(sn, cs, axis, p->ny[i], p->ny[i]);
			RotMatrix3(sn, cs, axis, p->nz[i], p->nz[i]);
		}
		p = p->next;
	}
}



void BendingMotion1D(CHAIN *CH)	// update syt's orientation (nx, ny, nz), only allowing bending about local z
{
	int i, n;
	double dotp, tauz[3];
	double tq, axis[3], da, sn, cs;
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		// don't ignore n=1 since there is still non-zero taum
		
		for(i=0; i<n; i++) {
			dotp=dotprod(p->nz[i], p->tauorient[i]);
			vecprod(p->nz[i], dotp, tauz);	// torque along local nz
			tq=norm(tauz);
			vecdiv(tauz, max2(tq,eps), axis);
			/*
			if(vecdivx(p->tauorient[i], max2(tq,eps), axis)==0) {
				printf("Error in BendingMotion: r=(%.4g, %.4g, %.4g), x=%.4g\n", p->tauorient[i][0], p->tauorient[i][1], p->tauorient[i][2], max2(tq,eps));
				exit(0);
			}
			*/
		#if FRICTN
			da=Murotdt*tq;
		#else
			da=Murot0dt*tq;
		#endif
			da=min2(da, tq/KBsytmax);
			sn=sin(da);
			cs=cos(da);
			
			RotMatrix3(sn, cs, axis, p->nx[i], p->nx[i]);
			RotMatrix3(sn, cs, axis, p->ny[i], p->ny[i]);
			//RotMatrix3(sn, cs, axis, p->nz[i], p->nz[i]);
		}
		p = p->next;
	}
}




void StericMotion(CHAIN *CH)
{
	int i;
	double dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		//vecprod(p->ftot, min2(Mu/(p->nseg), KVinv), drc);	// motion of the whole chain
		
		for(i=0; i<(p->nseg); i++) {	// motion of individual subunits
			//vecprod(p->f[i], KSeffinvmin, dr);
		#if FRICTN
			vecprod(p->f[i], Mudt, dr);	// new !!!!!!
		#else
			vecprod(p->f[i], Mu0dt, dr);	// new !!!!!!
		#endif
			limitvector(dr, DistMin);
			//vecadd(dr, drc, dr);	// net shift
			vecadd(p->r[i], dr, p->r[i]);
			//p->r[i][2] = max2(p->r[i][2], DLThalf);
		}
		p = p->next;
	}
}



void SytNoise(CHAIN *CH)
{
	int i, j;
	double drn, dan, dr[3], axis[3], da, sn, cs;
	CHAIN *p;
	
#if OPENMP
	int tid;
#endif

#if FRICTN
	drn=dRnoise;
	dan=dAnoise*sqrt(3);	// 3 degrees of rotational freedom
#else
	drn=dRnoise0;
	dan=dAnoise0*sqrt(3);
#endif
	
	//veczero(dr);
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i,j,dr,axis,da,sn,cs,tid)
#endif

		for(i=0; i<(p->nseg); i++) {	// motion of individual subunits
		#if OPENMP
			tid = omp_get_thread_num();
			for(j=0; j<3; j++) dr[j]=drn*(1-2*RNDMP(tid));
			getRndDirMP(tid, axis);
			da=dan*(1-2*RNDMP(tid));
		#else
			for(j=0; j<3; j++) dr[j]=drn*(1-2*RND());
			getRndDir(axis);
			da=dan*(1-2*RND());
		#endif
			
			vecadd(p->r[i], dr, p->r[i]);
			sn=sin(da);
			cs=cos(da);
			RotMatrix3(sn, cs, axis, p->nx[i], p->nx[i]);	// rotation about axis
			RotMatrix3(sn, cs, axis, p->ny[i], p->ny[i]);
			RotMatrix3(sn, cs, axis, p->nz[i], p->nz[i]);
		}
		
		p = p->next;
	}
}



void NormalizeChain(CHAIN *CH)
{
	int i;
	CHAIN *p;
	
	p = CH;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i)
#endif

		for(i=0; i<(p->nseg); i++) {
			normalize(p->nx[i], p->nx[i]);
			crossprod(p->nx[i], p->ny[i], p->nz[i]);	// use temp ny to get nz
			normalize(p->nz[i], p->nz[i]);
			crossprod(p->nz[i], p->nx[i], p->ny[i]);	// use nx, nz to get real ny
			normalize(p->ny[i], p->ny[i]);
		}
		
		p = p->next;
	}
}


void BoxChain(CHAIN *CH)
{
	int i, n, flg;
	double dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		flg=0;
		veczero(dr);
		if(fabs(p->rcenter[0]) > Lxhalfplus) {
			flg=1;
			if(p->rcenter[0] > 0) dr[0] = -Lx;
			else dr[0] = Lx;
		}
		if(fabs(p->rcenter[1]) > Lyhalfplus) {
			flg=1;
			if(p->rcenter[1] > 0) dr[1] = -Ly;
			else dr[1] = Ly;
		}
		
		if(flg==1) {	// need to move all segments in chain
			vecadd(p->rhead, dr, p->rhead);
			vecadd(p->rtail, dr, p->rtail);
			vecadd(p->rcenter, dr, p->rcenter);
			vecadd(p->rcenterprev, dr, p->rcenterprev);
			
			n = p->nseg;
			for(i=0; i<n; i++) vecadd(p->r[i], dr, p->r[i]);
			UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		}
		p = p->next;
	}
}



void CenterChain(CHAIN *CH)
{
	int i, n, nc;
	double dr[3];
	CHAIN *p;
	
	p = CH;	// update p->rcenter
	while(p) {
		n = p->nseg;
		
		veczero(dr);
		
		if(n>Nsytperring0 && !CLOSED) nc=Nsytperring0;	// get rcenter from the first nc syts
		else nc=max2(n,1);
			
		for(i=0; i<nc; i++) vecadd(dr, p->r[i], dr);
		vecdiv(dr, -1.0*nc, dr);	// dr=-rcenter
		dr[2] = 0;	// not changing the z value
		
		vecadd(p->rcenter, dr, p->rcenter);
		vecadd(p->rcenterprev, dr, p->rcenterprev);
		vecadd(p->rhead, dr, p->rhead);
		vecadd(p->rtail, dr, p->rtail);
		
		for(i=0; i<n; i++) vecadd(p->r[i], dr, p->r[i]);
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		p = p->next;
	}
}



void CheckSyt(CHAIN *CH)
{
	int i;
	double r;
	CHAIN *p;
	
	p = CH;
	while(p) {
		for(i=0; i < p->nseg; i++) {
			r=p->r[i][0]+p->r[i][1]+p->r[i][2];
			if(r<-1e3 || r>1e3) {	// problem
				printf("problem for i=%d\n", i);
				printf("r=(%.3g, %.3g, %.3g)\n", p->r[i][0], p->r[i][1], p->r[i][2]);
				printf("nx=(%.3g, %.3g, %.3g)\n", p->nx[i][0], p->nx[i][1], p->nx[i][2]);
				printf("ny=(%.3g, %.3g, %.3g)\n", p->ny[i][0], p->ny[i][1], p->ny[i][2]);
				printf("nz=(%.3g, %.3g, %.3g)\n", p->nz[i][0], p->nz[i][1], p->nz[i][2]);
				printf("f=(%.3g, %.3g, %.3g)\n", p->f[i][0], p->f[i][1], p->f[i][2]);
				printf("fb=(%.3g, %.3g, %.3g)\n", p->fb[i][0], p->fb[i][1], p->fb[i][2]);
				printf("taucc=(%.3g, %.3g, %.3g)\n", p->taust[i][0], p->taust[i][1], p->taust[i][2]);
				printf("taubd=(%.3g, %.3g, %.3g)\n", p->taubd[i][0], p->taubd[i][1], p->taubd[i][2]);
				printf("tauts=(%.3g, %.3g, %.3g)\n", p->tauts[i][0], p->tauts[i][1], p->tauts[i][2]);
				printf("tauorient=(%.3g, %.3g, %.3g)\n", p->tauorient[i][0], p->tauorient[i][1], p->tauorient[i][2]);
				printf("taum=(%.3g, %.3g, %.3g)\n", p->taum[i][0], p->taum[i][1], p->taum[i][2]);
				exit(0);
			}
		}
		p = p->next;
	}
}




void ShiftChain(CHAIN *CH)
{
	int i, n;
	double mudt_n, drtot[3], drexist[3], dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n>0) {	// include monomers
		#if FRICTN
			mudt_n = Mudt/n;
		#else
			mudt_n = Mu0dt/n;
		#endif
			vecprod(p->ftot, mudt_n, drtot);	// total force for each chain
			vecsub(p->rcenter, p->rcenterprev, drexist);	// existing displacement from forces
			vecsub(drtot, drexist, dr);	// net displacement to be adjusted
			limitvector(dr, drmemmin);
			for(i=0; i<n; i++) vecadd(p->r[i], dr, p->r[i]);
			UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		}
		p = p->next;
	}
}




void ShiftFixedChain(CHAIN *CH)	// no rcenterprev
{
	int i, n;
	double mudt_n, dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n>0) {	// include monomers
		#if FRICTN
			mudt_n = Mudt/n;
		#else
			mudt_n = Mu0dt/n;
		#endif
			vecprod(p->ftot, mudt_n, dr);	// total force for each chain
			limitvector(dr, drmemmin);
			for(i=0; i<n; i++) vecadd(p->r[i], dr, p->r[i]);
			UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		}
		p = p->next;
	}
}



void ShiftFixedChainZ(CHAIN *CH)	// not shifting x & y, no rcenterprev
{
	int i, n;
	double mudt_n, dr[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n>0) {	// include monomers
		#if FRICTN
			mudt_n = Mudt/n;
		#else
			mudt_n = Mu0dt/n;
		#endif
			vecprod(p->ftot, mudt_n, dr);
			limitvector(dr, drmemmin);
			dr[0] = dr[1] = 0;	// no x-/y- movement
			for(i=0; i<n; i++) vecadd(p->r[i], dr, p->r[i]);
			UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		}
		p = p->next;
	}
}



void getSytChainMu(CHAIN *p, double tau[3])	// rotational mobility of syt chain, see notes
{
	int i;
	double e[3], dr2, dre, dre2, sum;
	
	if(p->nseg==0) p->MuRotMemdt=0;
	else {
		normalize(tau, e);
		sum = 0;
		for(i=0; i<(p->nseg); i++) {	// tau = zeta1 * sum[r x (w x r)]
			dr2 = norm2(p->rc[i]);	// r x (w x r) = w ( r . r) - r (w . r)
			dre = dotprod(p->rc[i], e);
			dre2 = dre*dre;
			sum += dr2 - dre2;
		}
	#if FRICTN
		p->MuRotMemdt = Mudt/max2(sum,eps);
	#else
		p->MuRotMemdt = Mu0dt/max2(sum,eps);
	#endif
		p->Drdt = (p->MuRotMemdt)*kBT;
	}
}



void rotatesytring(double rc[3], double r[3], double sth, double cth, double axis[3])
{	// rotate r by th about axis which is through rc
	double dr[3];
	
	vecsub(r, rc, dr);
	RotMatrix3(sth, cth, axis, dr, dr);
	vecadd(rc, dr, r);
}



void RotateChain(CHAIN *CH)
{
	int i, j, n;
	double dr[3];
	double tau[3], tautot[3], omg[3], dth, axis[3];
	double sth, cth;
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n>0) {	// include monomers
			veczero(tautot);
			for(i=0; i<n; i++) {
				vecsub(p->rsm[i], p->rcenter, dr);	// torque from attractive force
				crossprod(dr, p->fsm[i], tau);
				vecadd(tautot, tau, tautot);
				
				for(j=0; j<4; j++) {
					crossprod(p->rc[i], p->frep[i][j], tau);
					vecadd(tautot, tau, tautot);
				}
			}
			
		#if NOSPIN
			tautot[2]=0;	// prevent spin about z
		#endif
			
			getSytChainMu(p, tautot);	// get p->MuRotMemdt
			
			vecprod(tautot, p->MuRotMemdt, omg);
			dth = norm(omg);
			dth = min2(dth, 0.01);	// limit angle !!!!!!!!!!!!!
			sth = sin(dth);
			cth = cos(dth);
			normalize(omg, axis);
			
			rotatesytring(p->rcenter, p->rhead, sth, cth, axis);	// rotation in center frame
			rotatesytring(p->rcenter, p->rtail, sth, cth, axis);
			for(i=0; i<n; i++) {
				// rotation about ring center
				rotatesytring(p->rcenter, p->r[i], sth, cth, axis);
				// rotation about itself
				RotMatrix3(sth, cth, axis, p->nx[i], p->nx[i]);
				RotMatrix3(sth, cth, axis, p->ny[i], p->ny[i]);
				RotMatrix3(sth, cth, axis, p->nz[i], p->nz[i]);
				RotMatrix3(sth, cth, axis, p->nnxt[i], p->nnxt[i]);
				RotMatrix3(sth, cth, axis, p->nxa[i], p->nxa[i]);
				RotMatrix3(sth, cth, axis, p->nya[i], p->nya[i]);
				RotMatrix3(sth, cth, axis, p->nza[i], p->nza[i]);
			}
		}
		
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		
		p = p->next;
	}
}


void SytChainNoise(CHAIN *CH)
{
	int i;
	double drn, dr[3], dan, da1, da2, da3;
	double sn1, sn2, sn3, cs1, cs2, cs3;
	CHAIN *p;
	
	veczero(dr);
	p = CH;
	while(p) {
	#if FRICTN
		drn=dRnoise/max2(p->nseg,1);
	#else
		drn=dRnoise0/max2(p->nseg,1);
	#endif
		dr[0]=drn*(1-2*RND());	// thermal fluctuation of the entire chain
		dr[1]=drn*(1-2*RND());
		dr[2]=drn*(1-2*RND());
		vecadd(p->rhead, dr, p->rhead);
		vecadd(p->rtail, dr, p->rtail);
		vecadd(p->rcenter, dr, p->rcenter);
		for(i=0; i<(p->nseg); i++) vecadd(p->r[i], dr, p->r[i]);	// translational motion
		
		dan=sqrt(6*(p->Drdt));
		da1=dan*(1-2*RND());
		da2=dan*(1-2*RND());
		da3=dan*(1-2*RND());
		sn1=sin(da1);	cs1=cos(da1);
		sn2=sin(da2);	cs2=cos(da2);
		sn3=sin(da3);	cs3=cos(da3);
		RotMatrix3(sn1, cs1, p->nx[i], p->ny[i], p->ny[i]);	// along x-axis
		RotMatrix3(sn1, cs1, p->nx[i], p->nz[i], p->nz[i]);
		RotMatrix3(sn2, cs2, p->ny[i], p->nx[i], p->nx[i]);	// along y-axis
		RotMatrix3(sn2, cs2, p->ny[i], p->nz[i], p->nz[i]);
		RotMatrix3(sn3, cs3, p->nz[i], p->nx[i], p->nx[i]);	// along z-axis
		RotMatrix3(sn3, cs3, p->nz[i], p->ny[i], p->ny[i]);
		
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		
		p = p->next;
	}
}




void MotionSyt(CHAIN *CH)
{
	PrepSytMotion(CH);
	
#if FIXSYT
	if(Vesicle) ShiftFixedChain(CH);
	else ShiftFixedChainZ(CH);
		
	#if !NO_ROT
		RotateChain(CH);
	#endif
	
	if(FIXSYT==1) {	// allow stretching and bending about local z
		UpdateSytStretchAll(CH);	// update dr1, dr2 (run this after UpdateSytSites())
		StretchMotion(CH);
		if(LP>0) BendingMotion1D(CH);
		StericMotion(CH);
	}
	if(noise) SytChainNoise(CH);
#else
	StretchMotion(CH);
	if(LP>0) BendingMotion(CH);	// update syt's orientation (nx, ny, nz)
	StericMotion(CH);
	ShiftChain(CH);	// shift entire chain
	if(noise) SytNoise(CH);
#endif

	NormalizeChain(CH);
	
#if OPENGL
	BoxChain(CH);
#endif
	
	//CheckSyt();
	if(Limit==0 && Vesicle==0 && Nchain==1) CenterChain(CH);	// test !!!!!!!!

	//UpdateChain(CH);	// update Nclosed, Nopen, rcenter, rmax, rave after motion
	UpdateSytSitesAll(CH);	// update syt-syt sites and syt-mnode sites
	UpdateSytStretchAll(CH);	// update dr1, dr2 (run this after UpdateSytSites())
}


