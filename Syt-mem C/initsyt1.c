// for Limit=0

void GetRndCoord(double nx[3], double ny[3], double nz[3])	// random directions
{
	int i;
	double r[3], r2, dotp;
	
	do {
		for(i=0; i<3; i++) r[i]=RND();	// 1st quadrant (x,y,z)>0
		r2=norm2(r);
	} while(r2>1 || r2==0);
	vecdiv(r, sqrt(r2), nx);
	
	do {
		for(i=0; i<3; i++) r[i]=RND();
		r2=norm2(r);
		if(r2>0) {
			vecdiv(r, sqrt(r2), r);
			/*
			if(vecdivx(r, sqrt(r2), r)==0) {
				printf("Error in GetRndCoord1b: r=(%.4g, %.4g, %.4g), x=%.4g\n", r[0], r[1], r[2], sqrt(r2));
				exit(0);
			}
			*/
			dotp=dotprod(r, nx);
		}
		else dotp=1;
	} while(r2>1 || r2==0 || fabs(dotp)==1);
	
	crossprod(nx, r, ny);
	normalize(ny, ny);
	
	crossprod(nx, ny, nz);
}



void GetRndFlatCoord(double nx[3], double ny[3], double nz[3])	// random directions with z pionting up
{
	int i;
	double r[3], r2, dotp;
	
	ny[0]=ny[1]=0;
	ny[2]=1;
	
	do {
		for(i=0; i<3; i++) r[i]=RND();	// 1st quadrant (x,y,z)>0
		r2=max2(norm2(r),eps);
		vecdiv(r, sqrt(r2), r);
		/*
		if(vecdivx(r, sqrt(r2), r)==0) {
			printf("Error in GetRndFlatCoord: r=(%.4g, %.4g, %.4g), x=%.4g\n", r[0], r[1], r[2], sqrt(r2));
			exit(0);
		}
		*/
		dotp=dotprod(r, ny);
	} while(r2>1 || r2==0 || fabs(dotp)==1);
	
	crossprod(ny, r, nx);
	
	crossprod(nx, ny, nz);
}



void GetRndSphCoord(double r[3], double nx[3], double ny[3], double nz[3])	// r is in vesicle's coord.
{
	int i;
	
	normalize(r, ny);	// normal direction
	
	do {
		for(i=0; i<3; i++) nz[i]=1-2*RND();
		crossprod(ny, nz, nx);
	} while(norm2(nx)==0);
	
	normalize(nx, nx);
	crossprod(nx, ny, nz);
}



void getLinkerpos(CHAIN *p)	// location of c2a-linker joint
{
	int i;
	double r, dr[3];
	
#if SYTCNT
	r=Rc2adb;
#else
	r=Rc2adb+Rc2b;
#endif
	
	for(i=0; i<(p->nseg); i++) {
		vecprod(p->ny[i], r, dr);
		vecadd(p->r[i], dr, p->rlinker[i]);
	}
}



void initTMDpos(CHAIN *p)	// map r to vesicle, run after p->rlinker is defined
{
	int i;
	double dr[3];
	
	for(i=0; i<(p->nseg); i++) {
		vecsub(p->rlinker[i], RCves, dr);	// dr in vesicle's center frame
		mapsphere(Rv, dr);
		vecadd(dr, RCves, p->rtmd[i]);
	}
}


// position of 4 beads in lab frame, run after p->r,nx,ny,nz are updated
void getSBeadSites(CHAIN *p)
{
	int i, j;
	
	for(i=0; i<(p->nseg); i++) {
		for(j=0; j<3; j++) {
			p->drC2AB[i][0][j] = RCc2a1[0]*(p->nx[i][j]) + RCc2a1[1]*(p->ny[i][j]) + RCc2a1[2]*(p->nz[i][j]);	// C2A1
			p->drC2AB[i][1][j] = RCc2a2[0]*(p->nx[i][j]) + RCc2a2[1]*(p->ny[i][j]) + RCc2a2[2]*(p->nz[i][j]);	// C2A2
			p->drC2AB[i][2][j] = RCc2b1[0]*(p->nx[i][j]) + RCc2b1[1]*(p->ny[i][j]) + RCc2b1[2]*(p->nz[i][j]);	// C2B1
			p->drC2AB[i][3][j] = RCc2b2[0]*(p->nx[i][j]) + RCc2b2[1]*(p->ny[i][j]) + RCc2b2[2]*(p->nz[i][j]);	// C2B2
			
			p->rC2AB[i][0][j] = p->r[i][j] + p->drC2AB[i][0][j];	// C2A1
			p->rC2AB[i][1][j] = p->r[i][j] + p->drC2AB[i][1][j];	// C2A2
			p->rC2AB[i][2][j] = p->r[i][j] + p->drC2AB[i][2][j];	// C2B1
			p->rC2AB[i][3][j] = p->r[i][j] + p->drC2AB[i][3][j];	// C2B2
		}
	}
}



void GetSSsites(CHAIN *p)	// syt-syt binding sites, run after p->r, p->nx and p->ny are updated
{
	int i;
	double dr[3], dr1[3], dr2[3];
	
	for(i=0; i<(p->nseg); i++) {
		veczero(p->dr1[i]);	// left stretching
		veczero(p->dr2[i]);	// right stretching
		
		vecprod(p->nx[i], -Lc2bhalf, dr1);	// left displacement
		vecprod(p->nx[i], Lc2bhalf, dr2);	// right displacement
		
	#if SYTCNT
		vecprod(p->ny[i], -Rc2b, dr);	// y-displacement for C2B
	#else
		veczero(dr);
	#endif
	
		vecadd(dr1, dr, dr1);
		vecadd(dr2, dr, dr2);
		
		vecadd(p->r[i], dr1, p->rc1[i]);	// left binding site
		vecadd(p->r[i], dr2, p->rc2[i]);	// right binding site
	}
}



void GetSMsites(CHAIN *p)	// after p->r, p->nx and p->ny are updated
{
	int i, j;
	
	for(i=0; i<(p->nseg); i++) {
		for(j=0; j<3; j++) {
			p->rsm[i][j] = p->r[i][j] + RCsm[0]*(p->nx[i][j]) + RCsm[1]*(p->ny[i][j]) + RCsm[2]*(p->nz[i][j]);
		}
	}
}



void UpdateSytSites(CHAIN *p)	// given p->r, update the rest positions excluding TMD position
{
	getSBeadSites(p);
	GetSSsites(p);
	GetSMsites(p);
	
	if(Vesicle) getLinkerpos(p);
}


void UpdateSytSitesAll(CHAIN *CH)
{
	CHAIN *p;
	
	p = CH;
	while(p) {
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		p = p->next;
	}
}


void NewChainProp(CHAIN *p)	// run after Nchain++
{
	double r[3], nr[3], a, rd;
	
	p->id = Nchain;
	p->nseg = 1;
#if CLOSED
	p->closed = 1;
#else
	p->closed = 0;
#endif
	p->Pclose = 0;
	
	if(NC0==1) {
		//veczero(r);
		r[0] = 0;
		r[1] = 0;	// 1/K0;
	}
	else {
		r[0] = 0.6*Lxhalf*(1-2*RND());
		r[1] = 0.6*Lyhalf*(1-2*RND());
	}
	
	if(Membrane && INIMEM==2 && NC0==1 && CLOSED) {
		readcontour();
		#if SYTCNT
			rd = Nsytperring*Lc2b/Pidb-Rc2b;
		#else
			rd = Nsytperring*Lc2b/Pidb;
		#endif
		r[2] = getcontour(rd)+Rc2b;
	}
	else r[2] = Rc2b;
	
	//r[2] = DLThalf;	//4*DLT;
	
	veccopy(r, p->r[0]);	// position of the 1st subunit
	
	veccopy(r, p->rhead);
	veccopy(r, p->rtail);
	veccopy(r, p->rcenter);
	p->rmax = p->rave = DLThalf;
	
	a = Pidb*RND();
	//a = 0;
	nr[0] = cos(a);	// vector from 0 to 1
	nr[1] = sin(a);
	nr[2] = 0;
	
	veccopy(nr, p->nx[0]);
	veccopy(nr, p->nnxt[0]);
	
	getnromalvector1(p->nx[0], p->ny[0]);	// ny is in the x-y plane (interaction surface facing side)
	//getnromalvector2(p->nx[0], p->ny[0]);	// ny is in the z-direction (interaction surface facing down)
	
	crossprod(p->nx[0], p->ny[0], p->nz[0]);
	
	UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
	initTMDpos(p);
	
	//p->next = NULL;
}



void AddNewSegments(CHAIN *p, int num)	// add (num) segments to *tail*, for init chains
{
	int i, n, nm1;
	double r[3], nx[3], ny[3], nz[3];
	
	if(num<=0) return;
	
	for(i=0; i<num; i++) {
		n = p->nseg;	// number of existing segments
		nm1 = n-1;	// tail seg-id
		
	#if CLOSED
		NewSegPosRing(p, nm1, 2, r, nx, ny, nz);	// requires p->(nx,ny,nz)[nm1], 2 means add to tail
	#else
		NewSegPos(p, nm1, 2, r, nx, ny, nz);	// requires p->(nx,ny,nz)[nm1], 2 means add to tail
		
		if(num>=Nsytperring0-1) {	// initial spiral
			if(Spiral==1) NewSegPos(p, nm1, 3, r, nx, ny, nz);	// planar spiral
			else if(Spiral==2) r[2]+=DLT/Nsytperring0;	// helical spiral
			
			/*
			if(Spiral==1 || Spiral==2) {
				vecsub(p->r[nm1], r, nx);	// adjust nx
				normalize(nx, nx);
			}
			
			if(Spiral==1) crossprod(nz, nx, ny);	// planar spiral: nz is unchanged, adjusting ny
			else if(Spiral==2) crossprod(nx, ny, nz);	// helical spiral: ny is unchanged, adjusting nz
			*/
		}
	#endif
	
		veccopy(r, p->rtail);
		veccopy(r, p->r[n]);	// new tail
		
		veccopy(nx, p->nx[n]);
		veccopy(ny, p->ny[n]);
		veccopy(nz, p->nz[n]);
		veccopy(p->nx[n], p->nnxt[n]);
		
		vecsub(p->r[n], p->r[nm1], r);
		normalize(r, p->nnxt[nm1]);
		
		(p->nseg)++;
	}
	
	
	UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
	initTMDpos(p);
	
	GetSytStretch(p);
}



CHAIN * AddNewChain(CHAIN *CH)	// append to chain, run after Nchain++, return the address of new chain
{
	CHAIN *p, *q;
	
	if(!CH) {
		p = (CHAIN *)malloc(sizeof(CHAIN));
		NewChainProp(p);	// run after Nchain++
		if(Limit==0) AddNewSegments(p, NS0-1);
		else AddNewSegments(p, NS1-1);
		p->next = NULL;
	}
	else {
		p = CH;
		while(p->next) p = p->next;	// p is the end of chain
		q = (CHAIN *)malloc(sizeof(CHAIN));
		NewChainProp(q);	// run after Nchain++
		if(Limit==0) AddNewSegments(q, NS0-1);
		else AddNewSegments(q, NS1-1);
		q->next = NULL;
		p->next = q;
		//q->prev = p;
	}
	
	return p;
}



void CheckChain(void)
{
	int i;
	CHAIN *p;
	
	printf("Chains:\n");
	p = chain;
	while(p) {
		printf("%d", p->nseg);
		//printf("%d: (%.2f, %.2f) -- (%.2f, %.2f)", p->id, p->rhead[0], p->rhead[1], p->rtail[0], p->rtail[1]);
		for(i=0; i<(p->nseg); i++) {
			printf("\t%d:\t", i);
			printf("%.2f, %.2f, %.2f\n", p->r[i][0], p->r[i][1], p->r[i][2]);
			//printf("%.2f, %.2f, %.2f\n", p->rtmd[i][0], p->rtmd[i][1], p->rtmd[i][2]);
			//printf("%.2f, %.2f, %.2f\n", p->nx[i][0], p->nx[i][1], p->nx[i][2]);
			//printf("%.2f, %.2f, %.2f\n", p->ny[i][0], p->ny[i][1], p->ny[i][2]);
			//printf("%.2f, %.2f, %.2f\n", p->nz[i][0], p->nz[i][1], p->nz[i][2]);
		}
		printf("\n");
		p = p->next;
	}
}



void getSphereUpPos(double r[3], double rn[3])
{
	double h, rinner, dr;

	/*
	do {
		r[0] = 1-2*RND();
		r[1] = 1-2*RND();
		r[2] = 0.2+0.8*RND();
	} while(norm2(r)>1);
	*/
	
	veczero(r);
	rn[0]=rn[1]=0;
	if(VesRing==1) rn[2]=1;	// ring on top of vesicle
	else rn[2]=-1;	// ring at bottom of vesicle
	
#if SYTCNT
	dr=Rc2bdb;
#else
	dr=Rc2b;
#endif
	
#if CLOSED
	rinner=Nsytperring*Lc2b/Pidb-dr;	// inner radius of closed ring
#else
	rinner=R0-dr;
#endif
	h=sqrt(max2(Rv*Rv-rinner*rinner,0));	// distance from vesicle center to ring plane
		
	if(VesRing==1) r[2]=Hv+Rv+h;	// ring on top
	else r[2]=Hv+Rv-h;	// ring at bottom
}



void RingOnVes(CHAIN *CH)
{
	int i, j;
	double r[3], rn[3], dr[3];
	CHAIN *p;
	
	UpdateChain(CH);	// get p->rcenter
	
	veczero(dr);
	
	p = CH;
	while(p) {
		if(VesRing!=0) {
			getSphereUpPos(r, rn);
			vecsub(r, p->rcenter, dr);
			vecadd(p->rhead, dr, p->rhead);
			vecadd(p->rtail, dr, p->rtail);
			vecadd(p->rcenter, dr, p->rcenter);
		}
		for(i=0; i<(p->nseg); i++) {
			vecadd(p->r[i], dr, p->r[i]);
			
			for(j=0; j<4; j++) vecadd(p->rC2AB[i][j], dr, p->rC2AB[i][j]);	// C2AB
			
			vecadd(p->rlinker[i], dr, p->rlinker[i]);
		
			vecadd(p->rc1[i], dr, p->rc1[i]);
			vecadd(p->rc2[i], dr, p->rc2[i]);
			vecadd(p->rsm[i], dr, p->rsm[i]);
		}
		
		initTMDpos(p);
		
		p = p->next;
	}
}



void InitSyt(void)
{
	int i;
	for(i=0; i<NC0; i++) {
		Nchain = i+1;
		chain = AddNewChain(chain);	// run after Nchain++
		
#if LGV_MC && BD_2ND
		chain_try = AddNewChain(chain_try);	// run after Nchain++
#endif
	}
	//CheckChain();
	
	if(Vesicle) RingOnVes(chain);
#if LGV_MC && BD_2ND
	if(Vesicle) RingOnVes(chain_try);
#endif
	
	//CheckChain();
	//exit(0);
}


void CleanSyt(void)
{
	CHAIN *p, *q;
	
	if(chain) {
		p = chain;
		while(p) {
			q = p->next;
			free(p);
			chain = p = q;
		}
		chain=NULL;
	}
	
#if LGV_MC && BD_2ND
	if(chain_try) {
		p = chain_try;
		while(p) {
			q = p->next;
			free(p);
			chain_try = p = q;
		}
		chain_try=NULL;
	}
#endif
}

