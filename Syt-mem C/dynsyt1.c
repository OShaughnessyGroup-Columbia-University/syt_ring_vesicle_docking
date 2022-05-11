// for Limit=0

void CopySegments(CHAIN *q, int j, CHAIN *p, int i)	// copy q-j to p-i
{
	int k;
	
	veccopy(q->r[j], p->r[i]);
	veccopy(q->rc[j], p->rc[i]);
	
	for(k=0; k<4; k++) {
		veccopy(q->rC2AB[j][k], p->rC2AB[i][k]);
		veccopy(q->drC2AB[j][k], p->drC2AB[i][k]);
	}
	
	veccopy(q->rtmd[j], p->rtmd[i]);
	veccopy(q->rlinker[j], p->rlinker[i]);
	veccopy(q->rc1[j], p->rc1[i]);
	veccopy(q->rc2[j], p->rc2[i]);
	
	veccopy(q->dr1[j], p->dr1[i]);
	veccopy(q->dr2[j], p->dr2[i]);
	
	veccopy(q->rsm[j], p->rsm[i]);
	
	p->E[i] = q->E[j];
	p->Poff[i] = q->Poff[j];
	veccopy(q->nx[j], p->nx[i]);
	veccopy(q->ny[j], p->ny[i]);
	veccopy(q->nz[j], p->nz[i]);
	veccopy(q->nxa[j], p->nxa[i]);
	veccopy(q->nya[j], p->nya[i]);
	veccopy(q->nza[j], p->nza[i]);
}



void Nucleation(void)
{
	if(RND()<Pnuc) {
		Nchain++;
		AddNewChain(chain);	// run after Nchain++
		
#if LGV_MC && BD_2ND
		AddNewChain(chain_try);	// run after Nchain++
#endif
	}
}



void UpdateClosedChain(CHAIN *p)	// update last segment after ring closure
{
	int i, im1;
	
	p->closed = 1;
	
	i = p->nseg - 1;	// last segment
	im1 = i-1;
	
	vecsub(p->rc1[0], p->rc2[i], p->dr2[i]);
	vecsub(p->rc2[i], p->rc1[0], p->dr1[0]);
	
	p->E[i] = 0;
	p->Poff[i] = Poff0;
}



// see if orientation is good for chain closing between p->i1 (tail) q->i2 (head) (CW)
int CheckCloseOrientation(CHAIN *p, CHAIN *q, int i1, int i2)
{
	double tng[3], a1, a2, a3, dr[3];
	
	if(dotprod(p->nz[i1], q->nz[i2]) < cTH) return 0;
	
	// bending angle
	if(K0==0) a1 = getanglePi(p->nx[i1], q->nx[i2]);	// naturally straight
	else {	// naturally curved
		RotMatrix3(-sTH, cTH, p->nz[i1], p->nx[i1], tng);	// tng = ideal angle for q->i2
		a1 = getanglePi(tng, q->nx[i2]);
	}
	if(a1>dTHclose) return 0;
	
	a2 = getanglePi(p->nz[i1], q->nz[i2]);	// torsional angle
	if(a2>dTHclose) return 0;
	
	vecsub(p->r[i1], q->r[i2], dr);	// dr is from head of q to tail of p
	a3 = getanglePi(dr, q->nx[i2]);
	if(a3>dTHclose) return 0;
	
	return 1;
}



void CheckClosure(CHAIN *p, int status)	// status=1: check if timing is right, otherwise no timing check
{
	int nm1, flg;
	double d2, dotp, r[3];
	
	flg = 0;
	nm1 = p->nseg - 1;
	vecsub(p->rhead, p->rtail, r);	// r is from tail to head
	d2 = norm2(r);
	if(d2<Distclose2) flg=1;	// 1: close enough
	if(flg==1) {
		dotp = dotprod(p->nx[0], p->nx[nm1]);
		if(dotp<=0) flg=0;	// 2: same direction
	}
	
	if(status==1) {	// check timing
		if(flg==1) {
			if(RND()>Pon0half) flg=0;	// 3: timing is good
		}
	}
	
	if(flg==1) {	// OK to close
		if(CheckCloseOrientation(p,p,nm1,0)==1) UpdateClosedChain(p);	// closed!
	}
}



void CloseChain(void)	// run this every dt to check closure
{
	CHAIN *p;
	p = chain;
	while(p) {
		if(p->nseg > 3) {
			//CheckClosure(p, 1);	// check closure with timing requirement
			CheckClosure(p, 0);	// check closure without timing requirement
		}
		p = p->next;
	}
}

/*
void CalcEnergy(void)
{
	int i, istart, iend, im1, nm1, nm2;
	double es, eb;
	CHAIN *p;
	
	p = chain;
	while(p) {
		nm1 = p->nseg - 1;
		if(p->closed == 0) {
			p->E[0] = KShalf*pow(p->l[0] - DLT, 2);	// head
			p->E[nm1] = 0;	// no energy for tail
		}
		if(p->nseg > 2) {
			nm2 = p->nseg - 2;
			if(p->closed == 0) {	// chain
				istart = 1;
				iend = nm2;
			}
			else {	// ring
				istart = 0;
				iend = nm1;
			}
			for(i=istart; i<=iend; i++) {
				im1 = i-1;
				if(im1<0) im1 = nm1;
				es = KShalf*(pow(p->l[im1] - DLT, 2)+pow(p->l[i] - DLT, 2));
				eb = KBsythalf*(pow(p->ax[i],2) + pow(p->ay[i],2) + pow(p->az[i],2));
			}
		}
		p = p->next;
	}
}
*/


void UpdateRates(void)
{
	int i;
	double e;
	CHAIN *p;
	
	//CalcEnergy();
	p = chain;
	while(p) {
		
#if OPENMP
#pragma omp parallel for private(i,e)
#endif

		for(i=0; i<(p->nseg); i++) {
			e=min2(p->E[i]/kBT, 100);
			p->Poff[i] = Poff0*exp(e);
		}
		
		p = p->next;
	}
}


void ShiftSegBack(CHAIN *p)	// shift segment values to the right by 1
{
	int i, n;
	n = p->nseg;
	for(i=n; i>0; i--) CopySegments(p, i-1, p, i);
}


void ShiftSegFront(CHAIN *p)	// shift segment values to the left by 1
{
	int i, nm1;
	
	nm1 = (p->nseg)-1;
	for(i=0; i<nm1; i++) CopySegments(p, i+1, p, i);
}



void NewSegPos(CHAIN *p, int tipid, int growpos, double r[3], double nx[3], double ny[3], double nz[3])
{
	double sgn, t1[3], t2[3], dr[3];
	double th, dth, sndth, csdth;
	
	if(growpos==1) sgn=1;	// add to head
	else sgn=-1;	// add to tail
	vecprod(p->nx[tipid], sgn, t1);	// t1 is from tip to new subunit (not tng)
	
	if(K0==0) {	// naturally straight
		veccopy(p->nx[tipid], nx);	// keep current direction
		veccopy(p->ny[tipid], ny);
		veccopy(p->nz[tipid], nz);
		vecprod(t1, Lc2b, dr);
	}
	else {	// naturally curved
		veccopy(p->nz[tipid], nz);
		
		if(growpos==1) {	// add to head
			RotMatrix3(-sTH, cTH, nz, t1, t2);	// curving cw
			crossprod(nz, t2, ny);	// cw: n = nz x t2
		}
		else if(growpos==2) {	// add to tail
			RotMatrix3(sTH, cTH, nz, t1, t2);	// curving ccw
			crossprod(t2, nz, ny);	// ccw: n = t2 x nz
		}
		else {	// planar spiral, add to tail
			th=Pidb*R0mDLThalf/DLT*(sqrt(1+tipid*DLT2/Pi/R0mDLThalf2)-1);
			dth=DLT/(R0mDLThalf+DLT*th/Pidb);
			sndth=sin(dth);
			csdth=cos(dth);
			RotMatrix3(sndth, csdth, nz, t1, t2);	// curving ccw
			crossprod(t2, nz, ny);	// ccw: n = t2 x nz
		}
		
		crossprod(ny, nz, nx);
		
		vecadd(t1, t2, dr);	// center-center direction is (t1+t2)/2
		normalize(dr, dr);
		vecprod(dr, Lc2b, dr);
	}
	
	//printf("\ndr=(%.3g, %.3g, %.3g) -> %.3g\n", dr[0],dr[1],dr[2],norm(dr));
	
	vecadd(p->r[tipid], dr, r);
}



void NewSegPosRing(CHAIN *p, int tipid, int growpos, double r[3], double nx[3], double ny[3], double nz[3])
{	// position of segments in a closed ring of any size (force to be close)
	double da, sn, cs;
	double sgn, t1[3], t2[3], dr[3];
	
	da=Pidb/Nsytperring;
	sn=sin(da);
	cs=cos(da);
	
	if(growpos==1) sgn=1;	// add to head
	else sgn=-1;	// add to tail
	vecprod(p->nx[tipid], sgn, t1);	// t1 is from tip to new subunit (not tng)
	veccopy(p->nz[tipid], nz);
	
	if(growpos==1) {	// add to head
		RotMatrix3(-sn, cs, nz, t1, t2);	// curving cw
		crossprod(nz, t2, ny);	// cw: ny = nz x t2
	}
	else {	// add to tail
		RotMatrix3(sn, cs, nz, t1, t2);	// curving ccw
		crossprod(t2, nz, ny);	// ccw: ny = t2 x nz
	}
	crossprod(ny, nz, nx);
	
	vecadd(t1, t2, dr);	// center-center direction is (t1+t2)/2
	normalize(dr, dr);
	vecprod(dr, Lc2b, dr);
	
	vecadd(p->r[tipid], dr, r);
}



void ChainGrowth(void)
{
	int n, nm1, nm2, flg, growpos;
	double r[3], dr[3], nx[3], ny[3], nz[3];
	CHAIN *p;
	
	p = chain;
	while(p) {
		if(p->closed == 0) {	// for open chains
			flg=0;
			
			if(RND()<Pon && (p->nseg) < NSegMax) flg=1;
				
		#if GRWGAP
			if(flg==1) {
				if((p->nseg)>2 && distance2(p->rhead, p->rtail) < DLT2) flg=0;
			}
		#endif
		
			if(flg==1) {	// OK to add one subunit
				
				if(GrwEnd==2) {	// grow from both ends
					if(RND()<0.5) growpos=1;	// add to head
					else growpos=2;	// add to tail
				}
				
				else if(GrwEnd==1) {	// grow from tip of tube
					if(p->rhead[2] > p->rtail[2]) growpos=1;	// add to head
					else growpos=2;	// add to tail
				}
				
				else {	// grow from bottom of tube
					if(p->rhead[2] < p->rtail[2]) growpos=1;	// add to head
					else growpos=2;	// add to tail
				}
				
				if(growpos==1) {	// add to head
					NewSegPos(p, 0, growpos, r, nx, ny, nz);
					
					veccopy(r, p->rhead);
					ShiftSegBack(p);	// shift segments to tail
					veccopy(r, p->r[0]);
					
					veccopy(nx, p->nx[0]);
					veccopy(ny, p->ny[0]);
					veccopy(nz, p->nz[0]);
					
					vecprod(p->ny[0], Rc2adb, dr);
					vecadd(p->r[0], dr, p->rlinker[0]);
					
					if(Vesicle) {
						vecsub(p->rlinker[0], RCves, dr);	// dr in vesicle's center frame
						mapsphere(Rv, dr);
						vecadd(dr, RCves, p->rtmd[0]);
					}
					
					p->E[0] = 0;
					p->Poff[0] = Poff0;
				}
				else {	// add to tail
					n = p->nseg;	// number of segments
					nm1 = n-1;	// tail seg-id
					nm2 = n-2;
					
					NewSegPos(p, nm1, growpos, r, nx, ny, nz);
					
					veccopy(r, p->rtail);
					veccopy(r, p->r[n]);	// new tail
					
					veccopy(nx, p->nx[nm1]);	// previous tail
					veccopy(ny, p->ny[nm1]);
					veccopy(nx, p->nx[n]);	// new tail
					veccopy(ny, p->ny[n]);
					veccopy(nz, p->nz[n]);
					
					vecprod(p->ny[n], Rc2adb, dr);
					vecadd(p->r[n], dr, p->rlinker[n]);
					
					if(Vesicle) {
						vecsub(p->rlinker[n], RCves, dr);	// dr in vesicle's center frame
						mapsphere(Rv, dr);
						vecadd(dr, RCves, p->rtmd[n]);
					}
					
					p->E[n] = 0;
					p->Poff[nm1] = Poff0;
				}
				(p->nseg)++;
				
				UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
				//GetSytStretch(p);
				
				if(p->nseg >3) CheckClosure(p, 0);	// check closure without timing requirement
			}	// end of flg=1
			
		}	// end of p->closed
		p = p->next;
	}
}



void MoveSeg2LastChain(CHAIN *p, int istart, int cnt)	// run after Nchain++
{
	int i, j, n;
	CHAIN *q, *r;
	
	// adding a new chain
	if(!chain) printf("chain error\n");
	q = chain;
	while(q->next) q = q->next;	// q is the end of chain
	
	r = (CHAIN *)malloc(sizeof(CHAIN));	// r is the last chain
	q->next = r;
	r->next = NULL;
	//r->prev = q;
	
	// copy segments to new chain
	r->id = Nchain;
	r->nseg = cnt;
	r->closed = 0;
	r->Pclose = 0;
	veccopy(p->r[istart], r->rhead);
	veccopy(p->rtail, r->rtail);
	
	n = p->nseg;
	j=0;
	for(i=istart; i<n; i++) {
		CopySegments(p, i, r, j);
		j++;
	}
	
	UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
}



void ReSortChain(CHAIN *p, int id)	// for breaking rings
{
	int i, j, n, nmid;
	CHAIN tempchain, *q;
	
	q = &tempchain;
	
	n = p->nseg;
	nmid = n - id;
	
	// swap i & j, make sure direction is the same (i to i+1)	
	for(i=0; i<n; i++) {	// i=new #
		if(i<nmid) j=i+id;	// j=old #
		else j=i-nmid;
		
		CopySegments(p, i, q, 0);	// backup i
		CopySegments(p, j, p, i);	// copy j to i
		CopySegments(q, 0, p, j);	// copy i to j
	}
	
	UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
}



void ChainBreak(void)
{
	int i, j, n, nm1, nm2, imax, flg, cnt;
	CHAIN *p;
	
	p = chain;
	while(p) {
		n = p->nseg;
		nm1 = n-1;
		nm2 = n-2;
		flg = 0;
		if(p->closed == 0) imax = nm1;
		else imax = n;
		for(i=0; i<imax && !flg; i++) {	// for each link between i & i+1
			if(RND() < p->Poff[i]) {	// time to break
			//if(RND() < p->Poff[i] || fabs(p->az[i])>1) {	// time to break
				flg=1;	// allow only 1 break per chain at each time
				
				if(p->closed == 0) {	// open chain
					if(i==0) {	// head
						ShiftSegFront(p);	// new i=0 is the old i=1
						veccopy(p->r[0], p->rhead);
						(p->nseg)--;
					}
					else if(i==nm2) {	// tail
						veccopy(p->r[nm2], p->rtail);	// new tail
						(p->nseg)--;
					}
					else {	// inside: break into 2 pieces
						j = i+1;
						cnt = n-j;	// number of segments linked to tail
						Nchain++;
						MoveSeg2LastChain(p, j, cnt);
						p->nseg -= cnt;
						veccopy(p->r[i], p->rtail);
					}
				}
				else {	// closed chain: break into an open chain
					ReSortChain(p, i);
					veccopy(p->r[0], p->rhead);
					veccopy(p->r[nm1], p->rtail);
					p->closed = 0;
				}
			}
		}
		p = p->next;
	}
}



void RemoveOrphans(void)	// remove chains with nseg<=1
{
	CHAIN *p, *q, *r;
	
	p = chain;
	while(p && p->nseg <= 1) {	// for head
		q = p->next;
		free(p);
		p = q;
		//p->prev = NULL;
		chain = p;
		Nchain--;
	}
	
	while(p) {	// for body
		q = p->next;
		while(q && q->nseg <= 1) {
			r = q->next;
			free(q);
			q = r;
			p->next = q;
			Nchain--;
		}
		//if(q) q->prev = p;
		p = p->next;
	}
}



void RemoveCoils(void)	// remove chains with nseg>40
{
	CHAIN *p, *q, *r;
	
	p = chain;
	while(p && p->closed == 0 && p->nseg > 30) {	// for head
		q = p->next;
		free(p);
		p = q;
		//p->prev = NULL;
		chain = p;
		Nchain--;
	}
	
	while(p) {	// for body
		q = p->next;
		while(q && q->closed == 0 && q->nseg > 40) {
			r = q->next;
			free(q);
			q = r;
			p->next = q;
			Nchain--;
		}
		//if(q) q->prev = p;
		p = p->next;
	}
}



void DynamicsSyt1(void)
{
	if(Knuc>0) Nucleation();
	CloseChain();
	ChainGrowth();
	
#if VRATES
	UpdateRates();
#endif
	
	//ChainBreak();
	//RemoveOrphans();
	//RemoveCoils();
}

