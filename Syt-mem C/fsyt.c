
void ResetSytForceTorque(CHAIN *CH)
{
	int i, j;
	CHAIN *p;
	
	p=CH;
	while(p) {
		veczero(p->ftot);
		veczero(p->tautot);
	
#if OPENMP
#pragma omp parallel for private(i,j)
#endif

		for(i=0; i<(p->nseg); i++) {
			veczero(p->f[i]);
			veczero(p->ftmd[i]);
			
			for(j=0; j<4; j++) veczero(p->frep[i][j]);
			veczero(p->fatt[i]);
			veczero(p->fext[i]);
			
			veczero(p->taust[i]);
			veczero(p->taurp[i]);
			veczero(p->taubd[i]);
			veczero(p->tauts[i]);
			veczero(p->taum[i]);
			veczero(p->tautmd[i]);
			veczero(p->tauorient[i]);
			veczero(p->fsm[i]);
			
			p->Eb[i] = p->Em[i] = p->Ev[i] = p->E[i] = 0;
		}
		
		p = p->next;
	}
}



void UpdateChain(CHAIN *CH)	// update Nclosed, Nopen, rcenter, rmax, rave after motion
{
	int i, n, nc;
	double r, rmax, rave;
	CHAIN *p;
	
	Nclosed = 0;
	p = CH;
	while(p) {
		n = p->nseg;
		if(n==0) {
			p = p->next;
			continue;
		}
		
		if(p->closed == 1) Nclosed++;
		
		veccopy(p->r[0], p->rhead);
		veccopy(p->r[n-1], p->rtail);
		
		if(n>Nsytperring0 && !CLOSED) nc=Nsytperring0;	// get rcenter from the first nc syts
		else nc=max2(n,1);
		
		veczero(p->rcenter);
		for(i=0; i<nc; i++) vecadd(p->rcenter, p->r[i], p->rcenter);
		vecdiv(p->rcenter, 1.0*nc, p->rcenter);
		
		rmax = rave = 0;
		for(i=0; i<n; i++) {
			vecsub(p->r[i], p->rcenter, p->rc[i]);
			r = norm(p->rc[i]);
			rmax = max2(rmax, r);
			rave += r; 
		}
		rave/=n;
		p->rmax = rmax;
		p->rave = rave;
		
		p = p->next;
	}
	
	Nopen = Nchain - Nclosed;
}



void GetSytStretch(CHAIN *p)	// update p->l, dl, dr. Run this after UpdateSytSitesAll()
{
	int i, n, nm1;
	int im1, ip1;
	
	n = p->nseg;
	
	if(n==1) {
		veczero(p->dr1[0]);
		veczero(p->dr2[0]);
		return;
	}
	
	nm1 = n-1;
	
#if OPENMP
#pragma omp parallel for private(i,im1,ip1)
#endif

	for(i=0; i<n; i++) {
		im1 = i-1;
		ip1 = i+1;
		if(im1<0) im1=nm1;
		if(ip1>=n) ip1=0;
		if(p->closed==0 && i==0) {	// head
			vecsub(p->rc2[ip1], p->rc1[i], p->dr1[i]);
			veczero(p->dr2[i]);
		}
		else if(p->closed==0 && i==nm1) {	// tail
			veczero(p->dr1[i]);
			vecsub(p->rc1[im1], p->rc2[i], p->dr2[i]);
		}
		else {
			vecsub(p->rc2[ip1], p->rc1[i], p->dr1[i]);	// vector from i to i+1
			vecsub(p->rc1[im1], p->rc2[i], p->dr2[i]);	// vector from i to i-1
		}
	}
	
}



void UpdateSytStretchAll(CHAIN *CH)	// get l, dl, dr
{
	CHAIN *p;
	
	p=CH;
	while(p) {
		GetSytStretch(p);
		p = p->next;
	}
}



void SelfAvoidF(CHAIN *CH)	// point-to-point
{
	int i, j, n, nm1, nm2;
	double dd, d, dr, r[3], f, fv[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n>2) {	// exclude immediate neighbours
			nm1 = n-1;
			nm2 = n-2;
			
#if OPENMP
#pragma omp parallel for private(i,j,dd,d,dr,r,f,fv)
#endif

			for(i=0; i<nm2; i++) {
				for(j=i+2; j<n; j++) {
					dd = distance2(p->r[i], p->r[j]);
					if(dd<RSSselfavoid2) {	// in contact
						d = sqrt(dd);
						dr = RSSselfavoid-d;
						vecsub(p->r[i], p->r[j], r);	// r is from j to i
						normalize(r, r);
						f = KV*dr;
						vecprod(r, f, fv);
#if OPENMP
#pragma omp critical(SelfAvoidF1)
{
#endif
						vecadd(p->f[i], fv, p->f[i]);
						vecsub(p->f[j], fv, p->f[j]);	// be careful, use omp critical !!!
#if OPENMP
}
#endif
					}
				}
			}
		}
		p = p->next;
	}
}



void RepulseF(CHAIN *CH)	// point-to-point
{
	int i, j;
	double d, d1, d2, dd, dr, r[3], f, fv[3];
	CHAIN *p, *q;
	
	p = CH;
	while(p->next) {
		q = p->next;
		while(q) {	// interactions between p & q
			d1 = distance(p->rcenter, q->rcenter);
			d2 = (p->rmax) + (q->rmax) + RSSselfavoid;
			if(d1<d2) {	// two chains are close enough
						
#if OPENMP
#pragma omp parallel for private(i,j,d,dd,dr,r,f,fv)
#endif

				for(i=0; i<(p->nseg); i++) {	// for each segment in p
					for(j=0; j<(q->nseg); j++) {	// for each segment in q
						dd = distance2(p->r[i], q->r[j]);
						if(dd<RSSselfavoid2) {	// in contact
							d=sqrt(dd);
							dr=RSSselfavoid-d;	// overlapping length
							vecsub(p->r[i], q->r[j], r);	// r is from q to p
							normalize(r, r);
							f=KV*dr;
							vecprod(r, f, fv);
#if OPENMP
#pragma omp critical(RepulseF)
{
#endif
							vecadd(p->f[i], fv, p->f[i]);
							vecsub(q->f[j], fv, q->f[j]);	// be careful, use omp critical !!!
#if OPENMP
}
#endif
						}
					}
				}
			}
			q = q->next;
		}
		p = p->next;
	}
}




void CarbonF(CHAIN *CH)
{
	int i, j;
	double z, f;
	CHAIN *p;
	
	p = CH;
	while(p) {
	
#if OPENMP
#pragma omp parallel for private(i,j,z,f)
#endif

		for(i=0; i<(p->nseg); i++) {
			for(j=0; j<4; j++) {
				z = p->rC2AB[i][j][2];
				if(z<0) {
					f = -KV*z;	// f>0
					p->f[i][2] += f;
					p->frep[i][j][2] += f;
				}
			}
		}
		
		p = p->next;
	}
}





void AllSytForce(CHAIN *CH)
{
	int i, j;
	CHAIN *p;
	
	p = CH;
	while(p) {
		veczero(p->ftot);
		for(i=0; i<(p->nseg); i++) {
			
			vecadd(p->fext[i], p->fatt[i], p->fext[i]);
			for(j=0; j<4; j++) vecadd(p->fext[i], p->frep[i][j], p->fext[i]);
			
			vecadd(p->ftot, p->fext[i], p->ftot);	// total force on chain
		}
		p = p->next;
	}
}



void StericForce(CHAIN *CH)
{
	SelfAvoidF(CH);
	RepulseF(CH);
	if(Carbon) CarbonF(CH);
	
	//AllSytForce(CH);
}



void ForceSyt(CHAIN *CH)
{
	//ResetSytForceTorque(CH);
	
	UpdateChain(CH);	// update Nclosed, rcenter, rmax, rave after motion
	UpdateSytStretchAll(CH);	// update p->l, dl, nnxt, nxa, nya, nza

	StericForce(CH);	// steric forces
}

