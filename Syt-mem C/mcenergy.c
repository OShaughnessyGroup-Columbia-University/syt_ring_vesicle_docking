
double MCEenergySYT(void)
{
	int i;
	double r, invr, e;
	CHAIN *p;

#if X_PIP2
	double d;
#endif
	
	r=Nsytperring*Lc2b/Pidb;
	invr=1.0/r;
	
	e = 0;
	p = chain;
	while(p) {
		if(p->nseg>0) {	// energy reduction due to oligomerization
			//if(p->closed) e -= (p->nseg)*ESS;	// additional drop due to closure
			//else e -= (p->nseg - 1)*ESS;
			
#if OPENMP
#pragma omp parallel for private(i) reduction(+:e)
#endif

			for(i=0; i<(p->nseg); i++) {
				p->E[i] = 0;
				if(FIXSYT!=2 || CLOSED) {
					if(CLOSED) p->Eb[i] = 0.5*LP*kBT*pow(invr-K0,2)*Lc2b;
					else p->Eb[i] = 0;	// need change !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					p->E[i] += p->Eb[i];
				}
				if(Membrane) {
					//p->Em[i] = ESMpersite - getEsm(p->rsm[i]);	// ground state at r=0
					p->Em[i] = - getEsm(p->rsm[i]);	// ground state at r=inf
					
				#if X_PIP2	// extra PIP2 molecules near syt
					d = getSMdistance(mnode, p->rsm[i]);	// distance to membrane
					p->Em[i] += -getEsmX_PIP2(d);	// extra energy
				#endif
				
					EsytES += p->Em[i];
					p->E[i] += p->Em[i];
				}
				
				if(Vesicle) {
					//p->Ev[i] = ESVpersite - getEsv(p->rsm[i]);	// ground state at r=0
					p->Ev[i] = - getEsv(p->rsm[i]);	// ground state at r=inf
					EsytES += p->Ev[i];
					p->E[i] += p->Ev[i];
				}
				e += p->E[i];
			}
			
		}
		
		p = p->next;
	}
	
	return e;
}



double MCEnergyMEM(void)
{
	int i, j;
	double e;
	MNODE *m;
	
	e = 0;
	
#if OPENMP
#pragma omp parallel for private(i,j,m) reduction(+:e)
#endif

	for(i=1; i<Nx; i++) {	// reset force
		for(j=1; j<Ny; j++) {
			m = &mnode[i][j];
			m->Eb = (KBmemhalf*pow(2*(m->H)-KMspon,2) + KGmem*(m->K)) * (m->area);	// bending with Gaussian modulus
			m->Ea = Gamma*(m->area - dAmemNode);	// areal
			m->E = m->Eb + m->Ea;
			//area += m->area;
			if(Carbon) {
				m->Ec = getEmc(m->r[2], m->area);	// lift from carbon adhesion dE>0
				m->E += m->Ec;
			}
			e += m->E;
		}
	}
	
	return e;
}



double MCEnergyVES(void)
{
	int i;
	double e;
	VNODE *p;
	
#if VSTRCH
	GammaVes = Gamma + KAmem*max2(0, (AreaVes-AreaVes0)/AreaVes0);	// Ka=A*dGamma/dA --> dGamma=Ka*dA/A
#endif

	Evesbend = Evestens = 0;

#if OPENMP
#pragma omp parallel for private(i,p) reduction(+:Evesbend,Evestens)
#endif

	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		p->Eb = KBmemhalf*pow(2*(p->H)-KMspon,2)*(p->area) - dEves0;	// change in bending, Gaussian term sums up to 0 for ves
		p->Ea = GammaVes*(p->area - dAvesNode);	// areal
		p->E = p->Eb + p->Ea;
		Evesbend += p->Eb;
		Evestens += p->Ea;
	}
	
	//Evestens = GammaVes*(AreaVes - AreaVes0);
	e = Evesbend + Evestens;
	return e;
}



double MCEnergySMrepN_ij(CHAIN *p, int id, int i, int j)	// repulsion between syt and mem nodes
{
	int k;
	double dr[3], d2, d;	//, sgn;
	double e;
	MNODE *m;
	
	e = 0;
	m = &mnode[i][j];	// membrane node-m
	
	for(k=2; k<4; k++) {	// for C2AB bead-k
		vecsub(p->rC2AB[id][k], m->r, dr);	// dr is from mnode grid to centers of C2AB beads
		d2 = norm2(dr);
		
		if(d2<drSMrepulsion2) {
			d = max2(sqrt(d2), eps);
			e += KVhalf*pow(drSMrepulsion - d, 2);

		}	// end of if
	}	// end of k
	
	return e;
}



double MCEnergySMrepF_ij(CHAIN *p, int id, int idf)	// repulsion between syt and center of mem faces
{
	int j;
	double dr[3], d2, d;	//, sgn;
	double e;
	MFACE *m;
	
	e = 0;
	m = &mface[idf];	// membrane face-m
	
	for(j=0; j<4; j++) {	// for C2AB bead-j
		vecsub(p->rC2AB[id][j], m->rc, dr);	// dr is from face center to syt center
		d2 = norm2(dr);
		
		if(d2<drSMrepulsion2) {
			d = max2(sqrt(d2), eps);
			e += KVhalf*pow(drSMrepulsion - d, 2);	// *(m->area);
			
		}	// end of if
	}	// end of k
	
	return e;
}



double MCEnergySMrep(void)
{
	int id, i, j, k, n;
	int imin, jmin, imax, jmax;
	int idx, idy, idf;
	double x, y, e;
	CHAIN *p;
	MZONE *z;
	
	e = 0;
	
	p = chain;
	while(p) {
		n = p->nseg;
		
#if OPENMP
#pragma omp parallel for private(i) reduction(+:e)
#endif
		
		for(id=0; id<n; id++) {
			x = p->r[id][0] + Lxhalf;	// 0<=x<=Lx
			y = p->r[id][1] + Lyhalf;	// 0<=y<=Ly
			
			imin = (int)(x/dX+1.5) - NSMnbrhalf;	// lower left corner of searching zones
			jmin = (int)(y/dY+1.5) - NSMnbrhalf;
			imax = imin + NSMnbr;
			jmax = jmin + NSMnbr;
			
			for(i=imin; i<=imax; i++) {
				for(j=jmin; j<=jmax; j++) {
					if(i>0 && i<=Nx && j>0 && j<=Ny) {	// zone[i][j]
						z = &mzone[i][j];
						
						for(k=0; k<(z->nnode); k++) {	// all nodes in zone
							idx = z->inode[k][0];
							idy = z->inode[k][1];
							e += MCEnergySMrepN_ij(p, id, idx, idy);
						}
						
						///*
						for(k=0; k<(z->nface); k++) {	// all faces in zone
							idf = z->iface[k];
							e += MCEnergySMrepF_ij(p, id, idf);
						}
						//*/
					}
				}	// end of j
			}	// end of i
			
		}	// end of id
		
		p = p->next;
		
	}	// end of while
	
	return e;
}

