
void MemStretch(MNODE **MN)
{
	int i, j, istart, iend, jstart, jend;
	double df, dfn[3], dfp[3], drn[3], drp[3], dr[3];
	
#if BNDCND
	istart=jstart=1;	iend=Nxm1;	jend=Nym1;
#else
	istart=jstart=2;	iend=Nxm1;	jend=Nym1;
#endif

#if OPENMP
#pragma omp parallel for private(i,j,df,dfn,dfp,drn,drp,dr)
#endif

	for(i=istart; i<=iend; i++) {
		for(j=jstart; j<=jend; j++) {
	
		#if CRCBND
			if(MN[i][j].circbnd == 1) continue;
		#endif
			//vecprod(mnode[i][j].f, Mumem0dt, dr);	// too fast, causing problems?
			//vecprod(mnode[i][j].f, 10*Mumemdt, dr);
			//vecprod(MN[i][j].f, MumemXdt, dr);
			df=dotprod(MN[i][j].f, MN[i][j].e3);
			vecprod(MN[i][j].e3, df, dfn);	// normal force
			vecsub(MN[i][j].f, dfn, dfp);	// parallel force
			
			vecprod(dfn, MumemPerpdt, drn);	// normal displacement
			vecprod(dfp, MumemParadt, drp);	// parallel displacement
			vecadd(drn, drp, dr);	// total displacement
			
			limitvector(dr, drmemmin);
		#if BNDCND
			if(i==istart || i==iend || j==jstart || j==jend) dr[0]=dr[1]=0;	// fix x & y of boundary nodes
		#endif
			vecadd(MN[i][j].r, dr, MN[i][j].r);
			
		#if SIHOLE
			if(norm2(MN[i][j].r)>Lxhalf2) MN[i][j].r[2]=0;
		#endif
		}
	}
}


void MemNoise(MNODE **MN)
{
	int i, j, istart, iend, jstart, jend;
	double rnd, dz;
#if OPENMP
	int tid;
#endif
	
#if BNDCND
	istart=jstart=1;	iend=Nxm1;	jend=Nym1;
#else
	istart=jstart=2;	iend=Nxm1;	jend=Nym1;
#endif

#if OPENMP
#pragma omp parallel for private(i,j,rnd,dz,tid)
#endif

	for(i=istart; i<=iend; i++) {
	#if OPENMP
		tid = omp_get_thread_num();
	#endif
		for(j=jstart; j<=jend; j++) {
		#if CRCBND
			if(MN[i][j].circbnd == 1) continue;
		#endif
		
		#if !OPENMP
			rnd = RND();
		#else
			rnd = RNDMP(tid);
		#endif
			dz = dlnoise*(1-2*rnd);
			MN[i][j].r[2] += dz;
			if(Carbon) MN[i][j].r[2] = max2(MN[i][j].r[2], 0);	// if has carbon support
			
			//for(k=0; k<3; k++) MN[i][j].r[k]+=dlnoise*(1-2*RND());
		}
	}
}



void MemCarbon(MNODE **MN)
{
	int istart, iend, jstart, jend;
	int i, j;
	MNODE *m;
	
#if BNDCND
	istart=jstart=1;	iend=Nxm1;	jend=Nym1;
#else
	istart=jstart=2;	iend=Nxm1;	jend=Nym1;
#endif


#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=istart; i<=iend; i++) {	// for each grid
		for(j=jstart; j<=jend; j++) {
			m = &MN[i][j];
		#if CRCBND
			if(m->circbnd == 1) continue;
		#endif
			if(m->r[2]<0) m->r[2]=0;
		}
	}
}



void MemNodeTether(MNODE **MN)	// tethering potential
{
	int i, j, k, inbr, jnbr;
	double dr[3], d2, d, dd;
	MNODE *p, *q;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,inbr,jnbr,dr,d2,d,dd,p,q)
#endif

	for(i=2; i<Nx; i++) {
		for(j=2; j<Ny; j++) {
			p = &MN[i][j];
			for(k=0; k<(p->nnbr); k++) {
				inbr = p->nbr[k][0];
				jnbr = p->nbr[k][1];
				
				q = &MN[inbr][jnbr];
				
				//if(i>=inbr) {	// to avoid double-counting, this will cause imbalance
					vecsub(q->r, p->r, dr);	// dr is from p to q
					
					d2 = norm2(dr);
					// repulsion
					///*
					if(d2 < LBmemmin2) {	// p & q are in contact
						d = max2(sqrt(d2), eps);
						if(i==1 || i==Nxm1 || j==1 || j==Nym1) dd = LBmemmin-d;
						else dd = (LBmemmin - d)/2.0;	// dd>0
						vecprod(dr, dd/d, dr);	// displacement for p
						vecsub(p->r, dr, p->r);	// push back
#if OPENMP	// be careful, use omp critical !!!
#pragma omp critical(MemNodeTether)	// share the same critical name as below
{
#endif
					if(inbr>1 && jnbr>1 && inbr<Nx && jnbr<Ny) vecadd(q->r, dr, q->r);
#if OPENMP
}
#endif
					}
					//*/
					
					// attraction
					if(d2>LBmemmax2) {
						d = sqrt(d2);
						if(i==1 || i==Nxm1 || j==1 || j==Nym1) dd = d-LBmemmax;
						else dd = (d - LBmemmax)/2.0;	// dd>0
						vecprod(dr, dd/d, dr);	// displacement for p
						vecadd(p->r, dr, p->r);
#if OPENMP	// be careful, use omp critical !!!
#pragma omp critical(MemNodeTether)	// share the same critical name as below
{
#endif
						if(inbr>1 && jnbr>1 && inbr<Nx && jnbr<Ny) vecsub(q->r, dr, q->r);
#if OPENMP
}
#endif
					}
				//}
			}
		}
	}
}






void UpdateMemFaceNodePos(MNODE **MN, MFACE *MF)
{
	int i, i1, j1, i2, j2, i3, j3;
	MFACE *f;
	MNODE *m1, *m2, *m3;
	
	for(i=0; i<NMface; i++) {
		f = &MF[i];
		i1 = f->inode[0][0];	j1 = f->inode[0][1];
		i2 = f->inode[1][0];	j2 = f->inode[1][1];
		i3 = f->inode[2][0];	j3 = f->inode[2][1];
		m1 = &MN[i1][j1];
		m2 = &MN[i2][j2];
		m3 = &MN[i3][j3];
		veccopy(m1->r, f->r1);
		veccopy(m2->r, f->r2);
		veccopy(m3->r, f->r3);
		vecave3(f->r1, f->r2, f->r3, f->rc);
	}
}



void MotionMem(MNODE **MN, MFACE *MF)
{
	MemStretch(MN);	// motion
	
	//if(noise) MemNoise(MN);
	
	GetMemNbrPositionsOnly(MN);	// update m->rnbr for RemapMem() after motion
	if(Carbon==1) MemCarbon(MN);
	
	GetMemNbrPositionsOnly(MN);
	MemNodeTether(MN);
	
#if FLPBND
	if(MN==mnode) FlipMem(MN);
#endif

	GetMemNbrPositionsOnly(MN);
	UpdateMemFaceNodePos(MN, MF);
}


