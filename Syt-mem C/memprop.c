
void CopyMemNode(MNODE **p, MNODE **q)	// copy from p to q
{
	int i, j, k;
	MNODE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,s,t)
#endif
	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			s = &p[i][j];
			t = &q[i][j];
			
			t->id[0] = s->id[0];
			t->id[1] = s->id[1];
			
			//printf("[i,j]=[%d, %d], id=(%d, %d)\n", i, j, t->id[0], t->id[1]);
			
			t->nnbr = s->nnbr;
			t->nface = s->nface;
			t->nfaceplus = s->nfaceplus;
			t->circbnd = s->circbnd;
			
			veccopy(s->r, t->r);
			t->area = s->area;
			t->Av = s->Av;
			t->C1 = s->C1;
			t->C2 = s->C1;
			t->H = s->H;
			t->K = s->K;
			t->LaplaceH = s->LaplaceH;
			
			veccopy(s->e1, t->e1);
			veccopy(s->e2, t->e2);
			veccopy(s->e3, t->e3);
			veccopy(s->f, t->f);
			veccopy(s->drb, t->drb);
			
			t->Eb = s->Eb;
			t->Ea = s->Ea;
			t->Ec = s->Ec;
			t->Es = s->Es;
			t->E = s->E;
			
			for(k=0; k<(s->nnbr); k++) {
				t->nbr[k][0] = s->nbr[k][0];
				t->nbr[k][1] = s->nbr[k][1];
				t->iface[k] = s->iface[k];
				veccopy(s->rnbr[k], t->rnbr[k]);
				t->lnbr[k] = s->lnbr[k];
				veccopy(s->drnbr[k], t->drnbr[k]);
				veccopy(s->dirnbr[k], t->dirnbr[k]);
				t->areanbr[k] = s->areanbr[k];
				veccopy(s->nrmnbr[k], t->nrmnbr[k]);
			}
			
			for(k=0; k<NnbrMaxPlus; k++) t->ifaceplus[k] = s->ifaceplus[k];
		}
	}
}



void CopyMemFace(MFACE *p, MFACE *q)
{
	int i, j;
	MFACE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,j,s,t)
#endif
	for(i=0; i<NMface; i++) {
		s = &p[i];
		t = &q[i];
		
		t->id = s->id;
		for(j=0; j<3; j++) {
			t->inode[j][0] = s->inode[j][0];
			t->inode[j][1] = s->inode[j][1];
		}
		t->nrmsign = s->nrmsign;
		veccopy(s->nrm, t->nrm);
		t->area = s->area;
		t->vol = s->vol;
		veccopy(s->rc, t->rc);
		veccopy(s->r1, t->r1);
		veccopy(s->r2, t->r2);
		veccopy(s->r3, t->r3);
	}
}



void CopyNodeProp(int i1, int j1, int i2, int j2)	// copy some properties of node(i1,j1) to node(i2,j2)
{
	veccopy(mnode[i1][j1].e3, mnode[i2][j2].e3);
	mnode[i2][j2].area = mnode[i1][j1].area;
	mnode[i2][j2].C1 = mnode[i1][j1].C1;
	mnode[i2][j2].C2 = mnode[i1][j1].C2;
	mnode[i2][j2].H = mnode[i1][j1].H;
	mnode[i2][j2].K = mnode[i1][j1].K;
}


void GetGhostGridPos(void)
{
	int i, j;
	double drx[3], dry[3], drxypp[3], drxypm[3];
	
	drx[0]=Lx;	drx[1]=drx[2]=0;
	dry[1]=Ly;	dry[0]=dry[2]=0;
	drxypp[0]=Lx;	drxypp[1]=Ly;	drxypp[2]=0;
	drxypm[0]=Lx;	drxypm[1]=-Ly;	drxypm[2]=0;
	
	// 4 corners
	vecsub(mnode[Nxm1][Nym1].r, drxypp, mnode[0][0].r);	// lower left
	CopyNodeProp(Nxm1,Nym1,0,0);
	
	vecadd(mnode[1][Nym1].r, drxypm, mnode[Nx][0].r);	// lower right
	CopyNodeProp(1,Nym1,Nx,0);
	
	vecsub(mnode[Nxm1][1].r, drxypm, mnode[0][Ny].r);	// upper left
	CopyNodeProp(Nxm1,1,0,Ny);
	
	vecadd(mnode[1][1].r, drxypp, mnode[Nx][Ny].r);	// upper right
	CopyNodeProp(1,1,Nx,Ny);
	
	j=0;	// bottom
	for(i=1; i<Nx; i++) {
		vecsub(mnode[i][Nym1].r, dry, mnode[i][j].r);
		CopyNodeProp(i,Nym1,i,j);
	}
	
	j=Ny;	// top
	for(i=1; i<Nx; i++) {
		vecadd(mnode[i][1].r, dry, mnode[i][j].r);
		CopyNodeProp(i,1,i,j);
	}
	
	i=0;	// left
	for(j=1; j<Ny; j++) {
		vecsub(mnode[Nxm1][j].r, drx, mnode[i][j].r);
		CopyNodeProp(Nxm1,j,i,j);
	}
	
	i=Nx;	// right
	for(j=1; j<Ny; j++) {
		vecadd(mnode[1][j].r, drx, mnode[i][j].r);
		CopyNodeProp(1,j,i,j);
	}
}




void GetMemNbrPositions(MNODE **MN)
{
	int i, j, k, ni, nj;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,ni,nj,m)
#endif

	for(i=1; i<Nx; i++) {	// all effective grids
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			
			for(k=0; k<(m->nnbr); k++) {	// neighbour-k
				ni = m->nbr[k][0];
				nj = m->nbr[k][1];
				veccopy(MN[ni][nj].r, m->rnbr[k]);
				
				vecsub(m->rnbr[k], m->r, m->drnbr[k]);	// displacement pointing to neighbors
				m->lnbr[k] = norm(m->drnbr[k]);	// distance to neighbors
				m->lnbr[k] = max2(m->lnbr[k], eps);
				vecdiv(m->drnbr[k], m->lnbr[k], m->dirnbr[k]);	// unit direction
				/*
				if(vecdivx(m->drnbr[k], m->lnbr[k], m->dirnbr[k])==0) {
					printf("Error in GetMemNbrPositions-a: r=(%.4g, %.4g, %.4g), x=%.4g\n", m->drnbr[k][0], m->drnbr[k][1], m->drnbr[k][2], m->lnbr[k]);
					exit(0);
				}
				*/
			}	// end of k
		}	// end of j
	}	// end of i

}


void GetMemNbrPositionsOnly(MNODE **MN)
{
	int i, j, k, ni, nj;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,ni,nj,m)
#endif

	for(i=1; i<Nx; i++) {	// all effective grids
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			
			for(k=0; k<(m->nnbr); k++) {	// neighbour-k
				ni = m->nbr[k][0];
				nj = m->nbr[k][1];
				veccopy(MN[ni][nj].r, m->rnbr[k]);
			}	// end of k
		}	// end of j
	}	// end of i

}


void GetMemArea(MNODE **MN, MFACE *MF, int status)
{
	int i, j, k;
	MNODE *n;
	MFACE *f;
	
	AreaMem=0;
	
#if OPENMP
#pragma omp parallel for private(i,f) reduction(+:AreaMem)
#endif

	for(i=0; i<NMface; i++) {
		f = &MF[i];
		f->area = getarea3(f->r1, f->r2, f->r3);	// triangle area
		AreaMem += f->area;
	}
	
	AreaMem-=AreaGhost;	// extra area from ghost nodes
	
	Have=0;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,n) reduction(+:Have)
#endif
	
	for(i=1; i<Nx; i++) {
		for(j=1; j<Ny; j++) {
			n = &MN[i][j];
			n->area = 0;
			for(k=0; k<(n->nface); k++) {
				n->areanbr[k] = (MF[n->iface[k]].area)/3.0;	// for area-weighted normal
				n->area += n->areanbr[k];
			}
			n->Av = n->area;	// temporary value
			Have += n->r[2];
		}
	}
	
	Have/=Nxm1Nym1;
	
	if(status==0) {	// initialize
		for(i=0; i<=Nx; i+=Nx) {	// left & right boundary
			for(j=0; j<=Ny; j++) MN[i][j].Av=dX*dY;	// Voronoi area of boundary nodes
		}
		for(i=1; i<Nx; i++) {	// top and bottom boundary
			for(j=0; j<=Ny; j+=Ny) MN[i][j].Av=dX*dY;	// Voronoi area of boundary nodes
		}
		//AreaMem0=AreaMem;
		AreaMem0=Lx*Ly;
		dAmemNode=AreaMem0/Nxm1Nym1;
		dAmemFace=AreaMem0/Nxm1Nym1;
	}
}



void GetMemNormal(MNODE **MN, MFACE *MF)	// run after area calculation, b/c using area-weighted average
{
	int i, j, k;
	double r12[3], r13[3], nrm[3];
	double nv[3];
	MNODE *m;
	MFACE *f;
	
#if OPENMP
#pragma omp parallel for private(i,r12,r13,nrm,f)
#endif

	for(i=0; i<NMface; i++) {	// normal of faces
		f = &MF[i];
		vecsub(f->r1, f->r2, r12);
		vecsub(f->r1, f->r3, r13);
		crossprod(r12, r13, nrm);
		vecprod(nrm, f->nrmsign, nrm);
		normalize(nrm, f->nrm);
	}
	
#if OPENMP
#pragma omp parallel for private(i,j,k,nv,nrm,m,f)
#endif

	for(i=1; i<Nx; i++) {	// normal of nodes
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			veczero(nrm);
			for(k=0; k<(m->nface); k++) {	// neighboring faces
				f = &MF[m->iface[k]];
				vecprod(f->nrm, m->areanbr[k], nv);	// area-weighted normal
				vecadd(nrm, nv, nrm);
			}
			normalize(nrm, m->e3);
		}
	}
}



void funcparaboloidmem(int xi, double r[8][3], double afunc[])	// note: 0<=xi<=3, afunc=dvector(1,ma);
{
	double x, y;
	
	x=r[xi-1][0];
	y=r[xi-1][1];
	
	afunc[1] = x*x;
	afunc[2] = x*y;
	afunc[3] = y*y;
}



void getMemCurv1(MNODE *m)	// find principal axes using paraboloid fitting
{
	int i, k;
	double ang, sn, cs, axis[3];
	double nz[3], rnbr[NnbrMax][3], a, b, c;
	double H, K, dlt;
	
	// fitting, don't use dynamic arrays b/c openmp will crash
	int ndata, ma, ndata1, ma1;
	int ix[NnbrMax+1];
	double z[NnbrMax+1], sig[NnbrMax+1], abc[4], v[4][4], w[4];
	double chisq;
	
	double u[NnbrMax+1][4];
	
	//ndata=NnbrMax;
	ndata = m->nnbr;
	ma=3;	// 3 parameters to fit
	
	ndata1=ndata+1;
	ma1=ma+1;
	
	nz[0]=nz[1]=0;	// nz in local coordinates
	nz[2]=1;
	
	for(i=1; i<=ndata; i++) ix[i] = i;	// shared index for vector (x,y)


	// rotate neighbors so that n->e3 matches nz
	crossprod(m->e3, nz, axis);	// rotate nrm to nz about (axis)
	normalize(axis, axis);
	
	ang = getanglePi(m->e3, nz);	// rotation angle
	sn = sin(ang);
	cs = cos(ang);
	
	// rotate all neighbors
	for(k=0; k<ndata; k++) RotMatrix3(sn, cs, axis, m->drnbr[k], rnbr[k]);
	
	// now fit neighbors with z=a*x^2+b*x*y+c*y^2 in lab coord.
	for(k=1; k<=ndata; k++) {
		z[k] = rnbr[k-1][2];
		sig[k] = drmax;
	}
	
	abc[1]=-0.01;	// initial values: abc[1]=a, abc[2]=b, abc[3]=c
	abc[2]=0;	// b=0 means no tilt
	abc[3]=-0.02;
	
	svdfit(ix,z,rnbr,sig,ndata,abc,ma,u,v,w,&chisq,funcparaboloidmem);
	a=abc[1];
	b=abc[2];
	c=abc[3];
	
	// get temporary curvatures
	
	H = -(a+c);	// H=a+c. Note: minus sign b/c convex surface (H>0) has a,c <0
	K = 4*a*c-b*b;	// K=4ac-b^2
	m->H = H;
	m->K = K;
	
	dlt = H*H - K;
	if(dlt>=0) {
		dlt = sqrt(dlt);
		m->C1 = H + dlt;
		m->C2 = H - dlt;
		
		// special case: fabs(c1)<fabs(c2) and c1>0, c2<0
		if(fabs(m->C1)<fabs(m->C2)) {	// principal direction has max fabs(curvature)!
			m->C1 = H - dlt;
			m->C2 = H + dlt;
		}
	}
	else m->C1 = m->C2 = 0;
}





void getMemCurv2(MNODE *m)	// Watanabe method
{
	int nnbr, nnbrm1, k, k1, k2;
	double phi[NnbrMax], phiave, kn[NnbrMax];
	double kphiave, intk1, intk2;
	double K, H, dlt;
	
	nnbr = m->nnbr;
	nnbrm1 = nnbr-1;
	
	for(k=0; k<nnbr; k++) {	// for each neighbor pair
		k2 = (k+1)%nnbr;	// ccw neighbor
		phi[k] = getanglePi(m->drnbr[k], m->drnbr[k2]);
		kn[k] = -dotprod(m->e3, m->drnbr[k]);	// kn = -2*(nrm.dr)/r^2, negative sign s.t. convex is positive
		kn[k] *= 2/pow(m->lnbr[k],2);
	}
	
	intk1 = intk2 = 0;
	for(k=0; k<nnbr; k++) {
		k1 = (k+nnbrm1)%nnbr;	// cw neighbor
		phiave = (phi[k]+phi[k1])/2.0;
		kphiave = kn[k]*phiave;
		intk1 += kphiave;
		intk2 += kn[k]*kphiave;
	}
	
	H = intk1/Pidb;
	K = 3*H*H - intk2/Pi;
	m->H = H;
	m->K = K;
	
	dlt = H*H-K;
	if(dlt>-1e-5) {
		dlt = sqrt(max2(dlt,0));
		m->C1 = H + dlt;	// principal curvatures
		m->C2 = H - dlt;
	}
	else {	// otherwise keep values of C1 & C2 from paraboloid fitting
		printf("Error in Membrane Curvature: K=%.3g, H=%.3g, dlt=%.3g\n", K, H, dlt);
		///*
		printf("r1=(%.4g, %.4g, %.4g)\n", m->rnbr[0][0], m->rnbr[0][1], m->rnbr[0][2]);
		printf("r2=(%.4g, %.4g, %.4g)\n", m->rnbr[1][0], m->rnbr[1][1], m->rnbr[1][2]);
		printf("r3=(%.4g, %.4g, %.4g)\n", m->rnbr[2][0], m->rnbr[2][1], m->rnbr[2][2]);
		printf("r4=(%.4g, %.4g, %.4g)\n", m->rnbr[3][0], m->rnbr[3][1], m->rnbr[3][2]);
		printf("r5=(%.4g, %.4g, %.4g)\n", m->rnbr[4][0], m->rnbr[4][1], m->rnbr[4][2]);
		printf("r6=(%.4g, %.4g, %.4g)\n", m->rnbr[5][0], m->rnbr[5][1], m->rnbr[5][2]);
		exit(0);
		//*/
		printf("ij=%d, %d\n", m->id[0], m->id[1]);
		printf("n=%d\n", nnbr);
		exit(0);
		
		m->C1 = m->C2 = 0;
	}
}





void getMemCurv3(MNODE **MN, MNODE *m)	// Meryer method
{
	int nnbr, nnbrm1;
	int k, k1, k2, ik, jk, ik1, jk1, ik2, jk2;
	double drki[NnbrMax][3], drkk1[NnbrMax][3], drkk2[NnbrMax][3];
	double r2ki[NnbrMax], r2kk2[NnbrMax];
	double cota[NnbrMax], cotb[NnbrMax];
	double sum, Amix, th, sumth;
	double K, H, dlt;
	MNODE *mk, *mk1, *mk2;
	
	nnbr = m->nnbr;
	nnbrm1 = nnbr-1;
	sum = sumth = Amix = 0;
	
	for(k=0; k<nnbr; k++) {	// for each neighbor k
		k1 = k-1;	// cw neighbor of k
		k2 = k+1;	// ccw neighbor of k
		if(k1<0) k1 = nnbrm1;
		if(k2>nnbrm1) k2 = 0;
		
		ik = m->nbr[k][0];	// k
		jk = m->nbr[k][1];
		
		ik1 = m->nbr[k1][0];	// k1
		jk1 = m->nbr[k1][1];
		
		ik2 = m->nbr[k2][0];	// k2
		jk2 = m->nbr[k2][1];
		
		mk = &MN[ik][jk];
		mk1 = &MN[ik1][jk1];
		mk2 = &MN[ik2][jk2];
		
		vecsub(m->r, mk->r, drki[k]);	// drki is from k to i
		vecsub(mk1->r, mk->r, drkk1[k]);	// drkk1 is from k to k1
		vecsub(mk2->r, mk->r, drkk2[k]);	// drkk2 is from k to k2
		
		r2ki[k] = norm2(drki[k]);
		r2kk2[k] = norm2(drkk2[k]);
		
		cota[k] = cotangent(drki[k], drkk1[k]);	// cota is between k-i and k-k1
		cotb[k] = cotangent(drki[k], drkk2[k]);	// cotb is between k-i and k-k2
	}
	
	for(k=0; k<nnbr; k++) {
		k1 = k-1;	// cw neighbor of k
		k2 = k+1;	// ccw neighbor of k
		if(k1<0) k1 = nnbrm1;
		if(k2>nnbrm1) k2 = 0;
		
		ik = m->nbr[k][0];	// k
		jk = m->nbr[k][1];
		mk = &MN[ik][jk];
		
		sum += dotprod(drki[k], m->e3) * (cotb[k1] + cota[k2]);	// alpha_k = b_k1, beta_k = a_k2, minus sign comes from laplace=-K in the paper
		
		Amix += getAmix_k(drki[k], drkk2[k], r2ki[k], r2kk2[k], r2ki[k2], cotb[k], cota[k2]);	// A_mix of triangle i-k-k2 (Meyer's original formula)
		//Amix += (cotb[k1] + cota[k2])*r2ki[k]/8.0;	// Meyer's formula w/o considering the shape of triangles
		
		th = getanglePi(drki[k], drki[k2]);
		sumth += th;
	}
	
	m->H = H = 0.25*sum/Amix;
	m->K = K = (Pidb - sumth)/Amix;
	
	dlt = H*H-K;
	if(dlt>-1e-5) {
		dlt = sqrt(max2(dlt,0));
		m->C1 = H + dlt;	// principal curvatures
		m->C2 = H - dlt;
	}
	else {
		m->C1 = m->C2 = 0;
	}
}




void getMemCurv(MNODE **MN, MNODE *m)
{
		if(CURVAT==0) getMemCurv1(m);	// paraboloid fitting
		else if(CURVAT==1) getMemCurv2(m);	// Watanabe
		else getMemCurv3(MN, m);	// Meyer
}



void GetMemCurvatures(MNODE **MN)
{
	int i, j;
	MNODE *m;

#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=1; i<Nx; i++) {	// for 1<=i<Nx & 1<=j<Ny (not including Nx & Ny !)
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			getMemCurv(MN, m);
		}
	}
}



void GetMemMaxMin(MNODE **MN)
{
	int i, j;
	MNODE *m;
	
	Curvmin=Zmin=1e10;
	Curvmax=Zmax=-1e10;
	
#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=1; i<Nx; i++) {
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			
#if OPENMP
#pragma omp critical(GetMemMaxMin)
{
#endif
			if(Zmin > m->r[2]) Zmin=m->r[2];
			if(Zmax < m->r[2]) Zmax=m->r[2];
			if(Curvmin > m->H) Curvmin=m->H;
			if(Curvmax < m->H) Curvmax=m->H;
			
#if OPENMP
}
#endif
		}
	}
	//printf("curv: %.3g\t%.3g\n", Curvmin, Curvmax);
	dHmem=Zmax-Zmin;
}



void GetMemZone(void)
{
	int i, j, idx, idy;
	MNODE *m;
	MFACE *f;
	MZONE *z;
	
#if OPENMP
#pragma omp parallel for private(i,j)
#endif

	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) mzone[i][j].nnode = mzone[i][j].nface = 0;
	}
	
#if OPENMP
#pragma omp parallel for private(i,j,idx,idy,m,z)
#endif

	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			m = &mnode[i][j];
			idx = (int)((m->r[0] + Lxhalf)/dX + 1.5);
			idy = (int)((m->r[1] + Lyhalf)/dY + 1.5);
			idx = min2(idx, Nx);
			idy = min2(idy, Ny);
			if(idx>0 && idy>0) {	// non-ghost grids
#if OPENMP
#pragma omp critical(GetMemZone)
{
#endif
				z = &mzone[idx][idy];
				z->inode[z->nnode][0] = i;
				z->inode[z->nnode][1] = j;
				if((z->nnode)<NodeZoneMaxm1) (z->nnode)++;
#if OPENMP
}
#endif
			}
		}
	}
	
	
#if OPENMP
#pragma omp parallel for private(i,idx,idy,f,z)
#endif

	for(i=0; i<NMface; i++) {
		f = &mface[i];
		idx = (int)((f->rc[0] + Lxhalf)/dX + 1.5);
		idy = (int)((f->rc[1] + Lyhalf)/dY + 1.5);
		idx = min2(idx, Nx);
		idy = min2(idy, Ny);
		if(idx>0 && idy>0) {	// non-ghost grids
#if OPENMP
#pragma omp critical(GetMemZone)
{
#endif
			z = &mzone[idx][idy];
			if(z->nface < NodeZoneMax) {
				z->iface[z->nface] = i;
				(z->nface)++;
			}
#if OPENMP
}
#endif
		}
	}	// end of i
}



void CheckMemZone(void)
{
	int i, j;
	MZONE *z;
	
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			z = &mzone[i][j];
			printf("zone[%d][%d]:", i, j);
			//for(k=0; k<(z->nnode); k++) printf("\t(%d, %d)", z->inode[k][0], z->inode[k][1]);
			//for(k=0; k<(z->nface); k++) printf("\t%d", z->iface[k]);
			printf("%d", z->nface);
			printf("\n");
		}
	}
}



void getMemNodeNbrFace_i(MNODE *m)	// get list of faces around memnode m, using zones
{
	int i, j, k;
	int imin, imax, jmin, jmax;
	double x, y;
	MZONE *z;
	
	x = m->r[0] + Lxhalf;	// 0<=x<=Lx
	y = m->r[1] + Lyhalf;	// 0<=y<=Ly
	
	imin = (int)(x/dX+1.5) - NMnbrzoneplushalf;	// lower left corner of searching zones
	jmin = (int)(y/dY+1.5) - NMnbrzoneplushalf;
	imax = imin + NMnbrzoneplus;
	jmax = jmin + NMnbrzoneplus;
	
	m->nfaceplus = 0;
	for(i=imin; i<=imax; i++) {
		for(j=jmin; j<=jmax; j++) {
			if(i>0 && i<=Nx && j>0 && j<=Ny) {	// zone[i][j]
				z = &mzone[i][j];
				for(k=0; k<(z->nface); k++) {
					if(m->nfaceplus < NnbrMaxPlus) {
						m->ifaceplus[m->nfaceplus] = z->iface[k];
						(m->nfaceplus)++;
					}
				}
			}	// end of if
		}	// end of j
	}	// end of i
}




void getMemNodeNbrFaces(void)	// get list of faces around each mem node
{
	int i, j;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=1; i<Nx; i++) {	// for 1<=i<Nx & 1<=j<Ny (not including Nx & Ny !)
		for(j=1; j<Ny; j++) {
			m = &mnode[i][j];
			getMemNodeNbrFace_i(m);
		}
	}
}




void getM_Av(MNODE **MN, MNODE *m)	// calculate the Voronoi area for node m
{
	int nnbr, nnbrm1;
	int k, k1, k2, ik, jk, ik1, jk1, ik2, jk2;
	double drki[NnbrMax][3], drkk1[NnbrMax][3], drkk2[NnbrMax][3];
	double r2ki[NnbrMax];
	double cota[NnbrMax], cotb[NnbrMax];
	double a;
	MNODE *mk, *mk1, *mk2;
	
	nnbr = m->nnbr;
	nnbrm1 = nnbr-1;
	a = 0;
	
	for(k=0; k<nnbr; k++) {	// for each neighbor k
		k1 = k-1;	// cw neighbor of k
		k2 = k+1;	// ccw neighbor of k
		if(k1<0) k1 = nnbrm1;
		if(k2>nnbrm1) k2 = 0;
		
		ik = m->nbr[k][0];	// k
		jk = m->nbr[k][1];
		
		ik1 = m->nbr[k1][0];	// k1
		jk1 = m->nbr[k1][1];
		
		ik2 = m->nbr[k2][0];	// k2
		jk2 = m->nbr[k2][1];
		
		mk = &MN[ik][jk];
		mk1 = &MN[ik1][jk1];
		mk2 = &MN[ik2][jk2];
		
		vecsub(m->r, mk->r, drki[k]);	// drki is from k to i
		vecsub(mk1->r, mk->r, drkk1[k]);	// drkk1 is from k to k1
		vecsub(mk2->r, mk->r, drkk2[k]);	// drkk2 is from k to k2
		
		r2ki[k] = norm2(drki[k]);
		
		cota[k] = cotangent(drki[k], drkk1[k]);	// cota is between k-i and k-k1
		cotb[k] = cotangent(drki[k], drkk2[k]);	// cotb is between k-i and k-k2
	}
	
	for(k=0; k<nnbr; k++) {
		k1 = k-1;	// cw neighbor of k
		k2 = k+1;	// ccw neighbor of k
		if(k1<0) k1 = nnbrm1;
		if(k2>nnbrm1) k2 = 0;
		a += (cotb[k1] + cota[k2]) * r2ki[k];	// Meyer's formula w/o considering the shape of triangles
	}
	m->Av = a/8.0;
}



void getMemNodeAv(MNODE **MN)
{
	int i, j;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=1; i<Nx; i++) {	// for 1<=i<Nx & 1<=j<Ny (not including Nx & Ny !)
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			getM_Av(MN, m);
		}
	}
}



void UpdateMem(MNODE **MN, MFACE *MF)
{
#if BNDCND
	GetGhostGridPos(MN);
#endif

	GetMemNbrPositions(MN);
	
	GetMemArea(MN, MF, 1);
	GetMemNormal(MN, MF);	// run after area calculation

	GetMemCurvatures(MN);
	
	GetMemZone();
	
#if !LAPLAC
	getMemNodeNbrFaces(MN);	// Belkin's formula
#else
	// LAPLAC==1 is for Meyer's formula
	if(LAPLAC==2) getMemNodeAv(MN);	// Vallet-Levy's formula
#endif
	
	//CheckMemZone();
	//exit(0);
}

