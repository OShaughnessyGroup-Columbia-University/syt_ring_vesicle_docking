
void CopyVesNode(VNODE *p, VNODE *q)	// copy from p to q
{
	int i, k;
	VNODE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,k,s,t)
#endif
	for(i=0; i<NVnode; i++) {
		s = &p[i];
		t = &q[i];
		
		t->id = s->id;
		t->nnbr = s->nnbr;
		t->nface = s->nface;
		t->nfaceplus = s->nfaceplus;
		t->nedge = s->nedge;
		
		veccopy(s->r0, t->r0);
		veccopy(s->r, t->r);
		veccopy(s->rc, t->rc);
		
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
		veccopy(s->fext, t->fext);
		t->Eb = s->Eb;
		t->Ea = s->Ea;
		t->Es = s->Es;
		t->E = s->E;
		
		for(k=0; k<(s->nnbr); k++) {
			t->nbr[k] = s->nbr[k];
			t->iface[k] = s->iface[k];
			t->iedge[k] = s->iedge[k];
			t->lnbr[k] = s->lnbr[k];
			veccopy(s->drnbr[k], t->drnbr[k]);
		}
		
		for(k=0; k<NnbrMaxPlus; k++) t->ifaceplus[k] = s->ifaceplus[k];
	}
}



void CopyVesEdge(VEDGE *p, VEDGE *q)
{
	int i;
	VEDGE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,s,t)
#endif
	for(i=0; i<NVedge; i++) {
		s = &p[i];
		t = &q[i];
		t->id = s->id;
		t->inode[0] = s->inode[0];
		t->inode[1] = s->inode[1];
		t->l0 = s->l0;
		t->l = s->l;
		t->dl = s->dl;
		veccopy(s->dir, t->dir);
	}
}



void CopyVesFace(VFACE *p, VFACE *q)
{
	int i;
	VFACE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,s,t)
#endif
	for(i=0; i<NVface; i++) {
		s = &p[i];
		t = &q[i];
		t->id = s->id;
		t->inode[0] = s->inode[0];
		t->inode[1] = s->inode[1];
		t->inode[2] = s->inode[2];
		t->nrmsign = s->nrmsign;
		veccopy(s->nrm, t->nrm);
		t->area = s->area;
		t->vol = s->vol;
		
		veccopy(s->rclab, t->rclab);
		veccopy(s->rcrc, t->rcrc);
		veccopy(s->r1, t->r1);
		veccopy(s->r2, t->r2);
		veccopy(s->r3, t->r3);
	}
}




void CalcVesAreaVolume(VNODE *VN, VFACE *VF, int status)	// get vface.area, vface.norm, vface.rclab, vface.rcrc, vface.vol, vnode.area
{
	int i, j;
	double r12[3], r13[3], crs[3], nrm, rc1[3];
	
	//----- volume -----
	
	AreaVes=VolVes=0;
	
#if OPENMP
#pragma omp parallel for private(i,j,r12,r13,crs,nrm,rc1) reduction(+:AreaVes,VolVes)
#endif

	for(i=0; i<NVface; i++) {
		vecsub(VF[i].r1, VF[i].r2, r12);
		vecsub(VF[i].r1, VF[i].r3, r13);
		crossprod(r12, r13, crs);
		nrm=max2(norm(crs),1e-9);
		VF[i].area=0.5*nrm;	// face area
		
		AreaVes+=VF[i].area;
		
		vecprod(crs, (VF[i].nrmsign)/nrm, VF[i].nrm);	// normal of face
		
		vecave3(VF[i].r1, VF[i].r2, VF[i].r3, VF[i].rclab);
		vecsub(VF[i].rclab, RCves, VF[i].rcrc);
		
		vecsub(VF[i].r1, RCves, rc1);
		VF[i].vol = fabs(dotprod(rc1, crs))/6.0;
		
		VolVes += VF[i].vol;
	}
	
	//----- area -----
#if OPENMP
#pragma omp parallel for private(i,j)
#endif
	
	for(i=0; i<NVnode; i++) {	// area of each vnode
		VN[i].area = 0;
		for(j=0; j<(VN[i].nface); j++) {
			VN[i].area += VF[VN[i].iface[j]].area;
		}
		VN[i].area /= 3.0;
	}
	
	
	if(status==0) {	// initialize
		VolVes0=VolVes;
		AreaVes0=AreaVes;
		dVolVes=0;
		dAvesNode=AreaVes0/NVnode;
		dAvesFace=AreaVes0/NVface;
		//for(i=0; i<NVnode; i++) vnode[i].area0 = vnode[i].area;
		//printf("dAvesNode0=%.4g,\tdAvesNode=%.4g,\tNv=%d\n", dAvesNode0, dAvesNode, NVnode);
	}
	else dVolVes=VolVes-VolVes0;
}



void GetVesGeometry(VNODE *VN, VFACE *VF)	// update rc, e3, lnbr, drnbr of each vnode
{
	int i, j, id;
	double nv[3], nvave[3];
	
	CalcVesAreaVolume(VN, VF, 1);
	
#if OPENMP
#pragma omp parallel for private(i,j,id,nv,nvave)
#endif

	for(i=0; i<NVnode; i++) {	// update normal at each vnode
		vecsub(VN[i].r, RCves, VN[i].rc);
		veczero(nvave);
		for(j=0; j<VN[i].nface; j++) {	// for each neighboring face
			id = VN[i].iface[j];	// face id
			//veccopy(vface[id].nrm, nv);	// unweighted
			vecprod(VF[id].nrm, VF[id].area, nv);	// area-weighted
			vecadd(nvave, nv, nvave);
		}
		normalize(nvave, VN[i].e3);	// normal of vnode
	}
	
#if OPENMP
#pragma omp parallel for private(i,j,id)
#endif

	for(i=0; i<NVnode; i++) {
		for(j=0; j<VN[i].nnbr; j++) {	// for each neighbor-id
			id = VN[i].nbr[j];
			vecsub(VN[id].r, VN[i].r, VN[i].drnbr[j]);
			VN[i].lnbr[j] = norm(VN[i].drnbr[j]);
			VN[i].lnbr[j] = max2(VN[i].lnbr[j], dRveseps);	// to avoid divide-by-zero
			//printf("%d\t%d\t%.4g\n", i,j,VN[i].lnbr[j]);
		}
	}
}



void funcparaboloidves(int xi, double r[7][3], double afunc[])	// note: 0<=xi<=3, afunc=dvector(1,ma);
{
	double x, y;
	
	x=r[xi-1][0];
	y=r[xi-1][1];
	
	afunc[1] = x*x;
	afunc[2] = x*y;
	afunc[3] = y*y;
}



void GetVesPrincipalDir(void)	// find principal axes using paraboloid fitting
{
	int i, j, id;
	double nz[3]={0,0,1};	// nz in lab
	double ang, sn, cs, axis[3];
	double dr[3], rnbr[NnbrMaxp1][3], a, b, c;
	double H, K, dlt;
	double th, sth, cth, sth2, sthcth, cth2;
	double l12inv, l22inv, e1[3];
	VNODE *p;
	
	// fitting, don't use dynamic arrays b/c openmp will crash
	int ndatamax, ndata, ma, ndata1, ma1;
	int ix[NnbrMaxp1];
	double z[NnbrMaxp1], sig[NnbrMaxp1], abc[4], v[4][4], w[4];
	double chisq;
	double u[NnbrMaxp1][4];
	
	ndatamax=NnbrMax;	// max # of neighbors
	ma=3;	// 3 parameters to fit
	ma1=ma+1;
	
	for(i=1; i<=ndatamax; i++) ix[i] = i;	// shared index for vector (x,y)
	
#if OPENMP
#pragma omp parallel for private(i,j,id,ang,sn,cs,axis,dr,rnbr,a,b,c,H,K,dlt,\
th,sth,cth,sth2,sthcth,cth2,l12inv,l22inv,e1,p,ndata,ndata1,ix,z,sig,abc,u,v,w,chisq)
#endif

	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		ndata = p->nnbr;
		ndata1 = ndata + 1;
		
		// rotate neighbors so that n->e3 matches nz
		crossprod(p->e3, nz, axis);	// rotate nrm to nz about (axis)
		normalize(axis, axis);
		
		ang = getanglePi(p->e3, nz);	// rotation angle
		sn = sin(ang);
		cs = cos(ang);
		
		// rotate all neighbors
		for(j=0; j<ndata; j++) {
			id = p->nbr[j];
			vecsub(vnode[id].r, p->r, dr);
			RotMatrix3(sn, cs, axis, dr, rnbr[j]);
		}
		
		// now fit neighbors with z=a*x^2+b*x*y+c*y^2 in lab coord.
		for(j=1; j<=ndata; j++) {
			z[j] = rnbr[j-1][2];
			sig[j] = dRveseps;
		}
		
		abc[1]=-0.01;	// initial values: abc[1]=a, abc[2]=b, abc[3]=c
		abc[2]=0;	// b=0 means no tilt
		abc[3]=-0.02;
		
		svdfit(ix,z,rnbr,sig,ndata,abc,ma,u,v,w,&chisq,funcparaboloidves);
		a=abc[1];
		b=abc[2];
		c=abc[3];
		
		/*
		get temporary curvatures, see "A comparison of Gaussian and mean curvature estimation methods 
		on triangular meshes of range image data" by Magid, Soldea and Rivlin
		*/
		
		H = -(a+c);	// H=a+c. Note: minus sign b/c convex surface (H>0) has a,c <0
		K = 4*a*c-b*b;	// K=4ac-b^2
		p->H = H;
		p->K = K;
		
		dlt = H*H - K;
		if(dlt>=0) {
			dlt = sqrt(dlt);
			p->C1 = H + dlt;
			p->C2 = H - dlt;
			
			// special case: fabs(c1)<fabs(c2) and c1>0, c2<0
			if(fabs(p->C1)<fabs(p->C2)) {	// principal direction has max fabs(curvature)!
				p->C1 = H - dlt;
				p->C2 = H + dlt;
			}
		}
		else p->C1 = p->C2 = 0;
		
		/*
		from (a,b,c), get tilt angle of principal axis in local coord.
		calc: lab (x,y)=RotMatrix*(x',y'), plug this into z=a*x^2+b*x*y+c*y^2, 
		let cross term x'*y'=0 and use identity cos^2(x)-sin^2(x)=cos(2x)
		-> (a-c)*sin(2*th)=b*cos(2*th)
		*/
		
		if(a!=c) th=0.5*atan(b/(a-c));	// tan(2*th)=b/(a-c)
		else th=Piquad;	// cos(2*th)=0
		sth = sin(th);
		cth = cos(th);
		sth2 = sth*sth;
		cth2 = cth*cth;
		sthcth = sth*cth;
		
		// determine which is e1
		l12inv = a*cth2 + b*sthcth + c*sth2;	// l12inv=1/(l1)^2, l22inv=1/(l2)^2
		l22inv = a*sth2 - b*sthcth + c*cth2;	// where l1, l2 are the lengths of semiaxes
		
		l12inv = fabs(l12inv);	// don't forget fabs()!!!!!!
		l22inv = fabs(l22inv);
		
		if(l12inv < l22inv) {	// if l1<l2, l2 is major axis. then theta+=pi/2
			th += Pihalf;
			cs = cth;	// do swap correctly
			cth = -sth;
			sth = cs;
		}
		e1[0] = cos(th);//cth;	// e1 is direction of principal axis for rotated neighbors
		e1[1] = sin(th);//sth;
		e1[2] = 0;
		
		// transform e1 back to un-rotated neighbors
		RotMatrix3(-sn, cs, axis, e1, p->e1);
		crossprod(p->e3, p->e1, p->e2);	// 2nd principal axis
		
		// extra guarantees, important!!!
		crossprod(p->e2, p->e3, p->e1);
		normalize(p->e1, p->e1);
		normalize(p->e2, p->e2);
	}
}




void getVesCurv(VNODE *p)	// Watanabe method: "Detection of salient curvature features on polygonal surfaces"
{
	int k, k1, k2, nnbr, nnbrm1;
	double phi[NnbrMax], phiave, kn[NnbrMax];
	double kphiave, intk1, intk2;
	double K, H, dlt;
	
	nnbr = p->nnbr;
	nnbrm1 = nnbr-1;
	
	for(k=0; k<nnbr; k++) {	// for each neighbor pair
		k2 = (k+1)%nnbr;	// ccw neighbor
		phi[k] = getanglePi(p->drnbr[k], p->drnbr[k2]);
		kn[k] = -dotprod(p->e3, p->drnbr[k]);	// kn = -2*(nrm.dr)/r^2, negative sign s.t. convex is positive
		kn[k] *= 2/pow(p->lnbr[k],2);
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
	
	dlt = H*H-K;
	if(dlt>-1e-5) {
		dlt = sqrt(max2(dlt,0));
		p->C1 = H + dlt;	// principal curvatures
		p->C2 = H - dlt;
	}
	else {	// otherwise keep values of C1 & C2 from paraboloid fitting
		/*
		printf("Error in Membrane Curvature: K=%.3g, H=%.3g, dlt=%.3g\n", K, H, dlt);
		printf("r1=(%.4g, %.4g, %.4g)\n", m->rnbr[0][0], m->rnbr[0][1], m->rnbr[0][2]);
		printf("r2=(%.4g, %.4g, %.4g)\n", m->rnbr[1][0], m->rnbr[1][1], m->rnbr[1][2]);
		printf("r3=(%.4g, %.4g, %.4g)\n", m->rnbr[2][0], m->rnbr[2][1], m->rnbr[2][2]);
		printf("r4=(%.4g, %.4g, %.4g)\n", m->rnbr[3][0], m->rnbr[3][1], m->rnbr[3][2]);
		exit(0);
		*/
	}
	
	p->H = (p->C1 + p->C2)/2.0;
	p->K = (p->C1) * (p->C2);
}



void GetVesCurvatures(VNODE *VN)
{
	int i;
	
#if OPENMP
#pragma omp parallel for private(i)
#endif
	for(i=0; i<NVnode; i++) getVesCurv(&VN[i]);
}




void GetVesZone(VNODE *VN, VFACE *VF)	// for ves-mem interactions
{
	int i, j, idx, idy;
	VNODE *n;
	VFACE *f;
	VZONE *z;
	
#if OPENMP
#pragma omp parallel for private(i,j)
#endif

	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			vzone[i][j].nnode = vzone[i][j].nface = 0;
			vzone[i][j].nnodedown = vzone[i][j].nfacedown = 0;
		}
	}
	
#if OPENMP
#pragma omp parallel for private(i,idx,idy,n,z)
#endif

	for(i=0; i<NVnode; i++) {
		n = &VN[i];
		idx = (int)((n->r[0] + Lxhalf)/dX + 1.5);
		idy = (int)((n->r[1] + Lyhalf)/dY + 1.5);
		idx = min2(idx, Nx);
		idy = min2(idy, Ny);
		if(idx>0 && idy>0) {	// non-ghost grids
#if OPENMP
#pragma omp critical(GetVesZone)
{
#endif
			z = &vzone[idx][idy];
			if(z->nnode < NodeZoneMax) {
				z->inode[z->nnode] = i;
				(z->nnode)++;
				if(n->e3[2]<0) {	// facing down
					z->inode[z->nnodedown] = i;
					(z->nnodedown)++;
				}
			}
#if OPENMP
}
#endif
		}	// end of if
	}	// end of i
	
	
#if OPENMP
#pragma omp parallel for private(i,idx,idy,f,z)
#endif

	for(i=0; i<NVface; i++) {
		f = &VF[i];
		idx = (int)((f->rclab[0] + Lxhalf)/dX + 1.5);
		idy = (int)((f->rclab[1] + Lyhalf)/dY + 1.5);
		idx = min2(idx, Nx);
		idy = min2(idy, Ny);
		if(idx>0 && idy>0) {	// non-ghost grids
#if OPENMP
#pragma omp critical(GetVesZone)
{
#endif
			z = &vzone[idx][idy];
			if((z->nface)<NodeZoneMax) {
				z->iface[z->nface] = i;
				(z->nface)++;
				if(n->e3[2]<0) {	// facing down
					z->iface[z->nfacedown] = i;
					(z->nfacedown)++;
				}
			}
#if OPENMP
}
#endif
		}
	}	// end of i
}



void CheckVesZone(void)
{
	int i, j, k;
	MZONE *z;
	
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			z = &mzone[i][j];
			printf("zone[%d][%d]:", i, j);
			for(k=0; k<(z->nnode); k++) printf("\t(%d, %d)", z->inode[k][0], z->inode[k][1]);
			printf("\n");
		}
	}
}



void getVesNodeNbrFace_i(VNODE *VN, VNODE *p)	// get list of faces around vesnode v, using zones
{
	int i, j, k;
	int imin, imax, jmin, jmax;
	double x, y;
	VZONE *z;
	VNODE *n;
	
	x = p->r[0] + Lxhalf;	// 0<=x<=Lx
	y = p->r[1] + Lyhalf;	// 0<=y<=Ly
	
	imin = (int)(x/dX+1.5) - NMnbrzoneplushalf;	// lower left corner of searching zones
	jmin = (int)(y/dY+1.5) - NMnbrzoneplushalf;
	imax = imin + NMnbrzoneplus;
	jmax = jmin + NMnbrzoneplus;
	
	p->nfaceplus = 0;
	for(i=imin; i<=imax; i++) {
		for(j=jmin; j<=jmax; j++) {
			if(i>0 && i<=Nx && j>0 && j<=Ny) {	// zone[i][j]
				z = &vzone[i][j];
				for(k=0; k<(z->nface); k++) {
					n = &VN[z->iface[k]];
					if(distance2(n->r, p->r)<h_belkincutoff2) {	// not in opposite hemispheres
						if(p->nfaceplus < NnbrMaxPlus) {
							p->ifaceplus[p->nfaceplus] = z->iface[k];
							(p->nfaceplus)++;
						}
					}
				}
			}	// end of if
		}	// end of j
	}	// end of i
}




void getVesNodeNbrFaces(VNODE *VN)	// get list of faces around each ves node
{
	int i;
	VNODE *p;
	
#if OPENMP
#pragma omp parallel for private(i,p)
#endif

	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		getVesNodeNbrFace_i(VN, p);
	}
}




void getV_Av(VNODE *p)	// calculate the Voronoi area for node p
{
	p->Av = 0;	// temp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}



void getVesNodeAv(VNODE *VN)
{
	int i;
	VNODE *p;
	
#if OPENMP
#pragma omp parallel for private(i,p)
#endif

	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		getV_Av(p);	// temp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
}



void UpdateVes(VNODE *VN, VFACE *VF)
{
	GetVesGeometry(VN, VF);
	GetVesCurvatures(VN);
	
	if(Membrane) GetVesZone(VN, VF);
	//CheckVesZone();
	
#if !LAPLAC
	getVesNodeNbrFaces(VN);	// Belkins' formula
#else
	// LAPLAC==1 is for Meyer's formula
	if(LAPLAC==2) getVesNodeAv(VN);	// Vallet-Levy's formula, not finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
}




