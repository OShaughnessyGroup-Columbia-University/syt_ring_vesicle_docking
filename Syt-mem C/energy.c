
double getEsmNear(double r, double cs, double radius)	// point-disc energy, near-field
{
	double radius2, radius3, radius5, radius7, exp_radius_L;
	double cs2, cs3, cs4, cs5, cs6, cs7, cs8;
	double c0, c1, c2, c3, c4, c5, c6, c7, c8;
	double p0, p1, p2, p3, p4, p5, p6, p7, p8;
	double r2, r3, r4, r5, r6, r7, r8;
	double energy;
	
	radius2 = radius*radius;	// radius = effective radius for membrane patch
	radius3 = radius2*radius;
	radius5 = radius3*radius2;
	radius7 = radius5*radius2;
	exp_radius_L = exp(-radius/LDebye);
		
	r2 = r*r;
	r3 = r2*r;
	r4 = r3*r;
	r5 = r4*r;
	r6 = r5*r;
	r7 = r6*r;
	r8 = r7*r;
	//r9 = r8*r;
	//r10 = r9*r;
	
	cs2 = cs*cs;
	cs3 = cs2*cs;
	cs4 = cs3*cs;
	cs5 = cs4*cs;
	cs6 = cs5*cs;
	cs7 = cs6*cs;
	cs8 = cs7*cs;
	//cs9 = cs8*cs;
	//cs10 = cs9*cs;
	//cs11 = cs10*cs;
	
	// coefficients of series expansion of potential: exp[-z/D]-exp[-sqrt(z^2+R^2)/D] = c0+c1*r+c2*r^2+...
	c0 = 1 - exp_radius_L;
	c1 = -1.0/LDebye;
	c2 = 0.5*(radius+exp_radius_L*LDebye)/radius/LDebye2;
	c3 = -1/6.0/LDebye3;
	c4 = (radius3 - 3*exp_radius_L*LDebye2*(radius+LDebye))/24.0/radius3/LDebye4;
	c5 = -1/120.0/LDebye5;
	c6 = (radius5+15*exp_radius_L*LDebye3*(radius2+3*radius*LDebye+3*LDebye2))/720.0/radius5/LDebye6;
	c7 = -1/5040.0/LDebye7;
	c8 = (radius7-105*exp_radius_L*LDebye4*(radius3+6*radius2*LDebye+15*radius*LDebye2+15*LDebye3))/40320.0/radius7/LDebye8;
	//c9 = -1/362880.0/LDebye9;
	
	p0 = 1;
	p1 = cs;	// Legendre polynomials
	p2 = (3*cs2-1)/2.0;
	p3 = (5*cs3-3*cs)/2.0;
	p4 = (35*cs4-30*cs2+3)/8.0;
	p5 = (63*cs5-70*cs3+15*cs)/8.0;
	p6 = (231*cs6-315*cs4+105*cs2-5)/16.0;
	p7 = (429*cs7-693*cs5+315*cs3-35*cs)/16.0;
	p8 = (6435*cs8-12012*cs6+6930*cs4-1260*cs2+35)/128.0;
	//p9 = (12155*cs9-25740*cs7+18018*cs5-4620*cs3+315*cs)/128.0;
	//p10 = (46189*cs10-109395*cs8+90090*cs6-30030*cs4+3465*cs2-63)/256.0;
	//p11 = (88179*cs11-230945*cs9+218790*cs7-90090*cs5+15015*cs3-693*cs)/256.0;
	
	//pot = c0*p0 + c1*p1*r + c2*p2*r2 + c3*p3*r3 + ...
	//pot *= Pidb*DebyeMemCoeff*LDebye;
	
	energy = c0*p0 + c1*p1*r + c2*p2*r2 + c3*p3*r3 + c4*p4*r4 + c5*p5*r5 + c6*p6*r6 + c7*p7*r7 + c8*p8*r8;
	energy *= Pidb*LDebye;	// without the DebyeMemCoeff (alpha) or DebyeVesCoeff factor
	
	return energy;	// energy >0
}



double getEsmFar(double r, double cs, double radius)	// point-disc energy, far-field
{	// from Rowan et al. Mol. Phys. 2000
	double radius2, radius4, radius6, exp_r_L;
	double cs2, cs4, sn2, sn4;
	double c1, c2, c3, c4, c5;
	double r2, r3, r4, r5, r6;
	double Z, Q, F;
	double energy;
	
	radius2 = radius*radius;	// radius = effective radius for membrane patch
	radius4 = radius2*radius2;
	radius6 = radius4*radius2;
	exp_r_L = exp(-r/LDebye);
		
	r2 = r*r;
	r3 = r2*r;
	r4 = r3*r;
	r5 = r4*r;
	r6 = r5*r;
	
	cs2 = cs*cs;
	cs4 = cs2*cs2;
	sn2 = 1-cs2;
	sn4 = sn2*sn2;
	
	Z = Pi*radius2;	// without the DebyeMemCoeff (alpha) term
	Q = -Pi/4*radius4;
	F = Pi/8*radius6;
	
	c1 = Z+Q/2/LDebye2*(cs2-1)+F/24/LDebye4*sn4;
	c2 = Q/2/LDebye*(3*cs2-1)-F/3/LDebye3*(sn2-5/32.0*sn4);
	c3 = Q/2*(3*cs2-1)+F/3/LDebye2*(1-6*sn2+45/8.0*sn4);
	c4 = F/8/LDebye*(35*cs4-30*cs2+3);
	c5 = F/8*(35*cs4-30*cs2+3);
	
	energy = c1/r + c2/r2 + c3/r3 + c4/r4 + c5/r5;
	energy *= exp_r_L;	// without the DebyeMemCoeff or DebyeVesCoeff factor
	return energy;	// energy>0
}



double getEsmPF(double r, double cs, double a)	// point-face attraction
{
	double radius;
	double dr, r0, r1, r2;
	double energy;
	
	radius = sqrt(a/Pi);
	r0 = 1.0*radius;
	
	dr = 0.1*r0;
	r1 = r0 - dr/2;
	r2 = r0 + dr/2;
	
	if(r<=r1) energy = max2(0, getEsmNear(r, cs, radius));
	else if(r>=r2) energy = max2(0, getEsmFar(r, cs, radius));
	else energy = max2(0, getEsmNear(r, cs, radius))*(1-(r-r1)/dr) + max2(0, getEsmFar(r, cs, radius))*(r-r1)/dr;	// linear transition
	
	// common prefactor for near & far is q*alpha=ESM0/LDebye0/2/pi
	return energy;	// energy>0 and decreases monotonically with r 
}



double getEsmPP(double r, double a)	// point-point attraction
{
	double rcut, ri, e;
	
	rcut = interpolate_func(r_cutoff_sm, a, amin_cutoff_sm, da_cutoff_sm, n_cutoff_smv);
	if(r < rcut) ri = rcut;
	else ri = r;
	
	e = exp(-r/LDebye) / r;
	
	// common prefactor is q*alpha=ESM0/LDebye0/2/pi
	return e;	// energy>0 and decreases monotonically with r 
}




double getEsm(double rsm[3])	// just one attraction site !!!!!
{
	int i, j, k, idx, idy;
	int imin, imax, jmin, jmax;
	double x, y, r[3], d2, d;
	double e, energy;
	MNODE *m;
	MZONE *z;
#if !ATTRCT
	double nr[3], cs;
#endif
	
	x = rsm[0] + Lxhalf;
	y = rsm[1] + Lyhalf;
	
	imin = (int)(x/dX+1.5) - NSMnbrhalf;	// lower left corner of searching zones
	jmin = (int)(y/dY+1.5) - NSMnbrhalf;
	imax = imin + NSMnbr;
	jmax = jmin + NSMnbr;
	
	energy=0;
	for(i=imin; i<=imax; i++) {
		for(j=jmin; j<=jmax; j++) {
			if(i>0 && i<=Nx && j>0 && j<=Ny) {	// zone[i][j]
				z = &mzone[i][j];
				for(k=0; k<(z->nnode); k++) {	// all nodes in zone
					idx = z->inode[k][0];
					idy = z->inode[k][1];
					m = &mnode[idx][idy];
					vecsub(m->r, rsm, r);	// r pointing from syt site to membrane grid
					d2 = norm2(r);
					if(d2<RSMmax2) {
						d = max2(sqrt(d2), eps);
					#if ATTRCT	// point-point attraction
						e = getEsmPP(d, m->area);	// e>0 and decreases monotonically with r
					#else	// point face attraction
						vecdiv(r, d, nr);	// nr pointing from syt site to membrane grid
						cs = fabs(dotprod(m->e3, nr));
						e = getEsmPF(d, cs, m->area);	// e>0 and decreases monotonically with r
					#endif
						energy += DebyeMemCoeff * e;	// energy>0 and decreases monotonically with r
					}
				}
			}
		}
	}
	return energy;	// energy>0 and decreases monotonically with r
}





double getSMdistance(MNODE **MN, double rsm[3])
{
	int i, j, k;
	int imin, jmin, imax, jmax, idx, idy;
	double x, y, d2, d2min, dn, dnx;
	double dr[3];
	MZONE *z;
	MNODE *m;
	MFACE *f;
	
	x = rsm[0] + Lxhalf;	// 0<=x<=Lx
	y = rsm[1] + Lyhalf;	// 0<=y<=Ly
	
	imin = (int)(x/dX+1.5) - NSMnbrhalf;	// lower left corner of searching zones
	jmin = (int)(y/dY+1.5) - NSMnbrhalf;
	imax = imin + NSMnbr;
	jmax = jmin + NSMnbr;
	
	imin = max2(imin, 1);
	imax = min2(imax, Nx);
	jmin = max2(jmin, 1);
	jmax = min2(jmax, Ny);
	
	d2min = inf;
	dnx = 0;
	
	for(i=imin; i<=imax; i++) {	// zone[i][j]
		for(j=jmin; j<=jmax; j++) {
			z = &mzone[i][j];
			
			for(k=0; k<(z->nnode); k++) {	// all nodes in zone
				idx = z->inode[k][0];
				idy = z->inode[k][1];
				m = &MN[idx][idy];
				vecsub(rsm, m->r, dr);	// from vertex
				d2 = norm2(dr);
				if(d2 < d2min) {
					dn = dotprod(dr, m->e3);
					if(dn>0) {	// site is outside membrane
						d2min = d2;
						dnx = dn;
					}
				}
			}
			
			for(k=0; k<(z->nface); k++) {	// all faces in zone
				idx = z->iface[k];
				f = &mface[idx];
				vecsub(rsm, f->rc, dr);	// from face center
				d2 = norm2(dr);
				if(d2 < d2min) {
					dn = dotprod(dr, f->nrm);
					if(dn>0) {	// site is outside membrane
						d2min = d2;
						dnx = dn;
					}
				}
			}
		}	// end of j-loop
	}	// end of i-loop
	
	return dnx;
}





void getSMdistanceplus(MNODE **MN, double rsm[3], double *dist, double nrm[3], int mi[3], int mj[3])	// normal is from membrane to KKKK
{
	int i, j, k;
	int imin, jmin, imax, jmax, idx, idy;
	double x, y, d2, d2min, dn;
	double dr[3];
	MZONE *z;
	MNODE *m;
	MFACE *f;
	
	x = rsm[0] + Lxhalf;	// 0<=x<=Lx
	y = rsm[1] + Lyhalf;	// 0<=y<=Ly
	
	imin = (int)(x/dX+1.5) - NSMnbrhalf;	// lower left corner of searching zones
	jmin = (int)(y/dY+1.5) - NSMnbrhalf;
	imax = imin + NSMnbr;
	jmax = jmin + NSMnbr;
	
	imin = max2(imin, 1);
	imax = min2(imax, Nx);
	jmin = max2(jmin, 1);
	jmax = min2(jmax, Ny);
	
	d2min = inf;
	*dist = 0;
	veczero(nrm);
	mi[0]=mi[1]=mi[2]=0;
	mj[0]=mj[1]=mj[2]=0;
	
	for(i=imin; i<=imax; i++) {	// zone[i][j]
		for(j=jmin; j<=jmax; j++) {
			z = &mzone[i][j];
			
			for(k=0; k<(z->nnode); k++) {	// all nodes in zone
				idx = z->inode[k][0];
				idy = z->inode[k][1];
				m = &MN[idx][idy];
				vecsub(rsm, m->r, dr);	// from vertex
				d2 = norm2(dr);
				if(d2 < d2min) {
					dn = dotprod(dr, m->e3);
					if(dn>0) {	// site is outside membrane
						d2min = d2;
						*dist = dn;
						vecdiv(dr, max2(sqrt(d2),1e-10), nrm);	// normal direction from mem to KKKK
						mi[0] = mi[1] = mi[2] = idx;	// just one mem node
						mj[0] = mj[1] = mj[2] = idy;
					}
				}
			}
			
			for(k=0; k<(z->nface); k++) {	// all faces in zone
				idx = z->iface[k];
				f = &mface[idx];
				vecsub(rsm, f->rc, dr);	// from face center
				d2 = norm2(dr);
				if(d2 < d2min) {
					dn = dotprod(dr, f->nrm);
					if(dn>0) {	// site is outside membrane
						d2min = d2;
						*dist = dn;
						vecdiv(dr, max2(sqrt(d2),1e-10), nrm);	// normal direction from mem to KKKK
						mi[0] = f->inode[0][0];	mj[0] = f->inode[0][1];	// force is shared among 3 vertices
						mi[1] = f->inode[1][0];	mj[1] = f->inode[1][1];
						mi[2] = f->inode[2][0];	mj[2] = f->inode[2][1];
					}
				}
			}
		}	// end of j-loop
	}	// end of i-loop
}



double getEsmX_PIP2(double r)	// point-point attraction between KKKK and an extra PIP2 molecule
{
	double e;
	
	if(r<=r0_PIP2) e = ESP;
	else e = QsytQpip2_4pie0e/r*exp(-r/LDebye);	// e>0 in pN*nm
	return e;	// energy>0 and decreases monotonically with r, 
}



double getEsv(double rsm[3])	// just one attraction site !!!!!
{
	int i;
	double r[3], d2, d;
	double e, energy;
	VNODE *n;
#if !ATTRCT
	double nr[3], cs;
#endif
	
	energy=0;
	for(i=0; i<NVnode; i++) {
		n = &vnode[i];
		vecsub(n->r, rsm, r);
		d2=norm2(r);
		if(d2<RSMmax2) {
			d=sqrt(d2);
		#if ATTRCT	// point-point attraction
			e = getEsmPP(d, n->area);	// e>0 and decreases monotonically with r
		#else	// point-face attraction
			vecdiv(r, d, nr);	// nr pointing from syt site to membrane grid
			cs=fabs(dotprod(n->e3, nr));
			e = getEsmPF(d, cs, n->area);	// e>0 and decreases monotonically with r
		#endif
			energy += DebyeVesCoeff * e;	// energy>0 and decreases monotonically with r
		}
	}
	return energy;	// energy>0 and decreases monotonically with r
}




void EnergySYT(void)
{
	int i;
	double r, invr;
	CHAIN *p;

#if X_PIP2
	double d;
#endif
	
	r=Nsytperring*Lc2b/Pidb;
	invr=1.0/r;
	
	Esyt = EsytES = 0;
	p = chain;
	while(p) {
		if(p->nseg>0) {	// energy reduction due to oligomerization
			if(p->closed) Esyt -= (p->nseg)*ESS;	// additional drop due to closure
			else Esyt -= (p->nseg - 1)*ESS;
			
#if OPENMP
#pragma omp parallel for private(i) reduction(+:Esyt,EsytES)
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
				Esyt += p->E[i];
			}
			
		}
		
		p = p->next;
	}
}



double getEmc(double z, double area)	// lift from adhesion dE>0
{
	double energy;

	if(z >0) energy = Eadh * (1 - exp(-z/Radh)) * area;
	else energy = 0;
	
	return energy;
}




void EnergyMEM(void)
{
	int i, j;
	MNODE *m;
	
	Emem = Emembend = Ememtens = 0;
	
#if OPENMP
#pragma omp parallel for private(i,j,m) reduction(+:Emem,Emembend,Ememtens)
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
			Emem += m->E;
			Emembend += m->Eb;
			Ememtens += m->Ea;
		}
	}
	
	//printf("A=%.4g, A0=%.4g, A00=%.4g\n", area, Nxm1Nym1*dAmemNode, Lx*Ly);
}



void EnergyVES(void)
{
	int i;
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
	Eves = Evesbend + Evestens;
}



void AllEnergy(void)
{
	Esys = 0;
	
	EnergySYT();
	Esys += Esyt;
	
	if(Membrane) {
		EnergyMEM();
		Esys += Emem;
	}
	
	if(Vesicle) {
		EnergyVES();
		Esys += Eves;
	}
}


