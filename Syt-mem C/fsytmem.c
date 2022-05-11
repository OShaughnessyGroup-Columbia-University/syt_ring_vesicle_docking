
double getFattrSytPIP2(double r)	// point-point attraction between syt and extra PIP2, f>0 for attraction
{
	double f;
	
	if(r<=r0_PIP2) f=0;
	else f=QsytQpip2_4pie0e/r*(LDebyeinv+1.0/r)*exp(-r/LDebye);
	
	return f;
}



double getFattrSytMemPP1(double r, double a)	// point-point attraction assuming a disc charge, f>0 for attraction
{
	double r02, sq, f;
	
	r02 = a/Pi;
	sq = sqrt(r02+r*r);
	/*
	the following comes from point charge being on the axis of a charged disc.
	potential = 2pi * alpha * L * [ exp(-z/L) - exp(-sqrt(z^2+R^2)/L) ]
	so force = -potential'
	since r!=z, this formula is incorrect!
	*/
	f = DebyeMemCoeff * Pidb * (exp(-r/LDebye) - r/sq*exp(-sq/LDebye));	// f>0 for attraction
	
	return f;
}



double getFattrSytMemPP(double r, double a)	// point-point attraction w/ cut-off distance, f>0 for attraction
{
	double rcut, invr, f;
	
	rcut = interpolate_func(r_cutoff_sm, a, amin_cutoff_sm, da_cutoff_sm, n_cutoff_smv);
	if(r <= rcut) f = 0;
	else {
		invr = 1.0/r;
		f = DebyeMemCoeff * invr * (invr + LDebyeinv) * exp(-r/LDebye);
	}
	
	return f;
}



double getEfieldNear(double r, double cs, double radius)	// point-disc force, near-field
{
	double radius2, radius3, radius5, radius7, exp_radius_L;
	double cs2, cs3, cs4, cs5, cs6, cs7, cs8;
	double c0, c1, c2, c3, c4, c5, c6, c7;
	double p1, p2, p3, p4, p5, p6, p7, p8;
	double r2, r3, r4, r5, r6, r7;
	double elec;
	
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
	//r8 = r7*r;
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
	
	// coefficients of series expansion of E-field: exp[-z/D]-z/sqrt(z^2+R^2)*exp[-sqrt(z^2+R^2)/D] = c0+c1*r+c2*r^2+...
	c0 = 1;
	c1 = -(radius + exp_radius_L*LDebye)/radius/LDebye;
	c2 = 0.5/LDebye2;
	c3 = (-radius3 + 3*exp_radius_L*LDebye2*(radius+LDebye))/6.0/radius3/LDebye3;
	c4 = -1/24.0/LDebye4;
	c5 = -(radius5 + 15*exp_radius_L*LDebye3*(radius2+3*radius*LDebye+3*LDebye2))/120.0/radius5/LDebye5;
	c6 = 1/720.0/LDebye6;
	c7 = (-radius7 + 105*exp_radius_L*LDebye4*(radius3+6*radius2*LDebye+15*radius*LDebye2+15*LDebye3))/5040.0/radius7/LDebye7;
	//c8 = 1.0/40320/LDebye8;
	//c9 = -(radius9 + 945*exp_radius_L*LDebye5*(radius4+10*radius3*LDebye+45*radius2*LDebye2+105*radius*LDbye3+105*LDebye4))/362880.0/radius9/LDebye9;
	//c10 = 1.0/3628800/LDebye10;
	
	//p0 = 1;
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
	
	//pot = c0'*p0 + c1'*p1*r + c2'*p2*r2 + c3'*p3*r3 + c4'*p4*r4 + c5'*p5*r5 + c6'*p6*r6;	// with different c1', c2'...
	//pot *= DebyeMemCoeff*Pidb*LDebye;
	// E ~ -d(Potential)/dr
	
	elec = c0*p1 + c1*r*p2 + c2*r2*p3 + c3*r3*p4 + c4*r4*p5 + c5*r5*p6 + c6*r6*p7 + c7*r7*p8;
	//elec += c8*r8*p9 + c9*r9*p10 + c10*r10*p11;
	
	elec *= Pidb;	// without the prefactor DebyeMemCoeff or DebyeVesCoeff
	//f = DebyeMemCoeff*elec;	// f>0 for attraction
	return elec;
}



double getEfieldFar(double r, double cs, double radius)	// point-disc force, far-field
{
	double radius2, radius4, radius6, exp_r_L;
	double cs2, cs4, sn2, sn4;
	double c1, c2, c3, c4, c5;
	double r2, r3, r4, r5, r6;
	double Z, Q, F;
	double elec, e1, e2;
	
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
	
	Z = Pi*radius2;	// without the DebyeMemCoeff term
	Q = -Pi/4*radius4;
	F = Pi/8*radius6;
	
	c1 = Z+Q/2/LDebye2*(cs2-1)+F/24/LDebye4*sn4;
	c2 = Q/2/LDebye*(3*cs2-1)-F/3/LDebye3*(sn2-5/32.0*sn4);
	c3 = Q/2*(3*cs2-1)+F/3/LDebye2*(1-6*sn2+45/8.0*sn4);
	c4 = F/8/LDebye*(35*cs4-30*cs2+3);
	c5 = F/8*(35*cs4-30*cs2+3);
	
	e1 = c1/r + c2/r2 + c3/r3 + c4/r4 + c5/r5;
	e2 = c1/r2 +2*c2/r3 + 3*c3/r4 + 4*c4/r5 + 5*c5/r6;	// derivative w/ respect to r
	elec = exp_r_L*(e1/LDebye + e2);	// E ~ -d(Potential)/dr
	//f = DebyeMemCoeff*elec;	// f>0 for attraction
	return elec;
}



double getFattrSytMemPF(double r, double cs, double a)	// point-face attraction
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
	
	f = DebyeMemCoeff*elec;	// f>0 for attraction
	return f;
}



void testfatt(void)
{
	int i, j;
	double a0, r0, rm;
	double dr, dth, r, th, cs, f;
	FILE *fp;
	
	a0=1.0;
	r0=sqrt(a0/Pi);
	rm=3*r0;
	dr=rm/60;
	dth=Pihalf/3.0;
	fp=fopen("fatt.dat", "w");
	fprintf(fp, "r\t0\t30\t60\t90\n");
	for(i=0; i<=60; i++) {	// for each r value
		r = dr*i;
		fprintf(fp, "%.4g", r);
		for(j=0; j<=3; j++) {
			th = j*dth;
			cs = cos(th);
			f = getFattrSytMemPF(r, cs, a0);
			fprintf(fp, "\t%.4g", f);
		}
		fprintf(fp, "\n");
	}
}




void SytMemAttract(MNODE **MN, CHAIN *p, int id, int i, int j)	// both i & j should be omp critical
{
	double r[3], d2, d, nr[3];
	double f, fv[3], dr[3], tau[3];
	MNODE *m;
	
#if !ATTRCT
	double cs;
#endif
	
	m = &MN[i][j];	// m should be omp critical
	
	vecsub(m->r, p->rsm[id], r);	// r pointing from syt site to membrane grid
	d2 = norm2(r);
	
	if(d2<RSMmax2) {
		d=max2(sqrt(d2), eps);
		vecdiv(r, d, nr);	// nr pointing from syt site to membrane grid
			
#if ATTRCT	// point-point attraction
		f=getFattrSytMemPP(d, m->area);	// f>0 for attraction
#else	// point-face attraction
		cs=fabs(dotprod(m->e3, nr));
		f=getFattrSytMemPF(d, cs, m->area);	// f>0 for attraction
#endif			

		vecprod(nr, f, fv);	// force on syt
		vecadd(p->f[id], fv, p->f[id]);	// attractive force
		vecadd(p->fsm[id], fv, p->fsm[id]);
		vecadd(p->fatt[id], fv, p->fatt[id]);
		
#if OPENMP
#pragma omp critical(getSMforce)	// share the same name with sytmemrepulse!
{
#endif
		vecsub(m->f, fv, m->f);	// be careful, use omp critical !!!
		TrackFsyt += fv[2];
		TrackFmem += -fv[2];
#if OPENMP
}
#endif
		
		vecsub(p->rsm[id], p->r[id], dr);
		crossprod(dr, fv, tau);
		vecadd(p->taum[id], tau, p->taum[id]);
	}
}



void SytMemRepulseN(MNODE **MN, CHAIN *p, int id, int i, int j)	// for each mnode[i][j]
{
	int k;
	double dr[3], d2, d;	//, sgn;
	double f, fv[3];
	MNODE *m;
	
	m = &MN[i][j];	// membrane node-m
	
	for(k=0; k<4; k++) {	// for C2AB bead-k
		vecsub(p->rC2AB[id][k], m->r, dr);	// dr is from mnode grid to centers of C2AB beads
		d2 = norm2(dr);
		
		if(d2<drSMrepulsion2) {
			d = max2(sqrt(d2), eps);
			//sgn = sign(dotprod(dr, m->e3));	// -f if bead is under membrane
			//f = sgn*KV*(drSMrepulsion - d);
			f = KV*(drSMrepulsion - d);
			vecprod(dr, f/d, fv);	// fv is the force on syt
			
			vecadd(p->f[id], fv, p->f[id]);
			vecadd(p->frep[id][k], fv, p->frep[id][k]);
		
#if OPENMP
#pragma omp critical(getSMforce)	// share the same name with sytmemattract!
{
#endif
			vecsub(m->f, fv, m->f);	// be careful, use omp critical !!!
			TrackFsyt += fv[2];
			TrackFmem -= fv[2];
#if OPENMP
}
#endif
		}	// end of if
	}	// end of k
}



void SytMemRepulseF(MNODE **MN, CHAIN *p, int id, int i)	// for each mface[i]
{
	int j, k, ix, jx;
	double dr[3], d2, d;	//, sgn;
	double f, fv[3], fv_3[3];
	MFACE *m;
	
	m = &mface[i];	// membrane face-m
	
	for(j=0; j<4; j++) {	// for C2AB bead-j
		vecsub(p->rC2AB[id][j], m->rc, dr);	// dr is from face center to syt center
		d2 = norm2(dr);
		
		if(d2<drSMrepulsion2) {
			d = max2(sqrt(d2), eps);
			//sgn = sign(dotprod(dr, m->nrm));	// -f if bead is under membrane
			//f = sgn*KV*(drSMrepulsion - d);	// drSMrepulsion=max2(Rc2abAVE, 0.5*max2(dX,dY)); the last part may cause problems!
			f = KV*(drSMrepulsion - d);	// *(m->area);
			vecprod(dr, f/d, fv);	// fv is the force on syt
			
			vecadd(p->f[id], fv, p->f[id]);
			vecadd(p->frep[id][j], fv, p->frep[id][j]);
	
			vecdiv(fv, 3.0, fv_3);
			
#if OPENMP
#pragma omp critical(getSMforce)	// share the same name with sytmemattract!
{
#endif
			for(k=0; k<3; k++) {	// for each vertex
				ix = m->inode[k][0];
				jx = m->inode[k][1];
				vecsub(MN[ix][jx].f, fv_3, MN[ix][jx].f);	// distribute to 3 mnode nodes
			}
			TrackFsyt += fv[2];
			TrackFmem -= fv[2];
#if OPENMP
}
#endif
		}	// end of if
	}	// end of k
}



void SytMemRepulseF2(MNODE **MN, CHAIN *p, int id, int i)	// for each mface[i]
{
	int j, k, ix, jx;
	double dr[3], d2, d, drn[3], cs;
	double f, fv[3], fv_3[3];
	MFACE *m;
	
	m = &mface[i];	// membrane face-m
	
	for(j=0; j<4; j++) {	// for C2AB bead-j
		vecsub(p->rC2AB[id][j], m->rc, dr);	// dr is from face center to syt center
		d2 = norm2(dr);
		
		if(d2<drSMrepulsion2) {
			d = max2(sqrt(d2), eps);
			
			normalize(dr, drn);
			cs = dotprod(drn, m->nrm);
			
			f = KV*(drSMrepulsion - d);	// drSMrepulsion=max2(Rc2abAVE, 0.5*max2(dX,dY)); the last part may cause problems!
			
			vecprod(dr, f/d*fabs(cs), fv);	// fv is the force on syt
			
			vecadd(p->f[id], fv, p->f[id]);
			vecadd(p->frep[id][j], fv, p->frep[id][j]);
			
			vecdiv(fv, 3.0, fv_3);
	
#if OPENMP
#pragma omp critical(getSMforce)	// share the same name with sytmemattract!
{
#endif
			for(k=0; k<3; k++) {	// for each vertex
				ix = m->inode[k][0];
				jx = m->inode[k][1];
				vecsub(MN[ix][jx].f, fv_3, MN[ix][jx].f);	// distribute to 3 mnode nodes
			}
			TrackFsyt += fv[2];
			TrackFmem -= fv[2];
#if OPENMP
}
#endif
		}	// end of if
	}	// end of k
}



void getSMforce(MNODE **MN, CHAIN *p, int id)
{
	int i, j, k;
	int imin, jmin, imax, jmax;
	int idx, idy, idf;
	double x, y;
	MZONE *z;
	
#if X_PIP2	// extra PIP2 molecules near syt
	int mi[3], mj[3];
	double d, nrm[3], f, fv[3], fv_3[3];
#endif
	
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
					
					if(attracton) SytMemAttract(MN, p, id, idx, idy);	// both ix & iy should be omp critical !!!
					
					SytMemRepulseN(MN, p, id, idx, idy);
				}
				
				///*
				for(k=0; k<(z->nface); k++) {	// all faces in zone
					idf = z->iface[k];
					SytMemRepulseF(MN, p, id, idf);
				}
				//*/
			}
		}	// end of j-loop
	}	// end of i-loop
	
#if X_PIP2	// extra PIP2 molecules near syt
	getSMdistanceplus(MN, p->rsm[id], &d, nrm, mi, mj);	// distance to membrane, normal (from mem to KKKK), membrane nodes ID
	f = getFattrSytPIP2(d);	// f>0 for attraction
	
	vecprod(nrm, -f, fv);	// force on syt
	
	vecadd(p->f[id], fv, p->f[id]);	// force on syt
	vecadd(p->fsm[id], fv, p->fsm[id]);
	vecadd(p->fatt[id], fv, p->fatt[id]);
	
	vecdiv(fv, -3, fv_3);	// force on 3 membrane nodes
	
	for(k=0; k<3; k++) {
		i = mi[k];
		j = mj[k];
		vecadd(MN[i][j].f, fv_3, MN[i][j].f);	// distribute force to 3 mnode nodes
	}
	
	TrackFsyt += fv[2];
	TrackFmem -= fv[2];
#endif
}



void FTSytMem(MNODE **MN, CHAIN *CH)	// OpenMP is inside each subroutine
{
	int i, n;
	CHAIN *p;

	p = CH;
	while(p) {
		n = p->nseg;
		
#if OPENMP
#pragma omp parallel for private(i) //reduction(+:TrackFsyt,TrackFmem)
#endif

		for(i=0; i<n; i++) getSMforce(MN, p, i);	// syt-membrane
		
		p = p->next;
	}
}

