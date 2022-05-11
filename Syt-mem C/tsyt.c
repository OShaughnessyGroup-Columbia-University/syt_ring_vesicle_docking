
void SytSytStretchTorque(CHAIN *CH)	// torque due to syt-syt stretching
{
	int i, n;
	double dr[3], tau[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n<=1) {
			p = p->next;
			continue;
		}
		
#if OPENMP
#pragma omp parallel for private(i,dr,tau)
#endif

		for(i=0; i<n; i++) {
			vecsub(p->dr2[i], p->dr1[i], dr);	// condense two forces
			//printf("i=%d,\tdr1=%.3g,\tdr2=%.3g\n", i, norm(p->dr1[i]), norm(p->dr2[i]));
			crossprod(p->nx[i], dr, tau);	// torque
			vecprod(tau, KSDLThalf, tau);
			vecadd(p->taust[i], tau, p->taust[i]);
			//printf("%d:\t(%.3g, %.3g, %.3g)\n", i, p->taust[i][0], p->taust[i][1], p->taust[i][2]);
			//printf("%d:\t(%.3g, %.3g, %.3g)\n", i, p->nx[i][0], p->nx[i][1], p->nx[i][2]);
			//printf("%d:\t(%.3g, %.3g, %.3g)\n", i, dr[0], dr[1], dr[2]);
		}
		
		p = p->next;
	}
}




void SytRepulsionTorque(CHAIN *CH)	// torque due to syt-mem & syt-ves repulsion
{
	int i, j, n;
	double tau[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		
#if OPENMP
#pragma omp parallel for private(i,j,tau)
#endif

		for(i=0; i<n; i++) {	// for each segment-i
			for(j=0; j<4; j++) {	// for each C2AB bead-j
				crossprod(p->drC2AB[i][j], p->frep[i][j], tau);
				vecadd(p->taurp[i], tau, p->taurp[i]);
			}
		}
		p = p->next;
	}
}




void SytSytOrientTorque(CHAIN *CH)	// torque due to syt-syt orientation
{
	int i, im1, n, nm1, istart;
	double vi[3], dotp, tmp[3];
	double vi1[3], vi2[3], a1, a2, a3, crs1[3], crs2[3], crs3[3];
	double nx[3], ny[3], nz[3], cs, sn, ny1[3];
	double tau1[3], tau2[3], tau12[3], tau3[3];
	CHAIN *p;
	
	p = CH;
	while(p) {
		n = p->nseg;
		if(n<=1) {
			p = p->next;
			continue;
		}
		nm1 = n-1;
		if(p->closed==0) istart=1;
		else istart=0;

#if OPENMP
#pragma omp parallel for private(i,im1,vi,dotp,tmp,vi1,vi2,\
a1,a2,a3,crs1,crs2,crs3,nx,ny,nz,cs,sn,ny1,tau1,tau2,tau12,tau3)
#endif

		for(i=istart; i<n; i++) {
			im1 = i-1;
			if(im1<0) im1=nm1;
			
			veccopy(p->nx[i], vi);
			
			veccopy(p->nx[im1], nx);	// based on (i-1)'s coordinates
			veccopy(p->ny[im1], ny);
			veccopy(p->nz[im1], nz);
			
			// bending
			
			if(K0==0) {	// intrinsically straight
				// bending angles will be a1=angle(v2, v1xy), a2=angle(v2, v1xz)
				dotp=dotprod(vi, nz);
				vecprod(nz, dotp, tmp);
				vecsub(vi, tmp, vi1);	// projection of current v1 to the x-y plane of i+1
				
				dotp=dotprod(vi, ny);	// ideal y-directions at i obtained from i+1
				vecprod(ny, dotp, tmp);
				vecsub(vi, tmp, vi2);	// projection of current v1 to the x-z plane of i+1
				
				a1=getanglePi(nx, vi1);	// 0<=a<=Pi
				a2=getanglePi(nx, vi2);
				crossprod(vi1, nx, crs1);	// cw chain
				crossprod(vi2, nx, crs2);
			}
			else {	// intrinsically curved
				RotMatrix3(sTH, cTH, nz, nx, nx);	// ideal x-direction at i obtained from i+1
				RotMatrix3(sTH, cTH, nz, ny, ny);	// ideal y-direction at i obtained from i+1
				
				// bending angles will be a1=angle(nx, vi1), a2=angle(nx, vi2)
				dotp=dotprod(vi, nz);
				vecprod(nz, dotp, tmp);
				vecsub(vi, tmp, vi1);	// projection of v1 to the x-y plane of i
				
				dotp=dotprod(vi, ny);
				vecprod(ny, dotp, tmp);
				vecsub(vi, tmp, vi2);	// projection of v1 to the x-z plane of i
				
				a1=getanglePi(nx, vi1);	// 0<=a<=Pi
				a2=getanglePi(nx, vi2);
				crossprod(vi1, nx, crs1);	// cw chain
				crossprod(vi2, nx, crs2);
			}
			
			normalize(crs1, crs1);
			normalize(crs2, crs2);
			vecprod(crs1, KBsyt*a1, tau1);	// radial bending
			vecprod(crs2, KBsytaxial*a2, tau2);	// bending in z (axis of helix)
			vecadd(tau1, tau2, tau12);
			
			vecadd(p->taubd[i], tau12, p->taubd[i]);	// total bending torque
			
#if OPENMP
#pragma omp critical(SytSytOrientTorque)	// share the same critical name as below
{
#endif
			vecsub(p->taubd[im1], tau12, p->taubd[im1]);	// be careful, use omp critical !!!
#if OPENMP
}
#endif
			
			// torsion
			cs=dotprod(vi, nx);
			sn=sqrt(1-min2(1,cs*cs));
			crossprod(vi, nx, crs3);	// rotate vi to nx
			normalize(crs3, crs3);
			RotMatrix3(sn, cs, crs3, p->ny[i], ny1);	// after aligning nx, ny[i] becomes ny1
			a3=getanglePi(ny1, ny);
			crossprod(ny1, ny, crs3);
			normalize(crs3, crs3);
			vecprod(crs3, Ktorsion*a3, tau3);
			vecadd(p->tauts[i], tau3, p->tauts[i]);
#if OPENMP
#pragma omp critical(SytSytOrientTorque)	// share the same critical name as below
{
#endif
			vecsub(p->tauts[im1], tau3, p->tauts[im1]);
#if OPENMP
}
#endif

		}
		
		p = p->next;
	}
}



void TotalSytTorque(CHAIN *CH)
{
	int i;
	CHAIN *p;
	
	p = CH;
	while(p) {
		
		veczero(p->tautot);
		
#if OPENMP
#pragma omp parallel for private(i)
#endif

		for(i=0; i < p->nseg; i++) {	// total self-orient torque
			if(Vesicle) vecadd6(p->taust[i], p->taurp[i], p->taubd[i], p->tauts[i], p->taum[i], p->tautmd[i], p->tauorient[i]);
			else vecadd5(p->taust[i], p->taurp[i], p->taubd[i], p->tauts[i], p->taum[i], p->tauorient[i]);
				
			vecadd(p->tautot, p->tauorient[i], p->tautot);	// total torque on chain
		}
		
		p = p->next;
	}
}



void TorqueSyt(CHAIN *CH)
{
	SytSytStretchTorque(CH);	// center-center torque
	SytRepulsionTorque(CH);	// repulsion torque
	SytSytOrientTorque(CH);	// orientational torque
	TotalSytTorque(CH);
}

