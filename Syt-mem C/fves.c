
void ResetVesForce(VNODE *VN)
{
	int i;
	VNODE *p;
	
	veczero(Fves);
	veczero(Tauves);
	
#if OPENMP
#pragma omp parallel for private(i,p)
#endif

	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		veczero(p->f);
		veczero(p->fext);
		p->Eb = p->Ea = p->E = 0;
	}
}



void VesLaplaceH0(VNODE *p)	// Belkin's Laplacian formula
{
	int i, j;
	double sum1, sum2, dr2;
	VFACE *vf;
	VNODE *vn[3];
	
	sum1 = 0;
	for(i=0; i<(p->nfaceplus); i++) {	// for each face around m
		vf = &vface[p->ifaceplus[i]];
		sum2 = 0;
		for(j=0; j<3; j++) {	// for each vertex of face-vf
			vn[j] = &vnode[vf->inode[j]];
			dr2 = distance2(p->r, vn[j]->r);
			sum2 += exp(-dr2/h_belkinx4)*(vn[j]->H - p->H);
		}
		sum1 += (vf->area)*sum2;
	}
	
	p->LaplaceH = sum1/Pih_belkin2x12;
}




void VesLaplaceH1(VNODE *p)	// Meyer's Laplacian formula
{
	int nnbr, nnbrm1;
	int k, k1, k2, ik, ik1, ik2;
	double drki[NnbrMax][3], drkk1[NnbrMax][3], drkk2[NnbrMax][3];
	double r2ki[NnbrMax], r2kk2[NnbrMax];
	double cota[NnbrMax], cotb[NnbrMax];
	double sum, Amix;
	VNODE *pk, *pk1, *pk2;
	
	nnbr = p->nnbr;
	nnbrm1 = nnbr-1;
	sum = Amix = 0;
	
	for(k=0; k<nnbr; k++) {	// for each neighbor k
		k1 = k-1;	// cw neighbor of k
		k2 = k+1;	// ccw neighbor of k
		if(k1<0) k1 = nnbrm1;
		if(k2>nnbrm1) k2 = 0;
		
		ik = p->nbr[k];	// k
		ik1 = p->nbr[k1];	// k1
		ik2 = p->nbr[k2];	// k2
		
		pk = &vnode[ik];
		pk1 = &vnode[ik1];
		pk2 = &vnode[ik2];
		
		vecsub(p->r, pk->r, drki[k]);	// drki is from k to i
		vecsub(pk1->r, pk->r, drkk1[k]);	// drkk1 is from k to k1
		vecsub(pk2->r, pk->r, drkk2[k]);	// drkk2 is from k to k2
		
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
		
		ik = p->nbr[k];	// k
		pk = &vnode[ik];
		
		sum += (pk->H - p->H) * (cotb[k1] + cota[k2]);	// alpha_k = b_k1, beta_k = a_k2, minus sign comes from laplace=-K in the paper
		
		Amix += getAmix_k(drki[k], drkk2[k], r2ki[k], r2kk2[k], r2ki[k2], cotb[k], cota[k2]);	// A_mix of triangle i-k-k2 (Meyer's original formula)
		//Amix += (cotb[k1] + cota[k2])*r2ki[k]/8.0;	// Meyer's formula w/o considering the shape of triangles
	}
	
	p->LaplaceH = 0.5*sum/Amix;
}




void VesLaplaceH2(VNODE *p)	// Vallet & Levy's Laplacian formula
{
	printf("Error: VesLaplaceH2 is not finished\n");
	exit(0);
}



void VesLaplaceH(VNODE *p)
{
#if !LAPLAC
	VesLaplaceH0(p);	// Belkin's formula
#else
	if(LAPLAC==1) VesLaplaceH1(p);	// Meyer's formula
	else VesLaplaceH2(p);	// Vallet-Levy's formula
#endif
}




void VesHelfrichForce(VNODE *VN)
{
	int i;
	double ptension, pcurv1, pcurv2, ptot, f, fv[3];
	VNODE *p;
	
#if VSTRCH
	GammaVes = Gamma + KAmem*max2(0, (AreaVes-AreaVes0)/AreaVes0);	// Ka=A*dGamma/dA --> dGamma=Ka*dA/A
#endif

	Pves=log(VolVes0/VolVes)/Beta;	// beta=-(dV/dP)/V, P=(1/beta)*ln(V0/V)
	
#if OPENMP
#pragma omp parallel for private(i,ptension,pcurv1,pcurv2,ptot,f,fv,p)
#endif

	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		VesLaplaceH(p);
		
		ptension = -2*GammaVes*(p->H);
		pcurv1 = KBvesx4*(p->H)*(pow(p->H,2)-(p->K));
		pcurv2 = KBvesdb*(p->LaplaceH);
		ptot = Pves + ptension + pcurv1 + pcurv2;
		f = ptot*(p->area);
		vecprod(p->e3, f, fv);
		vecadd(p->f, fv, p->f);
	}
}



void TotalVesForce(void)
{
	int i;
	for(i=0; i<NVnode; i++) vecadd(Fves, vnode[i].fext, Fves);	// total force on vesicle
}



void ForceVes(VNODE *VN)
{
	//ResetVesForce();
	
	VesHelfrichForce(VN);
	
	//TotalVesForce();	// get Fves
}


