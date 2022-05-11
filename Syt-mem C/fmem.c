
void ResetMemForce(MNODE **MN)
{
	int i, j;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,m)
#endif

	for(i=1; i<Nx; i++) {	// reset force
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			veczero(m->f);
			veczero(m->drb);
			m->Eb = m->Ea = m->Ec = m->E = 0;
		}
	}
	
}


void MemAdhesion(MNODE **MN)	// adhesion to carbon
{
	int i, j;
	double f;
	MNODE *m;
	
#if OPENMP
#pragma omp parallel for private(i,j,f,m)
#endif

	for(i=1; i<Nx; i++) {
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			if(m->r[2] > 0) {
				f = -Fadhmax*exp(-m->r[2]/Radh)*(m->area);
				m->f[2] += f;
			}
		}
	}
}



void MemLaplaceH0(MNODE **MN, MFACE *MF, MNODE *m)	// Belkin's Laplacian
{
	int i, j;
	double sum1, sum2, dr2;
	MFACE *mf;
	MNODE *mv;
	
	sum1 = 0;
	for(i=0; i<(m->nfaceplus); i++) {	// for each face around m
		mf = &MF[m->ifaceplus[i]];
		sum2 = 0;
		for(j=0; j<3; j++) {	// for each vertex of face-mf
			mv = &MN[mf->inode[j][0]][mf->inode[j][1]];
			dr2 = distance2(m->r, mv->r);
			sum2 += exp(-dr2/h_belkinx4)*(mv->H - m->H);
		}
		sum1 += (mf->area)*sum2;
	}
	
	m->LaplaceH = sum1/Pih_belkin2x12;	// minus sign? - NO.
}



double getAmix_k(double dr_ki[3], double dr_kk2[3], double r2_ki, double r2_kk2, double r2_k2i, double cotbk, double cotak2)
{
	double r2_sort[3], a;
	
	sort3(r2_ki, r2_kk2, r2_k2i, r2_sort);	// sort the first 3 numbers in descending order and save them in the 4th array
	
	if(r2_sort[0] <= r2_sort[1] + r2_sort[2]) a = (r2_ki*cotak2 + r2_k2i*cotbk) / 8.0;	// non-obtuse triangle
	else {	// obtuse triangle
		a = getarea2(dr_ki, dr_kk2);	// area of the entire triangle
		if(r2_kk2 > r2_ki + r2_k2i) a /= 2.0;	// angle at i is obtuse
		else a /= 4.0;	// angle at i is acute
	}	
	return a;
}



void MemLaplaceH1(MNODE **MN, MNODE *m)	// Meyer's Laplacian
{
	int nnbr, nnbrm1;
	int k, k1, k2, ik, jk, ik1, jk1, ik2, jk2;
	double drki[NnbrMax][3], drkk1[NnbrMax][3], drkk2[NnbrMax][3];
	double r2ki[NnbrMax], r2kk2[NnbrMax];
	double cota[NnbrMax], cotb[NnbrMax];
	double sum, Amix;
	MNODE *mk, *mk1, *mk2;
	
	nnbr = m->nnbr;
	nnbrm1 = nnbr-1;
	sum = Amix = 0;
	
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
		
		sum += (mk->H - m->H) * (cotb[k1] + cota[k2]);	// alpha_k = b_k1, beta_k = a_k2, minus sign comes from laplace=-K in the paper
		
		Amix += getAmix_k(drki[k], drkk2[k], r2ki[k], r2kk2[k], r2ki[k2], cotb[k], cota[k2]);	// A_mix of triangle i-k-k2 (Meyer's original formula)
		//Amix += (cotb[k1] + cota[k2])*r2ki[k]/8.0;	// Meyer's formula w/o considering the shape of triangles
	}
	
	m->LaplaceH = 0.5*sum/Amix;
}




void MemLaplaceH2(MNODE **MN, MNODE *m)	// 	// Vallet & Levy's Laplacian
{
	int nnbr, nnbrm1;
	int k, k1, k2, ik, jk, ik1, jk1, ik2, jk2;
	double drki[NnbrMax][3], drkk1[NnbrMax][3], drkk2[NnbrMax][3];
	double r2ki[NnbrMax], r2kk2[NnbrMax];
	double cota[NnbrMax], cotb[NnbrMax];
	double sum, Amix;
	MNODE *mk, *mk1, *mk2;
	
	nnbr = m->nnbr;
	nnbrm1 = nnbr-1;
	sum = Amix = 0;
	
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
		
		sum += (cotb[k1] + cota[k2]) * (mk->H - m->H) / sqrt(mk->Av);	// minus sign comes from laplace=-K in the paper
	}
	
	m->LaplaceH = 0.5*sum/sqrt(m->Av);
}



void MemLaplaceH(MNODE **MN, MFACE *MF, MNODE *m)
{
#if !LAPLAC
	MemLaplaceH0(MN, MF, m);	// Belkin's formula
#else
	if(LAPLAC==1) MemLaplaceH1(MN, m);	// Meyer's formula
	else MemLaplaceH2(MN, m);	// Vallet-Levy's formula
#endif
}



void HelfrichForceMem(MNODE **MN, MFACE *MF)
{
	int i, j;
	double Pext;
	double ptension, pcurv1, pcurv2, ptot, f, fv[3];
	MNODE *m;

#if BNDCND
	double Pbulk;
#endif
	
	if(Carbon) Pext=P0;
	else {
		#if BNDCND
			Pbulk=BKmod*Have/Hcell;	// pressure due to bulk modulus
			//Pbulk=log(Hcell/Have)/Beta;
			Pext=P0+Pbulk;
		#else
			//Pext=0;
			Pext=P0;
		#endif
	}
	
#if OPENMP
#pragma omp parallel for private(i,j,ptension,pcurv1,pcurv2,ptot,f,fv,m)
#endif

	for(i=1; i<Nx; i++) {
		for(j=1; j<Ny; j++) {
			m = &MN[i][j];
			MemLaplaceH(MN, MF, m);	// for 1<=i<Nx & 1<=j<Ny (not including Nx & Ny !)
			
			ptension = -2*Gamma*(m->H);
			pcurv1 = KBmemx4 * ((m->H)+KMsponhalf) * (pow(m->H,2) - KMsponhalf * (m->H)-(m->K));
			pcurv2 = KBmemdb * (m->LaplaceH);
			
			ptot = -Pext + ptension + pcurv1 + pcurv2;	// Pext has negative sign!!!
			f = ptot*(m->area);
			
			vecprod(m->e3, f, fv);
			vecadd(m->f, fv, m->f);
		}
	}
	
}



void ForceMem(MNODE **MN, MFACE *MF)
{
	//ResetMemForce();
	
	if(Carbon) MemAdhesion(MN);
	HelfrichForceMem(MN, MF);
}

