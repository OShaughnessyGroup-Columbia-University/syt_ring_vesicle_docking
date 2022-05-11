
double getEm1s1(double rm[3], double e3[3], double a, double rs[3])	// energy between membrane-rm and syt-rs
{
	double r[3], nr[3], d2, d, cs, energy;
	
	energy = 0;
	vecsub(rm, rs, r);
	if(fabs(r[0])<RSMmax && fabs(r[1])<RSMmax) {
		d2 = norm2(r);
		if(d2<RSMmax2) {
			d = max2(sqrt(d2), eps);
			vecdiv(r, d, nr);
			cs = fabs(dotprod(e3, nr));
			energy += DebyeMemCoeff*getEsmPF(d, cs, a);
		}
	}
	return energy;
}



double getEm1sall(MNODE *m)	// energy between membrane-m and all syts
{
	int i;
	double energy;
	CHAIN *p;
	
	energy = 0;
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			energy += getEm1s1(m->r, m->e3, m->area, p->rsm[i]);
		}
		p = p->next;
	}
	return -energy;	// -energy<0
}



double getMemNodeE4(MNODE **MN, MNODE *p, MNODE *q, MNODE *q1, MNODE *q2)	// calculate membrane energy at p+q+q1+q2
{	// need to recalculate e3, area and H for all 4 nodes
	int i, j, k, nnbr;
	double energy, a, nrm[3], nv[3];
	MNODE *m, *m1, *m2;
	
	energy = 0;
	
#if OPENMP
#pragma omp parallel for private(i,j,k,nnbr,a,nrm,nv,m,m1,m2) reduction(+:energy)
#endif

	for(k=0; k<4; k++) {	// update area and H for all 4 nodes
		if(k==0) m = p;
		else if(k==1) m = q;
		else if(k==2) m = q1;
		else m = q2;
		
		nnbr = m->nnbr;
		m->area = 0;
		veczero(nv);
		for(i=0; i<nnbr; i++) {
			j = (i+1)%nnbr;	// j is ccw neighbor of i
			m1 = &MN[m->nbr[i][0]][m->nbr[i][1]];
			m2 = &MN[m->nbr[j][0]][m->nbr[j][1]];
			getarea3normal(m->r, m1->r, m2->r, &a, nrm);	// get area & normal
			m->areanbr[i] = a/3.0;
			m->area += m->areanbr[i];
			
			vecprod(nrm, a*sign(dotprod(m->e3, nrm)), nrm);	// area-weighted outward normal
			vecadd(nv, nrm, nv);
		}
		normalize(nv, m->e3);	// new e3
		
		getMemCurv(MN, m);
		
		m->Eb = KBmemhalf*pow(2*(m->H)-KMspon,2)*(m->area);	// bending
		m->Ea = Gamma*(m->area - dAmemNode);	// areal
		m->Es = getEm1sall(m);	// Es<0
		m->E = m->Eb + m->Ea + m->Es;
		if(Carbon) {
			m->Ec = getEmc(m->r[2], m->area);	// carbon adhesion
			m->E += m->Ec;
		}
		energy += m->E;
	}
	
	return energy;
}




void ReconnectMemNode(int inbr, int nnbr, MNODE *p, MNODE *q, MNODE *q1, MNODE *q2, double d12)
{
	int i, j, k, k1, k2, n, flg;
	double dr12[3], dr21[3];
	
	vecsub(q2->r, q1->r, dr12);	// vector from q1 to q2
	vecprod(dr12, -1.0, dr21);	// vector from q2 to q1
	
	//-------- adjust p's neighboring nodes (removing q from p's list) --------
	
	for(i=inbr; i<(nnbr-1); i++) {	// shift the rest of the list forward by 1
		j = i+1;
		p->nbr[i][0] = p->nbr[j][0];
		p->nbr[i][1] = p->nbr[j][1];
		p->lnbr[i] = p->lnbr[j];
		veccopy(p->drnbr[j], p->drnbr[i]);
	}
	(p->nnbr)--;
	
	//-------- adjust q's neighboring nodes (removing p from q's list) --------
	
	k = -1;
	flg = 1;
	n = q->nnbr;
	for(i=0; i<n && flg; i++) {
		if(q->nbr[i][0] == p->id[0] && q->nbr[i][1] == p->id[1]) {
			k = i;	// p is q's k-th nbr
			flg = 0;
		}
	}
	
//	/*
	if(k<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectMemNode [q]: k=%d\n", k);
		exit(0);
	}
//	*/
	
	for(i=k; i<(n-1); i++) {
		j = i+1;
		q->nbr[i][0] = q->nbr[j][0];	// shift the rest of the list forward by 1
		q->nbr[i][1] = q->nbr[j][1];
		q->lnbr[i] = q->lnbr[j];
		veccopy(q->drnbr[j], q->drnbr[i]);
	}
	(q->nnbr)--;
	
	//-------- adjust q1's neighboring nodes (adding q2 to q1's list) --------
	
	k1 = k2 = -1;
	n = q1->nnbr;
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q1->nbr[i][0] == p->id[0] && q1->nbr[i][1] == p->id[1]) {
			k1 = i;	// p is q1's k1-th nbr
			flg = 0;
		}
	}
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q1->nbr[i][0] == q->id[0] && q1->nbr[i][1] == q->id[1]) {
			k2 = i;	// q is q1's k2-th nbr
			flg = 0;
		}
	}
	
//	/*
	if(k1<0 || k2<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectMemNode [q1]: k1=%d, k2=%d\n", k1, k2);
		printf("p->id = (%d, %d)\n", p->id[0], p->id[1]);
		printf("q->id = (%d, %d)\n", q->id[0], q->id[1]);
		printf("q1->id = (%d, %d)\n", q1->id[0], q1->id[1]);
		printf("p's neighbors = ");
		for(i=0; i<p->nnbr; i++) printf("(%d, %d), ", p->nbr[i][0], p->nbr[i][1]);
		printf("\n");
		printf("q's neighbors = ");
		for(i=0; i<q->nnbr; i++) printf("(%d, %d), ", q->nbr[i][0], q->nbr[i][1]);
		printf("\n");
		printf("q1's neighbors = ");
		for(i=0; i<q1->nnbr; i++) printf("(%d, %d), ", q1->nbr[i][0], q1->nbr[i][1]);
		printf("\n");
		exit(0);
	}
//	*/
	
	k = min2(k1, k2);
	k2 = max2(k1, k2);	// make sure k2>k1
	k1 = k;
	
	if(k2-k1>1) {	// k2=n-1 is the end, k1=0 is the start
		q1->nbr[n][0] = q2->id[0];	// append q2 to q1's nbr list
		q1->nbr[n][1] = q2->id[1];
		q1->lnbr[n] = d12;	// distance between q1 & q2
		veccopy(dr12, q1->drnbr[n]);	// pointing from q1 to q2
	}
	else {	// k2=k1+1
		for(i=n; i>k2; i--) {	// shift [k2,end] backwards by 1
			j = i-1;	// j is i's cw neighbor
			q1->nbr[i][0] = q1->nbr[j][0];
			q1->nbr[i][1] = q1->nbr[j][1];
			q1->lnbr[i] = q1->lnbr[j];
			veccopy(q1->drnbr[j], q1->drnbr[i]);
		}
		
		// q2 should be q1's new k2-th neighbor:
		q1->nbr[k2][0] = q2->id[0];
		q1->nbr[k2][1] = q2->id[1];
		q1->lnbr[k2] = d12;
		veccopy(dr12, q1->drnbr[k2]);
	}
	
	(q1->nnbr)++;
	
	//-------- adjust q2's neighboring nodes (adding q1 to q2's list) --------
	
	k1 = k2 = -1;
	n = q2->nnbr;
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q2->nbr[i][0] == p->id[0] && q2->nbr[i][1] == p->id[1]) {
			k1 = i;	// p is q2's k1-th nbr
			flg = 0;
		}
	}
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q2->nbr[i][0] == q->id[0] && q2->nbr[i][1] == q->id[1]) {
			k2 = i;	// q is q2's k2-th nbr
			flg = 0;
		}
	}
	
//	/*
	if(k1<0 || k2<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectMemNode [q2]: k1=%d, k2=%d\n", k1, k2);
		printf("p-id = (%d, %d), q2->id = (%d, %d)\n", p->id[0], p->id[1], q2->id[0], q2->id[1]);
		printf("p's neighbors = ");
		for(i=0; i<p->nnbr; i++) printf("(%d, %d), ", p->nbr[i][0], p->nbr[i][1]);
		printf("\n");
		printf("q2's neighbors = ");
		for(i=0; i<q2->nnbr; i++) printf("(%d, %d), ", q2->nbr[i][0], q2->nbr[i][1]);
		printf("\n");
		exit(0);
	}
//	*/
	
	k = min2(k1, k2);
	k2 = max2(k1, k2);	// make sure k2>k1
	k1 = k;
	
	if(k2-k1>1) {	// k2=n-1 is the end, k1=0 is the start
		q2->nbr[n][0] = q1->id[0];	// append q1 to q2's nbr list
		q2->nbr[n][1] = q1->id[1];
		q2->lnbr[n] = d12;	// distance between q1 & q2
		veccopy(dr21, q2->drnbr[n]);	// pointing from q2 to q1
	}
	else {	// k2=k1+1
		for(i=n; i>k2; i--) {	// shift [k2,end] backwards by 1
			j = i-1;	// j is i's cw neighbor
			q2->nbr[i][0] = q2->nbr[j][0];
			q2->nbr[i][1] = q2->nbr[j][1];
			q2->lnbr[i] = q2->lnbr[j];
			veccopy(q2->drnbr[j], q2->drnbr[i]);
		}
		
		// q1 should be q2's new k2-th neighbor:
		q2->nbr[k2][0] = q1->id[0];
		q2->nbr[k2][1] = q1->id[1];
		q2->lnbr[k2] = d12;
		veccopy(dr21, q2->drnbr[k2]);
	}
	
	(q2->nnbr)++;
}





void RedefineMemFace(MNODE *p, MNODE *q, MNODE *q1, MNODE *q2)
{
	int i, j, id1[3][2], id2[3][2];
	int pnb1, pnb2, fid1, fid2;
	int qnb1, qnb2, q1nb, q2nb;
	int nb1, nb2;
	int nf, nfm1;
	double r12[3], r13[3], nrm[3];
	MFACE *f;
	
	for(i=0; i<2; i++) {
		id1[0][i] = id2[0][i] = p->id[i];	// old face1 is p-q-q1, face2 is p-q-q2
		id1[1][i] = id2[1][i] = q->id[i];	// id1[] are nodes in face1
		id1[2][i] = q1->id[i];	// id2[] are nodes in face2
		id2[2][i] = q2->id[i];
	}
	
	//-------- find faces 1 & 2 in p's list --------
	fid1=fid2=-1;
	pnb1=pnb2=-1;
	for(i=0; i<(p->nface) && pnb1<0; i++) {	// for each neighboring face of p
		j = p->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id1,f)) {	// find face with vertices id1[3]
			pnb1=i;	// neighbor-id of face1 in p's list
			fid1=j;	// face-id of face1
		}
	}
	for(i=0; i<(p->nface) && pnb2<0; i++) {	// for each neighboring face of p
		j = p->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id2,f)) {	// find face with vertices id2[3]
			pnb2=i;	// neighbor-id of face2 in p's list
			fid2=j;	// face-id of face2
		}
	}
	
	//-------- find faces 1 & 2 in q's list --------
	qnb1=qnb2=-1;
	for(i=0; i<(q->nface) && qnb1<0; i++) {	// for each neighboring face of p
		j = q->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id1,f)) {	// find face with vertices id1[3]
			qnb1=i;	// neighbor-id of face1 in q's list
		}
	}
	for(i=0; i<(q->nface) && qnb2<0; i++) {	// for each neighboring face of p
		j = q->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id2,f)) {	// find face with vertices id2[3]
			qnb2=i;	// neighbor-id of face2 in q's list
		}
	}
	
	//-------- find face-1 in q1's list --------
	q1nb=-1;
	for(i=0; i<(q1->nface) && q1nb<0; i++) {	// for each neighboring face of q1
		j = q1->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id1,f)) {	// find face with vertices id1[3]
			q1nb=i;	// neighbor-id of face1 in q1's list
		}
	}
	
	//-------- find face-2 in q2's list --------
	q2nb=-1;
	for(i=0; i<(q2->nface) && q2nb<0; i++) {	// for each neighboring face of q2
		j = q2->iface[i];
		f = &mface[j];
		if(CompareMemNodeFace(id2,f)) {	// find face with vertices id2[3]
			q2nb=i;	// neighbor-id of face2 in q2's list
		}
	}
	
	if(fid1<0 || fid2<0 || pnb1<0 || pnb2<0 || q1nb<0 || q2nb<0) {
		printf("Error in RedefineMemFace: fid1=%d\tfid2=%d\tpnb1=%d\tpnb2=%d\tq1nb=%d\tq2nb=%d\n",
			fid1, fid2, pnb1, pnb2, q1nb, q2nb);
		
		printf("p=(%d, %d)\n", p->id[0], p->id[1]);
		printf("id1=p(%d, %d), q(%d, %d), q1(%d, %d)\n", id1[0][0], id1[0][1], id1[1][0], id1[1][1], id1[2][0], id1[2][1]);
		printf("id2=p(%d, %d), q(%d, %d), q2(%d, %d)\n", id2[0][0], id2[0][1], id2[1][0], id2[1][1], id2[2][0], id2[2][1]);
		
		printf("q1's faces:\n");
		for(i=0; i<(q1->nface); i++) printf("\t%d", q1->iface[i]);
		printf("\n");
		exit(0);
	}
	
	//-------- re-define face fid1 --------
	f = &mface[fid1];
	for(i=0; i<2; i++) {
		f->inode[0][i] = p->id[i];	// new face1 is p-q1-q2
		f->inode[1][i] = q1->id[i];
		f->inode[2][i] = q2->id[i];
	}
	
	veccopy(p->r, f->r1);
	veccopy(q1->r, f->r2);
	veccopy(q2->r, f->r3);
	vecave3(f->r1, f->r2, f->r3, f->rc);
	
	vecsub(f->r1, f->r2, r12);
	vecsub(f->r1, f->r3, r13);
	crossprod(r12, r13, nrm);
	if(nrm[2]>0) f->nrmsign=1;	// r12 x r13 is along outward normal
	else f->nrmsign=-1;	// -r12 x r13 is along outward normal
	vecprod(nrm, f->nrmsign, f->nrm);
	
	//-------- re-define face fid2 --------
	f = &mface[fid2];
	for(i=0; i<2; i++) {
		f->inode[0][i] = q1->id[i];	// new face2 is q1-q-q2
		f->inode[1][i] = q->id[i];
		f->inode[2][i] = q2->id[i];
	}
	
	veccopy(q1->r, f->r1);
	veccopy(q->r, f->r2);
	veccopy(q2->r, f->r3);
	vecave3(f->r1, f->r2, f->r3, f->rc);
	
	vecsub(f->r1, f->r2, r12);
	vecsub(f->r1, f->r3, r13);
	crossprod(r12, r13, nrm);
	if(nrm[2]>0) f->nrmsign=1;	// r12 x r13 is along outward normal
	else f->nrmsign=-1;	// -r12 x r13 is along outward normal
	vecprod(nrm, f->nrmsign, f->nrm);
	
	//-------- update pnb1, pnb2 from p's list --------
	nb1 = min2(pnb1, pnb2);	// nb1 is the low index
	nb2 = max2(pnb1, pnb2);	// nb2 is the high index
	
	p->iface[nb1] = fid1;	// define nb1 as new face1 (p-q1-q2)
	nfm1 = p->nface - 1;
	for(i=nb2; i<nfm1; i++) p->iface[i] = p->iface[i+1];	// remove nb2
	(p->nface)--;
	
	//-------- update qnb1, qnb2 from q's list --------
	nb1 = min2(qnb1, qnb2);	// nb1 is the low index
	nb2 = max2(qnb1, qnb2);	// nb2 is the high index
	
	q->iface[nb1] = fid2;	// define nb1 as new face2 (q1-q-q2)
	nfm1 = q->nface - 1;
	for(i=nb2; i<nfm1; i++) q->iface[i] = q->iface[i+1];	// remove qnb1
	(q->nface)--;
	
	//-------- update q1nb in q1's list --------
	q1->iface[q1nb] = fid2;	// define q1nb as new face2 (q1-q-q2)
	nf = q1->nface;
	for(i=nf; i>q1nb+1; i--) q1->iface[i] = q1->iface[i-1];	// make room for fid1
	q1->iface[q1nb+1] = fid1;	// q1nb+1 is new face1 (p-q1-q2)
	(q1->nface)++;
	
	//-------- update q2nb in q2's list --------
	q2->iface[q2nb] = fid1;	// define q2nb as new face1 (p-q1-q2)
	nf = q2->nface;
	for(i=nf; i>q2nb+1; i--) q2->iface[i] = q2->iface[i-1];	// make room for fid2
	q2->iface[q2nb+1] = fid2;	// q2nb+1 is new face2 (q1-q-q2)
	(q2->nface)++;
}



void MemFlipE(MNODE **MN, int k, int nnbr, int id2[2], double d1, double d2, MNODE *p, MNODE *q, MNODE *q1, MNODE *q2)
{
	int kx, l, nnbrx;
	int idx[2];
	double e1, e2, de, prob;
	
	e1 = getMemNodeE4(MN,p,q,q1,q2);	// energy before trial flip, openmp inside !!!!!!!!
	
	ReconnectMemNode(k,nnbr,p,q,q1,q2,d2);	// trial flip
	RedefineMemFace(p,q,q1,q2);
	
	e2 = getMemNodeE4(MN,p,q,q1,q2);	// energy after trial flip, openmp inside !!!!!!!!
	de = e2-e1;
	
	if(de>0) {	// Boltzmann factor
		prob = exp(-de/kBT);
		if(RND() > prob) {	// trial flip not accepted, make sure(RNDMP() is used if in OPENMP) !!!!!!!!!!!!!!
			nnbrx = q1->nnbr;
			kx = -1;
			for(l=0; l<nnbrx && kx<0; l++) {
				idx[0] = q1->nbr[l][0];
				idx[1] = q1->nbr[l][1];
				if(idx[0]==id2[0] && idx[1]==id2[1]) kx = l;	// index of q2 in q1's nbr list
			}
			if(kx<0) {
				printf("Error in MemBondFlip(): cannot find q2 in q1's nbrs.\n");
				exit(0);
			}
			
			ReconnectMemNode(kx,nnbrx,q1,q2,q,p,d1);	// flip back
			RedefineMemFace(q1,q2,q,p);
		}
	}
}



void MemBondFlip(MNODE **MN, int status)	// flip bonds according to a tethering potential (see Kroll & Gompper, Science 1992)
{	//	flip according to energy. don't use parallel !!!
	int i, j, k, k1, k2, nnbr, nnbrm1, flg;
	int id1[2], id2[2], id[2];
	double d1, d2;
	MNODE *p, *q, *q1, *q2;

/*	
#if OPENMP
#pragma omp parallel for private(i,j,k,k1,k2,nnbr,nnbrm1,flg,id1,id2,id,d1,d2,p,q,q1,q2)\
schedule(dynamic,2)
#endif
*/
	for(i=2; i<Nxm1; i++) {
		for(j=2; j<Nym1; j++) {
			p = &MN[i][j];
			nnbr = p->nnbr;
			nnbrm1 = nnbr-1;
			flg=1;
			if(nnbr > NnbrMin) {	// say, NnbrMin=4
				for(k=0; k<nnbr && flg; k++) {	// k-th neighbor of p
					id[0] = p->nbr[k][0];
					id[1] = p->nbr[k][1];
					q = &MN[id[0]][id[1]];
					
					if(i>=id[0] && q->nnbr > NnbrMin) {	// to avoid double-counting
					//if(q->nnbr > NnbrMin) {
						d1 = p->lnbr[k];
						k1 = (k+nnbrm1)%nnbr;	// cw neighbor of k
						k2 = (k+1)%nnbr;	// ccw neighbor of k

						id1[0] = p->nbr[k1][0];
						id1[1] = p->nbr[k1][1];
						id2[0] = p->nbr[k2][0];
						id2[1] = p->nbr[k2][1];
						q1 = &MN[id1[0]][id1[1]];
						q2 = &MN[id2[0]][id2[1]];
						
						d2 = distance(q1->r, q2->r);
						
						if(status==1) {	// flip with thermal
							if(d2>LBmemmin && d2<LBmemmax && max2(q1->nnbr, q2->nnbr)<NnbrMaxm1) {
								MemFlipE(MN, k, nnbr, id2, d1, d2, p, q, q1, q2);	// write this in a separate funcion, otherwise cause problems in OPENMP
								nnbr = p->nnbr;	// update nnbr
								flg=0;
							}
						}
						else {	// flip without thermal
							if(d2>LBmemmin && max2(q1->nnbr, q2->nnbr)<NnbrMaxm1 && d2<d1) {
								ReconnectMemNode(k,nnbr,p,q,q1,q2,d2);	// flip to short bond
								RedefineMemFace(p,q,q1,q2);
								nnbr = p->nnbr;	// update nnbr
								flg=0;
							}
						}
					}
				}
			}
		}
	}
}


void FlipMem(MNODE **MN)
{
	if(Cntflip_mem++ >= Nflip_mem) {
		GetMemNbrPositions(MN);
		if(noise) MemBondFlip(MN,1);	// requires p->lnbr
		else MemBondFlip(MN,0);
		Cntflip_mem = 0;
	}
}


