
double getEv1s1(double rv[3], double e3[3], double a, double rs[3])	// energy between ves-rv and syt-rs
{
	double r[3], nr[3], d2, d, cs, energy;
	
	energy = 0;
	vecsub(rv, rs, r);
	if(fabs(r[0])<RSMmax && fabs(r[1])<RSMmax) {
		d2 = norm2(r);
		if(d2<RSMmax2) {
			d = max2(sqrt(d2), eps);
			vecdiv(r, d, nr);
			cs = fabs(dotprod(e3, nr));
			energy += DebyeVesCoeff*getEsmPF(d, cs, a);
		}
	}
	return energy;
}



double getEv1sall(CHAIN *CH, VNODE *m)
{
	int i;
	double energy;
	CHAIN *p;
	
	energy = 0;
	p = CH;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			energy += getEv1s1(m->r, m->e3, m->area, p->rsm[i]);
		}
		p = p->next;
	}
	return -energy;	// -energy<0
}



double getVesNodeE4(CHAIN *CH, VNODE *VN, VNODE *p, VNODE *q, VNODE *q1, VNODE *q2)	// calculate membrane energy at p+q+q1+q2
{	// note: need to recalculate area and H for all 4 nodes
	int i, j, k, nnbr;
	double energy, a, nrm[3], nv[3];
	VNODE *m, *m1, *m2;
	
#if VSTRCH
	GammaVes = Gamma + KAmem*max2(0, (AreaVes-AreaVes0)/AreaVes0);	// Ka=A*dGamma/dA --> dGamma=Ka*dA/A
#endif

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
			m1 = &VN[m->nbr[i]];
			m2 = &VN[m->nbr[j]];
			getarea3normal(m->r, m1->r, m2->r, &a, nrm);
			m->area += a/3.0;
			
			vecprod(nrm, a*sign(dotprod(m->e3,nrm)), nrm);	// area-weighted outward normal
			vecadd(nv, nrm, nv);
		}
		normalize(nv, m->e3);	// new e3
		
		getVesCurv(m);
		
		m->Eb = KBmemhalf*pow(2*(m->H)-KMspon,2)*(m->area);	// bending
		m->Ea = GammaVes*(m->area - dAvesNode);	// areal
		m->Es = getEv1sall(CH, m);	// Es<0
		m->E = m->Eb + m->Ea + m->Es;
		energy += m->E;
	}
	
	return energy;
}




void ReconnectVesNode(int inbr, int nnbr, VNODE *p, VNODE *q, VNODE *q1, VNODE *q2, double d12)
{
	int i, j, k, k1, k2, n, flg;
	double dr12[3], dr21[3];
	
	vecsub(q2->r, q1->r, dr12);	// vector from q1 to q2
	vecprod(dr12, -1.0, dr21);	// vector from q2 to q1
	
	//-------- adjust p's neighboring nodes (removing q from p's list) --------
	
	for(i=inbr; i<(nnbr-1); i++) {	// shift the rest of the list forward by 1
		j = i+1;
		p->nbr[i] = p->nbr[j];
		p->lnbr[i] = p->lnbr[j];
		veccopy(p->drnbr[j], p->drnbr[i]);
	}
	(p->nnbr)--;
	
	//-------- adjust q's neighboring nodes (removing p from q's list) --------
	
	k = -1;
	flg = 1;
	n = q->nnbr;
	for(i=0; i<n && flg; i++) {
		if(q->nbr[i] == p->id) {
			k = i;	// p is q's k-th nbr
			flg = 0;
		}
	}
	
	/*
	if(k<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectVesNode [q]: k=%d\n", k);
		exit(0);
	}
	*/
	
	for(i=k; i<(n-1); i++) {
		j = i+1;
		q->nbr[i] = q->nbr[j];	// shift the rest of the list forward by 1
		q->lnbr[i] = q->lnbr[j];
		veccopy(q->drnbr[j], q->drnbr[i]);
	}
	(q->nnbr)--;
	
	//-------- adjust q1's neighboring nodes (adding q2 to q1's list) --------
	
	k = k1 = k2 = -1;
	n = q1->nnbr;
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q1->nbr[i] == p->id) {
			k1 = i;	// p is q1's k1-th nbr
			flg = 0;
		}
	}
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q1->nbr[i] == q->id) {
			k2 = i;	// q is q1's k2-th nbr
			flg = 0;
		}
	}
	
	/*
	if(k1<0 || k2<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectVesNode [q1]: k1=%d, k2=%d\n", k1, k2);
		printf("p-id = %d, q1->id = %d\n", p->id, q1->id);
		printf("p's neighbors = ");
		for(i=0; i<p->nnbr; i++) printf("%d, ", p->nbr[i]);
		printf("\n");
		printf("q1's neighbors = ");
		for(i=0; i<q1->nnbr; i++) printf("%d, ", q1->nbr[i]);
		printf("\n");
		exit(0);
	}
	*/
	
	k = min2(k1, k2);
	k2 = max2(k1, k2);	// make sure k2>k1
	k1 = k;
	
	if(k2-k1>1) {	// k2=n-1 is the end, k1=0 is the start
		q1->nbr[n] = q2->id;	// append q2 to q1's nbr list
		q1->lnbr[n] = d12;	// distance between q1 & q2
		veccopy(dr12, q1->drnbr[n]);	// pointing from q1 to q2
	}
	else {	// k2=k1+1
		for(i=n; i>k2; i--) {	// shift [k2,end] backwards by 1
			j = i-1;	// j is i's cw neighbor
			q1->nbr[i] = q1->nbr[j];
			q1->lnbr[i] = q1->lnbr[j];
			veccopy(q1->drnbr[j], q1->drnbr[i]);
		}
		
		// q2 should be q1's new k2-th neighbor:
		q1->nbr[k2] = q2->id;
		q1->lnbr[k2] = d12;
		veccopy(dr12, q1->drnbr[k2]);
	}
	
	(q1->nnbr)++;
	
	//-------- adjust q2's neighboring nodes (adding q1 to q2's list) --------
	
	k = k1 = k2 = 0;
	n = q2->nnbr;
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q2->nbr[i] == p->id) {
			k1 = i;	// p is q2's k1-th nbr
			flg = 0;
		}
	}
	
	flg = 1;
	for(i=0; i<n && flg; i++) {
		if(q2->nbr[i] == q->id) {
			k2 = i;	// q is q2's k2-th nbr
			flg = 0;
		}
	}
	
	/*
	if(k1<0 || k2<0) {	// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		printf("Error in ReconnectVesNode [q2]: k1=%d, k2=%d\n", k1, k2);
		printf("p-id = %d, q2->id = %d\n", p->id, q2->id);
		printf("p's neighbors = ");
		for(i=0; i<p->nnbr; i++) printf("%d, ", p->nbr[i]);
		printf("\n");
		printf("q2's neighbors = ");
		for(i=0; i<q2->nnbr; i++) printf("%d, ", q2->nbr[i]);
		printf("\n");
		exit(0);
	}
	*/
	
	k = min2(k1, k2);
	k2 = max2(k1, k2);	// make sure k2>k1
	k1 = k;
	
	if(k2-k1>1) {	// k2=n-1 is the end, k1=0 is the start
		q2->nbr[n] = q1->id;	// append q1 to q2's nbr list
		q2->lnbr[n] = d12;	// distance between q1 & q2
		veccopy(dr21, q2->drnbr[n]);	// pointing from q2 to q1
	}
	else {	// k2=k1+1
		for(i=n; i>k2; i--) {	// shift [k2,end] backwards by 1
			j = i-1;	// j is i's cw neighbor
			q2->nbr[i] = q2->nbr[j];
			q2->lnbr[i] = q2->lnbr[j];
			veccopy(q2->drnbr[j], q2->drnbr[i]);
		}
		
		// q1 should be q2's new k2-th neighbor:
		q2->nbr[k2] = q1->id;
		q2->lnbr[k2] = d12;
		veccopy(dr21, q2->drnbr[k2]);
	}
	
	(q2->nnbr)++;
}





void RedefineVesFace(VFACE *VF, VNODE *p, VNODE *q, VNODE *q1, VNODE *q2)
{
	int i, j, id1[3], id2[3];
	int pnb1, pnb2, fid1, fid2;
	int qnb1, qnb2, q1nb, q2nb;
	int nb1, nb2;
	int nf, nfm1;
	double r12[3], r13[3], nrm[3];
	VFACE *f;
	
	id1[0] = id2[0] = p->id;	// old face1 is p-q-q1, face2 is p-q-q2
	id1[1] = id2[1] = q->id;	// id1[] are nodes in face1
	id1[2] = q1->id;	// id2[] are nodes in face2
	id2[2] = q2->id;
	
	//-------- find faces 1 & 2 in p's list --------
	fid1=fid2=-1;
	pnb1=pnb2=-1;
	for(i=0; i<(p->nface) && pnb1<0; i++) {	// for each neighboring face of p
		j = p->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id1,f)) {	// find face with vertices id1[3]
			pnb1=i;	// neighbor-id of face1 in p's list
			fid1=j;	// face-id of face1
		}
	}
	for(i=0; i<(p->nface) && pnb2<0; i++) {	// for each neighboring face of p
		j = p->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id2,f)) {	// find face with vertices id2[3]
			pnb2=i;	// neighbor-id of face2 in p's list
			fid2=j;	// face-id of face2
		}
	}
	
	//-------- find faces 1 & 2 in q's list --------
	qnb1=qnb2=-1;
	for(i=0; i<(q->nface) && qnb1<0; i++) {	// for each neighboring face of p
		j = q->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id1,f)) {	// find face with vertices id1[3]
			qnb1=i;	// neighbor-id of face1 in q's list
		}
	}
	for(i=0; i<(q->nface) && qnb2<0; i++) {	// for each neighboring face of p
		j = q->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id2,f)) {	// find face with vertices id2[3]
			qnb2=i;	// neighbor-id of face2 in q's list
		}
	}
	
	//-------- find face-1 in q1's list --------
	q1nb=-1;
	for(i=0; i<(q1->nface) && q1nb<0; i++) {	// for each neighboring face of q1
		j = q1->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id1,f)) {	// find face with vertices id1[3]
			q1nb=i;	// neighbor-id of face1 in q1's list
		}
	}
	
	//-------- find face-2 in q2's list --------
	q2nb=-1;
	for(i=0; i<(q2->nface) && q2nb<0; i++) {	// for each neighboring face of q2
		j = q2->iface[i];
		f = &VF[j];
		if(CompareVesNodeFace(id2,f)) {	// find face with vertices id2[3]
			q2nb=i;	// neighbor-id of face2 in q2's list
		}
	}
	
	if(fid1<0 || fid2<0 || pnb1<0 || pnb2<0 || q1nb<0 || q2nb<0) {
		printf("Error in RedefineVesFace: fid1=%d\tfid2=%d\tpnb1=%d\tpnb2=%d\tq1nb=%d\tq2nb=%d\n",
			fid1, fid2, pnb1, pnb2, q1nb, q2nb);
		exit(0);
	}
	
	//-------- re-define face fid1 --------
	f = &VF[fid1];
	f->inode[0] = p->id;	// new face1 is p-q1-q2
	f->inode[1] = q1->id;
	f->inode[2] = q2->id;
	
	veccopy(p->r, f->r1);
	veccopy(q1->r, f->r2);
	veccopy(q2->r, f->r3);
	vecave3(f->r1, f->r2, f->r3, f->rclab);
	vecsub(f->rclab, RCves, f->rcrc);
	
	vecsub(f->r1, f->r2, r12);
	vecsub(f->r1, f->r3, r13);
	crossprod(r12, r13, nrm);
	if(dotprod(nrm, f->rcrc)>0) f->nrmsign=1;	// r12 x r13 is along outward normal
	else f->nrmsign=-1;	// -r12 x r13 is along outward normal
	vecprod(nrm, f->nrmsign, f->nrm);
	
	//-------- re-define face fid2 --------
	f = &VF[fid2];
	f->inode[0] = q1->id;	// new face2 is q1-q-q2
	f->inode[1] = q->id;
	f->inode[2] = q2->id;
	
	veccopy(q1->r, f->r1);
	veccopy(q->r, f->r2);
	veccopy(q2->r, f->r3);
	vecave3(f->r1, f->r2, f->r3, f->rclab);
	vecsub(f->rclab, RCves, f->rcrc);
	
	vecsub(f->r1, f->r2, r12);
	vecsub(f->r1, f->r3, r13);
	crossprod(r12, r13, nrm);
	if(dotprod(nrm, f->rcrc)>0) f->nrmsign=1;	// r12 x r13 is along outward normal
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




void VesFlipE(int j, int nnbr, int id2, double d1, double d2, CHAIN *CH, VNODE *VN, VFACE *VF, VNODE *p, VNODE *q, VNODE *q1, VNODE *q2)
{
	int nnbrx, jx, k, idx;
	double e1, e2, de, prob;
	
	e1 = getVesNodeE4(CH,VN,p,q,q1,q2);	// energy before trial flip, openmp inside !!!!!!!!
	ReconnectVesNode(j,nnbr,p,q,q1,q2,d2);
	RedefineVesFace(VF,p,q,q1,q2);
	
	e2 = getVesNodeE4(CH,VN,p,q,q1,q2);	// energy after trial flip, openmp inside !!!!!!!!
	de = e2-e1;
	
	if(de>0) {	// Boltzmann factor
		prob = exp(-de/kBT);
			if(RND() > prob) {	// trial flip not accepted
				nnbrx = q1->nnbr;
				jx = -1;
				for(k=0; k<nnbrx && jx<0; k++) {
					idx = q1->nbr[k];
					if(idx==id2) jx = k;	// index of q2 in q1's nbr list
				}
				if(jx<0) {
					printf("Error in VesBondFlip(): cannot find q2 in q1's nbrs.\n");
					exit(0);
				}
				
				ReconnectVesNode(jx,nnbrx,q1,q2,q,p,d1);	// flip back
				RedefineVesFace(VF,q1,q2,q,p);
			}
		}
}




void VesBondFlip(CHAIN *CH, VNODE *VN, VFACE *VF, int status)	// flip bonds according to a tethering potential (see Kroll & Gompper, Science 1992)
{	//	status=0: noise=0; status=1: noise=1, don't use parallel !!!
	int i, j, j1, j2, nnbr, nnbrm1, flg;
	int id1, id2, id;
	double d1, d2;
	VNODE *p, *q, *q1, *q2;

/*	
#if OPENMP
#pragma omp parallel for private(i,j,j1,j2,nnbr,flg,id1,id2,id,d1,d2,p,q,q1,q2)
#endif
*/
	for(i=0; i<NVnode; i++) {
		p = &VN[i];
		nnbr = p->nnbr;
		nnbrm1 = nnbr-1;
		flg=1;
		if(nnbr > NnbrMin) {	// say, NnbrMin=4
			for(j=0; j<nnbr && flg; j++) {	// j-th neighbor of p
				id = p->nbr[j];
				q = &VN[id];
				
				if(i>id && q->nnbr > NnbrMin) {	// to avoid double-counting
					d1 = p->lnbr[j];
					j1 = (j+nnbrm1)%nnbr;	// cw neighbor of j
					j2 = (j+1)%nnbr;	// ccw neighbor of j
					
					id1 = p->nbr[j1];
					id2 = p->nbr[j2];
					q1 = &VN[id1];
					q2 = &VN[id2];
					
					d2 = distance(q1->r, q2->r);
					
					if(status==1) {	// with noise
						if(d2>LBvesmin && d2<LBvesmax && max2(q1->nnbr, q2->nnbr)<NnbrMaxm1) {
							VesFlipE(j, nnbr, id2, d1, d2, CH, VN, VF, p, q, q1, q2);	// write this in a separate function, otherwise cause problems in OPENMP
							nnbr = p->nnbr;	// update nnbr
							flg=0;
						}
					}
					else {	// without noise
						if(d2>LBvesmin && max2(q1->nnbr, q2->nnbr)<NnbrMaxm1 && d2<d1) {
							ReconnectVesNode(j,nnbr,p,q,q1,q2,d2);
							RedefineVesFace(VF,p,q,q1,q2);
							nnbr = p->nnbr;	// update nnbr
							flg=0;
						}
					}
				}
			}
		}
	}
}

