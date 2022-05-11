
// for 2nd-order BD, copy forces & torques

void AveSytFT2bd(CHAIN *p, CHAIN *q)	// update p with averaged force & torque from p & q
{
	int i, j;
	CHAIN *s, *t;
	
	s = p;
	t = q;
	while(s) {
		
#if OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(i=0; i<(s->nseg); i++) {
			vecave2(s->dr1[i], t->dr1[i], s->dr1[i]);
			vecave2(s->dr2[i], t->dr2[i], s->dr2[i]);
			
			vecave2(s->f[i], t->f[i], s->f[i]);
			vecave2(s->fb[i], t->fb[i], s->fb[i]);
			vecave2(s->fsm[i], t->fsm[i], s->fsm[i]);
			vecave2(s->ftmd[i], t->ftmd[i], s->ftmd[i]);
			
			vecave2(s->taum[i], t->taum[i], s->taum[i]);
			vecave2(s->taust[i], t->taust[i], s->taust[i]);
			vecave2(s->taurp[i], t->taurp[i], s->taurp[i]);
			vecave2(s->taubd[i], t->taubd[i], s->taubd[i]);
			vecave2(s->tauts[i], t->tauts[i], s->tauts[i]);
			vecave2(s->tautmd[i], t->tautmd[i], s->tautmd[i]);
			vecave2(s->tauorient[i], t->tauorient[i], s->tauorient[i]);
			
			for(j=0; j<4; j++) vecave2(s->frep[i][j], t->frep[i][j], s->frep[i][j]);
			vecave2(s->fatt[i], t->fatt[i], s->fatt[i]);
			vecave2(s->fext[i], t->fext[i], s->fext[i]);
		}
		
		s = s->next;
		t = t->next;
	}
	
	vecave2(p->ftot, q->ftot, p->ftot);
	vecave2(p->tautot, q->tautot, p->tautot);
}




void AveMemFT2bd(MNODE **p, MNODE **q)	// update p with averaged force & torque from p & q
{
	int i, j;
	MNODE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,j,s,t)
#endif
	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			s = &p[i][j];
			t = &q[i][j];
			vecave2(s->f, t->f, s->f);
		}
	}
}


void AveVesFT2bd(VNODE *p, VNODE *q)	// update p with averaged force & torque from p & q
{
	int i;
	VNODE *s, *t;
	
#if OPENMP
#pragma omp parallel for private(i,s,t)
#endif
	for(i=0; i<NVnode; i++) {
		s = &p[i];
		t = &q[i];
		vecave2(s->f, t->f, s->f);
		vecave2(s->fext, t->fext, s->fext);
	}
	
	vecave2(Fves, Fves_old, Fves);
	vecave2(Tauves, Fves_old, Tauves);
}



void MotionSyt2ndBD(CHAIN *CH)	// for test-move
{
	PrepSytMotion(CH);
	
#if FIXSYT
	if(Vesicle) ShiftFixedChain(CH);
	else ShiftFixedChainZ(CH);
	RotateChain(CH);
	if(FIXSYT==1) {	// allow stretching and bending about local z
		UpdateSytStretchAll(CH);	// update dr1, dr2 (run this after UpdateSytSites())
		StretchMotion(CH);
		if(LP>0) BendingMotion1D(CH);
		StericMotion(CH);
	}
#else
	StretchMotion(CH);
	if(LP>0) BendingMotion(CH);	// update syt's orientation (nx, ny, nz)
	StericMotion(CH);
	ShiftChain(CH);	// shift entire chain
#endif

	NormalizeChain(CH);
	
	//UpdateChain(CH);	// update Nclosed, Nopen, rcenter, rmax, rave after motion
	UpdateSytSitesAll(CH);	// update syt-syt sites and syt-mnode sites
	UpdateSytStretchAll(CH);	// update dr1, dr2 (run this after UpdateSytSites())
}



void MotionMem2ndBD(MNODE **MN, MFACE *MF)
{
	MemStretch(MN);	// motion
	
	GetMemNbrPositionsOnly(MN);	// update m->rnbr for RemapMem() after motion
	if(Carbon==1) MemCarbon(MN);
	
	GetMemNbrPositionsOnly(MN);
	MemNodeTether(MN);

	GetMemNbrPositionsOnly(MN);
	UpdateMemFaceNodePos(MN, MF);
}




void MotionVes2ndBD(VNODE *VN, VFACE *VF, MNODE **MN, CHAIN *CH)
{
	VesTranslation(VN, VF, CH);
	if(Membrane) VesRotation(VN, CH);
	
	GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
	VesNodeTether(VN);
	
	GetVesGeometry(VN, VF);	// update e3, lnbr, drnbr of each vnode
	if(Membrane) {
		VesShift(VN, MN, CH, 2);	// fixing (x,y), shifting z according to force
		//VesShift(VN, MN, CH, 3);	// shifting (x,y,z) according to force
	}
	else {
		VesShift(VN, MN, CH, 5);	// fixing (x,y,z)
		//VesSytRotate(VN, CH);
	}
	
	UpdateVesFaceCoord(VN, VF);	// after moving nodes, update face coord.
}



void MotionTMD2ndBD(CHAIN *CH)
{
	TMDmotion(CH);
	MapTMD(CH);
}
