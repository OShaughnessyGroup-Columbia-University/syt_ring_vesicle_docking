// for Limit=1

void CopyChain(CHAIN *p, CHAIN *q)	// copy from p to q
{
	int i, j;
	CHAIN *s, *t;
	
	s = p;
	t = q;
	
	while(s) {
		t->id = s->id;
		t->nseg = s->nseg;
		t->closed = s->closed;
		t->Pclose = s->Pclose;
		
		veccopy(s->rhead, t->rhead);
		veccopy(s->rtail, t->rtail);
		veccopy(s->rcenter, t->rcenter);
		
		t->rmax = s->rmax;
		t->rave = s->rave;
		
#if OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(i=0; i<(s->nseg); i++) {
			veccopy(s->r[i], t->r[i]);
			veccopy(s->rc[i], t->rc[i]);
			
			for(j=0; j<4; j++) {
				veccopy(s->rC2AB[i][j], t->rC2AB[i][j]);
				veccopy(s->drC2AB[i][j], t->drC2AB[i][j]);
			}
			
			veccopy(s->rlinker[i], t->rlinker[i]);
			veccopy(s->rtmd[i], t->rtmd[i]);
			
			veccopy(s->rc1[i], t->rc1[i]);
			veccopy(s->rc2[i], t->rc2[i]);
			
			veccopy(s->dr1[i], t->dr1[i]);
			veccopy(s->dr2[i], t->dr2[i]);
			
			veccopy(s->rsm[i], t->rsm[i]);
			veccopy(s->fsm[i], t->fsm[i]);
			
			veccopy(s->nx[i], t->nx[i]);
			veccopy(s->ny[i], t->ny[i]);
			veccopy(s->nz[i], t->nz[i]);
			
			veccopy(s->nnxt[i], t->nnxt[i]);
			veccopy(s->nxa[i], t->nxa[i]);
			veccopy(s->nya[i], t->nya[i]);
			veccopy(s->nza[i], t->nza[i]);
			
			veccopy(s->nrmtmd[i], t->nrmtmd[i]);
			
			veccopy(s->f[i], t->f[i]);
			veccopy(s->fb[i], t->fb[i]);
			veccopy(s->ftmd[i], t->ftmd[i]);
			
			veccopy(s->taum[i], t->taum[i]);
			veccopy(s->taust[i], t->taust[i]);
			veccopy(s->taurp[i], t->taurp[i]);
			veccopy(s->taubd[i], t->taubd[i]);
			veccopy(s->tauts[i], t->tauts[i]);
			veccopy(s->tautmd[i], t->tautmd[i]);
			veccopy(s->tauorient[i], t->tauorient[i]);
			
			veccopy(s->fatt[i], t->fatt[i]);
			veccopy(s->fext[i], t->fext[i]);
			for(j=0; j<4; j++) veccopy(s->frep[i][j], t->frep[i][j]);
			
			t->Eb[i] = s->Eb[i];
			t->Em[i] = s->Em[i];
			t->Ev[i] = s->Ev[i];
			t->E[i] = s->E[i];
			t->Poff[i] = s->Poff[i];
		}
		
		s = s->next;
		t = t->next;
	}
	
	veccopy(p->ftot, q->ftot);
	veccopy(p->tautot, q->tautot);
}




void NewChainPropMono(CHAIN *p)	// for each monomer
{
	double r[3], r0, zmin;
	
	r0=Rv+DLThalf;
	zmin=sqrt(2)*DLThalf;
	
	p->id = Nchain;
	p->nseg = 1;
	p->closed = 0;
	p->Pclose = 0;

	if(Vesicle) {	// on vesicle
		do {
			do {	// within a unit sphere
				r[0] = 1-2*RND();
				r[1] = 1-2*RND();
				r[2] = 1-2*RND();
			} while (norm(r)>1);
			mapsphere(r0, r);
		} while(r[2]+Rv+Hv<=zmin);
		// r is in sphere's coord.
		vecadd(r, RCves, p->r[0]);	// position of the 1st subunit in lab coord
	}
	else {	// on target membrane
		r[0] = Lxhalf*(1-2*RND());
		r[1] = Lyhalf*(1-2*RND());
		r[2] = DLThalf;
		veccopy(r, p->r[0]);
	}
	
	veccopy(p->r[0], p->rhead);
	veccopy(p->r[0], p->rtail);
	veccopy(p->r[0], p->rcenter);
	p->rmax = p->rave = 0;
	
	veczero(p->ftot);
	
	p->E[0] = 0;
	p->Poff[0] = Poff0;
	
	if(Vesicle) {
		if(TMD) GetRndCoord(p->nx[0], p->ny[0], p->nz[0]);	// totally random directions
		else GetRndSphCoord(r, p->nx[0], p->ny[0], p->nz[0]);	// all syts bind to vesicle
	}
	else GetRndFlatCoord(p->nx[0], p->ny[0], p->nz[0]);	// random directions with z pionting up
	
	veccopy(p->nx[0], p->nnxt[0]);
	veccopy(p->nx[0], p->nxa[0]);
	veccopy(p->ny[0], p->nya[0]);
	veccopy(p->nz[0], p->nza[0]);
	
	if(Vesicle) {
		getLinkerpos(p);	// location of c2a-linker joint
		initTMDpos(p);	// location of tmd
	}
	
	p->next = NULL;
}


void AddNewChainMono(void)	// append to chain, run after Nchain++
{
	CHAIN *p, *q;
	
	if(!chain) {
		p = (CHAIN *)malloc(sizeof(CHAIN));
		NewChainPropMono(p);
		chain = p;
	}
	else {
		p = chain;
		while(p->next) p = p->next;	// p is the end of chain
		q = (CHAIN *)malloc(sizeof(CHAIN));
		NewChainPropMono(q);	// run after Nchain++
		p->next = q;
		//q->prev = p;
	}
}



void CheckChainMono(void)
{
	int i;
	CHAIN *p;
	
	printf("Chains:\n");
	p = chain;
	while(p) {
		printf("%d: (%.2f, %.2f)", p->id, p->r[0][0], p->r[0][1]);
		//printf("%d: ", p->id);
		for(i=0; i<(p->nseg); i++) printf("\t%d ", i+1);
		printf("\n");
		p = p->next;
	}
}


void InitSytMono(void)
{
	int i;
	for(i=0; i<NS1; i++) {
		Nchain=i+1;
		AddNewChainMono();	// run after Nchain++
	}
	//CheckChainMono();
}

