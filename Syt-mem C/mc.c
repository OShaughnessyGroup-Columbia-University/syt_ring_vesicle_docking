
void MCshiftchain(CHAIN *p, double dr[3])
{
	int i, j;
	CHAIN *q;
	
	q = p;
	while(q) {
		vecadd(q->rhead, dr, q->rhead);
		vecadd(q->rtail, dr, q->rtail);
		vecadd(q->rcenter, dr, q->rcenter);
		
#if OPENMP
#pragma omp parallel for private(i,j)
#endif
		for(i=0; i<(q->nseg); i++) {
			vecadd(q->r[i], dr, q->r[i]);
			vecadd(q->rc[i], dr, q->rc[i]);
			
			for(j=0; j<4; j++) {
				vecadd(q->rC2AB[i][j], dr, q->rC2AB[i][j]);
			}
			
			vecadd(q->rlinker[i], dr, q->rlinker[i]);
			
			vecadd(q->rc1[i], dr, q->rc1[i]);
			vecadd(q->rc2[i], dr, q->rc2[i]);
			
			vecadd(q->rsm[i], dr, q->rsm[i]);
		}
		q = q->next;
	}
}



void MCmove1(void)	// MC move syt
{
	double dr[3];
	double e1, e2, de;
	
	CopyChain(chain, chain_backup);	// make a copy
	e1 = MCEenergySYT() + MCEnergySMrep();	// old energy
	
	//for(i=0; i<3; i++) dr[i] = MCdr*(1-2*RND());
	dr[0] = dr[1] = 0;	dr[2] = MCdr*(1-2*RND());
	
	MCshiftchain(chain, dr);
	
	UpdateChain(chain);
	UpdateSytSitesAll(chain);
	
	e2 = MCEenergySYT() + MCEnergySMrep();	// new energy
	
	de = e2 - e1;
	//printf("%.4g\t", de);
	
	if(de>0) {
		if(RND() > exp(-de/temperature)) {	// reject move
			CopyChain(chain_backup, chain);	// revert change
			return;
		}
	}
	
	// accept new position
//	NormalizeChain(chain);
//	UpdateChain(chain);	// update Nclosed, Nopen, rcenter, rmax, rave after motion
//	UpdateSytSitesAll(chain);	// update syt-syt sites and syt-mnode sites
//	UpdateSytStretchAll(chain);	// update dr1, dr2 (run this after UpdateSytSites())
}



void MCmove2(void)	// MC scale membrane
{
	int i, j;
	double fct, dr[3];
	double e1, e2, de;
	MNODE *m;
	
	fct = 0.05*(1-2*RND());
	
	CopyChain(chain, chain_backup);	// make a copy
	CopyMemNode(mnode, mnode_backup);
	CopyMemFace(mface, mface_backup);
	
	e1 = MCEenergySYT() + MCEnergyMEM() + MCEnergySMrep();	// old energy
	
	// shift syt chain together
	dr[0] = dr[1] = 0;
	dr[2] = fct*(chain->rcenter[2]);
	
	MCshiftchain(chain, dr);
	
	UpdateChain(chain);
	UpdateSytSitesAll(chain);
	
	// shift membrane
#if OPENMP
#pragma omp parallel for private(i,j,m,dr)
#endif
	for(i=1; i<=Nxm1; i++) {
		for(j=1; j<=Nym1; j++) {
			m = &mnode[i][j];
		#if CRCBND
			if(m->circbnd == 1) continue;
		#endif
			
		#if MC_XYZ
			vecprod(m->e3, fct*(m->r[2]), dr);
		#else
			dr[0] = dr[1] = 0;
			dr[2] = fct*(m->r[2]);
		#endif
			vecadd(m->r, dr, m->r);	// update mnode.r
		}
	}
	
	GetMemNbrPositionsOnly(mnode);
	UpdateMemFaceNodePos(mnode, mface);
	UpdateMem(mnode, mface);
	
	e2 = MCEenergySYT() + MCEnergyMEM() + MCEnergySMrep();	// new energy
	
	de = e2 - e1;
	
	if(de>0) {
		if(RND() > exp(-de/temperature)) {	// reject move
			CopyChain(chain_backup, chain);	// revert change
			CopyMemNode(mnode_backup, mnode);
			CopyMemFace(mface_backup, mface);
			return;
		}
	}
	
	// accept new position
//	NormalizeChain(chain);
//	UpdateChain(chain);	// update Nclosed, Nopen, rcenter, rmax, rave after motion
//	UpdateSytSitesAll(chain);	// update syt-syt sites and syt-mnode sites
//	UpdateSytStretchAll(chain);	// update dr1, dr2 (run this after UpdateSytSites())
}



void MCmove3(void)	// MC Fourier modes
{
	int i, j, nmode;
	double ci, r, dr[3];
	double e1, e2, de;
	MNODE *m;
	
	nmode = (int)(MCnmode0*RND() + 0.99);	// pick a random mode
	ci = 0.2*(1-2*RND())/2;	// Fourier amplitude, 1/2 comes from 0.5*(1+cos(i*pi*r/rm))
	
	CopyChain(chain, chain_backup);	// make a copy
	CopyMemNode(mnode, mnode_backup);
	CopyMemFace(mface, mface_backup);
	
	e1 = MCEenergySYT() + MCEnergyMEM() + MCEnergySMrep();	// old energy
	
	// shift syt chain together
	r = chain->rave;
	if(r < Lxhalf) {
		dr[0] = dr[1] = 0;
		dr[2] = ci*(1+cos(nmode*Pi*r/Lxhalf));
		MCshiftchain(chain, dr);
		UpdateChain(chain);
		UpdateSytSitesAll(chain);
	}
	
	// shift membrane
#if OPENMP
#pragma omp parallel for private(i,j,m,dr)
#endif
	for(i=1; i<=Nxm1; i++) {
		for(j=1; j<=Nym1; j++) {
			m = &mnode[i][j];
		#if CRCBND
			if(m->circbnd == 1) continue;
		#endif
			
			r = sqrt(pow(m->r[0],2)+pow(m->r[1],2));
			
		#if MC_XYZ
			vecprod(m->e3, ci*(1+cos(nmode*Pi*r/Lxhalf)), dr);
		#else
			dr[0] = dr[1] = 0;
			dr[2] = ci*(1+cos(nmode*Pi*r/Lxhalf));
		#endif
		
			vecadd(m->r, dr, m->r);	// update mnode.r
		}
	}
	
	GetMemNbrPositionsOnly(mnode);
	UpdateMemFaceNodePos(mnode, mface);
	UpdateMem(mnode, mface);
	
	e2 = MCEenergySYT() + MCEnergyMEM() + MCEnergySMrep();	// new energy
	
	de = e2 - e1;
	
	if(de>0) {
		if(RND() > exp(-de/temperature)) {	// reject move
			CopyChain(chain_backup, chain);	// revert change
			CopyMemNode(mnode_backup, mnode);
			CopyMemFace(mface_backup, mface);
			return;
		}
	}
	
	// accept new position
//	NormalizeChain(chain);
//	UpdateChain(chain);	// update Nclosed, Nopen, rcenter, rmax, rave after motion
//	UpdateSytSitesAll(chain);	// update syt-syt sites and syt-mnode sites
//	UpdateSytStretchAll(chain);	// update dr1, dr2 (run this after UpdateSytSites())
}


