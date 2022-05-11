
void SaveData(void)
{
	int i, j;
	char fn[50];
	MNODE *m;
	MFACE *mf;
	CHAIN *p;
	VNODE *n;
	VFACE *f;
	FILE *fp;
	
	sprintf(fn, "P%.4u.dat", Nfile);
	fp=fopen(fn, "w");
	fprintf(fp, "%.4g\t%.4g\t%.4g\t", t, DLT, Rtmd);
	fprintf(fp, "// t DLT Rtmd\n\n");
	
	fprintf(fp, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t// Lc2a, Rc2a, Lc2b, Rc2b, Ac2a, Asm\n", 
		Lc2a, Rc2a, Lc2b, Rc2b, Ac2a, Asm);
	fprintf(fp, "%.4g\t%.4g\t%.4g\t// C2A1\n", RCc2a1[0], RCc2a1[1], RCc2a1[2]);
	fprintf(fp, "%.4g\t%.4g\t%.4g\t// C2A2\n", RCc2a2[0], RCc2a2[1], RCc2a2[2]);
	fprintf(fp, "%.4g\t%.4g\t%.4g\t// C2B1\n", RCc2b1[0], RCc2b1[1], RCc2b1[2]);
	fprintf(fp, "%.4g\t%.4g\t%.4g\t// C2B2\n", RCc2b2[0], RCc2b2[1], RCc2b2[2]);
	fprintf(fp, "%.4g\t%.4g\t%.4g\t// Binding site\n\n", RCsm[0], RCsm[1], RCsm[2]);
	
	// syt
	fprintf(fp, "%d\t%d\t%d\t// Syt: Nchain, Nopen, Nclosed\n\n", Nchain, Nopen, Nclosed);
	
	p = chain;
	while(p) {
		fprintf(fp, "%d\t%d\t%d", p->id, p->nseg, p->closed);
		fprintf(fp, "\t// SytID\tNseg\tClosed\n");
		for(i=0; i<p->nseg; i++) {
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->r[i][0], p->r[i][1], p->r[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->nx[i][0], p->nx[i][1], p->nx[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->nz[i][0], p->nz[i][1], p->nz[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->rlinker[i][0], p->rlinker[i][1], p->rlinker[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\n", p->rtmd[i][0], p->rtmd[i][1], p->rtmd[i][2]);
		}
		p = p->next;
	}
	
	// membrane
	if(Membrane) {
		GetMemMaxMin(mnode);
		
		fprintf(fp, "\n%.4g\t%.4g\t%d\t%d\t%d\t", Lx, Ly, Nx, Ny, NMface);
		fprintf(fp, "// Lx Ly Nx Ny NMface\n");
		
		fprintf(fp, "%.4g\t%.4g\t%.4g\t%.4g\t", Zmin, Zmax, Curvmin, Curvmax);
		fprintf(fp, "// Mem: Zmin, Zmax, Curvmin, Curvmax\n");
		
		for(i=0; i<=Nx; i++) {
			for(j=0; j<=Ny; j++) {
				m = &mnode[i][j];
				fprintf(fp, "%.4g\t%.4g\t%.4g\t%.4g\n", m->r[0], m->r[1], m->r[2], m->H);
			}
		}
		for(i=0; i<NMface; i++) {
			mf = &mface[i];
			fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, mf->inode[0][0], mf->inode[0][1], 
				mf->inode[1][0], mf->inode[1][1], mf->inode[2][0], mf->inode[2][1]);
		}
	}
	
	// vesicle
	if(Vesicle) {
		fprintf(fp, "\n%d\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t", NVnode, NVedge, NVface, Rv, RCves[0], RCves[1], RCves[2]);
		fprintf(fp, "// Ves: NVnode, NVedge, NVface, Rv, RC[3]\n");
		for(i=0; i<NVnode; i++) {
			n = &vnode[i];
			fprintf(fp, "%d\t%.4g\t%.4g\t%.4g\t%.4g\n", i, n->r[0], n->r[1], n->r[2], n->H);
		}
		for(i=0; i<NVface; i++) {
			f = &vface[i];
			fprintf(fp, "%d\t%d\t%d\t%d\n", i, f->inode[0], f->inode[1], f->inode[2]);
		}
	}
	
	fclose(fp);
	
	fp=fopen("id.dat", "w");	// record current file id
	fprintf(fp, "%d\n", Nfile);
	fclose(fp);
}



void SaveContour(void)	// save z(r)
{
	int i, j, k, n;
	char fn[50];
	double rmax, dr, rc[3];
	double *rbin, *zbin;
	double dri[3], ri, zi;
	MNODE *m;
	FILE *fp;
	
	if(NC0!=1) return;
	
	n = ContourBin;
	rmax = 0.5*min2(Lx, Ly);
	dr = rmax/(n-1);
	veccopy(chain->rcenter, rc);
	
	rbin = malloc(n*sizeof(double));
	zbin = malloc(n*sizeof(double));
		
	for(i=0; i<n; i++) rbin[i] = zbin[i] = 0;
	
	for(i=1; i<Nx; i++) {
		for(j=1; j<Ny; j++) {
			m = &mnode[i][j];
			vecsub(m->r, rc, dri);
			dri[2] = 0;
			ri = norm(dri);	// x-y radius
			zi = m->r[2];
			k = (int)(ri/dr+0.5);
			if(k<n) {
				rbin[k]+=1;
				zbin[k]+=zi;
			}
		}
	}
	
	for(i=0; i<n; i++) {
		if(rbin[i]>0) zbin[i]/=rbin[i];	// average
		if(zbin[i]<1e-4) zbin[i]=0;
	}
	
	sprintf(fn, "C%.4u.dat", Nfile);
	fp=fopen(fn, "w");
	fprintf(fp, "r\tz\n");
	for(i=0; i<n; i++) {
		if(rbin[i]>0) fprintf(fp, "%.4g\t%.4g\n", i*dr, zbin[i]);
	}
	fclose(fp);
	
	
	fp=fopen("contournew.ini", "w");
	fprintf(fp, "r\tz\n");
	for(i=0; i<n; i++) {
		if(rbin[i]>0) fprintf(fp, "%.4g\t%.4g\n", i*dr, zbin[i]);
	}
	fclose(fp);
	
	free(rbin);
	free(zbin);
}




void TestData(void)
{
	int i, j;
	char fn[50];
	MNODE *m;
	CHAIN *p;
	VNODE *n;
	VFACE *f;
	FILE *fp;
	
	sprintf(fn, "test%.4u.dat", Nfile);
	fp=fopen(fn, "w");
	
	fprintf(fp, "%.4g\t%.4g\t%.4g\t%d\t%d\t%.4g\t%.4g\t%.4g\n\n", t, Lx, Ly, Nx, Ny, DLT, Rv, Rtmd);
	
	fprintf(fp, "%d\t%d\t%d\n\n", Nchain, Nopen, Nclosed);
	
	p = chain;
	while(p) {
		fprintf(fp, "%d\t%d\t%d\n", p->id, p->nseg, p->closed);
		for(i=0; i<(p->nseg); i++) {	// for each segment on chain
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->r[i][0], p->r[i][1], p->r[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->nx[i][0], p->nx[i][1], p->nx[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->nz[i][0], p->nz[i][1], p->nz[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\t", p->rlinker[i][0], p->rlinker[i][1], p->rlinker[i][2]);
			fprintf(fp, "%.4g\t%.4g\t%.4g\n", p->rtmd[i][0], p->rtmd[i][1], p->rtmd[i][2]);
		}
		p = p->next;
	}
	
	fprintf(fp, "\n%.4g\t%.4g\t%.4g\t%.4g\n", Zmin, Zmax, Curvmin, Curvmax);
	for(i=0; i<Nx; i++) {
		for(j=0; j<Ny; j++) {
			m = &mnode[i][j];
			fprintf(fp, "%.4g\t%.4g\t%.4g\t%.4g\n", m->r[0], m->r[1], m->r[2], m->H);
		}
	}
	
	if(Vesicle) {
		fprintf(fp, "\n%d\t%d\t%d\t%.4g\t%.4g\t%.4g\n", NVnode, NVedge, NVface, RCves[0], RCves[1], RCves[2]);
		for(i=0; i<NVnode; i++) {
			n = &vnode[i];
			fprintf(fp, "%d\t%.4g\t%.4g\t%.4g\t%.4g\n", i, n->r[0], n->r[1], n->r[2], n->H);
		}
		for(i=0; i<NVface; i++) {
			f = &vface[i];
			fprintf(fp, "%d\t%d\t%d\t%d\n", i, f->inode[0], f->inode[1], f->inode[2]);
		}
	}
	
	fclose(fp);
}



void initsyt_resume(void)
{
	int i;
	CHAIN *p, *q;
	
	for(i=0; i<Nchain; i++) {	// for each chain
		if(!chain) {
			//printf("a%d\t", i);
			p = (CHAIN *)malloc(sizeof(CHAIN));
			p->next = NULL;
			chain = p;
		}
		else {
			//printf("b%d\t", i);
			p = chain;
			while(p->next) p = p->next;	// p is the end of chain
			q = (CHAIN *)malloc(sizeof(CHAIN));
			p->next = q;
			q->next = NULL;
		}
	}
	
	vnode=(VNODE *)malloc(NnodeMax*sizeof(VNODE));
	vedge=(VEDGE *)malloc(NedgeMax*sizeof(VEDGE));
	vface=(VFACE *)malloc(NfaceMax*sizeof(VFACE));
}


// to resume simulation, read data from existing files
int ReadData(void)
{
	int i, j, k;
	char fn[50];
	float tf;
	float x, y, z, h;
	MNODE *m;
	CHAIN *p;
	VNODE *n;
	VFACE *f;
	FILE *fp;
	
	fp=fopen("id.dat", "r");
	if(!fp) {
		printf("'id.dat' not found, cannot resume simulations.\n\n");
		exit(0);;
	}
	k = fscanf(fp, "%d", &Nfile);
	fclose(fp);
	
	sprintf(fn, "P%.4u.dat", Nfile);
	fp=fopen(fn, "r");
	if(!fp) { printf("Cannot read data from '%s'. Use default initial conditions.\n", fn); return 0; }
	
	k = fscanf(fp, "%f %f %f %d %d %f %f %f", &tf, &Lx, &Ly, &Nx, &Ny, &DLT, &Rv, &Rtmd);
	skip(fp);
	
	t=(double)(tf);
	printf("Resume simulation from t=%.4g ms\n", t*1e3);
	
	// --- for syt ---
	k = fscanf(fp, "%d %d %d", &Nchain, &Nopen, &Nclosed);
	skip(fp);
	
	if(Nchain==0) { printf("Nchain = 0\n"); exit(0); }
	if(Nx*Ny==0) { printf("Nx*Ny = 0\n"); exit(0); }
	
	initsyt_resume();	// initialize chains after getting Nchain
	
	p = chain;
	while(p) {
		k = fscanf(fp, "%d %d %d", &(p->id), &(p->nseg), &(p->closed));
		skip(fp);
		for(i=0; i<(p->nseg); i++) {	// for each segment on chain
			k = fscanf(fp, "%f %f %f", &x, &y, &z);	// read to float type
			p->r[i][0]=x;	p->r[i][1]=y;	p->r[i][2]=z;	// double type
			k = fscanf(fp, "%f %f %f", &x, &y, &z);
			p->nx[i][0]=x;	p->nx[i][1]=y;	p->nx[i][2]=z;
			k = fscanf(fp, "%f %f %f", &x, &y, &z);
			p->nz[i][0]=x;	p->nz[i][1]=y;	p->nz[i][2]=z;
			k = fscanf(fp, "%f %f %f", &x, &y, &z);
			p->rlinker[i][0]=x;	p->rlinker[i][1]=y;	p->rlinker[i][2]=z;
			k = fscanf(fp, "%f %f %f", &x, &y, &z);
			p->rtmd[i][0]=x;	p->rtmd[i][1]=y;	p->rtmd[i][2]=z;
			
			crossprod(p->nz[i], p->nx[i], p->ny[i]);
		}
		
		UpdateSytSites(p);	// given p->r, update the rest positions excluding TMD position
		GetSytStretch(p);
		
		p = p->next;
	}
	
	//CheckChain();
	
	// --- for membrane ---
	if(Membrane) {
		k = fscanf(fp, "%f %f %f %f", &x, &y, &z, &h);
		skip(fp);
		Zmin=x;	Zmax=y;	Curvmin=z;	Curvmax=h;
		
		
		for(i=0; i<Nx; i++) {
			for(j=0; j<Ny; j++) {
				m = &mnode[i][j];
				k = fscanf(fp, "%f %f %f %f", &x, &y, &z, &h);
				m->r[0]=x;	m->r[1]=y;	m->r[2]=z;	m->H=h;
			}
		}
	}
	//testGrid();
	
	// --- for vesicle ---
	if(Vesicle) {
		k = fscanf(fp, "%d %d %d %f %f %f", &NVnode, &NVedge, &NVface, &x, &y, &z);
		RCves[0]=x;	RCves[1]=y;	RCves[2]=z;
		for(i=0; i<NVnode; i++) {
			n = &vnode[i];
			k = fscanf(fp, "%d %f %f %f %f", &(n->id), &x, &y, &z, &h);
			n->r[0] = x;
			n->r[1] = y;
			n->r[2] = z;
			n->H = h;
		}
		for(i=0; i<NVface; i++) {
			f = &vface[i];
			k = fscanf(fp, "%d %d %d %d", &(f->id), &(f->inode[0]), &(f->inode[1]), &(f->inode[2]));
			for(j=0; j<3; j++) {
				f->r1[j] = vnode[f->inode[0]].r[j];
				f->r2[j] = vnode[f->inode[1]].r[j];
				f->r3[j] = vnode[f->inode[2]].r[j];
			}
		}
	}
	
	fclose(fp);
	
	//TestData();
	return 1;
}




double SytMemDist(void)
{
	int i, n, cnt;
	double d, dave;
	CHAIN *p;

	p = chain;
	cnt = 0;
	dave = 0;
	while(p) {
		n = p->nseg;
		d = 0;
		if(n>0) {
			for(i=0; i<n; i++) {
				d += getSMdistance(mnode, p->rsm[i]);
				//printf("%.3g\t", getSMdistance(mnode, p->rsm[i]));
			}
			d = d/n;
			dave += d;
		}
		cnt++;
		p = p->next;
	}
	cnt = max2(cnt, 1);
	dave /= cnt;
	return dave;
}





void SaveEnergy(void)
{
	int n;
	double zs, zrel, dave;
	FILE *fp;
	
#if LGV_MC
	if(t==0) {
		fp = fopen("energy.dat", "w");
		fprintf(fp, "t\t");
#else
	if(Nfile==0) {
		fp = fopen("energy.dat", "w");
		fprintf(fp, "Loop\t");
#endif
		fprintf(fp, "H\tHrel\tdist\tEsyt\tEes\tEmem\tEves\tEtot\tEes/N\tEmem/N\tEves/N\tEesmv/N\n");
	}
	else fp = fopen("energy.dat", "a");
	
	if(Limit==0) n = NS0;
	else n = NS1;
	
	if(Membrane) {
		GetMemMaxMin(mnode);	// get Zmax
		dave=SytMemDist();
		zs = chain->rcenter[2];
		zs -= Rc2b;	// bottom of syt ring
		zrel = Zmax - zs;
		if(Zmax<1e-6) Zmax=0;
		if(zrel<1e-6) zrel=0;
	}
	else {
		Zmax = 0;
		zrel = 0;
		dave = 0;
	}
	
	
#if LGV_MC
	fprintf(fp, "%.4g\t", t);
#else
	fprintf(fp, "%d\t", MCloopcnt);
#endif
	
	fprintf(fp, "%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n", 
		Zmax, zrel, dave, Esyt/kBT, EsytES/kBT, Emem/kBT, Eves/kBT, Esys/kBT, 
		EsytES/kBT/n, Emem/kBT/n, Eves/kBT/n, (EsytES+Emem+Eves)/kBT/n);
	
	fclose(fp);
}



