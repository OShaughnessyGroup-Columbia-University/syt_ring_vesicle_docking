
void InitMemVars(void)
{
	int i;
	
	mnode=(MNODE **)calloc(Nxp1, sizeof(MNODE *));
	for(i=0; i<=Nx; i++) {
		mnode[i]=(MNODE *)calloc(Nyp1, sizeof(MNODE));
		if(!mnode[i]) { printf("Cannot initialize mnode\n"); exit(0); }
	}
	
	mface=(MFACE *)calloc(NMface, sizeof(MFACE));
	
	mzone=(MZONE **)calloc(Nxp1, sizeof(MZONE *));
	vzone=(VZONE **)calloc(Nxp1, sizeof(VZONE *));
	for(i=0; i<=Nx; i++) {
		mzone[i]=(MZONE *)calloc(Nyp1, sizeof(MZONE));
		vzone[i]=(VZONE *)calloc(Nyp1, sizeof(VZONE));
		if(!mzone[i] || !vzone[i]) { printf("Cannot initialize mzone or vzone\n"); exit(0); }
	}
	
	
#if LGV_MC && BD_2ND
	mnode_try=(MNODE **)calloc(Nxp1, sizeof(MNODE *));
	for(i=0; i<=Nx; i++) {
		mnode_try[i]=(MNODE *)calloc(Nyp1, sizeof(MNODE));
		if(!mnode_try[i]) { printf("Cannot initialize mnode_try\n"); exit(0); }
	}
	
	mface_try=(MFACE *)calloc(NMface, sizeof(MFACE));
#endif
	
}


double getbuclkedmemz(double r, double r0)
{
	double r1, dr, h, z;
	
	r1 = r0 - Rc2b;	// radius through KKKK
	
	if(r>=r0) z = 0;	// flat
	else if(r>=r1) {	// touching C2B
		dr = r0 - r;
		h = sqrt(max2(Rc2b2 - pow(dr, 2),0));
		z = Rc2b - h;
	}
	else {	// in air
		h = sqrt(max2(r1*r1-r*r,0));
		z = Rc2b + h;
	}
	return z;
}



void readcontour(void)
{
	int i, n;
	char ch, fn[50], a[50];
	FILE *fp;
	
	if(Resume==0) sprintf(fn, "contour.ini");
	else {
		sprintf(fn, "contournew.ini");
		fp=fopen(fn, "r");
		if(!fp) sprintf(fn, "contour.ini");
		else fclose(fp);
	}
	
	fp=fopen(fn, "r");
	if(!fp) {
		printf("Cannot open '%s'.\n", fn);
		exit(0);
	}
	
	Ncontour=0;
	do {
		ch=fgetc(fp);
		if(ch=='\n') Ncontour++;
	} while(ch!=EOF);
	
	fclose(fp);
	
	Ncontour--;
	
	Rcontour=malloc(Ncontour*sizeof(float));
	Zcontour=malloc(Ncontour*sizeof(float));
	
	fp=fopen(fn, "r");
	n=fscanf(fp, "%s %s", a, a);
	//n=fscanf(fp, "%f %f %f %f", &Rcontour[0], &Zcontour[0], &tmp, &tmp);
	for(i=0; i<Ncontour; i++) n=fscanf(fp, "%f %f", &Rcontour[i], &Zcontour[i]);
	fclose(fp);
	
	//for(i=0; i<Ncontour; i++) printf("%d\t%.3f\t%.3f\n", i, Rcontour[i], Zcontour[i]);
	//printf("\n");
}



double getcontour(double r)
{
	int i, flg;
	double r1, r2, k, z;
	
	if(r<=Rcontour[0]) z=Zcontour[0];
	else if(r>=Rcontour[Ncontour-1]) z=0;
	else {
		i=flg=0;
		do {
			i++;
			if(r<Rcontour[i]) flg=1;
		} while(i<Ncontour-1 && !flg);
		
		r1=Rcontour[i-1];
		r2=Rcontour[i];
		k=(r-r1)/(r2-r1);
		z=Zcontour[i-1]+(Zcontour[i]-Zcontour[i-1])*k;
	}
	
	return z;
}


void testcontour(void)
{
	int i, n;
	double r, rm, dr;
	FILE *fp;
	
	n=30;
	rm=20.0;
	dr=rm/(n-1);
	fp=fopen("ct.dat", "w");
	for(i=0; i<n; i++) {
		r=i*dr;
		fprintf(fp, "%f\t%f\n", r, getcontour(r));
	}
	fclose(fp);
}


void InitMemNode(void)
{
	int i, j, ip1, im1, jp1, jm1, flg, num;
	int k, s;
	double x, y, r2, r0, r02, z;
	MNODE *m;
	
	if(INIMEM==2) {
		readcontour();
		//testcontour();
	}
	
	for(i=0; i<=Nx; i++) {	// (0,y)=(Nx-1,y) & (x,0)=(x,Ny-1) are ghost grids
		ip1 = i+1;
		im1 = i-1;
		for(j=0; j<=Ny; j++) {	// (Nx,y)=(1,y) & (x,Ny)=(x,1) are duplicates
			jp1 = j+1;
			jm1 = j-1;
			
			m = &mnode[i][j];
			
			m->id[0] = i;
			m->id[1] = j;
			
			m->nnbr = 6;
			m->nface = 6;
			
			// WARNING: m->nbr indices can go beyond range !!!!!!!!!!!!!!!!!!!!!!!!
			if(j%2==0) {
				m->nbr[0][0] = ip1;		m->nbr[0][1] = j;	// right
				m->nbr[1][0] = ip1;		m->nbr[1][1] = jp1;	// upper right
				m->nbr[2][0] = i;			m->nbr[2][1] = jp1;	// upper
				m->nbr[3][0] = im1;		m->nbr[3][1] = j;	// left
				m->nbr[4][0] = i;			m->nbr[4][1] = jm1;	// lower
				m->nbr[5][0] = ip1;		m->nbr[5][1] = jm1;	// lower right
			}
			else {
				m->nbr[0][0] = ip1;		m->nbr[0][1] = j;	// right
				m->nbr[1][0] = i;			m->nbr[1][1] = jp1;	// upper
				m->nbr[2][0] = im1;		m->nbr[2][1] = jp1;	// upper left
				m->nbr[3][0] = im1;		m->nbr[3][1] = j;	// left
				m->nbr[4][0] = im1;		m->nbr[4][1] = jm1;	// lower left
				m->nbr[5][0] = i;			m->nbr[5][1] = jm1;	// lower
			}
		}
	}
	
	// remove out-of-range neighbors
	for(i=0; i<=Nx; i+=Nx) {	// left & right boundaries
		for(j=0; j<=Ny; j++) {
			m = &mnode[i][j];
			k=0;
			
			while(k < m->nnbr) {
				flg=0;
				if(m->nbr[k][0]<0 || m->nbr[k][0]>Nx || m->nbr[k][1]<0 || m->nbr[k][1]>Ny) flg=1;
				if(flg) {	// neighbor is out of boundary, remove it!
					for(s=k; s<(m->nnbr)-1; s++) {	// shift the rest forward
						m->nbr[s][0] = m->nbr[s+1][0];
						m->nbr[s][1] = m->nbr[s+1][1];
					}
					(m->nnbr)--;
				}
				k++;
			}	// end of while
			
		}	// end of j
	}	// end of i
	
#if CLOSED
	r0 = Nsytperring*Lc2b/Pidb;	// inner radius through center of C2B
#else
	r0 = R0;	// related to spontaneous curvature
#endif

#if !SYTCNT
	r0 -= Rc2b;
#endif

	r02 = r0*r0;
	
	if(Limit==1) num=NS1;
	else num=NS0;

	for(i=0; i<=Nx; i++) {
		x=(i-1)*dX-Lxhalf;
		for(j=0; j<=Ny; j++) {
			y=(j-1)*dY-Lyhalf;
			m = &mnode[i][j];
			z=0;
			
			//if(Membrane && (Limit==0 || (Vesicle && VesRing==2)) && NC0==1 && (CLOSED || NS0>=17)) {
			if(Membrane && !Vesicle && Limit==0 && NC0==1 && (CLOSED || NS0>=17)) {
				r2=x*x+y*y;
				flg=0;
				if(Vesicle) {
					if(r2<r02 && r2>0.4*r02) flg=1;	// with a dent at center
				}
				else {
					if(INIMEM==2) flg=1;
					else if(r2<r02) flg=1;
				}
				
				if(SYT && flg) {
					if(Vesicle) z=Lc2a/2.0;	// step
					else {
						if(Spiral==2 && num>Nsytperring0 && !CLOSED) {
							z=Lc2a*(0.5+(int)(0.99*num/Nsytperring0));	// step
						}
						else {
							if(INIMEM==1) {	// curved
								z=h_mem_ini*(1+cos(Pi*sqrt(r2)/r0))/2.0;	// cosine
								//z=h_mem_ini;	// step
								//z=sqrt(r02-x*x-y*y);	// sphere
								//z=20*(1-sqrt(x*x+y*y)/(Lx/sqrt(10)));	// cone
								//z = getbuclkedmemz(sqrt(r2),r0);
							}
							else if(INIMEM==2) z=getcontour(sqrt(r2));	// read from contour_ini.dat
							else z=0;	// flat
						}
					}
				}
				else z=0;
				//z=0;	// test !!!!!
			}
			
			/*
			if(fabs(y)<Lx/4) z=sqrt(pow(Lx/4,2)-y*y);	// half-cylinder along x
			//if(fabs(x)<Lx/4) z=sqrt(pow(Lx/4,2)-x*x);	// half-cylinder along y
			else z=0;
			*/
			
		#if CRCBND
			if(Membrane) {
				if(x*x+y*y<=Lxhalf*Lxhalf) m->circbnd = 0;	// assuming Lx=Ly
				else {
					m->circbnd = 1;
					z=0;
				}
			}
		#endif
			
			m->r[0]=x;
			m->r[1]=y;
			m->r[2]=z;
			
			m->C1 = m->C2 = 0;
			veczero(m->e1);
			veczero(m->e2);
			veczero(m->e3);
		}
	}
	
	if(INIMEM==2) {
		free(Rcontour);
		free(Zcontour);
	}
}





int CompareMemNodeFace(int id[3][2], MFACE *f)	// return 1 if nodes-id[3] are vertices of mface-f
{
	int i, j, flg, z;
	
	flg = 0;
	for(i=0; i<3 && flg<3; i++) {
		for(j=0; j<3 && flg<3; j++) {
			if(id[i][0] == f->inode[j][0] && id[i][1] == f->inode[j][1]) flg++;
		}
	}
	if(flg==3) z=1;
	else z=0;
	
	return z;
}




void InitMemFace(void)
{
	int i, j, k, l, m, flg;
	int ip1, jp1, id;
	int i1, i2, i3, j1, j2, j3;
	int nnbr, nbrid[3][2];
	double r12[3], r13[3], nrm[3];
	MNODE *n;
	MFACE *f;
	
	id = 0;
	for(i=0; i<Nx; i++) {
		ip1 = i+1;
		for(j=0; j<Ny; j++) {
			jp1 = j+1;
			f = &mface[id];	// lower-right triangle
			if(j%2==0) {
				f->id = id;
				f->inode[0][0] = i;		f->inode[0][1] = j;
				f->inode[1][0] = ip1;	f->inode[1][1] = j;
				f->inode[2][0] = ip1;	f->inode[2][1] = jp1;
				id++;
				
				f = &mface[id];	// upper-left triangle
				f->id = id;
				f->inode[0][0] = i;		f->inode[0][1] = j;
				f->inode[1][0] = ip1;	f->inode[1][1] = jp1;
				f->inode[2][0] = i;		f->inode[2][1] = jp1;
				id++;
			}
			else {	// lower-left triangle
				f->id = id;
				f->inode[0][0] = i;		f->inode[0][1] = j;
				f->inode[1][0] = ip1;	f->inode[1][1] = j;
				f->inode[2][0] = i;		f->inode[2][1] = jp1;
				id++;
				
				f = &mface[id];	// upper-right triangle
				f->id = id;
				f->inode[0][0] = ip1;		f->inode[0][1] = j;
				f->inode[1][0] = ip1;	f->inode[1][1] = jp1;
				f->inode[2][0] = i;		f->inode[2][1] = jp1;
				id++;
			}
		}
	}
	
	NMface = id;	// actual # of faces
	
	for(i=0; i<NMface; i++) {
		f = &mface[i];
		i1 = f->inode[0][0];	j1 = f->inode[0][1];
		i2 = f->inode[1][0];	j2 = f->inode[1][1];
		i3 = f->inode[2][0];	j3 = f->inode[2][1];
		
		veccopy(mnode[i1][j1].r, f->r1);
		veccopy(mnode[i2][j2].r, f->r2);
		veccopy(mnode[i3][j3].r, f->r3);
		vecave3(f->r1, f->r2, f->r3, f->rc);
		
		vecsub(f->r1, f->r2, r12);
		vecsub(f->r1, f->r3, r13);
		crossprod(r12, r13, nrm);
		if(nrm[2]>0) f->nrmsign = 1;	// r12 x r13 is along outward normal
		else f->nrmsign = -1;	// -r12 x r13 is along outward normal
	}
	
	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			n = &mnode[i][j];
			nnbr = n->nnbr;
			
			for(k=0; k<nnbr; k++) {	// triangle (i,j)-k-l
				l = (k+1)%nnbr;
				nbrid[0][0] = i;	nbrid[0][1]=j;
				nbrid[1][0] = n->nbr[k][0];	nbrid[1][1] = n->nbr[k][1];
				nbrid[2][0] = n->nbr[l][0];	nbrid[2][1] = n->nbr[l][1];
				
				flg = 1;
				for(m=0; m<NMface && flg; m++) {
					f = &mface[m];
					if(CompareMemNodeFace(nbrid, f)) {
						n->iface[k] = f->id;
						flg = 0;
					}
				}
			}
		}
	}
}



void InitMemZone()
{
	int i, j;
	
	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			mzone[i][j].nnode = mzone[i][j].nface = 0;
			vzone[i][j].nnode = vzone[i][j].nface = 0;
		}
	}
}



void testGrid(void)
{
	int i, j;
	FILE *fp;
	
	fp=fopen("mnode.dat", "w");
	fprintf(fp, "i\tj\tx\ty\n");
	for(i=0; i<=Nx; i++) {
		for(j=0; j<=Ny; j++) {
			fprintf(fp, "%d\t%d\t%.4g\t%.4g\n", i, j, mnode[i][j].r[0], mnode[i][j].r[1]);
		}
	}
	fclose(fp);
	
	for(i=0; i<mnode[1][1].nnbr; i++) printf("%d: %d\t%d\n", i, mnode[1][1].nbr[i][0], mnode[1][1].nbr[i][1]);
	printf("r00=%.3g, %.3g\n", mnode[0][0].r[0], mnode[0][0].r[1]);
}



void InitMem(void)
{
	InitMemVars();
	InitMemNode();
	InitMemFace();
	InitMemZone();
	GetMemArea(mnode, mface, 0);
	//testGrid();
}



void CleanMem(void)
{
	int i;
	
	if(mnode) {
		for(i=0; i<=Nx; i++) free(mnode[i]);
		free(mnode);
	}
	mnode=NULL;
	
	if(mface) { free(mface); mface=NULL; }
	
	if(mzone) {
		for(i=0; i<=Nx; i++) free(mzone[i]);
		free(mzone);
	}
	mzone=NULL;
	
	if(vzone) {
		for(i=0; i<=Nx; i++) free(vzone[i]);
		free(vzone);
	}
	vzone=NULL;
	
#if LGV_MC && BD_2ND
	if(mnode_try) {
		for(i=0; i<=Nx; i++) free(mnode_try[i]);
		free(mnode_try);
	}
	mnode_try=NULL;
	
	if(mface_try) { free(mface_try); mface_try=NULL; }
#endif
}

