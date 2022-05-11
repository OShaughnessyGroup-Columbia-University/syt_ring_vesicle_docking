
void PreProcessVes(void)
{
	int i;
	
	NVface=NVedge=NVnode=0;
	
	NnodeMax=12;
	NedgeMax=30;
	NfaceMax=20;
	for(i=0; i<Ndiv; i++) {
		NnodeMax=NnodeMax+3*NfaceMax/2;
		NedgeMax=2*NedgeMax+3*NfaceMax;
		NfaceMax*=4;
	}
	//printf("NnodeMax=%d,\tNedgeMax=%d,\tNfaceMax=%d\n", NnodeMax, NedgeMax, NfaceMax);

	vnode=(VNODE *)malloc(NnodeMax*sizeof(VNODE));
	vedge=(VEDGE *)malloc(NedgeMax*sizeof(VEDGE));
	vface=(VFACE *)malloc(NfaceMax*sizeof(VFACE));

	for(i=0; i<NnodeMax; i++) {
		vnode[i].nnbr = vnode[i].nface = vnode[i].nedge = 0;
	}
	for(i=0; i<NfaceMax; i++) {
		vface[i].inode[0] = vface[i].inode[1] = vface[i].inode[2] = -1;
	}
	
#if LGV_MC && BD_2ND
	vnode_try=(VNODE *)malloc(NnodeMax*sizeof(VNODE));
	vedge_try=(VEDGE *)malloc(NedgeMax*sizeof(VEDGE));
	vface_try=(VFACE *)malloc(NfaceMax*sizeof(VFACE));

	for(i=0; i<NnodeMax; i++) {
		vnode_try[i].nnbr = vnode_try[i].nface = vnode_try[i].nedge = 0;
	}
	for(i=0; i<NfaceMax; i++) {
		vface_try[i].inode[0] = vface_try[i].inode[1] = vface_try[i].inode[2] = -1;
	}
#endif
}



void addface(double v1[], double v2[], double v3[])	// vertices are stored in vface.r1 - vface.r3
{
	veccopy(v1, vface[NVface].r1);
	veccopy(v2, vface[NVface].r2);
	veccopy(v3, vface[NVface].r3);
	NVface++;
	if(NVface>NfaceMax) {printf("Error: NVface =%d > %d!\n", NVface, NfaceMax); exit(0); }
}



void griddiv(double v1[], double v2[], double v3[], int depth)	// taken from redbook pg.59
{
	int i;
	double v12[3], v23[3], v31[3];

	if (depth == 0) {
		addface(v1, v2, v3);
		return;
	}
	for (i = 0; i < 3; i++) {
	v12[i] = v1[i]+v2[i];
	v23[i] = v2[i]+v3[i];
	v31[i] = v3[i]+v1[i];
	}
	normalize(v12, v12);	// on unit sphere
	normalize(v23, v23);
	normalize(v31, v31);

	griddiv(v1, v12, v31, depth-1);
	griddiv(v2, v23, v12, depth-1);
	griddiv(v3, v31, v23, depth-1);
	griddiv(v12, v23, v31, depth-1);
}



void AdjustVesFaceGrid(void)	// map to a sphere of radius Rv
{
	int i;
	for(i=0; i<NVface; i++) {
		vecprod(vface[i].r1, 1.0*Rv, vface[i].r1);	// map to sphere
		vecprod(vface[i].r2, 1.0*Rv, vface[i].r2);
		vecprod(vface[i].r3, 1.0*Rv, vface[i].r3);
		vecave3(vface[i].r1, vface[i].r2, vface[i].r3, vface[i].rcrc);
		vecadd(vface[i].r1, RCves, vface[i].r1);	// in lab frame
		vecadd(vface[i].r2, RCves, vface[i].r2);
		vecadd(vface[i].r3, RCves, vface[i].r3);
		vecave3(vface[i].r1, vface[i].r2, vface[i].r3, vface[i].rclab);
	}
}



int CompareVesNodeFace(int id[3], VFACE *f)	// return 1 if nodes-id[3] are vertices of vface-f
{
	int i, j, flg, z;
	
	flg = 0;
	for(i=0; i<3 && flg<3; i++) {
		for(j=0; j<3 && flg<3; j++) {
			if(id[i] == f->inode[j]) flg++;
		}
	}
	if(flg==3) z=1;
	else z=0;
	
	return z;
}



void GetVesNodes(void)	// find connections and copy vface.r to vnode.r
{
	int i, j, k, m, flg, nf, i1, i2, ix;
	double torr, vertex[3][3];
	VFACE *p;

	torr=1.0e-3*sqrt(4*Pi*Rv*Rv/NnodeMax);	// nodes within this distance is regarded as identical

	for(i=0; i<NVface; i++) {	// scan through faces for nodes
		p = &vface[i];
		veccopy(p->r1, vertex[0]);
		veccopy(p->r2, vertex[1]);
		veccopy(p->r3, vertex[2]);

		for(j=0; j<3; j++) {	// for each vertex in the face
			flg=1;
			for(k=0; k<NVnode && flg; k++) {	// search through existing nodes to see if there's a duplicate
				if(vectorequal(vertex[j], vnode[k].r, torr)) {	// vertex[j] already exists in vnode
					flg=0;
					m=k;	// don't use k outside b/c k=k+1 after this loop
				}
			}
			if(flg) {	// no duplicates
				vface[i].inode[j] = NVnode;	// update vface
				vnode[NVnode].id = NVnode;	// update vnode
				
				veccopy(vertex[j], vnode[NVnode].r);
				vecsub(vnode[NVnode].r, RCves, vnode[NVnode].rc);	// vnode.rc: vnode's coord in RC's frame
				veccopy(vnode[NVnode].rc, vnode[NVnode].r0);	// reference coord in RC's frame (for remapping)
				
				nf=vnode[NVnode].nface;
				vnode[NVnode].iface[nf] = i;
				(vnode[NVnode].nface)++;
				NVnode++;
				if(NVnode>NnodeMax) {printf("Error: NVnode =%d > %d!\n", NVnode, NnodeMax); exit(0); }
			}
			else {
				vface[i].inode[j] = m;	// update with the existing vnode
				nf=vnode[m].nface;
				vnode[m].iface[nf] = i;
				(vnode[m].nface)++;
			}
		}
	}
	
	for(i=0; i<NVface; i++) {	// scan through faces for edges
		p = &vface[i];
		
		for(j=0; j<3; j++) {	// go through all 3 edges
			flg=1;
			k=(j+1)%3;
			i1 = p->inode[j];
			i2 = p->inode[k];
			if(i1>i2) { ix=i1; i1=i2; i2=ix; }	// swap s.t. edge pair is always small->big

			for(k=0; k<NVedge && flg; k++) {	// search through existing edges to see if there's a duplicate
				if(vedge[k].inode[0]==i1 && vedge[k].inode[1]==i2) flg=0;
			}
			if(flg) {	// no duplicates
				vnode[i1].iedge[vnode[i1].nedge] = NVedge;	// update vnode
				vnode[i2].iedge[vnode[i2].nedge] = NVedge;
				(vnode[i1].nedge)++;
				(vnode[i2].nedge)++;

				vnode[i1].nbr[vnode[i1].nnbr] = i2;
				vnode[i2].nbr[vnode[i2].nnbr] = i1;
				(vnode[i1].nnbr)++;
				(vnode[i2].nnbr)++;

				vedge[NVedge].id = NVedge;	// update vedge
				vedge[NVedge].inode[0]=i1;
				vedge[NVedge].inode[1]=i2;
				NVedge++;
				if(NVedge>NedgeMax) {printf("Error: NVedge =%d > %d!\n", NVedge, NedgeMax); exit(0); }
			}
		}
	}
	// vedge is not sorted
}




void sortvesnbrid(int id)	// after getting vnode.rc
{
	int i, j, k, l, nnbr, nface, flg;
	int nbr[NnbrMax], fc[NnbrMax], vertex[3];
	double dr[NnbrMax][3], nr[3], crs[3];
	double a[NnbrMax], ax;
	VNODE *p;
	VFACE *f;
	
	p = &vnode[id];
	nnbr = p->nnbr;
	nface = p->nface;
	vecsub(p->r, RCves, nr);	// roughly outward normal
	
	for(i=0; i<nnbr; i++) {
		nbr[i] = p->nbr[i];	// original nbr list
		vecsub(vnode[p->nbr[i]].r, p->r, dr[i]);	// vector from vnode-id to nbr[i]
	}
	
	a[0]=0;
	for(i=1; i<nnbr; i++) {	// keep nbr[0] intact
		crossprod(dr[0], dr[i], crs);
		if(dotprod(crs, nr)>0) a[i] = getanglePi(dr[0], dr[i]);	// ccw neighbors
		else a[i] = Pidb - getanglePi(dr[0], dr[i]);	// cw neighbors
	}
	
	// now, sort nbr[] according to a[]
	for(i=1; i<nnbr-1; i++) {	// keep nbr[0] intact
		for(j=i+1; j<nnbr; j++) {
			if(a[j]<a[i]) {	// sort a[i] from small to large
				ax=a[i];	// swap a[i] a[j]
				a[i]=a[j];
				a[j]=ax;
				k=nbr[i];	// swap nbr[i] nbr[j] accordingly
				nbr[i]=nbr[j];
				nbr[j]=k;
			}
		}
	}
	
	// redefine neighbors according to list[]
	for(i=1; i<nnbr; i++) p->nbr[i] = nbr[i];
	
	// then, find corresponding faces
	for(i=0; i<nface; i++) fc[i] = p->iface[i];	// keep a copy of current faces for future use
	
	for(i=0; i<nnbr; i++) {	// for each triangle id-i-(i+1), define it to be vface-i
		j=(i+1)%nnbr;
		vertex[0] = id;
		vertex[1] = p->nbr[i];
		vertex[2] = p->nbr[j];
		
		flg=0;
		for(j=0; j<nface && flg==0; j++) {	// check if vertex[] matches existing neighboring vface-j
			f = &vface[fc[j]];
			flg=CompareVesNodeFace(vertex,f);
			if(flg==1) l=j;
		}
		p->iface[i] = fc[l];	// update the id of vface-i
		
	}	// end of i
}



void SortVesNbr(void)	// sort neighboring nodes such that vnode.nbr[i] -> vnode.nbr[i+1] is ccw
{
	int i;
	for(i=0; i<NVnode; i++) sortvesnbrid(i);
}



void PostProcessVes(void)
{
	int i;
	double r12[3], r13[3], nrm[3];
	
	// get signs for normal directions
	/*
	// not needed because all neighbors are ccw, see sortnbr()
	for(i=0; i<NVnode; i++) {	// for each vnode
		i1 = vnode[i].nbr[0];
		i2 = vnode[i].nbr[1];
		vecsub(vnode[i1].r, vnode[i].r, r12);
		vecsub(vnode[i2].r, vnode[i].r, r13);
		crossprod(r12, r13, nrm);
		if(dotprod(nrm, vnode[i].r0)>0) vnode[i].nrmsign=1;
		else vnode[i].nrmsign=-1;
	}
	*/
	
	for(i=0; i<NVface; i++) {	// for each face
		vecsub(vface[i].r1, vface[i].r2, r12);
		vecsub(vface[i].r1, vface[i].r3, r13);
		crossprod(r12, r13, nrm);
		if(dotprod(nrm, vface[i].rcrc)>0) vface[i].nrmsign=1;	// r12 x r13 is along outward normal
		else vface[i].nrmsign=-1;	// -r12 x r13 is along outward normal
		
		//printf("%d\t%d\n", i, vface[i].nrmsign);
	}
}


/*
void GetFID(void)	// determine which vnode the external force applies to
{
	int i;
	double x, xmin, xmax;
	
	xmin=Rv;
	xmax=-Rv;
	
	for(i=0; i<NVnode; i++) {
		x=vnode[i].r[0];
		if(x>xmax) { xmax=x; FID[0]=i; }
		if(x<xmin) { xmin=x; FID[1]=i; }
	}
}
*/


void testvesmesh(void)
{
	int i, j;
	
	for(i=0; i<NVnode; i++) {
		printf("i=%d:", i);
		for(j=0; j<vnode[i].nface; j++) printf("\tf%d", vnode[i].iface[j]);
		printf("\n");
	}
}


void InitVes(void)
{
	int i;
	
	PreProcessVes();
	
	// divide each of the 20 faces, vface.r are calculated
	for(i=0; i<20; i++) griddiv(Vico[Fico[i][0]], Vico[Fico[i][1]], Vico[Fico[i][2]], Ndiv);
	AdjustVesFaceGrid();	// adjust vface.r according to Rv and RC
	
	GetVesNodes();	// find connections and define vnode.r
	SortVesNbr();	// sort neighboring nodes such that vnode.nbr[i] -> vnode.nbr[i+1] is ccw
	
	PostProcessVes();	// update coordinates
	//GetFID();
	
	CalcVesAreaVolume(vnode,vface,0);	// get AreaVes0 & VolVes0
	
	//testvesmesh();
	//printf("NVnode=%d,\tNedge=%d,\tNface=%d\n", NVnode, NVedge, NVface);
	//printf("NVnode=%d, dA=%.3g\n", NVnode, dAvesNode0);
}



void CleanVes(void)
{
	if(vnode) { free(vnode); vnode=NULL; }
	if(vedge) { free(vedge); vedge=NULL; }
	if(vface) { free(vface); vface=NULL; }
		
#if LGV_MC && BD_2ND
	if(vnode_try) { free(vnode_try); vnode_try=NULL; }
	if(vedge_try) { free(vedge_try); vedge_try=NULL; }
	if(vface_try) { free(vface_try); vface_try=NULL; }
#endif
}
