
void Rep_VesN_MemN(VNODE *p, MNODE *m)	// ves-node -  mem-node repulsion
{
	double dr[3], r, r2, f, fv[3];
	double tau[3];
	
	vecsub(p->r, m->r, dr);	// dr from membrane node to vesicle node
	r2 = norm2(dr);
	
	if(r2<drVMrepulsion2) {	// in contact
		r = max2(sqrt(r2), eps);
		f = KV*(drVMrepulsion-r);	// f>0
		vecprod(dr, f/r, fv);	// force on vesicle vnode
		vecadd(p->f, fv, p->f);
		vecadd(p->fext, fv, p->fext);
		
#if OPENMP
#pragma omp critical(MemVesRepulsion1)
{
#endif
		vecsub(m->f, fv, m->f);	// be careful, use omp critical !!!
		TrackFves += fv[2];
		TrackFmem -= fv[2];
	//vecadd(Fves, fv, Fves);
		
		crossprod(p->rc, fv, tau);	// torque on vesicle
		vecadd(Tauves, tau, Tauves);
#if OPENMP
}
#endif

	}
}


void Rep_VesN_MemF(MNODE **MN, VNODE *p, MFACE *mf)	// ves-node - mem-face repulsion
{
	int k, i[3], j[3];
	double dr[3], r, r2, f, fv[3], fv_3[3];
	double tau[3];
	MNODE *m1, *m2, *m3;
	
	for(k=0; k<3; k++) {	// for each vertex of face-mf
		i[k] = mf->inode[k][0];
		j[k] = mf->inode[k][1];
	}
	m1 = &MN[i[0]][j[0]];	// 3 vertices of membrane face
	m2 = &MN[i[1]][j[1]];
	m3 = &MN[i[2]][j[2]];
	
	vecsub(p->r, mf->rc, dr);	// dr from mem-center to vesicle node
	r2 = norm2(dr);
	
	if(r2<drVMrepulsion2) {	// in contact
		r = max2(sqrt(r2), eps);
		f = KV*(drVMrepulsion-r);	// f>0
		vecprod(dr, f/r, fv);	// force on vesicle vnode
		vecdiv(fv, 3.0, fv_3);
		vecadd(p->f, fv, p->f);
		vecadd(p->fext, fv, p->fext);
		
#if OPENMP
#pragma omp critical(MemVesRepulsion1)
{
#endif
		vecsub(m1->f, fv_3, m1->f);	// be careful, use omp critical !!!
		vecsub(m2->f, fv_3, m2->f);
		vecsub(m3->f, fv_3, m3->f);
		
		TrackFves += fv[2];
		TrackFmem -= fv[2];
	//vecadd(Fves, fv, Fves);
		
		crossprod(p->rc, fv, tau);	// torque on vesicle
		vecadd(Tauves, tau, Tauves);
#if OPENMP
}
#endif

	}
}


void Rep_VesF_MemN(VNODE *VN, VFACE *vf, MNODE *m)	// ves-face -  mem-node repulsion
{
	int i1, i2, i3;
	double dr[3], r, r2, f, fv[3], fv_3[3];
	double tau[3];
	VNODE *n1, *n2, *n3;
	
	i1 = vf->inode[0];	// 3 vertices of ves-face
	i2 = vf->inode[1];
	i3 = vf->inode[2];
	n1 = &VN[i1];
	n2 = &VN[i2];
	n3 = &VN[i3];
	
	vecsub(vf->rclab, m->r, dr);	// dr from membrane node to vesicle face center
	r2 = norm2(dr);
	
	if(r2<drVMrepulsion2) {	// in contact
		r = max2(sqrt(r2), eps);
		f = KV*(drVMrepulsion-r);	// f>0
		vecprod(dr, f/r, fv);
		vecdiv(fv, 3.0, fv_3);	// force on each vesicle vnode
		
		vecadd(n1->f, fv_3, n1->f);	// force on each vertex of ves face
		vecadd(n2->f, fv_3, n2->f);
		vecadd(n3->f, fv_3, n3->f);
		
		vecadd(n1->fext, fv_3, n1->fext);
		vecadd(n2->fext, fv_3, n2->fext);
		vecadd(n3->fext, fv_3, n3->fext);
		
#if OPENMP
#pragma omp critical(MemVesRepulsion2)
{
#endif
		vecsub(m->f, fv, m->f);	// be careful, use omp critical !!!
		TrackFves += fv[2];
		TrackFmem -= fv[2];
	//vecadd(Fves, fv, Fves);
		
		crossprod(vf->rcrc, fv, tau);	// torque on vesicle
		vecadd(Tauves, tau, Tauves);
#if OPENMP
}
#endif

	}
}


void Rep_VesF_MemF(VNODE *VN, VFACE *vf, MNODE **MN, MFACE *mf)	// ves-face - mem-face repulsion
{
	int i1, i2, i3;
	int k, i[3], j[3];
	double dr[3], r, r2, f, fv[3], fv_3[3];
	double tau[3];
	MNODE *m1, *m2, *m3;
	VNODE *n1, *n2, *n3;
	
	i1 = vf->inode[0];	// 3 vertices of ves-face
	i2 = vf->inode[1];
	i3 = vf->inode[2];
	n1 = &VN[i1];
	n2 = &VN[i2];
	n3 = &VN[i3];
	
	for(k=0; k<3; k++) {	// for each vertex of mem face-mf
		i[k] = mf->inode[k][0];
		j[k] = mf->inode[k][1];
	}
	m1 = &MN[i[0]][j[0]];	// 3 vertices of membrane face
	m2 = &MN[i[1]][j[1]];
	m3 = &MN[i[2]][j[2]];
	
	vecsub(vf->rclab, mf->rc, dr);	// dr from mem-center to vesicle face center
	r2 = norm2(dr);
	
	if(r2<drVMrepulsion2) {	// in contact
		r = max2(sqrt(r2), eps);
		f = KV*(drVMrepulsion-r);	// f>0
		vecprod(dr, f/r, fv);	// force on vesicle vnode
		vecdiv(fv, 3.0, fv_3);
		
		vecadd(n1->f, fv_3, n1->f);	// force on each vertex of ves face
		vecadd(n2->f, fv_3, n2->f);
		vecadd(n3->f, fv_3, n3->f);
		
		vecadd(n1->fext, fv_3, n1->fext);
		vecadd(n2->fext, fv_3, n2->fext);
		vecadd(n3->fext, fv_3, n3->fext);
		
#if OPENMP
#pragma omp critical(MemVesRepulsion2)
{
#endif
		vecsub(m1->f, fv_3, m1->f);	// be careful, use omp critical !!!
		vecsub(m2->f, fv_3, m2->f);
		vecsub(m3->f, fv_3, m3->f);
		
		TrackFves += fv[2];
		TrackFmem -= fv[2];
	//vecadd(Fves, fv, Fves);
		
		crossprod(vf->rcrc, fv, tau);	// torque on vesicle
		vecadd(Tauves, tau, Tauves);
#if OPENMP
}
#endif

	}
}



void VesMemRepulsion(VNODE *VN, VFACE *VF, MNODE **MN, MFACE *MF)	// for each ves node, find repulsion from membrane
{
	int i, j, k, l, imin, imax, jmin, jmax;
	int idx, idy;
	double hmax, x, y;
	VNODE *vn;
	VFACE *vf;
	MNODE *mn;
	MFACE *mf;
	MZONE *z;
	
	GetMemMaxMin(MN);
	hmax = Zmax + drVMrepulsion;	// critical height that vnode will interact w/ membrane
	
#if OPENMP
#pragma omp parallel for private(i,j,k,l,imin,imax,jmin,jmax,idx,idy,x,y,vn,mn,mf,z) \
schedule(dynamic,2)
#endif

	for(k=0; k<NVnode; k++) {	// for each vesicle node
		vn = &VN[k];
		if(vn->r[2] <= hmax) {	// has interaction
			x = vn->r[0] + Lxhalf;
			y = vn->r[1] + Lyhalf;
			imin = (int)(x/dX+1.5) - NVMnbrhalf;	// lower left corner of searching zones
			jmin = (int)(y/dY+1.5) - NVMnbrhalf;
			imax = imin + NVMnbr;
			jmax = jmin + NVMnbr;
			
			for(i=imin; i<imax; i++) {
				for(j=jmin; j<=jmax; j++) {
					z = &mzone[i][j];	// zone [i][j]
					
					for(l=0; l<(z->nnode); l++) {	// for each m-node in zone
						idx = z->inode[l][0];
						idy = z->inode[l][1];
						mn = &MN[idx][idy];
						Rep_VesN_MemN(vn,mn);	// ves-node - mem-node repulsion
					}
					
					for(l=0; l<(z->nface); l++) {	// for each m-face-center in zone
						idx = z->iface[l];
						mf = &MF[idx];
						Rep_VesN_MemF(MN,vn,mf);	// ves-node - mem-face repulsion
					}
				}
			}
		}
	}
	
#if OPENMP
#pragma omp parallel for private(i,j,k,l,imin,imax,jmin,jmax,idx,idy,x,y,vf,mn,mf,z) \
schedule(dynamic,2)
#endif

	for(k=0; k<NVface; k++) {	// for each vesicle face center
		vf = &VF[k];
		if(vf->rclab[2] <= hmax) {	// has interaction
			x = vf->rclab[0] + Lxhalf;
			y = vf->rclab[1] + Lyhalf;
			imin = (int)(x/dX+1.5) - NVMnbrhalf;	// lower left corner of searching zones
			jmin = (int)(y/dY+1.5) - NVMnbrhalf;
			imax = imin + NVMnbr;
			jmax = jmin + NVMnbr;
			
			for(i=imin; i<imax; i++) {
				for(j=jmin; j<=jmax; j++) {
					z = &mzone[i][j];	// zone [i][j]
					
					for(l=0; l<(z->nnode); l++) {	// for each m-node in zone
						idx = z->inode[l][0];
						idy = z->inode[l][1];
						mn = &MN[idx][idy];
						Rep_VesF_MemN(VN,vf,mn);	// ves-face - mem-node repulsion
					}
					
					for(l=0; l<(z->nface); l++) {	// for each m-face-center in zone
						idx = z->iface[l];
						mf = &MF[idx];
						Rep_VesF_MemF(VN,vf,MN,mf);	// ves-face - mem-face repulsion
					}
				}
			}
		}
	}
}




void FTMemVes(VNODE *VN, VFACE *VF, MNODE **MN, MFACE *MF)
{
	VesMemRepulsion(VN, VF, MN, MF);
}

