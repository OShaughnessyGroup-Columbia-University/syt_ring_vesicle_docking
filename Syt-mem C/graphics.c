
void DrawBottom(void)
{
	int i, j;
	float x, y, z;
	
	z=-0.2;
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_QUADS);
		i=0;	j=0;
		x=(float)(mnode[i][j].r[0]);
		y=(float)(mnode[i][j].r[1]);
		glVertex3f(x,z,-y);
		i=Nxm1;	j=0;
		x=(float)(mnode[i][j].r[0]);
		y=(float)(mnode[i][j].r[1]);
		glVertex3f(x,z,-y);
		i=Nxm1;	j=Nym1;
		x=(float)(mnode[i][j].r[0]);
		y=(float)(mnode[i][j].r[1]);
		glVertex3f(x,z,-y);
		i=0;	j=Nym1;
		x=(float)(mnode[i][j].r[0]);
		y=(float)(mnode[i][j].r[1]);
		glVertex3f(x,z,-y);
	glEnd();
}


void getcolor(double x, double xmin, double xmax, float rgb[3])
{// from http://paulbourke.net/texture_colour/colourspace/
	double dx;
	
	rgb[0]=rgb[1]=rgb[2]=1;
	
	x=max2(min2(x, xmax), xmin);
	dx=xmax-xmin;
	
	if(x<xmin+0.25*dx) {
		rgb[0]=0;
		rgb[1]=4.0*(x-xmin)/dx;
	}
	else if(x<xmin+0.5*dx) {
		rgb[0]=0;
		rgb[2]=1+4.0*(xmin+0.25*dx-x)/dx;
	}
	else if(x<xmin+0.75*dx) {
		rgb[0]=4.0*(x-xmin-0.5*dx)/dx;
		rgb[2]=0;
	}
	else {
		rgb[1]=1+4.0*(xmin+0.75*dx-x)/dx;
		rgb[2]=0;
	}
	/*
	rgb[0]+=0.1;
	rgb[1]+=0.1;
	rgb[2]+=0.1;
	*/
}


void drawlineclr1(double v1[], double v2[], double x, double xmin, double xmax)
{
	int i;
	float v1f[3], v2f[3], rgb[3];
	
	for(i=0; i<3; i++) {
		v1f[i]=(float) (v1[i]);
		v2f[i]=(float) (v2[i]);
	}
	getcolor(x, xmin, xmax, rgb);
	glColor3f(rgb[0], rgb[1], rgb[2]);
	
	glVertex3f(v1f[0], v1f[2], -v1f[1]);
	glVertex3f(v2f[0], v2f[2], -v2f[1]);
}


void drawlineclr(double v1[], double v2[], double x1, double x2, double xmin, double xmax)
{
	int i;
	float v1f[3], v2f[3], rgb[3];
	
	for(i=0; i<3; i++) {
		v1f[i]=(float) (v1[i]);
		v2f[i]=(float) (v2[i]);
	}
	getcolor(x1, xmin, xmax, rgb);
	glColor3f(rgb[0], rgb[1], rgb[2]);
	glVertex3f(v1f[0], v1f[2], -v1f[1]);	// GL(x,y,z)=Lab(x,z,-y)
	
	getcolor(x2, xmin, xmax, rgb);
	glColor3f(rgb[0], rgb[1], rgb[2]);
	glVertex3f(v2f[0], v2f[2], -v2f[1]);
}

// ---------- vesicle ------------
/*
void drawvesicle(void)
{
	float x, y, z;
	
	x=(float)(RCves[0]);
	y=(float)(RCves[1]);
	z=(float)(RCves[2]);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glCallList(vesicle);
	glTranslatef(-x, -z, y);	// GL(x,y,z)=Lab(x,z,-y)
}
*/


void drawvesnodes(void)
{
	int i;
	float xc, yc, zc;
	float x, y, z;
	VNODE *p;
	
	xc=(float)(RCves[0]);
	yc=(float)(RCves[1]);
	zc=(float)(RCves[2]);
	glTranslatef(xc, zc, -yc);
	
	glDisable(GL_LIGHTING);
	glPointSize(3.0);
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POINTS);
	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		x=(float)(1.002*p->rc[0]);	// slightly above surface to be noticable
		y=(float)(1.002*p->rc[1]);
		z=(float)(1.002*p->rc[2]);
		glVertex3f(x,z,-y);	// GL(x,y,z)=Lab(x,z,-y)
	}
	glEnd();
	glTranslatef(-xc, -zc, yc);
}



void drawvescenter(void)
{
	float xc, yc, zc;
	
	xc=(float)(RCves[0]);
	yc=(float)(RCves[1]);
	zc=(float)(RCves[2]);
	
	glDisable(GL_LIGHTING);
	glPointSize(30.0);
	glColor3f(1.0, 1.0, 1.0);
	
	glBegin(GL_POINTS);
	glVertex3f(xc,zc,-yc);	// GL(x,y,z)=Lab(x,z,-y)
	glEnd();
}



void DrawVesMesh(void)
{
	int i, j, i1, i2;
	double d2, rmin, rmax, rave, VesCurvmin, VesCurvmax;
	//double r1, r2;
	VNODE *p;
	
	rmin=inf;
	rmax=-inf;
	VesCurvmin=inf;
	VesCurvmax=-inf;
	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		if(showcolor==1) {	// based on vnode-center distance
			d2 = norm2(p->rc);
			if(d2<rmin) rmin=d2;
			if(d2>rmax) rmax=d2;
		}
		else if(showcolor==2) {
			if(p->H < VesCurvmin) VesCurvmin = p->H;
			if(p->H > VesCurvmax) VesCurvmax = p->H;
		}
	}
	if(showcolor==1) {
		rmin = sqrt(rmin);
		rmax = sqrt(rmax);
		if(rmax-rmin < 0.04*Rv) {
			rave = (rmin+rmax)/2.0;
			rmin = rave - 0.02*Rv;
			rmax = rave + 0.02*Rv;
		}
	}
	
	glDisable(GL_LIGHTING);
	glColor3f(0.7, 0.7, 0.7);
	glLineWidth(2.0);
	glBegin(GL_LINES);
	/*
	for(i=0; i<NVedge; i++) {
		i1=vedge[i].inode[0];
		i2=vedge[i].inode[1];
		drawline(vnode[i1].r, vnode[i2].r);
	}
	*/
	for(i=0; i<NVface; i++) {	// for each face
		for(j=0; j<3; j++) {	// 3 vertices
			i1 = vface[i].inode[j];
			i2 = vface[i].inode[((j+1)%3)];
			
			drawline(vnode[i1].r, vnode[i2].r);
			/*
			if(showcolor==0) drawline(vnode[i1].r, vnode[i2].r);
			else if(showcolor==1) {
				r1 = norm(vnode[i1].rc);
				r2 = norm(vnode[i2].rc);
				drawlineclr(vnode[i1].r, vnode[i2].r, r1, r2, rmin, rmax);
			}
			else drawlineclr(vnode[i1].r, vnode[i2].r, vnode[i1].H, vnode[i2].H, VesCurvmin, VesCurvmax);
			*/
		}
	}
	glEnd();
}



void DrawVesFace(void)
{
	int i, j, i1, i2, i3;
	float xc, yc, zc;
	float r1[3], r2[3], r3[3];
	VFACE *p;
	
	xc=(float)(RCves[0]);
	yc=(float)(RCves[1]);
	zc=(float)(RCves[2]);
	glTranslatef(xc, zc, -yc);	// GL(x,y,z)=Lab(x,z,-y)
	
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);
	for(i=0; i<NVface; i++) {
		p = &vface[i];
		i1 = p->inode[0];
		i2 = p->inode[1];
		i3 = p->inode[2];
		for(j=0; j<3; j++) {	// for x, y, z
			r1[j]=(float)(0.998*vnode[i1].rc[j]);	// offset to avoid aliasing
			r2[j]=(float)(0.998*vnode[i2].rc[j]);
			r3[j]=(float)(0.998*vnode[i3].rc[j]);
		}
		glBegin(GL_POLYGON);
			glVertex3f(r1[0], r1[2], -r1[1]);	// GL(x,y,z)=Lab(x,z,-y)
			glVertex3f(r2[0], r2[2], -r2[1]);
			glVertex3f(r3[0], r3[2], -r3[1]);
		glEnd();
	}
	glTranslatef(-xc, -zc, yc);
}



void DrawVesNormal(int status)
{
	int i;
	double normfct, r[3], len, lenhalf, r1[3], r2[3];
	
	normfct=0.08*Rv;
	len = 0.3*dRves;
	lenhalf = 0.5*len;
	
	glLineWidth(2);
	
	if(status==1) {	// face normal
		glColor3f(1, 1, 0);
		glBegin(GL_LINES);
		for(i=0; i<NVface; i++) {
			vecprod(vface[i].nrm, normfct, r);
			vecadd(vface[i].rclab, r, r);
			drawline(vface[i].rclab, r);
		}
		glEnd();
	}
	else if(status==2) {	// node normal
		glColor3f(1, 1, 0);
		glBegin(GL_LINES);
		for(i=0; i<NVnode; i++) {
			vecprod(vnode[i].e3, normfct, r);
			vecadd(vnode[i].r, r, r);
			drawline(vnode[i].r, r);
		}
		glEnd();
	}
	else {	// vnode e1-e2-e3
		glColor3f(1, 0, 0);	// e1: red
		glBegin(GL_LINES);
		for(i=0; i<NVnode; i++) {
			vecprod(vnode[i].e1, len, r);
			vecsub(vnode[i].r, r, r1);
			vecadd(vnode[i].r, r, r2);
			drawline(r1, r2);
		}
		glEnd();
		
		glColor3f(0, 1, 0);	// e2: green
		glBegin(GL_LINES);
		for(i=0; i<NVnode; i++) {
			vecprod(vnode[i].e2, len, r);
			vecsub(vnode[i].r, r, r1);
			vecadd(vnode[i].r, r, r2);
			drawline(r1, r2);
		}
		glEnd();
		
		glColor3f(0.5, 0.5, 1);	// e3: blue
		glBegin(GL_LINES);
		for(i=0; i<NVnode; i++) {
			vecprod(vnode[i].e3, len, r);
			vecsub(vnode[i].r, r, r1);
			vecadd(vnode[i].r, r, r2);
			drawline(r1, r2);
		}
		glEnd();
	}
}



void DrawVesForce(void)
{
	int i;
	double fct, r1[3], r2[3], dr[3];
	VNODE *p;
	
	fct=0.1;
	
	glDisable(GL_LIGHTING);
	glColor3f(1, 0, 0);
	glLineWidth(2.0);
	
	glBegin(GL_LINES);
	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		veccopy(p->r, r1);
		vecprod(p->f, fct, dr);
		vecadd(r1, dr, r2);
		drawline(r1, r2);
	}
	glEnd();
}



void MarkVes(void)
{
	int i, j, k, len;
	float xc, yc, zc;
	float x, y, z;
	char str[100];
	VNODE *p;
	
	xc=(float)(RCves[0]);
	yc=(float)(RCves[1]);
	zc=(float)(RCves[2]);
	glTranslatef(xc, zc, -yc);
	
	glColor3f(0.8, 0.8, 0.8);
	for(i=0; i<NVnode; i++) {
		p = &vnode[i];
		//sprintf(str, "(%.2g, %.2g)", p->C1, p->C2);
		//sprintf(str, "(%.3g, %.3g, %.5g)", p->r[0], p->r[1], p->r[2]);
		//sprintf(str, "%.2g", p->r[2]);
		sprintf(str, "%d", p->id);
		len=(int) strlen(str);
		x=(float)(1.01*(p->rc[0]));
		y=(float)(1.01*(p->rc[1]));
		z=(float)(1.01*(p->rc[2]));
		glRasterPos3f(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
		for(j=0; j<len; j++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[j]);
	}
	
	glColor3f(1, 0, 1);
	for(i=0; i<NVface; i++) {
		sprintf(str, "%d", i);
		len=(int) strlen(str);
		x=(float)(1.1*(vface[i].rcrc[0]));
		y=(float)(1.1*(vface[i].rcrc[1]));
		z=(float)(1.1*(vface[i].rcrc[2]));
		glRasterPos3f(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
		for(k=0; k<len; k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[k]);
	}
	
	/*
	i=0;	//FID[0];
	
	glColor3f(0.8, 0.8, 0);
	p=&vnode[i];
	for(j=0; j<p->nnbr; j++) {
		sprintf(str, "n%d", j);
		len=(int) strlen(str);
		x=(float)(1.05*(vnode[p->nbr[j]].rc[0]));
		y=(float)(1.05*(vnode[p->nbr[j]].rc[1]));
		z=(float)(1.05*(vnode[p->nbr[j]].rc[2]));
		glRasterPos3f(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
		for(k=0; k<len; k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[k]);
	}
	
	glColor3f(0.4, 0.4, 1);
	p=&vnode[i];
	for(j=0; j<p->nface; j++) {
		sprintf(str, "f%d", j);
		len=(int) strlen(str);
		x=(float)(1.05*(vface[p->iface[j]].rc[0]));
		y=(float)(1.05*(vface[p->iface[j]].rc[1]));
		z=(float)(1.05*(vface[p->iface[j]].rc[2]));
		glRasterPos3f(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
		for(k=0; k<len; k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[k]);
	}
	*/
	
	glTranslatef(-xc, -zc, yc);
}



// ---------- membrane ------------

void drawmemnodes(void)
{
	int i, j;
	float x, y, z;
	
	glPointSize(3.0);
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_POINTS);
		for(i=1; i<=Nx; i++) {
			for(j=1; j<=Ny; j++) {
				x=(float)(1.002*mnode[i][j].r[0]);	// slightly above surface to be noticable
				y=(float)(1.002*mnode[i][j].r[1]);
				z=(float)(1.002*mnode[i][j].r[2]);
				glVertex3f(x,z,-y);
			}
		}
	glEnd();
}


void DrawMemMesh(void)
{
	int i, j, k, i1, j1;
	double r1[3], r2[3];
	double zm;
	MNODE *m, *m1;
	
	if(!showcolor) glColor3f(0.4, 0.4, 0.4);
//	glDisable(GL_LINE_SMOOTH);
	glLineWidth(2.0);
	glBegin(GL_LINES);
	
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			m = &mnode[i][j];
			for(k=0; k<(m->nnbr); k++) {
				i1 = m->nbr[k][0];
				j1 = m->nbr[k][1];
				//if(i1>=1 && i1<=Nx && j1>=1 && j1<=Ny && i1>=i && j1>=j) {	// avoid double-counting
				if(i1>=1 && i1<=Nx && j1>=1 && j1<=Ny) {
					m1 = &mnode[i1][j1];
					veccopy(m->r, r1);
					veccopy(m1->r, r2);
					
					if(showcolor==0) drawline(r1, r2);	// no color
					else if(showcolor==1) {
						zm=max2(Zmax, Zmin+DLT);
						drawlineclr(r1, r2, r1[2], r2[2], Zmin, zm);	// show z-height
					}
					else drawlineclr(r1, r2, m->H, m1->H, Curvmin, Curvmax);	// show x-curvature
				}
			}
		}
	}
	
	glEnd();
//	glEnable(GL_LINE_SMOOTH);
}



void DrawMemFace(int color)	// background color 0=mono, 1=color
{
	int i, i1, i2, i3, j1, j2, j3;
	float zeps, r1[3], r2[3], r3[3];
	float rgb1[3], rgb2[3], rgb3[3];
	MFACE *f;
	
	zeps=0.1;	// +DLThalf;
	
	glDisable(GL_LIGHTING);
	if(color==0) {	// mono background
		//glColor3f(0.1, 0.1, 0.1);
		glColor3f(0, 0, 0);
	}
	
	for(i=0; i<NMface; i++) {
		f = &mface[i];
		i1 = f->inode[0][0];	j1 = f->inode[0][1];
		i2 = f->inode[1][0];	j2 = f->inode[1][1];
		i3 = f->inode[2][0];	j3 = f->inode[2][1];
		
		if(i1>0 && i2>0 && i3>0 && j1>0 && j2>0 && j3>0) {
			r1[0]=(float)(mnode[i1][j1].r[0]);
			r1[1]=(float)(mnode[i1][j1].r[1]);
			r1[2]=(float)(mnode[i1][j1].r[2]-zeps);

			r2[0]=(float)(mnode[i2][j2].r[0]);
			r2[1]=(float)(mnode[i2][j2].r[1]);
			r2[2]=(float)(mnode[i2][j2].r[2]-zeps);

			r3[0]=(float)(mnode[i3][j3].r[0]);
			r3[1]=(float)(mnode[i3][j3].r[1]);
			r3[2]=(float)(mnode[i3][j3].r[2]-zeps);
			
			if(color) {
				getcolor(mnode[i1][j1].r[2], Zmin, max2(Zmax,0.1), rgb1);
				getcolor(mnode[i2][j2].r[2], Zmin, max2(Zmax,0.1), rgb2);
				getcolor(mnode[i3][j3].r[2], Zmin, max2(Zmax,0.1), rgb3);
			}
			
			glBegin(GL_POLYGON);
				if(color) glColor3f(rgb1[0], rgb1[1], rgb1[2]);
				glVertex3f(r1[0], r1[2], -r1[1]);	// GL(x,y,z)=Lab(x,z,-y)
				
				if(color) glColor3f(rgb2[0], rgb2[1], rgb2[2]);
				glVertex3f(r2[0], r2[2], -r2[1]);
				
				if(color) glColor3f(rgb3[0], rgb3[1], rgb3[2]);
				glVertex3f(r3[0], r3[2], -r3[1]);
			glEnd();
		}
	}
}



void DrawMemDirection(void)
{
	int i, j;
	double len, lenhalf, r1[3], r2[3], dr[3];
	
	len = 0.5*dX;
	lenhalf = 0.5*len;
	
	glDisable(GL_LIGHTING);
	glLineWidth(2);
	
	glColor3f(1, 0, 0);	// e1: red
	glBegin(GL_LINES);
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			vecprod(mnode[i][j].e1, len, dr);
			vecsub(mnode[i][j].r, dr, r1);
			vecadd(mnode[i][j].r, dr, r2);
			drawline(r1, r2);
		}
	}
	glEnd();
	
	glColor3f(0, 1, 0);	// e2: green
	glBegin(GL_LINES);
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			vecprod(mnode[i][j].e2, len, dr);
			vecsub(mnode[i][j].r, dr, r1);
			vecadd(mnode[i][j].r, dr, r2);
			drawline(r1, r2);
		}
	}
	glEnd();
	
	glColor3f(0.6, 0.6, 1);	// e3: blue
	glBegin(GL_LINES);
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			vecprod(mnode[i][j].e3, len, dr);
			veccopy(mnode[i][j].r, r1);
			vecadd(r1, dr, r2);
			drawline(r1, r2);
		}
	}
	glEnd();
}



void DrawMemForce(void)
{
	int i, j;
	double fct, r1[3], r2[3], dr[3];
	MNODE *m;
	
	fct=0.1;
	glDisable(GL_LIGHTING);
	glLineWidth(2);
	glColor3f(1, 1, 1);
	glBegin(GL_LINES);
	
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			m = &mnode[i][j];
			veccopy(m->r, r1);
			vecprod(m->f, fct, dr);
			vecadd(r1, dr, r2);
			drawline(r1, r2);
		}
	}
	glEnd();
}


void MarkMemNode(void)
{
	int i, j, k, len;
	float x, y, z;
	char str[50];
	MNODE *m;
	
	glColor3f(0.5, 0.5, 0.5);
	for(i=1; i<=Nx; i++) {
		for(j=1; j<=Ny; j++) {
			m = &mnode[i][j];
			//sprintf(str, "%.2g", m->H);
			//sprintf(str, "%.2g", m->LaplaceH);
			//sprintf(str, "(%.3g, %.3g, %.5g)", mnode[i][j].r[0], mnode[i][j].r[1], mnode[i][j].r[2]);
			//sprintf(str, "%.2g", mnode[i][j].r[2]);
			//sprintf(str, "%.4g", mnode[i][j].area);
			sprintf(str, "(%d, %d)", i, j);
			//sprintf(str, "%d", m->nnbr);
			len=(int) strlen(str);
			x=(float)(mnode[i][j].r[0]);
			y=(float)(mnode[i][j].r[1]);
			z=(float)(mnode[i][j].r[2]);
			glRasterPos3f(x, z+1, -y);	// GL(x,y,z)=Lab(x,z,-y)
			for(k=0; k<len; k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[k]);
		}
	}
}


void MarkMemFace(void)
{
	int i, j, len;
	float x, y, z;
	char str[50];
	MFACE *f;
	
	glColor3f(0.5, 0.5, 0.5);
	for(i=0; i<NMface; i++) {
		f = &mface[i];
		sprintf(str, "%.4g", f->area);
		//sprintf(str, "%d", f->id);
		len=(int) strlen(str);
		x=(float)(f->rc[0]);
		y=(float)(f->rc[1]);
		z=(float)(f->rc[2]);
		glRasterPos3f(x, z+1, -y);	// GL(x,y,z)=Lab(x,z,-y)
		for(j=0; j<len; j++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[j]);
	}
}


void MarkMem(void)
{
	MarkMemNode();
	//MarkMemFace();
}



//--------------- tmd & linkers below -------------
void DrawTMD(void)
{
	int i;
	float x, y, z;
	//double r[3];
	CHAIN *p;
	
	glColor3f(0.8, 0.8, 0);
	glLineWidth(3);
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) drawline(p->rlinker[i], p->rtmd[i]);
		p = p->next;
	}
	glEnd();
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {	// GL(x,y,z)=Lab(x,z,-y)
			x=(float)(p->rtmd[i][0]);
			y=(float)(p->rtmd[i][1]);
			z=(float)(p->rtmd[i][2]);
			//glVertex3f(x,z,-y);
			glTranslatef(x,z,-y);
			glCallList(tmd);
			glTranslatef(-x,-z,y);
		}
		p = p->next;
	}
	
	/*
	glColor3f(1, 0.2, 1);
	glLineWidth(3);
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			vecprod(p->nrmtmd[i], 5.0, r);
			vecadd(p->rtmd[i], r, r);
			drawline(p->rtmd[i], r);
		}
		p = p->next;
	}
	glEnd();
	*/
}




//--------------- syt below --------------------

void alignaxis(double n1[3], double n2[3], double axis[3], double *a)	// to rotate n1 to n2: about axis by angle a
{
	double r[3];
	*a = getanglePi(n1, n2);	// a in rad
	crossprod(n1, n2, r);
	normalize(r, axis);	// this also prevents r=0
}


void DrawSYT(void)
{
	int i, j;
	float x, y, z;
	float thf, phif;
	float axisf[3];
	double nx[3], ny[3], nz[3];
	double th, phi, sn, cs, axis[3], r1[3], r2[3];
	CHAIN *p;
	
	nx[1] = nx[2] = 0;	nx[0] = 1;
	ny[0] = ny[2] = 0;	ny[1] = 1;
	nz[0] = nz[1] = 0;	nz[2] = 1;
	
	//glLineWidth(1);
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			x = (float)(p->r[i][0]);
			y = (float)(p->r[i][1]);
			z = (float)(p->r[i][2]);
			
			alignaxis(nz, p->nz[i], axis, &th);	// rotate lab-z directly to p->nz, around axis by th in rad
			for(j=0; j<3; j++) axisf[j] = (float)(axis[j]);
			
			sn=sin(th);
			cs=cos(th);
			RotMatrix3(sn, cs, axis, nx, r1);	// rotate lab-x in the same way
			phi=getanglePi(r1, p->nx[i]);	// phi is the angle to be rotated about z
			crossprod(r1, p->nx[i], r2);	// r2 is the axis about which r1 can be rotated to p->nx
			if(dotprod(r2, p->nz[i])<0) phi*=-1;
			
			thf=(float)(th*r2d);	// angle in degrees
			phif=(float)(phi*r2d);
			
			glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
			
			//glRotatef(phif, p->nz[i][0], p->nz[i][2], -p->nz[i][1]);	// either this, or 2 lines below will do the rotation
			glRotatef(thf, axisf[0], axisf[2], -axisf[1]);	// GL(x,y,z)=Lab(x,z,-y)
			glRotatef(phif, 0, 1, 0);
			
			if(p->closed==0) glCallList(C2AB);
			else glCallList(C2ABclose);
			/*
			if(p->closed == 0) {
				if(showboxsphere) glCallList(monomer1A);
				else glCallList(monomer1B);
			}
			else {
				if(showboxsphere) glCallList(monomer2A);
				else glCallList(monomer2B);
			}
			*/
			
			glRotatef(-phif, 0, 1, 0);
			glRotatef(-thf, axisf[0], axisf[2], -axisf[1]);
			//glRotatef(-phif, p->nz[i][0], p->nz[i][2], -p->nz[i][1]);
				
			glTranslatef(-x, -z, y);
		}
		p = p->next;
	}
}



void DrawSytBeads(void)
{
	int i, j, k, len;
	char str[50];
	float x, y, z;
	CHAIN *p;
	
	glColor3f(1, 1, 1);
	glPointSize(5);
	glBegin(GL_POINTS);
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {	// for each segment
			for(j=0; j<4; j++) {	// for each C2AB bead
				x=(float)(p->rC2AB[i][j][0]);
				y=(float)(p->rC2AB[i][j][1]);
				z=(float)(p->rC2AB[i][j][2]);
				glVertex3f(x,z,-y);	// GL(x,y,z)=Lab(x,z,-y)
			}
			
			x=(float)(p->rsm[i][0]);
			y=(float)(p->rsm[i][1]);
			z=(float)(p->rsm[i][2]);
			glVertex3f(x,z,-y);	// GL(x,y,z)=Lab(x,z,-y)
		}
		p = p->next;
	}
	
	glEnd();
	
	
	
	glColor3f(1, 1, 1);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {	// for each segment
			for(j=0; j<4; j++) {	// for each C2AB bead
				sprintf(str, "%d", j);
				len=(int) strlen(str);
				x=(float)(p->rC2AB[i][j][0]);
				y=(float)(p->rC2AB[i][j][1]);
				z=(float)(p->rC2AB[i][j][2]);
				glRasterPos3f(x, z+0.2, -y);	// GL(x,y,z)=Lab(x,z,-y)
				for(k=0; k<len; k++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[k]);
			}
		}
		p = p->next;
	}
}



void DrawSYTlinks(void)
{
	int i, ip1, imax;
	CHAIN *p;
	
	//glEnable(GL_LINE_SMOOTH);
	glColor3f(0.5, 0.5, 1);
	glLineWidth(4);
	glBegin(GL_LINES);
	
	p = chain;
	while(p) {
		if(p->closed == 0) imax=(p->nseg)-1;
		else imax=(p->nseg);
		for(i=0; i<imax; i++) {
			ip1 = i+1;
			if(ip1>=(p->nseg)) ip1=0;
			drawline(p->r[i], p->r[ip1]);
		}
		p = p->next;
	}
	glEnd();
}



void DrawSYTnrm(int s)
{
	int i;
	double r[3];
	CHAIN *p;
	
	glLineWidth(2);
	
	glColor3f(1, 0, 0);	// nx: red
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			if(s==1) vecprod(p->nx[i], DLT, r);
			else vecprod(p->nxa[i], DLT, r);
			vecadd(p->r[i], r, r);
			drawline(p->r[i], r);
		}
		p = p->next;
	}
	glEnd();
	
	glColor3f(0, 1, 0);	// ny: green
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			if(s==1) vecprod(p->ny[i], DLT, r);
			else vecprod(p->nya[i], DLT, r);
			vecadd(p->r[i], r, r);
			drawline(p->r[i], r);
		}
		p = p->next;
	}
	glEnd();
	
	glColor3f(0.5, 0.5, 1);	// nz: blue
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			if(s==1) vecprod(p->nz[i], DLT, r);
			else vecprod(p->nza[i], DLT, r);
			vecadd(p->r[i], r, r);
			drawline(p->r[i], r);
		}
		p = p->next;
	}
	glEnd();
}



void DrawSYTnnxt(void)
{
	int i;
	double r[3];
	CHAIN *p;
	
	glLineWidth(2);
	
	glColor3f(1, 0, 0);	// nx: red
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			vecprod(p->nnxt[i], DLT, r);
			vecadd(p->r[i], r, r);
			drawline(p->r[i], r);
		}
		p = p->next;
	}
	glEnd();
}




void DrawSYTforce(void)
{
	int i, j;
	double fct, r1[3], r2[3], dr[3];
	float x, y, z;
	CHAIN *p;
	
	fct=0.1;
	glLineWidth(2);
	
	glColor3f(1, 1, 1);	// white: attraction force
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			veccopy(p->rsm[i], r1);
			vecprod(p->fatt[i], fct, dr);
			//printf("%d-%d\t%.3g\t%.3g\t%.3g\n", i, s, dr[0], dr[1], dr[2]);
			vecadd(r1, dr, r2);
			drawline(r1, r2);
		}
		p = p->next;
	}
	glEnd();
				
	glColor3f(0.5, 0.5, 1);	// blue: repulsion force
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			for(j=0; j<4; j++) {
				veccopy(p->rC2AB[i][j], r1);
				vecprod(p->frep[i][j], fct, dr);
				//printf("%d-%d\t%.3g\t%.3g\t%.3g\n", i, s, dr[0], dr[1], dr[2]);
				vecadd(r1, dr, r2);
				drawline(r1, r2);
			}
		}
		p = p->next;
	}
	glEnd();
				
	glColor3f(1, 0.2, 0.2);	// red: total external force
	glBegin(GL_LINES);
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			veccopy(p->r[i], r1);
			vecprod(p->fext[i], fct, dr);
			//printf("%.3g\t", norm(p->fext[i]));
			vecadd(r1, dr, r2);
			drawline(r1, r2);
		}
		p = p->next;
	}
	glEnd();
	
	
	glColor3f(1, 1, 0);
	glPointSize(5);
	glBegin(GL_POINTS);
	
	p = chain;	// mark site
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			x=(float)(p->rsm[i][0]);	// GL(x,y,z)=Lab(x,z,-y)
			y=(float)(p->rsm[i][1]);
			z=(float)(p->rsm[i][2]);
			glVertex3f(x,z,-y);
		}
		p = p->next;
	}
	glEnd();
}



void DrawSYTstretchForce(void)
{
	int i;
	double fct, r1[3], r2[3], dr[3];
	float x, y, z;
	CHAIN *p;
	
	fct=10.1;
	
	//glColor3f(1, 1, 1);
	glLineWidth(1);
	glBegin(GL_LINES);
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			veccopy(p->rc1[i], r1);
			vecprod(p->dr1[i], fct, dr);
			vecadd(r1, dr, r2);
			glColor3f(1, 0, 0);
			drawline(r1, r2);
			
			veccopy(p->rc2[i], r1);
			vecprod(p->dr2[i], fct, dr);
			vecadd(r1, dr, r2);
			glColor3f(0, 0, 1);
			drawline(r1, r2);
		}
		p = p->next;
	}
	glEnd();
	
	glColor3f(1, 1, 0);
	glPointSize(5);
	glBegin(GL_POINTS);
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			x=(float)(p->rc1[i][0]);
			y=(float)(p->rc1[i][1]);
			z=(float)(p->rc1[i][2]);
			glVertex3f(x,z,-y);
			x=(float)(p->rc2[i][0]);
			y=(float)(p->rc2[i][1]);
			z=(float)(p->rc2[i][2]);
			glVertex3f(x,z,-y);
		}
		p = p->next;
	}
	glEnd();
}




void DrawSYTtorque(void)
{
	int i;
	double fct, r1[3], r2[3], dr[3];
	CHAIN *p;
	
	fct=0.02;
	glColor3f(1, 1, 1);
	glLineWidth(1);
	glBegin(GL_LINES);
	
	p = chain;
	while(p) {
		for(i=0; i<(p->nseg); i++) {
			veccopy(p->r[i], r1);
			vecprod(p->tauorient[i], fct, dr);
			//vecprod(p->taust[i], fct, dr);
			vecadd(r1, dr, r2);
			drawline(r1, r2);
		}
		p = p->next;
	}
	glEnd();
}



void MarkSYTchain(void)
{
	float r[3], color[3];
	char str[20];
	CHAIN *p;
	
	color[0]=color[2]=0; color[1]=1;
	p = chain;
	while(p) {	// GL(x,y,z)=Lab(x,z,-y)
		r[0] = (float)(p->rcenter[0]);
		r[1] = (float)(p->rcenter[2]);
		r[2] = -(float)(p->rcenter[1]);
		sprintf(str, "%d", p->id);
		//sprintf(str, "%d", p->closed);
		drawstring3(str, r, color, font1);
		p = p->next;
	}
}


void MarkSYTseg(void)
{
	int i, n;
	float r[3], color[3];
	char str[100];
	CHAIN *p;
	
	color[0]=color[1]=color[2]=0.8;
	p = chain;
	while(p) {
		n = p->nseg;
		for(i=0; i<n; i++) {	// GL(x,y,z)=Lab(x,z,-y)
			r[0] = (float)(p->r[i][0] + DLT*(p->ny[i][0]));
			r[1] = (float)(p->r[i][2] + DLT*(p->ny[i][2]));
			r[2] = -(float)(p->r[i][1] + DLT*(p->ny[i][1]));
			//sprintf(str, "%d", i+1);
			sprintf(str, "%.1f", p->Em[i]/kBT);
			//sprintf(str, "(%.1f, %.1f, %.2f)", p->r[i][0],p->r[i][1],p->r[i][2]);
			//sprintf(str, "(%.3g, %.3g, %.3g)", p->taust[i][0], p->taust[i][1], p->taust[i][2]);
			//sprintf(str, "%d:(%.3g, %.3g, %.3g)", i, p->tauorient[i][0], p->tauorient[i][1], p->tauorient[i][2]);
			drawstring3(str, r, color, font2);
		}
		p = p->next;
	}
}

//-------- general ----------------

void DrawCoord(void)
{
	float r;
	r=(float)(0.6*max2(Lx, Ly));
	glLineWidth(5.0);
	glBegin(GL_LINES);	// GL(x,y,z)=Lab(x,z,-y), Lab(x,y,z)=GL(x,-z,y)
		glColor3f(1,0,0);	// lab x-axis
		glVertex3f(0, 0, 0);
		glVertex3f(r, 0, 0);
		
		glColor3f(0,1,0);	// lab y-axis
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, -r);
		
		glColor3f(0,0,1);	// lab z-axis
		glVertex3f(0, 0, 0);
		glVertex3f(0, r, 0);
	glEnd();
}


void drawall(void)
{
	glScalef(zoom, zoom, zoom);
	glTranslatef(PanX, PanY, 0.0);
	glRotatef(rottheta,1,0,0);
	glRotatef(rotphi,0,1,0);
	
#if BNDCND
	//DrawBottom();
#endif
	
	if(showsyt) {	// syt
		DrawSYT();
		//DrawSytBeads();
	}
	
	if(showtmd) DrawTMD();
	
	if(Membrane && showmem) {
		DrawMemMesh();	// membrane mesh
		if(showden) DrawMemFace(1);	// membrane face
		else DrawMemFace(0);
	}
	
	if(Vesicle && showvesicle) {	// vesicle
		DrawVesFace();
		DrawVesMesh();
		//drawvescenter();
	}
	
	if(shownormal==1) {
		DrawSYTlinks();
		DrawSYTnrm(1);
		//DrawSYTnrm(2);
		//DrawSYTnnxt();
	}
	else if(shownormal==2) {
		if(Membrane) DrawMemDirection();
	}
	else if(shownormal==3) {
		//DrawVesNormal(1);	// face normal
		if(Vesicle) DrawVesNormal(2);	// node normal
	}
	
	if(showforce) {
		/*
		if(showforce==1) DrawSYTforce();
		else if(showforce==2 && Membrane) DrawMemForce();
		else if(showforce==3 && Vesicle) DrawVesForce();
		*/
		DrawSYTtorque();
		DrawSYTstretchForce();
	}
	
	if(showmark==1) {
		MarkSYTchain();
		MarkSYTseg();
	}
	else if(showmark==2 && Membrane) MarkMem();
	else if(showmark==3 && Vesicle) MarkVes();
	
	if(showcoord) DrawCoord();
	
	glLoadIdentity();
}



void showtext1(void)
{
	float pos[3], color[3];
	char str[100];
	
	color[0]=color[1]=color[2]=0.8;
	pos[0]=-1.1*Lx;
	pos[1]=0.72*Ly;
	pos[2]=0;
	if(t<1e-6) sprintf(str, "t = %.3g ns,  H = %.2f nm", t*1e9, Zmax);
	else if(t<1e-3) sprintf(str, "t = %.3g us,  H = %.2f nm", t*1e6, Zmax);
	else if(t<1) sprintf(str, "t = %.3g ms,  H = %.2f nm", t*1e3, Zmax);
	else sprintf(str, "t = %.3g s,  H = %.2f nm", t, Zmax);
	drawstring3(str, pos, color, font);
	
	pos[0]=0.65*Lx;
	pos[1]=0.72*Ly;
	pos[2]=0;
	//sprintf(str, "Chain = %d :  closed = %d,  open = %d", Nchain, Nclosed, Nopen);
	sprintf(str, "SYT:  closed = %d,  open = %d", Nclosed, Nopen);
	drawstring3(str, pos, color, font);
}



void showtext(void)
{
	//int i;
	float pos[3], color[3];
	char str[200];
	//double fext;
	//CHAIN *p;
	
	color[0]=color[1]=color[2]=0.8;
	pos[0]=-1.05*Lx;
	pos[1]=0.8*Ly;
	pos[2]=0;
	if(t<1e-6) sprintf(str, "t = %.3g ns,", t*1e9);
	else if(t<1e-3) sprintf(str, "t = %.3g us,", t*1e6);
	else if(t<1) sprintf(str, "t = %.3g ms,", t*1e3);
	else sprintf(str, "t = %.3g s,", t);
	drawstring3(str, pos, color, font2);
	
	pos[0]=-0.85*Lx;
	pos[1]=0.8*Ly;
	pos[2]=0;
	sprintf(str, "dH = %.2f nm,  Zsyt = %.3g nm,", Zmax-Zmin, chain->rcenter[2]);
	drawstring3(str, pos, color, font2);
	
	pos[0]=-0.35*Lx;
	pos[1]=0.8*Ly;
	pos[2]=0;
	if(Membrane && Vesicle) sprintf(str, "Esyt = %.3g kT,  Ees = %.3g kT,  Emem = %.3g kT, Eves = %.3g kT, Etot = %.3g kT", 
		Esyt/kBT, EsytES/kBT, Emem/kBT, Eves/kBT, (EsytES+Emem+Eves)/kBT);
	else if(Membrane) sprintf(str, "Esyt = %.3g kT,  Ees = %.3g kT,  Emb = %.3g kT,  Emt = %.3g kT,  Emem = %.3g kT,  Etot = %.3g kT", 
		Esyt/kBT, EsytES/kBT, Emembend/kBT, Ememtens/kBT, Emem/kBT, (EsytES+Emem+Eves)/kBT);
	else if(Vesicle) sprintf(str, "Esyt = %.3g kT,  Ees = %.3g kT,  Evb = %.3g kT,  Evt = %.3g kT,  Eves = %.3g kT,  Etot = %.3g kT", 
		Esyt/kBT, EsytES/kBT, Evesbend/kBT, Evestens/kBT, Eves/kBT, (EsytES+Emem+Eves)/kBT);
	else sprintf(str, "Esyt = %.3g kT,  Ees = %.3g kT,  Etot = %.3g kT", 
		Esyt/kBT, EsytES/kBT, (EsytES+Emem+Eves)/kBT);
		
	drawstring3(str, pos, color, font2);
	
	pos[0]=-1.05*Lx;
	pos[1]=-0.85*Ly;
	pos[2]=0;
	
	if(Vesicle) sprintf(str, "Av00 = %.4g,  Av0 = %.4g,  Av = %.4g", AreaVes0est, AreaVes0, AreaVes);
	else sprintf(str, "A0 = %.4g,  A = %.4g", AreaMem0, AreaMem);
		
	drawstring3(str, pos, color, font2);
	
	/*
	color[0]=color[1]=color[2]=0.8;
	pos[0]=0.0*Lx;
	pos[1]=0.7*Ly;
	pos[2]=0;
	fext=0;
	
	p = chain;
	for(i=0; i<(p->nseg); i++) fext += p->fext[i][2];
	
	//sprintf(str, "Fsyt = %.3g, Fves = %.3g, dF = %.3g", fext, Fves[2], fext+Fves[2]);
	if(Vesicle) sprintf(str, "Fsyt = %.3g, Fmem = %.3g, Fves = %.3g, dF = %.3g", 
		TrackFsyt, TrackFmem, TrackFves, TrackFsyt+TrackFmem+TrackFves);
	else sprintf(str, "Fsyt = %.3g, Fmem = %.3g, dF = %.3g", TrackFsyt, TrackFmem, TrackFsyt+TrackFmem);
	drawstring3(str, pos, color, font);
*/
	
}




void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	
	//rotphi+=0.2;
	
	if(Membrane) GetMemMaxMin(mnode);
	
	AllEnergy();
	
	//showtext1();
	showtext();
	
	drawall();

	//pause=1;	// test !!!!!!!!!!!!!!!!!!!

//-------- snapshots --------
#if SAVPIC
	if(img_cnt%img_cntm==0) {
		nframe++;
		printf("frame = %d\n", nframe);
		sprintf(img_name, "frames/p%.4u.png", nframe);
		img_data=(GLubyte*)malloc((Width)*(Height)*3*sizeof(GLubyte));
		ilTexImage(Width, Height, 1, 3, GL_RGB, GL_UNSIGNED_BYTE, img_data);
	
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0, 0, Width, Height, GL_RGB, GL_UNSIGNED_BYTE, img_data);

		ilSetData(img_data);
		ilSave(IL_PNG, img_name);
		free(img_data);
	}
	img_cnt++;
#endif
//---------------------------

	glPopMatrix();
	glutSwapBuffers();
}



void GLmain(void)
{
#if LGV_MC
	int i;
#endif
	
	if(Membrane) UpdateMem(mnode, mface);
	if(Vesicle) {
		UpdateVes(vnode, vface);
		UpdateVesFaceCoord(vnode, vface);
	}
	
	display();
	
	if(pause==0) {
	#if LGV_MC
		for(i=0; i<showskip; i++) DoSimu1();	// Langevin
	#else
		DoSimu2();	// Monte-Carlo
	#endif
	}
}


void RunGL(int *argc, char **argv)
{
	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(Width,Height);
	glutCreateWindow("SYT + Membrane");
	
	init_GL();
#if SAVPIC
	ilInit();	// initializing IL
#endif
	display_list();
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	glutMotionFunc(drag);
	glutReshapeFunc(changeSize);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(arrow_keys);
	glutDisplayFunc(GLmain);
	glutIdleFunc(GLmain);
	glEnable(GL_DEPTH_TEST);
	glutMainLoop();
}


