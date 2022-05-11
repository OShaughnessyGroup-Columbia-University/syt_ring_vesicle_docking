
void drawline(double v1[], double v2[])
{
	int i;
	float v1f[3], v2f[3];
	
	for(i=0; i<3; i++) {
		v1f[i]=(float) (v1[i]);
		v2f[i]=(float) (v2[i]);
	}
	glVertex3f(v1f[0], v1f[2], -v1f[1]);	// GL(x,y,z)=Lab(x,z,-y)
	glVertex3f(v2f[0], v2f[2], -v2f[1]);
}



void drawstring2(const char *str, float pos[2], float color[3], void *font)
{
	glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
	glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
	//glDisable(GL_TEXTURE_2D);

	glColor3fv(color);          // set text color
	glRasterPos2fv(pos);        // place text position

	// loop all characters in the string
	while(*str) {
		glutBitmapCharacter(font, *str);
		++str;
		}

	// glEnable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glPopAttrib();
}



void drawstring3(const char *str, float pos[3], float color[3], void *font)
{
	glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
	glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
	//glDisable(GL_TEXTURE_2D);

	glColor3fv(color);          // set text color
	glRasterPos3fv(pos);        // place text position

	// loop all characters in the string
	while(*str) {
		glutBitmapCharacter(font, *str);
		++str;
		}

	// glEnable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glPopAttrib();
}



void display_syt_cube(void)
{
	float r;
	float a, b, r1, r2;
	
	r=(float)(1*DLT);
	a=0.25;	// thickness of blue region
	b=1-a;
	r1=a/2*r;	// center of non-interaction region
	r2=(a-1)/2*r;	// center of interaction region
	
	
	monomersideA=glGenLists(1);	// interaction side: blue
	glNewList(monomersideA, GL_COMPILE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	//glEnable(GL_LIGHT3);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blueDiffuseMaterial);
	glScalef(1, 1, a);
	//glutSolidSphere(r, 30, 30);
	glutSolidCube(r);
	//glutWireCube(r);
	glScalef(1, 1, 1/a);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHTING);
	glEndList();
	
	
	monomer1A=glGenLists(1);	// molecules in chain: green
	glNewList(monomer1A, GL_COMPILE);
	glTranslatef(0, 0, -r1);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glScalef(1, 1, b);
	//glutSolidSphere(r, 30, 30);
	glutSolidCube(r);
	glScalef(1, 1, 1/b);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(0, 0, r1);
	
	glTranslatef(0, 0, -r2);
	glCallList(monomersideA);
	glTranslatef(0, 0, r2);
	glEndList();
	
	
	monomer2A=glGenLists(1);	// molecules in ring: red
	glNewList(monomer2A, GL_COMPILE);
	glTranslatef(0, 0, -r1);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	//glEnable(GL_LIGHT2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glScalef(1, 1, b);
	//glutSolidSphere(r, 30, 30);
	glutSolidCube(r);
	glScalef(1, 1, 1/b);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHTING);
	glTranslatef(0, 0, r1);
	
	glTranslatef(0, 0, -r2);
	glCallList(monomersideA);
	glTranslatef(0, 0, r2);
	glEndList();
}



void display_syt_sphere(void)
{
	float r1, r2, z;
	
	r1=(float)(0.15*DLT);	// radius of membrane-binding site
	r2=(float)(DLThalf);	// radius of syt
	
	monomersideB=glGenLists(1);	// interaction side: blue
	glNewList(monomersideB, GL_COMPILE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	//glEnable(GL_LIGHT3);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blueDiffuseMaterial);
	glutSolidSphere(r1, 30, 30);
	//glutWireCube(r1);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHTING);
	glEndList();
	
	z = (float)(DLThalf);
	monomer1B=glGenLists(1);	// molecules in chain: green
	glNewList(monomer1B, GL_COMPILE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glutSolidSphere(r2, 30, 30);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	
	glTranslatef(0, 0, z);	// GL(x,y,z)=Lab(x,z,-y)
	glCallList(monomersideB);
	glTranslatef(0, 0, -z);
	glEndList();
	
	monomer2B=glGenLists(1);	// molecules in ring: red
	glNewList(monomer2B, GL_COMPILE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	//glEnable(GL_LIGHT2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glutSolidSphere(r2, 30, 30);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHTING);
	
	glTranslatef(0, 0, z);
	glCallList(monomersideB);
	glTranslatef(0, 0, -z);
	glEndList();
}



void display_syt_cylinders(void)
{
	float rsite, ra, rb;
	float x, y, z;
	GLUquadricObj *qobj;
	
	qobj = gluNewQuadric();
	gluQuadricDrawStyle(qobj, GLU_FILL); /* flat shaded */
	gluQuadricNormals(qobj, GLU_FLAT);
	
	rsite=(float)(0.4*Rc2abAVE);	// radius of membrane-binding site
	ra=(float)(Rc2a);	// radius of syt beads
	rb=(float)(Rc2b);	// radius of syt beads
	
	// lysine patch: blue
	x=(float)(RCsm[0]);
	y=(float)(RCsm[1]);
	z=(float)(RCsm[2]);
	
	LYS=glGenLists(1);
	glNewList(LYS, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blueDiffuseMaterial);
	glutSolidSphere(rsite, 20, 20);
	//glutWireSphere(rsite, 10, 10);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-1 open
	x=(float)(RCc2b1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b1[1]);
	z=(float)(RCc2b1[2]);
	
	C2B1=glGenLists(1);
	glNewList(C2B1, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glutSolidSphere(rb, 30, 30);
	//glutWireSphere(rb, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-1 close
	x=(float)(RCc2b1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b1[1]);
	z=(float)(RCc2b1[2]);
	
	C2B1close=glGenLists(1);
	glNewList(C2B1close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glutSolidSphere(rb, 30, 30);
	//glutWireSphere(rb, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-2 open
	x=(float)(RCc2b2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b2[1]);
	z=(float)(RCc2b2[2]);
	
	C2B2=glGenLists(1);
	glNewList(C2B2, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glutSolidSphere(rb, 30, 30);
	//glutWireSphere(rb, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-2 close
	x=(float)(RCc2b2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b2[1]);
	z=(float)(RCc2b2[2]);
	
	C2B2close=glGenLists(1);
	glNewList(C2B2close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glutSolidSphere(rb, 30, 30);
	//glutWireSphere(rb, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-3 open: cylinder
	x=(float)(RCc2b1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b1[1]);
	z=(float)(RCc2b1[2]);
	//x=(float)(-dc2bhalf);	// GL(x,y,z)=Lab(x,z,-y)
	//y=(float)(-Rc2b);
	//z=0;
	
	C2B3=glGenLists(1);	// cylinder of C2B
	glNewList(C2B3, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glRotatef(90, 0, 1, 0);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	gluCylinder(qobj, rb, rb, 2*dc2bhalf, 30, 2);	// gluCylinder(qobj, Rbase, Rtop, Height, Slice, Stack) along z!
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glRotatef(-90, 0, 1, 0);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2B-3 close: cylinder
	x=(float)(RCc2b1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2b1[1]);
	z=(float)(RCc2b1[2]);
	//x=(float)(-dc2bhalf);	// GL(x,y,z)=Lab(x,z,-y)
	//y=(float)(-Rc2b);
	//z=0;
	
	C2B3close=glGenLists(1);	// cylinder of C2B
	glNewList(C2B3close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glRotatef(90, 0, 1, 0);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	gluCylinder(qobj, rb, rb, 2*dc2bhalf, 30, 2);	// gluCylinder(qobj, Rbase, Rtop, Height, Slice, Stack) along z!
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glRotatef(-90, 0, 1, 0);
	glTranslatef(-x, -z, y);
	glEndList();
	
	
	// C2A-1 open
	x=(float)(RCc2a1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a1[1]);
	z=(float)(RCc2a1[2]);
	
	C2A1=glGenLists(1);
	glNewList(C2A1, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glutSolidSphere(ra, 30, 30);
	//glutWireSphere(ra, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	
	// C2A-1 close
	x=(float)(RCc2a1[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a1[1]);
	z=(float)(RCc2a1[2]);
	
	C2A1close=glGenLists(1);
	glNewList(C2A1close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glutSolidSphere(ra, 30, 30);
	//glutWireSphere(ra, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	
	// C2A-2 open
	x=(float)(RCc2a2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a2[1]);
	z=(float)(RCc2a2[2]);
	
	C2A2=glGenLists(1);
	glNewList(C2A2, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	glutSolidSphere(ra, 30, 30);
	//glutWireSphere(ra, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	
	// C2A-2 close
	x=(float)(RCc2a2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a2[1]);
	z=(float)(RCc2a2[2]);
	
	C2A2close=glGenLists(1);
	glNewList(C2A2close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	glutSolidSphere(ra, 30, 30);
	//glutWireSphere(ra, 10, 10);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glTranslatef(-x, -z, y);
	glEndList();
	
	
	// C2A-3 open: cylinder
	x=(float)(RCc2a2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a2[1]);
	z=(float)(RCc2a2[2]);
	
	C2A3=glGenLists(1);
	glNewList(C2A3, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glRotatef(90, -1, 0, 0);
	glRotatef(Ac2a, 0, 1, 0);	// old GL-z is now GL-y
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, greenDiffuseMaterial);
	gluCylinder(qobj, ra, ra, 2*dc2ahalf, 30, 2);	// gluCylinder(qobj, Rbase, Rtop, Height, Slice, Stack) along z!
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glRotatef(-Ac2a, 0, 1, 0);
	glRotatef(-90, -1, 0, 0);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2A-3 close: cylinder
	x=(float)(RCc2a2[0]);	// GL(x,y,z)=Lab(x,z,-y)
	y=(float)(RCc2a2[1]);
	z=(float)(RCc2a2[2]);
	
	C2A3close=glGenLists(1);
	glNewList(C2A3close, GL_COMPILE);
	glTranslatef(x, z, -y);	// GL(x,y,z)=Lab(x,z,-y)
	glRotatef(90, -1, 0, 0);
	glRotatef(Ac2a, 0, 1, 0);	// old GL-z is now GL-y
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	//glEnable(GL_LIGHT1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, redDiffuseMaterial);
	gluCylinder(qobj, ra, ra, 2*dc2ahalf, 30, 2);	// gluCylinder(qobj, Rbase, Rtop, Height, Slice, Stack) along z!
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHTING);
	glRotatef(-Ac2a, 0, 1, 0);
	glRotatef(-90, -1, 0, 0);
	glTranslatef(-x, -z, y);
	glEndList();
	
	// C2AB open
	C2AB=glGenLists(1);	// molecules in ring: red
	glNewList(C2AB, GL_COMPILE);
	glCallList(LYS);
	
	glCallList(C2B1);
	glCallList(C2B2);
	glCallList(C2B3);
	
	glCallList(C2A1);
	glCallList(C2A2);
	glCallList(C2A3);
	glEndList();
	
	
	// C2AB close
	C2ABclose=glGenLists(1);	// molecules in ring: red
	glNewList(C2ABclose, GL_COMPILE);
	glCallList(LYS);
	
	glCallList(C2B1close);
	glCallList(C2B2close);
	glCallList(C2B3close);
	
	glCallList(C2A1close);
	glCallList(C2A2close);
	glCallList(C2A3close);
	glEndList();
}





void display_syt(void)
{
	//display_syt_cube();
	//display_syt_sphere();
	display_syt_cylinders();
}



void display_vesicle(void)
{
	vesicle=glGenLists(1);	// vesicle
	glNewList(vesicle, GL_COMPILE);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glColor4f(0.4, 0.4, 1, 0.8);
	//glColor4f(1, 1, 1, 0.8);
	glutSolidSphere(Rv, 60, 40);
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glEndList();
}


void display_tmd(void)
{
	tmd=glGenLists(1);	// tmd
	glNewList(tmd, GL_COMPILE);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, yellowDiffuseMaterial);
	glutSolidSphere(Rtmd, 10, 10);
	/*
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glColor4f(1, 1, 0, 0.8);
	//glColor4f(1, 1, 0, 0.8);
	glutSolidSphere(Rtmd, 10, 10);
	glDisable(GL_BLEND);
	*/
	glDisable(GL_LIGHTING);
	glEndList();
}



void display_list(void)
{
	display_syt();
	//display_vesicle();
	display_tmd();
}


