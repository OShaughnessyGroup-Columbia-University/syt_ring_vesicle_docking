
void init_GL(void)
{
	float x=-4*Lx, y=8*Lx, z=8*Lx;
	
	GLfloat light_ambient[]={0.3, 0.3, 0.3, 0.5};
	GLfloat light_diffuse[]={0.5, 0.5, 0.5, 0.5};			// 4th value is transparency
	GLfloat light_specular[] = {0, 0, 0, 0};		// low alpha means transparent
	GLfloat light_shininess[] = {100};
	GLfloat light_position[] = {x, y, z, 1};
	
	GLfloat light1_ambient[]={0, 1, 0, 0.5};	// green
	GLfloat light1_diffuse[]={0, 1, 0, 0.5};			// 4th value is transparency
	GLfloat light1_specular[] = {0, 0, 0, 1};		// low alpha means transparent
	GLfloat light1_shininess[] = {100};
	GLfloat light1_position[] = {x, y, z, 1};
	
	GLfloat light2_ambient[]={1, 0, 0, 0.5};	// red
	GLfloat light2_diffuse[]={1, 0, 0, 0.5};			// 4th value is transparency
	GLfloat light2_specular[] = {0, 0, 0, 1};		// low alpha means transparent
	GLfloat light2_shininess[] = {100};
	GLfloat light2_position[] = {x, y, z, 1};
	
	GLfloat light3_ambient[]={0.5, 0.5, 1, 0.5};	// blue
	GLfloat light3_diffuse[]={0.5, 0.5, 1, 0.5};			// 4th value is transparency
	GLfloat light3_specular[] = {0, 0, 0, 1};		// low alpha means transparent
	GLfloat light3_shininess[] = {100};
	GLfloat light3_position[] = {x, y, z, 1};
	
	GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat mat_ambient[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat mat_ambient_color[] = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat mat_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat no_shininess[] = { 0.0 };
	GLfloat low_shininess[] = { 5.0 };
	GLfloat high_shininess[] = { 100.0 };
	GLfloat mat_emission[] = {0.3, 0.2, 0.2, 0.0};
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);	// default
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);	// light 1
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	
	glLightfv(GL_LIGHT2, GL_AMBIENT, light2_ambient);	// light 2
	glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
	glLightfv(GL_LIGHT2, GL_SPECULAR, light2_specular);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	
	glLightfv(GL_LIGHT3, GL_AMBIENT, light3_ambient);	// light 3
	glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_diffuse);
	glLightfv(GL_LIGHT3, GL_SPECULAR, light3_specular);
	glLightfv(GL_LIGHT3, GL_POSITION, light3_position);
	
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, light_shininess);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, light_diffuse);
	
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	

#if AALIAS
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glHint(GL_POINT_SMOOTH_HINT, GL_DONT_CARE);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
#endif

//	glClearColor(1.0, 1.0, 1.0, 1.0);	// white background
	glClearColor(0.0, 0.0, 0.0, 0.0);	// black background
//	glClearColor(0.4, 0.4, 0.4, 1.0);
	
	glShadeModel(GL_SMOOTH);

	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	glEnable(GL_DEPTH_TEST);
	
}


void changeSize(int w, int h)
{
	float ratio, dist;
	if(h == 0) h = 1;
	ratio = 1.0* w / h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(20,ratio,max2(1,0.1*Ly),10*Ly);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	dist=(float)(2.0*(Lx+Ly+0.8*Rv+Hv));
	gluLookAt(0, 0, dist, 
						0, 0, 0,
						0, 1, 0);
}


void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
		case 27:		// When Escape Is Pressed...
			terminate();	// Exit The Program
			break;		// Ready For Next Case
		case '+':
		case '=':
			zoom*=1.1;
			break;
			break;
		case '-':
			zoom/=1.1;
			break;
		case ' ':	// reset
			initview();
			break;
		case 'a':
		case 'A':
			attracton=(attracton+1)%2;
			break;
		case 'b':
		case 'B':
			showboxsphere=(showboxsphere+1)%2;
			break;
		case 'c':
		case 'C':
			showcolor=(showcolor+1)%3;
			break;
		case 'd':
		case 'D':
			showden=(showden+1)%2;
			break;
		case 'f':
		case 'F':
			if(Vesicle) showforce=(showforce+1)%4;
			else showforce=(showforce+1)%3;
			break;
		case 'k':
		case 'K':
			if(Vesicle) showmark=(showmark+1)%4;
			else showmark=(showmark+1)%3;
			break;
		case 'm':
		case 'M':
			showmem=(showmem+1)%2;
			break;
		case 'n':
		case 'N':
			noise=(noise+1)%2;
			break;
		case 'p':
		case 'P':
			pause=(pause+1)%2;
			break;
		case 'r':
		case 'R':
			readpara();
			showpara();
			getconst();
			reset();
			break;
		case 's':
		case 'S':
			showsyt=(showsyt+1)%2;
			break;
		case 't':
		case 'T':
			showtmd=(showtmd+1)%2;
			break;
		case 'x':
		case 'X':
			showcoord=(showcoord+1)%2;
			break;
		case 'v':
		case 'V':
			showvesicle=(showvesicle+1)%2;
			break;
		case 'z':
		case 'Z':
			if(Vesicle) shownormal=(shownormal+1)%4;
			else shownormal=(shownormal+1)%3;
			break;
		default:		// Now Wrap It Up
			break;
	}
}


void arrow_keys(int a_keys, int x, int y)  // Create Special Function (required for arrow keys)
{
	switch (a_keys) {
	case GLUT_KEY_UP:	// When Up Arrow Is Pressed...
		PanY+=0.05*Lx/zoom;
		break;
	case GLUT_KEY_DOWN:	// When Down Arrow Is Pressed...
		PanY-=0.05*Lx/zoom;
		break;
	case GLUT_KEY_RIGHT:
		PanX+=0.05*Lx/zoom;
		break;
	case GLUT_KEY_LEFT:
		PanX-=0.05*Lx/zoom;
		break;
	default:
		break;
  }
}


void mouse(int button, int state, int x, int y)
{
	if(button==GLUT_LEFT_BUTTON) {
		if(state==GLUT_DOWN) {
			mouse_left=1;
			mouse_x=x;
			mouse_y=y;
		}
		else mouse_left=0;
	}

	if(button==GLUT_RIGHT_BUTTON) {
		if(state==GLUT_DOWN) {
			mouse_right=1;
			mouse_x=x;
			mouse_y=y;
		}
		else mouse_right=0;
	}
}



void mousewheel(int button, int direction, int x, int y)
{
	if(direction>0) zoom*=1.1;	// zoom in
	else zoom/=1.1;	// zoom out
}


void drag(int x, int y) 
{ 
	if(mouse_left==1) {
		rotphi+=0.5*(x-mouse_x);
		rotphi+=0.5*(y-mouse_y);
		mouse_x=x;
		mouse_y=y; 
		glutPostRedisplay(); 
	}

	if(mouse_right==1) {
		rottheta+=0.5*(x-mouse_x);
		rottheta+=0.5*(y-mouse_y);
		mouse_x=x;
		mouse_y=y; 
		glutPostRedisplay(); 
	} 
}


