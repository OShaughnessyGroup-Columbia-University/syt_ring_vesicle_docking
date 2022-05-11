
//------- NR utils ---------------

#define NR_END 1
#define FREE_ARG char*
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


// -------------------------

#if INLINE
inline 
#endif
double norm(double r[3])
{	return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);	}


#if INLINE
inline 
#endif
double norm2(double r[3])
{	return r[0]*r[0]+r[1]*r[1]+r[2]*r[2];	}


#if INLINE
inline 
#endif
void normalize(double r[], double rn[])
{
	double z;
	z=norm(r);
	if(z>1e-10) { rn[0]=r[0]/z; rn[1]=r[1]/z; rn[2]=r[2]/z; }
	else { rn[0]=rn[1]=0.0; rn[2]=1.0; }	// !!!!!!!!!!!!! correct r to rn !!!!!!!!!!!!!
}


#if INLINE
inline 
#endif
double dotprod(double r1[], double r2[])
{
	return r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2];
}


#if INLINE
inline 
#endif
void crossprod(double r1[], double r2[], double r3[])
{
	r3[0]=r1[1]*r2[2]-r1[2]*r2[1];
	r3[1]=r1[2]*r2[0]-r1[0]*r2[2];
	r3[2]=r1[0]*r2[1]-r1[1]*r2[0];
}

#if INLINE
inline 
#endif
void getnromalvector1(double tng[3], double nrm[3])	// nrm = nz x tng, nrm is in the x-y plane
{
	nrm[0] = -tng[1];
	nrm[1] = tng[0];
	nrm[2] = 0;
}

#if INLINE
inline 
#endif
void getnromalvector2(double tng[3], double nrm[3])	// nrm is in the z-direction
{
	if(tng[0]!=0 && tng[2]!=0) {	// nrm = tng x ny
		nrm[0] = -tng[2];
		nrm[1] = 0;
		nrm[2] = tng[0];
	}
	else {	// nrm = nz
		nrm[0] = nrm[1] = 0;
		nrm[2] = 1;
	}
}

#if INLINE
inline 
#endif
void veczero(double r[])
{
	r[0]=r[1]=r[2]=0.0;
}

#if INLINE
inline 
#endif
void veccopy(double r1[], double r2[])	// r[3]
{
	r2[0]=r1[0];
	r2[1]=r1[1];
	r2[2]=r1[2];
}


#if INLINE
inline 
#endif
void vecadd(double r1[], double r2[], double r3[])	// r[3]
{
	r3[0]=r1[0]+r2[0];
	r3[1]=r1[1]+r2[1];
	r3[2]=r1[2]+r2[2];
}


#if INLINE
inline 
#endif
void vecadd3(double r1[], double r2[], double r3[], double r4[])	// r[3]
{
	r4[0]=r1[0]+r2[0]+r3[0];
	r4[1]=r1[1]+r2[1]+r3[1];
	r4[2]=r1[2]+r2[2]+r3[2]; 
}


#if INLINE
inline 
#endif
void vecadd4(double r1[], double r2[], double r3[], double r4[], double r5[])	// r[3]
{
	r5[0]=r1[0]+r2[0]+r3[0]+r4[0];
	r5[1]=r1[1]+r2[1]+r3[1]+r4[1];
	r5[2]=r1[2]+r2[2]+r3[2]+r4[2];
}


#if INLINE
inline 
#endif
void vecadd5(double r1[], double r2[], double r3[], double r4[], double r5[], double r6[])	// r[3]
{
	r6[0]=r1[0]+r2[0]+r3[0]+r4[0]+r5[0];
	r6[1]=r1[1]+r2[1]+r3[1]+r4[1]+r5[1];
	r6[2]=r1[2]+r2[2]+r3[2]+r4[2]+r5[2];
}


#if INLINE
inline 
#endif
void vecadd6(double r1[], double r2[], double r3[], double r4[], double r5[], double r6[], double r7[])	// r[3]
{
	r7[0]=r1[0]+r2[0]+r3[0]+r4[0]+r5[0]+r6[0];
	r7[1]=r1[1]+r2[1]+r3[1]+r4[1]+r5[1]+r6[1];
	r7[2]=r1[2]+r2[2]+r3[2]+r4[2]+r5[2]+r6[2];
}


#if INLINE
inline 
#endif
void vecsub(double r1[], double r2[], double r3[])	// r[3]
{
	r3[0]=r1[0]-r2[0];
	r3[1]=r1[1]-r2[1];
	r3[2]=r1[2]-r2[2];
}


#if INLINE
inline 
#endif
void vecprod(double r1[], double r, double r2[])	// r[3]
{
	r2[0]=r*r1[0];
	r2[1]=r*r1[1];
	r2[2]=r*r1[2];
}


#if INLINE
inline 
#endif
void vecdiv(double r1[], double r, double r2[])	// r[3]
{
	if(r==0.0) {printf("vecdiv: divided by zero\n"); exit(0);}
	else {
		r2[0]=r1[0]/r;
		r2[1]=r1[1]/r;
		r2[2]=r1[2]/r;
	}
}


#if INLINE
inline 
#endif
int vecdivx(double r1[], double r, double r2[])	// r[3]
{
	if(r==0.0) {printf("vecdiv: divided by zero\n"); return 0;}
	else {
		r2[0]=r1[0]/r;
		r2[1]=r1[1]/r;
		r2[2]=r1[2]/r;
		return 1;
	}
}


#if INLINE
inline 
#endif
void vecave2(double r1[], double r2[], double r3[])	// r[3]
{
	r3[0]=(r1[0]+r2[0])/2.0;
	r3[1]=(r1[1]+r2[1])/2.0;
	r3[2]=(r1[2]+r2[2])/2.0;
}


#if INLINE
inline 
#endif
void vecave3(double r1[], double r2[], double r3[], double r4[])	// r[3]
{
	r4[0]=(r1[0]+r2[0]+r3[0])/3.0;
	r4[1]=(r1[1]+r2[1]+r3[1])/3.0;
	r4[2]=(r1[2]+r2[2]+r3[2])/3.0;
}


#if INLINE
inline 
#endif
void vecave4(double r1[], double r2[], double r3[], double r4[], double r5[])	// r[3]
{
	r5[0]=(r1[0]+r2[0]+r3[0]+r4[0])/4.0;
	r5[1]=(r1[1]+r2[1]+r3[1]+r4[1])/4.0;
	r5[2]=(r1[2]+r2[2]+r3[2]+r4[2])/4.0;
}


int vectorequal(double r1[], double r2[], double torr)
{
	double dr[3];
	vecsub(r1, r2, dr);
	if(norm(dr)<torr) return 1;
	else return 0;
}


void limitvector(double r[], double s)	// make |r[3]|<=s
{
	int i;
	double rr;
	rr=norm(r);
	if(rr>s) {
		rr=s/rr;
		for(i=0; i<3; i++) r[i]*=rr;
	}
}

#if INLINE
inline 
#endif
double randsign(void)
{
	return ran3(&mseed)>0.5?1.0:-1.0;
}


#if INLINE
inline 
#endif
double distance2(double r1[], double r2[])
{
	return pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2)+pow(r1[2]-r2[2],2);
}


#if INLINE
inline 
#endif
double distance(double r1[], double r2[])
{
	return sqrt(distance2(r1,r2));
}


#if INLINE
inline 
#endif
double avoidzero(double x)
{
	if(fabs(x)<1.0e-9) return sign(x)*1.0e-9;
	else return x;
}


double cotangent(double r1[3], double r2[3])	// returns cot of angle between vectors r1 & r2
{	// cot=cos/sin, cos=(a.b)/(|a|.|b|)
	double dotp, denom;
	
	dotp=dotprod(r1,r2);
	denom=norm2(r1)*norm2(r2)-pow(dotp,2);
	denom=sqrt(max2(denom,eps));
	return dotp/denom;
}


#if INLINE
inline 
#endif
double srchLambertW(double w, double val)	// val = w*exp(w)
{
	return w*exp(w) - val;
}


double LambertW(double x)	// Lambert W function, w-root of x=w*exp(w)
{
	double L1, L2, L12, L13, L14, L15, L22, L23, L24, w;
	
	if(x==0) w=0;
	else if(x<1) w=x-pow(x,2)+1.5*pow(x,3)-8/3.0*pow(x,4)+125/24.0*pow(x,5)-54/5.0*pow(x,6)+16807/720.0*pow(x,7);
	else if(x>=3) {
		L1=log(x);
		L2=log(L1);
		L12=L1*L1;	L13=L12*L1;	L14=L13*L1;	L15=L14*L1;
		L22=L2*L2;	L23=L22*L2;	L24=L23*L2;
		w = L1 - L2 + L2/L1 + L2*(L2-2)/2.0/L12 + L2*(2*L22-9*L2+6)/6.0/L13 + L2*(3*L23-22*L22+36*L2-12)/12.0/L14 + 
			L2*(12*L24-125*L23+350*L22-300*L2+60)/60.0/L15;
	}
	else {
		w=zbrent2(srchLambertW, x, 0.0, 10.0, eps);
	}
	
	return w;
}



void RotMatrix2(double sth, double cth, double r1[2], double r2[2])	// rotate 2D vector ccw by th
{
	double rr[2];
	rr[0]=cth*r1[0]-sth*r1[1];
	rr[1]=sth*r1[0]+cth*r1[1];
	r2[0]=rr[0];
	r2[1]=rr[1];
}


void RotMatrix3z(double sth, double cth, double r1[3], double r2[3])	// rotate 3D vector ccw about z by th
{
	double rr[3];
	rr[0]=cth*r1[0]-sth*r1[1];
	rr[1]=sth*r1[0]+cth*r1[1];
	rr[2]=r1[2];
	veccopy(rr,r2);
}


// rotate r1 ccw by th about unit vector (axis) to r2
void RotMatrix3(double sth, double cth, double axis[3], double r1[3], double r2[3])
{	// Rodrigues' rotation formula
	int i;
	double crossp[3], dotp, dotp1mcs, rr[3];
	
	crossprod(axis, r1, crossp);
	dotp=dotprod(axis, r1);
	dotp1mcs=dotp*(1-cth);
	for(i=0; i<3; i++) rr[i] = r1[i]*cth + crossp[i]*sth + axis[i]*dotp1mcs;
	veccopy(rr, r2);
}


#if INLINE
inline 
#endif
int getCurvatureSign(double dz1, double dz2)	// convex is 1, concave is -1
{
	if(dz1>=dz2) return 1;
	else return -1;
}


double getanglez(double dz1, double dz2, double dx)	// 2D angle: convex is positive, concave is negative
{
	int i;
	double r1[2], r2[2], d1, d2, cs, a;
	
	r1[0]=r2[0]=dx;
	r1[1]=dz1;
	r2[1]=dz2;
	d1=sqrt(pow(r1[0],2)+pow(r1[1],2));
	d2=sqrt(pow(r2[0],2)+pow(r2[1],2));
	for(i=0; i<2; i++) {
		r1[i]/=d1;
		r2[i]/=d2;
	}
	cs=dotprod(r1,r2);
	cs=min2(-1, max2(1, cs));
	a=acos(cs);
	a*=getCurvatureSign(dz1, dz2);	// convex is positive, concave is negative
	return a;
}


double getanglePi(double r1[3], double r2[3])	// 0<=angle<pi
{
	double dotp, denom;
	denom=norm(r1)*norm(r2);
	if(denom==0) return 0;
	
	dotp=dotprod(r1,r2);
	if(dotp>=denom) return 0;
	else if(dotp<=-denom) return Pi;
	else return acos(dotp/denom);
}


double getanglepmPi2D(double x, double y)	// get ccw angle starting from 3 O'clock (-pi, pi]
{
	double r;

	r=sqrt(x*x+y*y);
	if(r==0) return 0;
	if(y>=0) return acos(x/r);
	else return -acos(x/r);
}


double getarea2(double r1[3], double r2[3])	// area of a triangle 0-r1-r2
{
	double rc[3];
	crossprod(r1, r2, rc);
	return 0.5*norm(rc);
}


double getarea3(double r1[3], double r2[3], double r3[3])	// area of a triangle r1-r2-r3
{
	double ra[3], rb[3], rc[3];
	
	vecsub(r1, r3, ra);
	vecsub(r2, r3, rb);
	crossprod(ra, rb, rc);
	return 0.5*norm(rc);
}



void getarea3normal(double r1[3], double r2[3], double r3[3], double *a, double nrm[3])	// area of a triangle r1-r2-r3
{
	double ra[3], rb[3], rc[3], d;
	
	vecsub(r1, r3, ra);
	vecsub(r2, r3, rb);
	crossprod(ra, rb, rc);
	d=norm(rc);
	*a=0.5*d;	// area
	vecdiv(rc, d, nrm);	// normal
}



double getarea3db(double r1[3], double r2[3], double r3[3])	// 2x area of a triangle r1-r2-r3
{
	double ra[3], rb[3], rc[3];
	
	vecsub(r1, r3, ra);
	vecsub(r2, r3, rb);
	crossprod(ra, rb, rc);
	return norm(rc);
}


double getarea4(double z1, double z2, double z3, double z4)	// r1-r4 is CW from lower left corner
{	// area of trapizoid is averaged over two folding configurations
	double r12[3], r23[3], r43[3], r14[3];
	double s1[3], s2[3], c[3], a;
	r12[0]=0;	r12[1]=dY;	r12[2]=z2-z1;
	r23[0]=dX;	r23[1]=0;	r23[2]=z3-z2;
	r14[0]=dX;	r14[1]=0;	r14[2]=z4-z1;
	r43[0]=0;	r43[1]=dY;	r43[2]=z4-z3;
	vecadd(r12,r43,s1);
	vecadd(r14,r23,s2);
	crossprod(s1,s2,c);
	a=norm(c)/4.0;
	return a;
}



#if INLINE
inline 
#endif
void sort3(double a, double b, double c, double d[3])	// sort (a,b,c) in descending order into d
{
	if(a>=b) {
		if(c>a) { d[0] = c; d[1] = a; d[2] = b; }	// c>a>=b
		else {
			d[0] = a;
			if(b>=c) { d[1] = b; d[2] = c; }	// a>=b>=c
			else { d[1] = c; d[2] = b; }	// a>=c>=b
		}
	}
	else {	// b>a
		if(a>=c) { d[0] = b; d[1] = a; d[2] = c; }	// b>a>=c
		else {
			d[2] = a;
			if(c>=b) { d[0] = c; d[1] = b; }	// c>=b>a
			else { d[0] = b; d[1] = c; }	// b>c>a
		}
	}
}



void mapsphere(double r0, double r[3])	// sphere centered at origin
{
	int i;
	double rn;
	rn=norm(r)/r0;
	if(rn>0) {
		for(i=0; i<3; i++) r[i]/=rn;
	}
	else {
		r[0]=r[1]=0;
		r[2]=r0;
	}
}


void getRndDir(double axis[3])	// don't put this in OPENMP!
{
	int i;
	double r[3], r2;
	
	do {
		for(i=0; i<3; i++) r[i]=1-2*RND();
		r2=norm2(r);
	} while(r2>1 || r2<=0);	// within a unit sphere
	vecdiv(r, sqrt(r2), axis);	// random axis
}


void getRndDirMP(int tid, double axis[3])	// for OPENMP
{
	int i;
	double r[3], r2;
	
	do {
		for(i=0; i<3; i++) r[i]=1-2*RNDMP(tid);
		r2=norm2(r);
	} while(r2>1 || r2==0);	// within a unit sphere
	vecdiv(r, sqrt(r2), axis);	// random axis
}


double interpolate_func(double *func, double x, double xmin, double dx, int n)
{
	int i;
	double x1, frac, y;
	
	i = (int)floor((x-xmin)/dx);
	if(i<=1) y = func[0];
	else if(i >= n-1) y = func[n-1];
	else {
		x1 = xmin + i*dx;
		frac = (x-x1)/dx;
		y = func[i] + (func[i+1]-func[i])*frac;
	}
	return y;
}


//----------------------------------------------------

double **NewArray2d(int nx, int ny)	// use: array = NewArray2d(nx, ny)
{
	int i;
	double **a;
	
	a=(double **)calloc(nx, sizeof(double *));
	for(i=0; i<nx; i++) {
		a[i]=(double *)calloc(ny, sizeof(double));
		if(!a[i]) {
			printf("Array2i: cannot allocate memory\n");
			exit(0);
		}
	}
	return a;
}

/*
void InitArray(double ***a, int nx, int ny)	// use: InitArray(&array, nx, ny)
{
	int i;
	
	*a=(double **)calloc(nx, sizeof(double *));	// either calloc or malloc
	//*a=malloc(nx*sizeof(double *));
	for(i=0; i<nx; i++) {
		(*a)[i]=(double *)calloc(ny, sizeof(double));	// either calloc or malloc
		//(*a)[i]=malloc(ny*sizeof(double));
		if(!(*a)[i]) {
			printf("Cannot allocate memory\n");
			exit(0);
		}
	}
}
*/

void CleanArray2d(double ***a, int nx)	// use: CleanArray2d(&array, nx)
{
	int i;
	
	if(*a) {
		for(i=0; i<nx; i++) free((*a)[i]);
		free(*a);
	}
}


void MatrixProd4(double m[4][4], double x[4], double y[4])	// y=m*x
{
	int i, j;
	
	for(i=0; i<4; i++) {
		y[i]=0;
		for(j=0; j<4; j++) y[i]+=m[i][j]*x[j];
	}
}


//----------------------------------------------------

// see migration 3d - v19

#if INLINE
inline 
#endif
void vecsub2D(double r1[], double r2[], double r3[])	// r[2]
{
	r3[0]=r1[0]-r2[0];	r3[1]=r1[1]-r2[1];
}

double crossprod2D(double r1[2], double r2[2])
{
	return r1[0]*r2[1]-r2[0]*r1[1];
}


double triangleareasigndb2D(double r1[2], double r2[2], double r3[2])
{	// calc. area of triangle r1-r2-r3
	double r12[2], r13[2], z;
	
	vecsub2D(r2, r1, r12);
	vecsub2D(r3, r1, r13);
	z=crossprod2D(r12, r13);
	return z;
}


// segment-segment intersection from "Computational Geometry in C" by O'Rourke
#if INLINE
inline 
#endif
int SegSegLeft2D(double r1[2], double r2[2], double r3[2])
{
	//if(triangleareasigndb2D(r1, r2, r3) > 0) return 1;	// r3 is on the left of r1->r2
	if(triangleareasigndb2D(r1, r2, r3) >= 0) return 1;	// r3 is on the left or on r1->r2
	else return 0;
}


#if INLINE
inline 
#endif
int SegSegCollinear2D(double r1[2], double r2[2], double r3[2])
{
	if(triangleareasigndb2D(r1, r2, r3) == 0) return 1;
	else return 0;
}


#if INLINE
inline 
#endif
int Xor(int x, int y)	// exclusive or: =1 iff exactly one argument is true
{
	return !x ^ !y;
}


#if INLINE
inline 
#endif
int SegIntersect2D(double r1[2], double r2[2], double r3[2], double r4[2])
{	// return 1 if r1-r2 & r3-r4 intersect, 0 if not
	
/*	full version
	if(SegSegCollinear2D(r1,r2,r3) || SegSegCollinear2D(r1,r2,r4) || 
		SegSegCollinear2D(r1,r3,r4) || SegSegCollinear2D(r2,r3,r4)) return 0;
	else return Xor(SegSegLeft2D(r1,r2,r3), SegSegLeft2D(r1,r2,r4)) && 
		Xor(SegSegLeft2D(r1,r3,r4), SegSegLeft2D(r2,r3,r4));
*/
///* simplified version
	return Xor(SegSegLeft2D(r1,r2,r3), SegSegLeft2D(r1,r2,r4)) && 
		Xor(SegSegLeft2D(r1,r3,r4), SegSegLeft2D(r2,r3,r4));
//*/
}


//----------------------------------------------------

// find point rp where line 0-r0 intersects with plane containing triangle r1-r2-r3, see "Computational Geometry in C" pg. 226
void trisegintersect3D(double r0[3], double r1[3], double r2[3], double r3[3], double nrm[3], double rp[3])	// all in RC frame!
{
	double rx[3], num, denom, r;
	
	vecprod(r0, 10.0, rx);	// extend r0 to rx
	num=dotprod(r1, nrm);
	denom=dotprod(rx, nrm);
	r=fabs(num/avoidzero(denom));	// 0<r<1
	vecprod(rx, r, rp);	// intersection point in RC frame
}


// project 3D tirangle r1-r2-r3 and point rp to one of the x-y, y-z, and z-x planes. 2D coordinates are in v's
void projection32(double r1[3], double r2[3], double r3[3], double nrm[3], double rp[3], 
			double v1[2], double v2[2], double v3[2], double vp[2])	// all in RC frame!
{
	int i, j, plane;
	double nx, ny, nz;
	
	nx=fabs(nrm[0]);
	ny=fabs(nrm[1]);
	nz=fabs(nrm[2]);
	
	if(nz>=max2(nx,ny)) plane=0;	// project to x-y plane
	else if(nx>=max2(ny,nz)) plane=1;	// project to y-z plane
	else plane=2;	// project to z-x plane
	
	for(i=0; i<2; i++) {
		j=(i+plane)%3;
		v1[i]=r1[j];
		v2[i]=r2[j];
		v3[i]=r3[j];
		vp[i]=rp[j];
	}
}


int pointinside2D(double v1[2], double v2[2], double v3[2], double vp[2])	// check if vp is inside triangle v1-v2-v3
{
	int cnt;
	cnt=0;
	if(SegIntersect2D(vp, Rinf2d, v1, v2)==1) cnt++;
	if(SegIntersect2D(vp, Rinf2d, v2, v3)==1) cnt++;
	if(SegIntersect2D(vp, Rinf2d, v3, v1)==1) cnt++;
	cnt=cnt%2;	// point is inside if cnt is odd, outside if cnt is even
	return cnt;	// 0=outside, 1=inside
}


int checkface3D(double r0[3], double r1[3], double r2[3], double r3[3], double nrm[3], double rp[3])	// all in RC frame
{	// check if r0 is inside triangle r1-r2-r3 (with nrm being the normal direction)
	int z;
	double v1[2], v2[2], v3[2], vp[2];
	
	trisegintersect3D(r0, r1, r2, r3, nrm, rp);	// rp is the line-plane intersection
	projection32(r1, r2, r3, nrm, rp, v1, v2, v3, vp);	// project 3D r's to 2D v's
	z=pointinside2D(v1, v2, v3, vp);
	return z;	// 0=outside, 1=inside
}


//----------------------------------------------------
// don't use the following randx's for RNG !!!!!!!!!!!!!!
/*
double randx1(int *seed)	// NR ran0()?
{
	double r;

	*seed = ( *seed % 65536 );
	*seed = ( ( 3125 * *seed ) % 65536 );
	r = ( double ) ( *seed ) / 65536.0;

	return r;
}


double randx_a(unsigned long *seed)	// use NR numbers for quick & dirty RNG
{
	double r;

	*seed = ( *seed % 1013904223L );
	*seed = ( ( 1664525L * *seed ) % 1013904223L );
	r = ( double ) ( *seed ) / 1013904223.0;

	return r;
}


double randx(unsigned long *seed)	// use NR numbers for quick & dirty RNG (even dirtier)
{
	unsigned long itemp;
	static unsigned long jflone = 0x3f800000;	// don't use on 64-bit machines!
	static unsigned long jflmsk = 0x007fffff;
	double r;

	*seed = 1664525L*(*seed) + 1013904223L;
	itemp = jflone | (jflmsk & (*seed));
	r = (*(float *)&itemp)-1.0;

	return r;
}
*/


//----- RNG: 64-bit ----
#if RNG_64

Ullong rand_int64(void) // Ranq1.int64() in NR(3e)
{
	U64 ^= U64 >> 21;
	U64 ^= U64 << 35;
	U64 ^= U64 >> 4;
	return U64*2685821657736338717LL;
}

void randSeed64(void) // based on Ranq1 constructor in NR(3e)
{
	Ullong seed = (int) (mseed); // random-ish seed
	U64 ^= seed;
	U64 = rand_int64();
}

double rand64(void) // Ranq1.doub() in NR(3e)
{
	return 5.42101086242752217e-20*rand_int64();
}

Ullong rand_int64mp(int i) // Ranq1.int64() in NR(3e)
{
	U64MP[i] ^= U64MP[i] >> 21;
	U64MP[i] ^= U64MP[i] << 35;
	U64MP[i] ^= U64MP[i] >> 4;
	return U64MP[i]*2685821657736338717LL;
}

void randSeed64mp(int i) // based on Ranq1 constructor in NR(3e)
{
	Ullong seed = (int) (mseedarray[i]); // random-ish seed
	U64MP[i] ^= seed;
	U64MP[i] = rand_int64mp(i);
}

double rand64mp(int i) // Ranq1.doub() in NR(3e)
{
	return 5.42101086242752217e-20*rand_int64mp(i);
}

#endif

//----- RNG: 32-bit ----
#if !RNG_64

Uint rand_int32(void)
{
	Uint x, y;
	U32 = U32 * 2891336453U + 1640531513U;
	V32 ^= V32 >> 13;
	V32 ^= V32 << 17;
	V32 ^= V32 >> 5;
	W132 = 33378 * (W132 & 0xffff) + (W132 >> 16);
	W232 = 57225 * (W232 & 0xffff) + (W232 >> 16);
	x = U32 ^ (U32 << 9);
	x ^= x >> 17;
	x ^= x << 6;
	y = W132 ^ (W132 << 17);
	y ^= y >> 15;
	y ^= y << 5;
	return (x + V32) ^ (y + W232);
}

void randSeed32(void)
{
	Uint seed = (int) (mseed); // random-ish seed
	U32 = seed ^ V32;
	rand_int32();
	V32 = U32;
	rand_int32();
}

double rand32(void)
{
	return 2.32830643653869629E-10 * rand_int32();
}

Uint rand_int32mp(int i)
{
	Uint x, y;
	U32MP[i] = U32MP[i] * 2891336453U + 1640531513U;
	V32MP[i] ^= V32MP[i] >> 13;
	V32MP[i] ^= V32MP[i] << 17;
	V32MP[i] ^= V32MP[i] >> 5;
	W132MP[i] = 33378 * (W132MP[i] & 0xffff) + (W132MP[i] >> 16);
	W232MP[i] = 57225 * (W232MP[i] & 0xffff) + (W232MP[i] >> 16);
	x = U32MP[i] ^ (U32MP[i] << 9);
	x ^= x >> 17;
	x ^= x << 6;
	y = W132MP[i] ^ (W132MP[i] << 17);
	y ^= y >> 15;
	y ^= y << 5;
	return (x + V32MP[i]) ^ (y + W232MP[i]);
}

void randSeed32mp(int i)
{
	Uint seed = (int) (mseedarray[i]); // random-ish seed
	//printf("seed-%d = %d\n", i, seed);
	U32MP[i] = seed ^ V32MP[i];
	rand_int32mp(i);
	V32MP[i] = U32MP[i];
	rand_int32mp(i);
}

double rand32mp(int i)
{
	return 2.32830643653869629E-10 * rand_int32mp(i);
}

#endif

//----- RNG: ran2 -----

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


//----- RNG: ran3 -----

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


//------ RNG: MT19937p ------

/* A C-program for MT19937: Real number version                */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long seed, mt19937p *config)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
	config->mti = N+1;
	config->mag01[0] = 0x0;
	config->mag01[1] = MATRIX_A;
    config->mt[0]= seed & 0xffffffff;
    for (config->mti=1; config->mti<N; config->mti++)
        config->mt[config->mti] = (69069 * config->mt[config->mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double genrand(mt19937p *config)
{
    unsigned long y;
//    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (config->mti >= N) { /* generate N words at one time */
        int kk;

/*        if (config->mti == N+1)*/   /* if sgenrand() has not been called, */
/*            sgenrand(4357);*/ /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (config->mt[kk]&UPPER_MASK)|(config->mt[kk+1]&LOWER_MASK);
            config->mt[kk] = config->mt[kk+M] ^ (y >> 1) ^ config->mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (config->mt[kk]&UPPER_MASK)|(config->mt[kk+1]&LOWER_MASK);
            config->mt[kk] = config->mt[kk+(M-N)] ^ (y >> 1) ^ config->mag01[y & 0x1];
        }
        y = (config->mt[N-1]&UPPER_MASK)|(config->mt[0]&LOWER_MASK);
        config->mt[N-1] = config->mt[M-1] ^ (y >> 1) ^ config->mag01[y & 0x1];

        config->mti = 0;
    }
  
    y = config->mt[config->mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}


//------------ SVD least square fitting ------------------

void svbksb(double u[5][4], double w[], double v[4][4], int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s,*tmp;

	//tmp=dvector(1,n);
	tmp=malloc((n+1)*sizeof(double));
	for (j=1;j<=n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	//free_dvector(tmp,1,n);
	free(tmp);
}



double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}



//void svdcmp(double **a, int m, int n, double w[], double **v)
void svdcmp(double a[5][4], int m, int n, double w[], double v[4][4])
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	//rv1=dvector(1,n);
	rv1=malloc((n+1)*sizeof(double));
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=max2(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=min2(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				//if(nm==0) printf("nm=%d\n", nm);
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	//free_dvector(rv1,1,n);
	free(rv1);
}


// modified svdfit

//#define TOL 1.0e-5	// default value
#define TOL 1.0e-4

void svdfit(int x[], double y[], double r[8][3], double sig[], int ndata, double a[], int ma,
//	double **u, double **v, double w[], double *chisq, void (*funcs)(int, double [8][3], double []))
	double u[5][4], double v[4][4], double w[], double *chisq, void (*funcs)(int, double [8][3], double []))
{
	//void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
	//void svdcmp(double **a, int m, int n, double w[], double **v);
	void svbksb(double u[5][4], double w[], double v[4][4], int m, int n, double b[], double x[]);
	void svdcmp(double a[5][4], int m, int n, double w[], double v[4][4]);
	int j,i;
	double wmax,tmp,thresh,sum,*b,*afunc;

	//b=dvector(1,ndata);
	//afunc=dvector(1,ma);
	b=malloc((ndata+1)*sizeof(double));
	afunc=malloc((ma+1)*sizeof(double));
	
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],r,afunc);
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
		b[i]=y[i]*tmp;
	}
	svdcmp(u,ndata,ma,w,v);
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=TOL*wmax;
	for (j=1;j<=ma;j++)
		if (w[j] < thresh) w[j]=0.0;
	svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],r,afunc);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
	}
	//free_dvector(afunc,1,ma);
	//free_dvector(b,1,ndata);
	free(afunc);
	free(b);
}
#undef TOL



//------------ exponential integral function En(x) --------------

#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS 1.0e-7

double expint(int n, double x)
{
	void nrerror(char error_text[]);
	int i,ii,nm1;
	double a,b,c,d,del,fact,h,psi,ans;

	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
	nrerror("bad arguments in expint");
	else {
		if (n == 0) ans=exp(-x)/x;
		else {
			if (x == 0.0) ans=1.0/nm1;

			else {
				if (x > 1.0) {
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans=h*exp(-x);
							return ans;
						}
					}
					nrerror("continued fraction failed in expint");
				} else {
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					nrerror("series failed in expint");
				}
			}
		}
	}
	return ans;
}
#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER



//---------------- Root Finding: Brent's --------------------

#define NRANSI
#define ITMAX 100
#define EPS 3.0e-8

double zbrent(double (*func)(double), double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double f_a=(*func)(a),f_b=(*func)(b),f_c,p,q,r,s,tol1,xm;

	if ((f_a > 0.0 && f_b > 0.0) || (f_a < 0.0 && f_b < 0.0))
	//	nrerror("Root must be bracketed in zbrent");
		return 0.0;
	f_c=f_b;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((f_b > 0.0 && f_c > 0.0) || (f_b < 0.0 && f_c < 0.0)) {
			c=a;
			f_c=f_a;
			e=d=b-a;
		}
		if (fabs(f_c) < fabs(f_b)) {
			a=b;
			b=c;
			c=a;
			f_a=f_b;
			f_b=f_c;
			f_c=f_a;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || f_b == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(f_a) > fabs(f_b)) {
			s=f_b/f_a;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=f_a/f_c;
				r=f_b/f_c;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		f_a=f_b;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		f_b=(*func)(b);
	}
	//nrerror("Maximum number of iterations exceeded in zbrent");
	return b;
}
#undef ITMAX
#undef EPS
#undef NRANSI

//--------------- Root Finding: Brent's w/ 2 indices --------------

#define NRANSI
#define ITMAX 100
#define EPS 3.0e-8

double zbrent2(double (*func)(double, double), double val, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double f_a=(*func)(a,val),f_b=(*func)(b,val),f_c,p,q,r,s,tol1,xm;

	if ((f_a > 0.0 && f_b > 0.0) || (f_a < 0.0 && f_b < 0.0))
	//	nrerror("Root must be bracketed in zbrent");
		return 0.0;
	f_c=f_b;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((f_b > 0.0 && f_c > 0.0) || (f_b < 0.0 && f_c < 0.0)) {
			c=a;
			f_c=f_a;
			e=d=b-a;
		}
		if (fabs(f_c) < fabs(f_b)) {
			a=b;
			b=c;
			c=a;
			f_a=f_b;
			f_b=f_c;
			f_c=f_a;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || f_b == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(f_a) > fabs(f_b)) {
			s=f_b/f_a;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=f_a/f_c;
				r=f_b/f_c;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		f_a=f_b;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		f_b=(*func)(b,val);
	}
	//nrerror("Maximum number of iterations exceeded in zbrent");
	return b;
}
#undef ITMAX
#undef EPS
#undef NRANSI



//--------------- Root Finding: Brent's w/ 3 indices --------------

#define NRANSI
#define ITMAX 100
#define EPS 3.0e-8

double zbrent3(double (*func)(double, double, double), double v1, double v2, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double f_a=(*func)(a,v1,v2),f_b=(*func)(b,v1,v2),f_c,p,q,r,s,tol1,xm;

	if ((f_a > 0.0 && f_b > 0.0) || (f_a < 0.0 && f_b < 0.0))
	//	nrerror("Root must be bracketed in zbrent");
		return 0.0;
	f_c=f_b;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((f_b > 0.0 && f_c > 0.0) || (f_b < 0.0 && f_c < 0.0)) {
			c=a;
			f_c=f_a;
			e=d=b-a;
		}
		if (fabs(f_c) < fabs(f_b)) {
			a=b;
			b=c;
			c=a;
			f_a=f_b;
			f_b=f_c;
			f_c=f_a;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || f_b == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(f_a) > fabs(f_b)) {
			s=f_b/f_a;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=f_a/f_c;
				r=f_b/f_c;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		f_a=f_b;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		f_b=(*func)(b,v1,v2);
	}
	//nrerror("Maximum number of iterations exceeded in zbrent");
	return b;
}
#undef ITMAX
#undef EPS
#undef NRANSI












