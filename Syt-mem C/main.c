#define INLINE 0	// use inline function or not (=1 if C99)

#define OPENMP 0	// use openmp or not
#define SAVDAT 1	// save simulation snapshots

#define LGV_MC 1	// use Langevin=1 or MC=0 (note: MC for ves is not finished!)
#define BD_2ND 0	// using 2nd order BD or not

#define OPENGL 0	// show graphics
#define SAVPIC 0	// save pics
#define AALIAS 0	// anti-alias

#define NOISE_ 0	// with thermal noise

// syt
#define FIXSYT 2	// fixing syt: 0=flexible, 1=allow stretching and bending about local z-axis, 2=totally rigid
#define CLOSED 1	// initial oligomer is a closed ring, regardless of # of segments

#define NOSPIN 1	// no spin of chain about z (for FIXSYT=1)
#define NO_ROT 1	// no rotation of syt ring of any kind (including x,y,z) (for FIXSYT=1)
#define FRICTN 0	// friction: 1=from membrane, 0=from cytosol

#define ATTRCT 0	// syt-mem/ves attraction: 0=point-face, 1=point-point
#define GRWGAP 0	// growth requires sufficient gap btwn head & tail

#define SYTCNT 1	// syt center: 1=between C2A&C2B, 0=on C2B
#define VRATES 0	// variable off-rates or not

// mem
#define LAPLAC 1	// laplace formula: 0=Belkin, 1=Meyer, 2=Vallet-Levy (2 not finished for Ves)
#define X_PIP2 0	// extra PIP2 molecule on membrane near each KKKK !!!

#define INIMEM 2	// initial membrane shape: 0=flat, 1=dome, 2=read from contour.ini
#define CRCBND 1	// circular boundary: 0=no, 1=yes

#define FLPBND 1	// flip mem bond or not
#define MC_XYZ 0	// move membrane nodes in all (x,y,z) direction (=1) or just z (=0)

#define CURVAT 1	// method for curvature calculation: 0=paraboloid fitting, 1=Watanabe, 2=Meyer
#define BNDCND 0	// boundary condition: 0=fixed (no pressure), 1=periodic (w/ pressure)
#define SIHOLE 0	// silicon hole of diameter Lx

// ves
#define VSTRCH 1	// consider stretching modulus for vesicle tension
#define XREPUL 0	// extra repulsion distance (5nm) to syt or not

// RNG
#define RNG_64 0	// 64-bit (ranq1 of NR3) or 32-bit (MT19937p faster!) RNG


#include "include.h"


void init(void)
{
	readpara();
	showpara();
	getconst();

	initseedarray();	// run after Nthread is read

	noise=NOISE_;
	RunCnt=1;
	attracton=1;

	// GL
	pause=0;
	showcolor=1;	// 0=no color, 1=height, 2=curvature
	showmem=1;
	showsyt=1;
	showvesicle=1;
	showsyt=1;
	showden=0;
	showforce=0;
	shownormal=0;
	showcoord=0;
	showmark=0;
	showdiagonal=1;
	showtmd=Vesicle*TMD;
	showboxsphere=0;

#if OPENGL && SAVPIC
	img_cnt=0;
	nframe=0;
#endif
	
#if OPENMP
	omp_set_num_threads(Nthread);
#endif

	//testRNG();
}



void initview(void)
{
	zoom=ZOOM;
	if(Vesicle) zoom=(float)(0.7*zoom);
	zoominv=1.0f/zoom;

	PanX=0.0;
	//PanY=-5;
	if(Vesicle) {
		PanY=-0.85f*Rv-max2(5,Hvmax);
		rottheta=3.0;
	}
	else {
		PanY=-4;
		rottheta=20.0;
	}
	rotphi=0;

	//rotphi=180;
	//PanY-=6;

	mouse_x=mouse_y=mouse_left=mouse_right=0;
}


void reset(void)
{
	t=0;
	Nchain=Nclosed=Nopen=0;
	NM=NM0;
	Pon=Pon0;
	Pnuc=Pnuc0;
	
	Nfile=SimCnt1=SimCnt2=0;
	Cntflip_mem=Cntflip_ves=0;

	RCves[0]=RCves[1]=0;
	RCves[2]=Rv+Hv;
	
	dHmemprev=-1;
	dHmem=0;
	
	if(SYT) CleanSyt();
	if(Membrane) CleanMem();
	if(Vesicle) CleanVes();
	
	if(SYT) {
		if(Vesicle && VesRing) InitSyt();
		else {
			if(Limit==0) InitSyt();
			else InitSytMono();
		}
	}
	
	if(Membrane) InitMem();
	if(Vesicle) InitVes();
	
	/*
	if(Resume) {
		if(SYT) CleanSyt();
		ReadData();
	}
	*/
	
	if(SYT) {
		UpdateChain(chain);
		//printf("\nrc=%.3g\n", chain->rcenter[2]);
		//printf("ok\n");
		if(Limit==0 && NC0==1) {
			CenterChain(chain);
			if(Vesicle) RingOnVes(chain);	// update syt sites after centering
		}

		UpdateSytSitesAll(chain);	// update syt-syt sites and syt-mem sites
		UpdateSytStretchAll(chain);
	}

	if(Membrane) UpdateMem(mnode, mface);
	if(Vesicle) UpdateVes(vnode, vface);

#if !LGV_MC
	temperature = 0.1;
	MCcyclecnt=MCloopcnt=0;

	if(SYT) {
		CleanSytBackup();
		InitChainBackup();
	}

	if(Membrane) {
		CleanMemBackup();
		InitMemBackup();
	}
#endif
	
	time(&start_time);
}



void SaveMemH(void)
{
	FILE *fp;
	
	fp=fopen("Hmem.dat", "w");
	fprintf(fp, "%.4g\t%.4g\t%.4g\n", ESM/kBT, Salt, Zmax-Zmin);
	fclose(fp);
}



void SimTime(void)
{
	int days, hours, minutes, seconds;
	double dtime, xtime;
	char tsim[200];
	FILE *fp;

	time(&end_time);
	
	dtime = difftime(end_time, start_time);

	xtime = dtime;
	days = (int)(xtime/(3600*24));
	xtime -= days*(3600*24);
	hours = (int)(xtime/3600);
	xtime -= hours*3600;
	minutes = (int)(xtime/60);
	xtime -= minutes*60;
	seconds = (int)(xtime);

	if(days!=0) sprintf(tsim, "%dd %dh %dm %ds", days, hours, minutes, seconds);
	else if(hours!=0) sprintf(tsim, "%dh %dm %ds", hours, minutes, seconds);
	else if(minutes!=0) sprintf(tsim, "%dm %ds", minutes, seconds);
	else sprintf(tsim, "%ds", seconds);

	printf("Simulation time = %s\n", tsim);

	fp=fopen("SimTime.dat", "w");
	fprintf(fp, "%s\n", tsim);
	fclose(fp);
}


void CleanUp(void)
{
	if(SYT) CleanSyt();
	if(Membrane) CleanMem();
	if(Vesicle) CleanVes();

	if(mseedarray) free(mseedarray);
	if(mtseeds) free(mtseeds);
	
#if RNG_64
	if(U64MP) free(U64MP);
#else
	/*
	if(U32MP) free(U32MP);
	if(V32MP) free(V32MP);
	if(W132MP) free(W132MP);
	if(W232MP) free(W232MP);
	*/
#endif

#if !LGV_MC
	if(SYT) CleanSytBackup();
	if(Membrane) CleanMemBackup();
#endif
}


void terminate(void)
{
#if OPENGL
	//SaveRings();
#endif
	//FlushData();
	//DoStat();	// run after FlushData()
	
	//testfatt();

	if(Membrane) SaveMemH();
	
	SimTime();
	CleanUp();
	
	exit(0);
}



int main(int argc, char **argv)
{
#if OPENGL
	getmseed();
#else
	if(argc==1) getmseed();	// no input
	else mseed=-abs(atoi(argv[1]));	// 1st input: random seed
#endif

	init();
	initview();
	reset();

#if OPENGL
	RunGL(&argc, argv);
#else
	#if LGV_MC
	while(t<=Tmax) DoSimu1();	// Langevin
	#else
	while(MCcyclecnt<MCcycle) DoSimu2();	// Monte-Carlo
	#endif
	terminate();
#endif

	return 0;
}

