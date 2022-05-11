
void skip(FILE *fp)
{
	char ch;
	do {
		ch=fgetc(fp);
	} while(ch!='\n' && ch!=EOF);
}


void readpara(void)
{
	int n;
	char var[50];
	FILE *fp;

	fp=fopen("paras.ini", "r");
	if(!fp) {
		printf("Error: cannot open file 'para.ini'.\n\n");
		//printf("Press Enter to exit...");
		//(void) getchar();
		exit(0);
	}
	
	n = fscanf(fp, "%s %d", var, &Nthread);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Resume);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &SYT);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Membrane);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Vesicle);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &TMD);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &VesRing);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Carbon);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Lx);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Ly);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Lz);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &Nx);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Ny);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &MeshLevel);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Rv);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Hv);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Hvmax);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &KAmem);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &KBmem);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Gauss);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Gamma);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Beta);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &KMspon);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &P0);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Hcell);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &h_belkin_fct);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &Limit);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &C0);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &NC0);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &NS0);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &NS1);
	skip(fp);
	
	
	n = fscanf(fp, "%s %d", var, &Spiral);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &GrwEnd);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Knuc);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Kon_C);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Koff0);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &K0);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &LP);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Bndaxial);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Torsn);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &KS);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &KV);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &ESS);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &ESM0);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &ESP0);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &PS_M);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &PIP2_M);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &PS_V);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &PIP2_V);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &LDebye0);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Salt);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Eadh);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Radh);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Lc2a);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Rc2a);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Lc2b);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Rc2b);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Ac2a);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Asm);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Ltmd);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &Rtmd);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &eta0);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &eta);
	skip(fp);
	
	n = fscanf(fp, "%s %f", var, &Tmax);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Nsave);
	skip(fp);
	
	n = fscanf(fp, "%s %f, %f, %f", var, &MCprob[0], &MCprob[1], &MCprob[2]);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &MCcycle);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &MCloop);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &MCperloop);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &Width);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &Height);
	skip(fp);
	n = fscanf(fp, "%s %f", var, &ZOOM);
	skip(fp);
	
	n = fscanf(fp, "%s %d", var, &meshskip);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &showskip);
	skip(fp);
	n = fscanf(fp, "%s %d", var, &img_cntm);	// skip time = dt*showskip*img_cntm
	skip(fp);
	
	fclose(fp);
}



void showpara(void)
{
	printf("================\n");
	printf("Nthread = %d\n", Nthread);
	printf("Resume = %d\n", Resume);
	
	printf("SYT = %d\n", SYT);
	printf("Membrane = %d\n", Membrane);
	printf("Vesicle = %d\n", Vesicle);
	printf("TMD = %d\n", TMD);
	
	printf("VesRing = %d\n", VesRing);
	printf("Carbon = %d\n", Carbon);
	
	printf("Lx = %.4g\n", Lx);
	printf("Ly = %.4g\n", Ly);
	printf("Lz = %.4g\n", Lz);
	printf("Nx = %d\n", Nx);
	printf("Ny = %d\n", Ny);
	
	printf("Rv = %.4g\n", Rv);
	printf("Hv = %.4g\n", Hv);
	printf("Hvmax = %.4g\n", Hvmax);
	
	printf("KAmem = %.4g\n", KAmem);
	printf("KBmem = %.4g\n", KBmem);
	printf("Gauss = %.4g\n", Gauss);
	printf("Gamma = %.4g\n", Gamma);
	printf("Beta = %.4g\n", Beta);
	
	printf("KMspon = %.4g\n", KMspon);
	printf("P0 = %.4g\n", P0);
	printf("Hcell = %.4g\n", Hcell);
	printf("h_belkin_fct = %.4g\n", h_belkin_fct);
	
	printf("Limit = %d\n", Limit);
	printf("C0 = %.4g\n", C0);
	printf("NC0 = %d\n", NC0);
	printf("NS0 = %d\n", NS0);
	printf("NS1 = %d\n", NS1);
	
	printf("Spiral = %d\n", Spiral);
	printf("GrwEnd = %d\n", GrwEnd);
	
	printf("Knuc = %.4g\n", Knuc);
	printf("Kon_C = %.4g\n", Kon_C);
	printf("Koff0 = %.4g\n", Koff0);
	
	printf("K0 = %.4g\n", K0);
	printf("LP = %.4g\n", LP);
	printf("Bndaxial = %.4g\n", Bndaxial);
	printf("Torsn = %.4g\n", Torsn);
	printf("KS = %.4g\n", KS);
	printf("KV = %.4g\n", KV);
	
	printf("ESS = %.4g\n", ESS);
	printf("ESM0 = %.4g\n", ESM0);
	printf("ESP0 = %.4g\n", ESP0);
	
	printf("PS_M = %.4g\n", PS_M);
	printf("PIP2_M = %.4g\n", PIP2_M);
	
	printf("PS_V = %.4g\n", PS_V);
	printf("PIP2_V = %.4g\n", PIP2_V);
	
	printf("LDebye0 = %.4g\n", LDebye0);
	printf("Salt = %.4g\n", Salt);
	
	printf("Eadh = %.4g\n", Eadh);
	printf("Radh = %.4g\n", Radh);
	
	printf("Lc2a = %.4g\n", Lc2a);
	printf("Rc2a = %.4g\n", Rc2a);
	printf("Lc2b = %.4g\n", Lc2b);
	printf("Rc2b = %.4g\n", Rc2b);
	printf("Ac2a = %.4g\n", Ac2a);
	
	printf("Ltmd = %.4g\n", Ltmd);
	printf("Rtmd = %.4g\n", Rtmd);
	
	printf("eta0 = %.4g\n", eta0);
	printf("eta = %.4g\n", eta);
	
	printf("Tmax = %.4g\n", Tmax);
	printf("Nsave = %d\n", Nsave);
	
	printf("MCprob = %.3g, %.3g, %.3g\n", MCprob[0], MCprob[1], MCprob[2]);
	printf("MCcycle = %d\n", MCcycle);
	printf("MCloop = %d\n", MCloop);
	printf("MCperloop = %d\n", MCperloop);
	
	printf("Width = %d\n", Width);
	printf("Height = %d\n", Height);
	printf("ZOOM = %.4g\n", ZOOM);
	
	printf("meshskip = %d\n", meshskip);
	printf("showskip = %d\n", showskip);
	printf("saveskip = %d\n", img_cntm);
}



void showconst(void)
{
	printf("----------------\n");
	printf("dt = %.4g ns\n", dt*1e9);
	printf("ESMpersite = %.4g kT\n", ESMpersite/kBT);
	printf("ESVpersite = %.4g kT\n", ESVpersite/kBT);
	if(Membrane==1) printf("dA_mem = %.4g\n", Lx*Ly/Nx/Ny);
	if(Vesicle==1) printf("dA_ves = %.4g\n", dAvesNode0);
	printf("Mumem = %.4g\n", Mumem);
	printf("KSMudt = %.4g\n", KSMudt);
	printf("Nflip_mem = %d\n", Nflip_mem);
	printf("Nflip_ves = %d\n", Nflip_ves);
	printf("================\n");
}




void getconstmem(void)
{
	double dri;
	
	ESM0*=(float)(kBT);	// energy in pN*nm
	ESP0*=(float)(kBT);
	Eadh*=(float)(kBT);
	
	ESM = (PS_M+4*PIP2_M)/25.0*ESM0;	// 25%PS is the standard value for membrane
	ESV = (PS_V+4*PIP2_V)/25.0*ESM0;	// PS=1e, PIP2=4e, 140 mM salt: 15%PS is the standard value for vesicle
	
	KBmem*=(float)(kBT);	// now kBmem is in unit of pN*nm !!!
	KGmem=Gauss*KBmem;	// Gaussian modulus
	BKmod=1.0/Beta;	// bulk modulus of water
	KMsponhalf=0.5*KMspon;
	
	Nxp1=Nx+1;
	Nyp1=Ny+1;
	Nxm1=Nx-1;
	Nym1=Ny-1;
	Nxm2=Nx-2;
	Nym2=Ny-2;
	NxNy=Nx*Ny;
	Nxm1Nym1=Nxm1*Nym1;
	
	NMface = 2*NxNy;
	NodeZoneMaxm1 = NodeZoneMax-1;
	
	Lxhalf=0.5*Lx;
	Lyhalf=0.5*Ly;
	Lxhalf2=Lxhalf*Lxhalf;
	Lyhalf2=Lyhalf*Lyhalf;
	dX=Lx/Nxm1;
	dY=Ly/Nym1;
	dX2=dX*dX;
	dY2=dY*dY;
	dXdb=2.0*dX;
	dYdb=2.0*dY;
	
	dArea=dX*dY;
	AreaGhost=(Nx+Ny-1)*dArea;
	Curv0=1.0/(2*min2(dX,dY));
	
	if(SYT && Membrane && max2(dX,dY)>3*Lc2b/2.0) {
		printf("============================================\n");
		printf("Mem grid is too sparse\n");
		printf("Make sure Nx & Ny > %d\n", (int)(max2(Lx,Ly)/(3*Lc2b/2.0)+0.5));
		printf("============================================\n");
		exit(0);
	}
	
	// each dArea has 2 triangles
	//LBmemmin=sqrt(dArea)/1.32;	// see Kroll & Gompper (Science 1992) and Fosnaric et al. (JCP 2009)
	LBmemmin=sqrt(2/sqrt(3)*dArea)/1.32;	// corresponds to the diameter of hard spheres
	LBmemmax=1.68*LBmemmin;
	Lflip_mem=0.15*LBmemmin;
	
	LBmemmin2=LBmemmin*LBmemmin;
	LBmemmax2=LBmemmax*LBmemmax;
	
	KBmemhalf=0.5*KBmem;
	KBmemdX=KBmem*dX;
	KBmemdY=KBmem*dY;
	KBmemdXYave=KBmem*(dX+dY)/2.0;
	
	KBmemdb=2.0*KBmem;
	KBmemx4=4.0*KBmem;
	
	Kbx=KBmem*dY/dX;	// derived from a bent ribbon: bending of a strip of dy along x-direction
	Kby=KBmem*dX/dY;	// bending of a strip of dx along y-direction
	
	Ksx=KAmem*Nx/Ny;	// tension=KA*dA/A=KA*(2*dL/L)
	Ksy=KAmem*Ny/Nx;
	Ksxydb=2.0*(Ksx+Ksy);
	
	Gammahalf=Gamma/2.0;
	Gamma_sqrt2=Gamma/sqrt(2.0);
	
	F0x=Gamma*dY;	// tension in each x-spring
	F0y=Gamma*dX;	// tension in each y-spring
	
	// mnode-carbon adhesion
	Fadhmax = Eadh/Radh;	// Fadh=Eadh/Radh: force density (per area)
	
	// mnode-syt adhesion: electrostatic w/ Debye screening phi=a/r*exp[-(r-rmin)/lambda]
	LDebye=LDebye0*sqrt(140.0/Salt);	// physiological salt ~ 140 mM
	LDebyeinv=1.0/LDebye;
	LDebye2=LDebye*LDebye;
	LDebye3=LDebye2*LDebye;
	LDebye4=LDebye3*LDebye;
	LDebye5=LDebye4*LDebye;
	LDebye6=LDebye5*LDebye;
	LDebye7=LDebye6*LDebye;
	LDebye8=LDebye7*LDebye;
	
	ESMpersite0=ESM;	// syt-t_membrane energy in pN*nm per syt site under physiological salt
	ESVpersite0=ESV;	// syt-vesicle energy in pN*nm per syt site under physiological salt
	
	ESMpersite=ESMpersite0*LDebye/LDebye0;	// binding energy in current salt
	ESVpersite=ESVpersite0*LDebye/LDebye0;
	
	DebyeMemCoeff=ESMpersite0/Pidb/LDebye0;	// coefficient for syt-mnode attraction = q1*q2/4pi*e0*er
	DebyeVesCoeff=ESVpersite0/Pidb/LDebye0;	// coefficient for syt-ves attraction
	
	Inv4Pie0e = 1.0/(4*Pi*8.85e-12*80);	// 1/4*pi*e_0*e in SI unit
	ESP=ESP0*LDebye/LDebye0;	// Syt-PIP2 1-to-1 binding affinity under current salt, in pN*nm
	
	// ESP = [q1*q2/(4*pi*e0*e)/r0]*exp(-r0/LDebye) = [a/r0]*exp(-r0/LDebye) --> r0 = LDebye*LambertW(a/ESP/LDebye)
	QsytQpip2_4pie0eSI = Q_KKKK*Q_PIP2*pow(e_charge,2)*Inv4Pie0e;	// in SI unit: J*m
	QsytQpip2_4pie0e = QsytQpip2_4pie0eSI * 1e30;	// in pN*nm*nm
	
	r0_PIP2 = LDebye * LambertW(QsytQpip2_4pie0e/ESP/LDebye);	// cut-off distance for KKKK-PIP2 interaction in nm (E=ESP at this distance)
	
	
	RSMmax=4*LDebye+sqrt(dX2+dY2);	//3*LDebye;	// max range for syt-membrane attraction
	RSMmax2=RSMmax*RSMmax;
	
	if(SYT && Membrane && RSMmax<=max2(dX,dY)) {
		printf("============================================\n");
		printf("Syt-Mem force range %.3g < grid size %.3g\n", RSMmax, min2(dX,dY));
		printf("Make sure Nx & Ny > %d\n", (int)(max2(Lx,Ly)/RSMmax+0.5));
		printf("============================================\n");
		exit(0);
	}
	
	NSMnbr=(int)(3*(RSMmax+DLThalf)/dX+0.5);	// number of grids to search in x or y: (-NSMnbr/2, NSMnbr/2)
	if(NSMnbr%2==0) NSMnbr++;	// make sure NSMnbr is odd
	NSMnbr=max2(NSMnbr,3);
	NSMnbr2=NSMnbr*NSMnbr;
	NSMnbrhalf=(NSMnbr-1)/2;
	//printf("%d\t%d\n", NSMnbr, NSMnbrhalf);
	
	h_belkin=h_belkin_fct*Lx;
	h_belkinx4 = h_belkin*4.0;
	h_belkincutoff = 3.0*h_belkin;
	h_belkincutoff2 = h_belkincutoff*h_belkincutoff;
	Pih_belkin2x12 = 12.0*Pi*pow(h_belkin,2);
	NMnbrzoneplushalf=(int)(h_belkincutoff/dX+0.5);	// number of zones to search in x or y (for Belkin's formula)
	NMnbrzoneplushalf=max2(NMnbrzoneplushalf,1);
	NMnbrzoneplus=2*NMnbrzoneplushalf;
	
	if(LAPLAC==0) printf("h_belkin=%.4g,\tNMnbrzoneplushalf=%d\n", h_belkin, NMnbrzoneplushalf);
	
	Rinf2d[0] = 3*(Lx+Ly+Rv);
	Rinf2d[1] = 2*(Lx+Ly+Rv);
	
	LxpdX=Lx+dX;
	LypdY=Ly+dY;
	drmemmin=0.2*min2(dX,dY);
	
	dri=sqrt(dArea/Pi);	// radius of each membrane patch
	
	//d=5.0;	// membrane thickness
	//Mumem = 1.0/(6*Pi*eta*sqrt(dX*d/2)/2);	// rough estimate of mobility constant (5nm is the thickness of membrane)
	Mumem = 1.0/(4*Pi*eta*dX);	// assuming it's a stick of length dX and thickness d (Berg)
	Mumem0 = 1.0/(16*eta0*dri);	// a patch of membrane moving in solution
	
	Mumemr0 = 3.0/(32*eta0*pow(dri,3));	// rotational
	
	MumemPerp = 1.0/(16*eta0*dri);
	MumemPara = log(eta/eta0/dri)/(4*Pi*eta);	// in-plane mobility: zeta=4*pi*eta/log(eta/eta0/ri)
	
	Dmem = kBT*MumemPara;	// diffusion of membrane patch in membrane
	Dmem0 = kBT*Mumem0;	// diffusion of membrane patch in solution
}




void getconstvesicle(void)
{
	double dri;
	
	Ndiv=max2(MeshLevel-1,0);
	
	NVnode=(int)(10*pow(4,Ndiv)+2);
	NVface=(int)(20*pow(4,Ndiv));
	NVedge=(int)(30*pow(4,Ndiv));
	
	KBves=KBmem;	// in pN*nm, after defining KBmem=KBmem*kBT
	KBvesdb=2.0*KBves;
	KBvesx4=4.0*KBves;
	VolVes0est=Pi*4.0/3*pow(Rv,3);
	
	AreaVes0est=4.0*Pi*Rv*Rv;
	dAvesNode0=AreaVes0est/NVnode;
	dAvesFace0=AreaVes0est/NVface;
	dRves=sqrt(4*dAvesFace0/sqrt(3));	// side length of a equilateral triangle: area=sqrt(3)/4*dr^2
	dRveseps=0.05*dRves;
	dri=sqrt(dAvesNode0/Pi);	// radius of each patch
	
	Eves0=8*Pi*KBves;	// bending energy
	dEves0=Eves0/NVnode;

#if !VSTRCH
	GammaVes = Gamma;
#endif

	drSVrepulsion=max2((Rc2a+Rc2b)/2.0, 0.6*dRves);
#if XREPUL
	drSVrepulsion=max2(drSVrepulsion, 5.0);	// extra repulsion due to proteins on vesicle
#endif
	
	drVMrepulsion=0.2*max2(max2(dX, dY), dRves);
	drVMrepulsion=max2(drVMrepulsion, 1.0);	// at least 1.0 nm for hydration repulsion
	
	drSVrepulsion2=drSVrepulsion*drSVrepulsion;
	drVMrepulsion2=drVMrepulsion*drVMrepulsion;
	
	NVMnbr=(int)(2.5*drVMrepulsion/dX+0.5);	// number of grids to search in x or y: (-NVMnbr/2, NVMnbr/2)
	if(NVMnbr%2==0) NVMnbr++;	// make sure NVMnbr is odd
	NVMnbr=max2(NVMnbr,3);
	NVMnbr2=NVMnbr*NVMnbr;
	NVMnbrhalf=(NVMnbr-1)/2;
	
	LBvesmin=dRves/1.32;	// see Kroll & Gompper (Science 1992) and Fosnaric et al. (JCP 2009)
	LBvesmax=1.68*LBvesmin;	// ~sqrt(2.8)
	Lflip_ves=0.15*LBvesmin;
	
	LBvesmin2=LBvesmin*LBvesmin;
	LBvesmax2=LBvesmax*LBvesmax;
	
	KBvesdr=KBves*dRves;
	
	ZetaVesT=16*eta0*(dRves/2.0);	// translational friction of a disc of diameter dRves
	ZetaVesR=32.0/3*eta0*pow(dRves/2.0, 3);	// rotational friction of a disc through its center
	
	MuVesT=1.0/ZetaVesT;
	MuVesR=1.0/ZetaVesR;
	
	MuVesTotT=1.0/(6*Pi*eta0*Rv);
	MuVesTotR=1.0/(8*Pi*eta0*pow(Rv,3));
	
	/*
	effective rotational friction of a disc through an axis of 2R=drmem away from its center
	tau=zetareff*w=2R*F+zetar*w, F=zetat*v=zetat*w*2R. so, zetareff=4R^2*zetat+zetar
	from Berg's book, zetat=16*eta*R, zetar=32/3*eta*R^3. so, zetareff=224/3*eta*R^3
	*/
	ZetaVesReff=224.0/3*eta0*pow(dRves/2.0, 3);
	
	MuVesReff=1.0/ZetaVesReff;
	MuVesReff_T=MuVesReff/MuVesT;
	//printf("MuVesR=%.3g, MuVesReff=%.3g, MuVesT=%.3g, MuVesReff_T=%.3g\n", MuVesR, MuVesReff, MuVesT, MuVesReff_T);
	
	Dves0 = Dves = kBT*MuVesT;
	DvestotT = kBT*MuVesTotT;
	DvestotR = kBT*MuVesTotR;
	
	
	MuvesPerp=1.0/(16*eta0*dri);
	MuvesPara=log(eta/eta0/dri)/(4*Pi*eta);	// in-plane mobility: zeta=4*pi*eta/log(eta/eta0/ri)
	
	Dves = kBT*MuvesPara;	// diffusion of membrane patch in membrane
}



void getconstsyt(void)
{
	ESS*=(float)(kBT);
	
	Ac2arad=Ac2a*d2r;	// in radian
	Asmrad=Asm*d2r;
	
	SinAc2a=sin(Ac2arad);
	CosAc2a=cos(Ac2arad);
	SinAsm=sin(Asmrad);
	CosAsm=cos(Asmrad);
	
	Rc2adb=2*Rc2a;
	Rc2bdb=2*Rc2b;
	Lc2ahalf=0.5*Lc2a;
	Lc2bhalf=0.5*Lc2b;
	dc2ahalf=Lc2ahalf-Rc2a;	// half length of cylindrical part
	dc2bhalf=Lc2bhalf-Rc2b;
	
	Rc2b2=Rc2b*Rc2b;
	Rc2abAVE=(Rc2a+Rc2b)/2.0;
	Rc2abAVE2=Rc2abAVE*Rc2abAVE;
	
#if SYTCNT
	// cetner between C2A & C2B
	RCc2a1[0] = dc2ahalf*SinAc2a;	// top C2A bead
	RCc2a1[1] = Rc2a;
	RCc2a1[2] = dc2ahalf*CosAc2a;
	
	RCc2a2[0] = -RCc2a1[0];	// bottom C2A bead
	RCc2a2[1] = RCc2a1[1];
	RCc2a2[2] = -RCc2a1[2];
	
	RCc2b1[0] = -dc2bhalf;	// left C2B bead
	RCc2b1[1] = -Rc2b;
	RCc2b1[2] = 0;
	
	RCc2b2[0] = dc2bhalf;	// right C2B bead
	RCc2b2[1] = -Rc2b;
	RCc2b2[2] = 0;
	
	RCsm[0] = 0;	// mem binding site
	RCsm[1] = -Rc2b*(1+CosAsm);
	RCsm[2] = -Rc2b*SinAsm;
#else
	// cetner in C2B
	RCc2a1[0] = dc2ahalf*SinAc2a;	// top C2A bead
	RCc2a1[1] = Rc2a+Rc2b;
	RCc2a1[2] = dc2ahalf*CosAc2a;
	
	RCc2a2[0] = -RCc2a1[0];	// bottom C2A bead
	RCc2a2[1] = RCc2a1[1];
	RCc2a2[2] = -RCc2a1[2];
	
	RCc2b1[0] = -dc2bhalf;	// left C2B bead
	RCc2b1[1] = 0;
	RCc2b1[2] = 0;
	
	RCc2b2[0] = dc2bhalf;	// right C2B bead
	RCc2b2[1] = 0;
	RCc2b2[2] = 0;
	
	RCsm[0] = 0;	// mem binding site
	RCsm[1] = -Rc2b*CosAsm;
	RCsm[2] = -Rc2b*SinAsm;
#endif
	
	//RSSselfavoid=sqrt(Rc2b*Rc2b+Lc2bhalf*Lc2bhalf);
	RSSselfavoid=2*Rc2b;
	RSSselfavoid2=RSSselfavoid*RSSselfavoid;
	
	Nsytperring0=(int)(Pidb/max2(K0,1e-6)/Lc2b+0.5);	// # of syt per ring
#if CLOSED
	if(Limit==0) Nsytperring=NS0;
	else Nsytperring=NS1;
#else
	Nsytperring=Nsytperring0;
#endif
	
	NSegMaxm1=NSegMax-1;
	
	Lxhalf=Lx/2.0;
	Lyhalf=Ly/2.0;
	Area=Lx*Ly;
	
	DLT=Lc2b;	// temp !
	DLThalf=DLT/2.0;
	DLThalf2=DLThalf*DLThalf;
	DLT2=DLT*DLT;
	
	DLThalfmHSM=DLThalf;
	DLThalfmHSM2=DLThalfmHSM*DLThalf;
	
	Lxhalfplus=1.05*Lxhalf;	// for display purpose
	Lyhalfplus=1.05*Lyhalf;
	R0=1.0/max2(K0,eps);
	R02=R0*R0;
	
	R0mDLThalf=R0-DLThalf;
	R0mDLThalf2=R0mDLThalf*R0mDLThalf;
	
	DistMin=0.2*DLT;
	
	dTHclose=30.0*d2r;	// max angle allowed for chain closure
	csdTHclose=cos(dTHclose);
	sndTHclose=sin(dTHclose);
	Distclose=1.5*Lc2b;
	Distclose2=Distclose*Distclose;
	
	TH=DLT*K0;
	sTH=sin(TH);
	cTH=cos(TH);
	
	NM0=(int)(Lx*Ly*Lz*(C0*602));	// initial number of syt monomers
	Kon0=Kon_C*C0;
	//if(NC0==1) Tmax=min2(Tmax, NSegMax/(Kon0*2));	// max simulation time for 1 chain
	
	KVinv=1.0/KV;
	KVhalf=0.5*KV;
	KBsyt=kBT*LP/Lc2b;	// bending stiffness between adjacent subunits on ring
	KBsytaxial=KBsyt*Bndaxial;
	Ktorsion=Torsn*KBsyt;
	KBsytmax=KBsyt+KBsytaxial+Ktorsion;
	
	Zeta=10*eta*DLThalf;	// a disc moving in membrane
	Zeta0=6*Pi*eta0*DLThalf;	// a sphere moving in solution
	
	D=kBT/Zeta;	// syt diffusion constant in solution
	D0=kBT/Zeta0;	// syt diffusion constant in solution
	
	Mu=1.0/Zeta;
	Mu0=1.0/Zeta0;
	
	// rotational diffusion constant of a syt
	DR=kBT/(10*eta*pow(DLThalf,3));	// a disc moving in membrane
	DR0=kBT/(8*Pi*eta0*pow(DLThalf,3));	// a sphere moving in solution
	
	Murot=DR/kBT;
	Murot0=DR0/kBT;
	
	KShalf=KS/2.0;
	KBsythalf=KBsyt/2.0;
	KBsytaxialhalf=KBsytaxial/2.0;
	KSeff=2*KS+4*KBsyt*(1+Bndaxial+Torsn)/3.0/DLT2;	// effective total spring constant calc'd from energy
	KSeffinv=1.0/KSeff;
	KSeffinvmin=min2(KSeffinv, Mu);
	KSDLThalf=KS*DLThalf;
	
	kssinv=min2(1.0/(2*KS+2*KV), Mu);
	kbsinv=min2(1.0/(KBsyt*(1+Bndaxial+Torsn)/3.0/DLT2), Mu);
}



void getconstlinkertmd(void)
{
	Llinker0=Nlinker*dLlinker;	// contour length of linker
	Flinker0=kBT/0.6;	// force prefactor in Marko-Siggia WLC model
	
	rlinkercutoff=0.95;
	Flinkermax=Flinker0*(0.25/pow(1-rlinkercutoff,2)+rlinkercutoff-0.25);
	
	ZetaTMD=4*Pi*eta*Ltmd/(log(Ltmd/Rtmd)+0.5);
	DTMD=kBT/ZetaTMD;
	MuTMD=1.0/ZetaTMD;
	printf("Flinker0 = %.3g pN\n", Flinker0);
	printf("DTMD = %.3g um^2/s\n", DTMD*1e-6);
}



double func_rcutoff(double x, double a, double d)
{	// x=distance from center, a=area of vertex, d=mesh size
	double d2, x2, y1, y2, y3, y4, y5, y6, rin, y;
	
	d2 = d*d;
	x2 = x*x;
	
	y1 = exp(-x/LDebye)/x;	// center node
	
	y2 = exp(-sqrt(x2+d2-sqrt(3)*x*d)/LDebye)/sqrt(x2+d2-sqrt(3)*x*d);	// 1st ring neighbors in the same plane
	y3 = exp(-sqrt(x2+d2)/LDebye)/sqrt(x2+d2);
	y4 = exp(-sqrt(x2+d2+sqrt(3)*x*d)/LDebye)/sqrt(x2+d2+sqrt(3)*x*d);
	
	y5 = exp(-sqrt(x2+3*d2)/LDebye)/sqrt(x2+3*d2);	// 2nd ring neighbors
	y6 = exp(-sqrt(x2+4*d2)/LDebye)/sqrt(x2+4*d2);
	
	// use continuum approx. to estimate potential of nodes outside the 2 rings
	rin = sqrt(19*a/Pi);	// radius of the 2 inner rings
	
	y = y1 + 2*y2 + 2*y3 + 2*y4 + 6*y5 + 6*y6 - LDebye*Pidb/a*(1-exp(-rin/LDebye));
	
	return y;
}



double solve_rcutoff(double a, double d)
{
	double rcut, x1, x2;
	x1 = 1e-3;
	x2 = 10;
	
	rcut = zbrent3(func_rcutoff, a, d, x1, x2, 1e-4);	// a=vertex area, d=mesh size, x1-x2=search range
	return rcut;
}



void getconstSytMemVes(void)
{
	int i;
	double am, av;
	
	amin_cutoff_sm = 0.2*(2*sqrt(3)/4.0*LBmemmin2);	// factor 2 comes from 6 triangles w/ 1/3 area each for vertex's area
	amax_cutoff_sm = 2.0*(2*sqrt(3)/4.0*LBmemmax2);
	
	amin_cutoff_sv = 0.2*(2*sqrt(3)/4.0*LBvesmin2);	// factor 2 comes from 6 triangles w/ 1/3 area each for vertex's area
	amax_cutoff_sv = 2.0*(2*sqrt(3)/4.0*LBvesmax2);
	
	da_cutoff_sm = (amax_cutoff_sm - amin_cutoff_sm) / (n_cutoff_smv - 1);
	da_cutoff_sv = (amax_cutoff_sv - amin_cutoff_sv) / (n_cutoff_smv - 1);
	
	for(i=0; i<n_cutoff_smv; i++) {
		am = amin_cutoff_sm + i*da_cutoff_sm;
		r_cutoff_sm[i] = solve_rcutoff(am, LBmemmin);
		
		av = amin_cutoff_sv + i*da_cutoff_sv;
		r_cutoff_sv[i] = solve_rcutoff(av, LBvesmin);
		
		//printf("%.3g\t%.3g\t%.3g\t%.3g\n", am, r_cutoff_sm[i], av, r_cutoff_sv[i]);
	}
}



void getconstMC(void)
{
	int i;
	double psum;
	
	psum = 0;
	for(i=0; i<3; i++) psum += MCprob[i];
	
	MCp[0] = MCprob[0]/psum; 
	for(i=1; i<3; i++) MCp[i] = MCprob[i]/psum + MCp[i-1];
	
	//for(i=0; i<3; i++) printf("MCp[%d]=%.4g\n", i, MCp[i]);
	
	MCdr = 0.06*dX;
	
	MCnmode0 = (int)(1.0*Lxhalf/dX);
}



void getconst(void)
{
	double dft, dfr, dl1, dl2, dl;
	double cnt;
	
	getconstmem();
	getconstsyt();
	getconstvesicle();	// run this after getconstmem()
	getconstlinkertmd();
	
#if ATTRCT
	getconstSytMemVes();
#endif
	
#if !LGV_MC
	getconstMC();
#endif

	//drSMrepulsion=max2(Rc2abAVE, 0.5*max2(dX,dY));
	drSMrepulsion=Rc2abAVE;
	drSMrepulsion2=drSMrepulsion*drSMrepulsion;
	
//----------- dt -------------
	dt = 1.0;
	
#if FRICTN	// friction from membrane
	dft=D;
	dfr=DR;
#else	// friction from cytosol
	dft=D0;
	dfr=DR0;
#endif
	
	// SYT
	if(SYT && FIXSYT!=2) {
		dt = min2(dt, pow(0.1*DLT,2)/2/dft);	// diffusion of SYT in dt is less than 0.1*DLT
		dt = min2(dt, 0.1/(eps+max2(Knuc*Area, max2(Kon0, Koff0))));	// limited by reaction rates
		
		if(FIXSYT==0) {
			dt = min2(dt, 0.4*kBT*dft/2.0/(KS+(KBsyt+KBsytaxial+Ktorsion)/DLT2));	// stretching & bending due to translational motion
			dt = min2(dt, 0.4*kBT/dfr/2.0/(KS*DLT2+KBsyt+KBsytaxial+Ktorsion));	// bending due to rotational motion
		}
		else if(FIXSYT==1) {	// only allow stretching and bending about local z
			dt = min2(dt, 0.4*kBT*dft/2.0/(KS+KBsyt/DLT2));	// stretching & bending due to translational motion
			dt = min2(dt, 0.4*kBT/dfr/2.0/(KS*DLT2+KBsyt));	// bending due to rotational motion
		}
		dt = min2(dt, 0.4*pow(min2(DLT,R0),4)/(LP*DLT*dft));	// stability of 4th order diffusion equation (biharmonic) for syt ring
	}
	
	if(Membrane) {
		dt = min2(dt, 0.4/Mumemr0/(Kbx+Kby));	// from membrane bending: dth=mur*kb*th*dt < th on each grid
		
		//dl1 = -LDebye*log(1-min2(kBT/max2(ESM,kBT),0.9999));	// for each binding site: E(h)=Esm*exp(-h/lambda), dl1 satisfies E(0)-E(dl1)=1kT for inf. plane
		dl1 = drmemmin;
		dl2 = drmemmin;
		dl = min2(dl1, dl2);
		dt = min2(dt, 0.4*dl*dl/Dmem0);	// limited by diffusion-induced energy change
		
		// Helfrich pressure: dP ~ 2KB*H" ~ 2*KB*(1/Rc2b)/dlt^2, so force f ~ dP*dlt^2 ~ 2*KB/R2cb, so dt<dl/mu/f
		dt = min2(dt, 0.4*dl/Mumem0/(2*KBmem/Rc2b));
	}
	
	if(Vesicle) {
		dt = min2(dt, 0.4/MuVesR/(Kbx+Kby));	// from ves bending: dth=mur*kb*th*dt < th on each grid
		
		dt = min2(dt, pow(0.2*dX,2)/2/(kBT*MuVesTotT)); // translational diffusion of whole vesicle
		dl1 = -LDebye*log(1-min2(kBT/max2(ESV,kBT),0.9999));	// for each binding site: E(h)=Esv*exp(-h/lambda), dl1 satisfies E(0)-E(dl1)=1kT for inf. plane
		dl2 = dRveseps;
		dl = min2(dl1, dl2);
		dt = min2(dt, 0.4*dl*dl/Dves0);	// limited by diffusion-induced energy change
		// Helfrich pressure: dP ~ 2KB*H" ~ 2*KB*(1/Rc2b)/dlt^2, so force f ~ dP*dlt^2 ~ 2*KB/R2cb, so dt<dl/mu/f
		dt = min2(dt, 0.4*dl/MuVesT/(2*KBves/Rc2b));
	}
	
	if(TMD) {
		dt = min2(dt, pow(0.2*dX,2)/2/DTMD);	// diffusion of TMD on vesicle
	}
	
#if BD_2ND
	dt *= 0.8;	// test !!!
#else
	dt *= 0.8;	// test !!!
#endif

//----------- end of dt -------------
	
	// below requires dt
	dt_dX=dt/dX;
	dt_dY=dt/dY;
	dt_dXdb=dt/dXdb;
	dt_dYdb=dt/dYdb;
	
	dl1=sqrt(kBT/Ksxydb);	// from membrane stretching
	dl2=sqrt(2*kBT/(Kbx/dX2+Kby/dY2));	// from membrane bending
	dlnoise=sqrt(3.0)*(dl1*dl2)/(dl1+dl2);
	//dlnoise=sqrt(2*kBT/(Kbx/dX2+Kby/dY2));	// from membrane bending
	//dlnoise*=0.2;	// test !!!!!!!!!!!!
	
	Mumemdt=Mumem*dt;
	Mumem0dt=Mumem0*dt;
	
	MumemPerpdt=MumemPerp*dt;
	MumemParadt=MumemPara*dt;
	
	//MumemXdt=Mumemdt;
	//MumemXdt=Mumem0dt;
	//MumemXdt=sqrt(Mumemdt*Mumem0dt);
	
	Tflip_mem=0.5*pow(Lflip_mem,2)/(kBT*Mumem0);
	Nflip_mem=(int)(Tflip_mem/dt+0.5);	// waiting steps for bond flipping
	Nflip_mem=max2(Nflip_mem,1);
	
	Tflip_ves=0.5*pow(Lflip_ves,2)*(ZetaVesT/kBT);
	Nflip_ves=(int)(Tflip_ves/dt+0.5);	// waiting steps for bond flipping
	Nflip_ves=max2(Nflip_ves,1);
	
	drmax=0.2*min2(dX,dY);	// for curvature fitting
	
	Pon0 = 2*Kon0*dt;	// b/c of 2 open ends
	Pon0half=Pon0/2.0;
	Poff0 = Koff0*dt;
	Pnuc0 = Knuc*Area*dt;
	
	Mudt=Mu*dt;
	Mu0dt=Mu0*dt;
	
	Murotdt=Murot*dt;
	Murot0dt=Murot0*dt;
	
	KSMudt=KS*Mudt;
	KSMu0dt=KS*Mu0dt;
	
	LPDLTDdt=LP*DLT*D*dt;
	
	MuTMDdt=MuTMD*dt;
	
	MuVesTotTdt=MuVesTotT*dt;
	MuVesTotRdt=MuVesTotR*dt;
	MuVesTdt=MuVesT*dt;
	
	MuvesPerpdt=MuvesPerp*dt;
	MuvesParadt=MuvesPara*dt;
	
	dRnoise=sqrt(3.0)*sqrt(2*D*dt);	// sqrt(3) comes from RND=U(-1,1)
	dAnoise=sqrt(3.0)*sqrt(2*DR*dt);	// in membrane
	
	dRnoise0=sqrt(3.0)*sqrt(2*D0*dt);	// sqrt(3) comes from RND=U(-1,1)
	dAnoise0=sqrt(3.0)*sqrt(2*DR0*dt);	// in solution
	
	dTMDnoise=sqrt(3.0)*sqrt(2*DTMD*dt);
	
	dVRnoise=sqrt(3.0)*sqrt(2*DvestotT*dt);
	dVAnoise=sqrt(3.0)*sqrt(2*DvestotR*dt);
	
	//SaveCntMax=Tmax/dt/Nsave+0.1;	// this many simulations between saving data (this could be a huge number)
	cnt=Tmax/dt/Nsave;	// this could be a very large number
	
	// use n=n1*1e6+n2 to represent large number
	SaveCntMax1=(long)(cnt/1e6);	// floor of double
	SaveCntMax2=(long)(cnt-1e6*SaveCntMax1);
	//printf("SaveCntMax=%d,\t%d\n", SaveCntMax1, SaveCntMax2);
	
	HWratio=(float)(1.0f*Height/Width);
	
	showconst();
}

