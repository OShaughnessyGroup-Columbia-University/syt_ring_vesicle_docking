
// SYT
#define Q_KKKK 4	// # of charges in KKKK
#define Q_PIP2 4	// # of charges in PIP2

#define Nlinker 60	// # of residues in syt's linker region
#define dLlinker 0.365	// length per residue in linker in nm

#define NSegMax 80	// max number of segments in each ring
#define StatBatch 20

// membrane
#define h_mem_ini 2.0	// initial membrane height if INIMEM=1
#define NodeZoneMax 30	// max # of nodes per zone
#define ContourBin 41	// # of bins for membrane contour

#define NnbrMax 12	// max # of neighbors per mnode or vnode
#define NnbrMin 4	// min # of neighbors per mnode or vnode
#define NnbrMaxm1 (NnbrMax-1)
#define NnbrMaxp1 (NnbrMax+1)

#define NnbrMaxPlus 80	// max # of nbr faces per node, for Belkin's formula

// syt-mem interaction
#define n_cutoff_smv 201	// # of bins for r_cutoff interpolation (syt-mem/-ves PP interaction)

// constants
#define kBT 4.1	// 1kT in pN*nm
#define e_charge 1.6e-19	// e in Coulomb

#define Pi 3.14159265
#define Pihalf (Pi/2.0)
#define Piquad (Pi/4.0)
#define Pidb (2.0*Pi)
#define r2d (180.0/Pi)
#define d2r (Pi/180.0)
#define eps 1.0e-10
#define inf 1.0e10

// system
int Nthread, Resume, SYT, Membrane, Vesicle, TMD;
int VesRing, Carbon;

// membrane constants
int Nx, Ny, MeshLevel;
float Lx, Ly, Lz;
float Rv, Hv, Hvmax;
float KAmem, KBmem, Gauss, Gamma, Beta;
float KMspon, P0, Hcell;
float h_belkin_fct;

// syt constants
int Limit, NC0, NS0, NS1, Spiral, GrwEnd;
float C0, Knuc, Kon_C, Koff0;
float K0, LP, Bndaxial, Torsn, KS, KV;
float ESS, ESM0, ESP0, PS_M, PIP2_M, PS_V, PIP2_V, LDebye0, Salt;
float Eadh, Radh;
float Lc2a, Rc2a, Lc2b, Rc2b, Ac2a, Asm;
float Ltmd, Rtmd;
float eta0, eta;

float Tmax;
int Nsave;

float MCprob[3];
int MCcycle, MCloop, MCperloop;

int Width, Height;
float ZOOM;
int meshskip, showskip, saveskip;

// general variables
int noise;
int RunCnt;
double t, dt;
double LDebye, ESV, ESM, ESP;
double Eves0, dEves0;
double Inv4Pie0e, QsytQpip2_4pie0eSI, QsytQpip2_4pie0e, ESPcutoff;
double Esyt, Emem, Eves, Esys;
double EsytES;	// electrostatic energy between SYT & Mem+Ves
double Evesbend, Evestens;
double Emembend, Ememtens;
double Lflip_mem, Lflip_ves;	// value "s" in Kroll & Gomper
double Tflip_mem, Tflip_ves;
int Nflip_mem, Nflip_ves;
int Cntflip_mem, Cntflip_ves;

int NMnbrzoneplus, NMnbrzoneplushalf;	// for Belkin's formula

// RNG
long mseed, *mseedarray;

#if RNG_64
	typedef unsigned long long int Ullong;
	static Ullong U64 = 4101842887655102017LL; // don't use this value as the seed!
	static Ullong *U64MP;
#else
	typedef unsigned int Uint;
	static Uint U32, V32=2244614371U, W132=521288629U, W232=362436069U;
	static Uint *U32MP, *V32MP, *W132MP, *W232MP;
#endif

// membrane variables
int NMface, NodeZoneMaxm1;
int Nxp1, Nyp1, Nxm1, Nym1, Nxm2, Nym2, NxNy, Nxm1Nym1;
double dX, dY, dX2, dY2, dXdb, dYdb;
double Lxhalf, Lyhalf, Lxhalfplus, Lyhalfplus, Lxhalf2, Lyhalf2;
double LxpdX, LypdY;
double KGmem, KBmemhalf;
double KBmemdX, KBmemdY, KBmemdXYave;
double KBmemdb, KBmemx4;
double BKmod;
double drSMrepulsion, drSMrepulsion2;
double dlnoise, drmax, Zmin, Zmax;
double Curvmax, Curvmin, Curv0;
double Fadhmax;
double Gammahalf, Gamma_sqrt2;
double KMsponhalf;

double Area, Area0, dArea, Have, AreaGhost;
double AreaMem, AreaMem0, dAmemNode0, dAmemFace0, dAmemNode, dAmemFace;
double Kbx, Kby, Ksx, Ksy, Ksxydb;
double dt_dX, dt_dY, dt_dXdb, dt_dYdb;
double Mumem, Mumemdt, Mumem0, Mumem0dt, MumemXdt;
double Mumemr0, Mumemr0dt;
double Dmem, Dmem0;
double MumemPerp, MumemPara;
double MumemPerpdt, MumemParadt;

double Peff, Reff, T0, F0x, F0y;
double dlzmax, Ebsim;
double Rinf2d[2], drmemmin;

double LBmemmax, LBmemmin, LBmemmax2, LBmemmin2;
double dHmem, dHmemprev;

double h_belkin, h_belkinx4, Pih_belkin2x12, h_belkincutoff, h_belkincutoff2;

int Ncontour;
float *Rcontour, *Zcontour;

// vesicle variables
int Ndiv, NVnode, NVface, NVedge;
int NnodeMax, NfaceMax, NedgeMax;
int NVMnbr, NVMnbrhalf, NVMnbr2;
double RCves[3], RCvesprev[3], dRCves[3];
double Fves[3], Tauves[3];
double KBves, KBvesdb, KBvesx4, KBvesdr;
double AreaVes, AreaVes0, AreaVes0est;
double GammaVes;
double VolVes, VolVes0, VolVes0est, dVolVes;
double Pves, Pves_3;
double dRves, dRveseps;
double dAvesFace0, dAvesFace, dAvesNode0, dAvesNode;
double drSVrepulsion, drSVrepulsion2;
double drVMrepulsion, drVMrepulsion2;
double LBvesmax, LBvesmin, LBvesmax2, LBvesmin2;
double ZetaVesT, ZetaVesR, ZetaVesReff;
double MuVesT, MuVesR, MuVesReff, MuVesReff_T;
double Dves, Dves0, DvestotT, DvestotR;
double MuVesTdt;
double MuVesTotT, MuVesTotTdt;
double MuVesTotR, MuVesTotRdt;
double dVRnoise, dVAnoise;
double MuvesPerp, MuvesPara;
double MuvesPerpdt, MuvesParadt;


// syt variables
int NM0, NM, Nchain, Nopen, Nclosed;
int NSegMaxm1;
int NSMnbr, NSMnbrhalf, NSMnbr2;
int Nsytperring0, Nsytperring;
float DLT;
double Ac2arad, Asmrad, SinAc2a, CosAc2a, SinAsm, CosAsm;
double Rc2adb, Rc2bdb, Rc2abAVE, Rc2b2, Rc2abAVE2;
double Lc2ahalf, Lc2bhalf, dc2ahalf, dc2bhalf;
double RCc2a1[3], RCc2a2[3], RCc2b1[3], RCc2b2[3], RCsm[3];	// local coordinates of 4 beads
double RSSselfavoid, RSSselfavoid2;
double DLThalf, DLThalf2, DLThalfmHSM, DLThalfmHSM2, DLT2, Ntot0, Ntot;
double LXhalfplus, LYhalfplus;
double KBsyt, KBsythalf, KBsytaxial, KBsytaxialhalf, Ktorsion, KBsytmax;
double R0, R02, KVinv, KVhalf;
double R0mDLThalf, R0mDLThalf2;
double D, DR, Zeta, Mu, Murot;
double D0, DR0, Zeta0, Mu0, Murot0;
double KSeff, KShalf, KSeffinv, KSeffinvmin, DistMin, kssinv, kbsinv;
double TH, sTH, cTH, dTHclose, csdTHclose, sndTHclose;
double dRnoise, dRnoise0, dAnoise, dAnoise0;
double Kon0, Pon0, Pon, Pon0half, Poff0, Pnuc0, Pnuc;
double Distclose, Distclose2;
double Mudt, KSMudt, Murotdt, LPDLTDdt, KSDLThalf;
double Mu0dt, KSMu0dt, Murot0dt;

double RSMmax, RSMmax2, ESMpersite0, ESMpersite, ESVpersite0, ESVpersite;
double DebyeMemCoeff, DebyeVesCoeff;
double LDebyeinv, LDebye2, LDebye3, LDebye4, LDebye5, LDebye6, LDebye7, LDebye8;

double r0_PIP2;
double FsmCoeff, FsvCoeff;

double Llinker0, Flinker0, rlinkercutoff, Flinkermax;
double DTMD, ZetaTMD, MuTMD, MuTMDdt, dTMDnoise;

// syt-mem interaction
double r_cutoff_sm[n_cutoff_smv], r_cutoff_sv[n_cutoff_smv];
double amin_cutoff_sm, amax_cutoff_sm, da_cutoff_sm;
double amin_cutoff_sv, amax_cutoff_sv, da_cutoff_sv;


// system
double TrackFves, TrackFsyt, TrackFmem;

// 2nd-order BD
double Fves_old[3], Tauves_old[3];

//long SimCnt, SaveCntMax;
long SimCnt1, SimCnt2, SaveCntMax1, SaveCntMax2;	// use 2 parts to handle very large numbers

int Nfile;
time_t start_time, end_time;

// MC
int MCcyclecnt, MCloopcnt;
double MCdr, MCp[3];
double temperature;
int MCnmode0;


//----- GL -----
float HWratio;

int savecnt;
int attracton, showdiagonal;
int pause, showcolor, showden, showmem, shownormal;
int showvesicle, showsyt, showcoord, showmark, showforce, showtmd;
int showboxsphere;
int mouse_x, mouse_y, mouse_left, mouse_right;

float zoom, zoominv, PanX, PanY;
float rotphi, rottheta;

int nframe;
int img_cnt, img_cntm;
char img_name[20];

#if OPENGL
void *font = GLUT_BITMAP_HELVETICA_18;
void *font1 = GLUT_BITMAP_HELVETICA_18;
void *font2 = GLUT_BITMAP_HELVETICA_12;
void *font3 = GLUT_BITMAP_HELVETICA_10;
GLuint box, vesicle, tmd;
GLuint monomersideA, monomersideB, monomer1A, monomer2A, monomer1B, monomer2B;
GLuint LYS;
GLuint C2A1, C2A2, C2A3, C2B1, C2B2, C2B3;
GLuint C2A, C2B, C2AB;
GLuint C2A1close, C2A2close, C2A3close, C2B1close, C2B2close, C2B3close;
GLuint C2Aclose, C2Bclose, C2ABclose;

GLfloat redDiffuseMaterial[] = {1.0, 0.5, 0.4};
GLfloat greenDiffuseMaterial[] = {0.5, 1, 0.5};
GLfloat blueDiffuseMaterial[] = {0.6, 0.6, 1.0};
GLfloat yellowDiffuseMaterial[] = {1, 1, 0.2};

#if SAVPIC
static GLubyte* img_data=0;
#endif
#endif


//----- SYT chain -----
typedef struct CHAIN CHAIN;

struct CHAIN {	// chains
	int id, nseg, closed;	// segments are numbered cw
	double MuRotMemdt, Drdt;	// rotational mu*dt for the entire chain on membrane
	double Pclose;	// probability of chain close
	double rhead[3], rtail[3], rcenter[3], rcenterprev[3];	// positions of head, tail & center
	double rmax, rave;	// max radius of gyration, average radius
	
	double r[NSegMax][3];	// position of each segment center
	double rc[NSegMax][3];	// position in rcenter frame
	
	double rC2AB[NSegMax][4][3];	// centers of C2AB (diam=30 len=50 tube incl end caps)
	double drC2AB[NSegMax][4][3];	// displacement from syt center to C2AB beads
	
	double rlinker[NSegMax][3];	// position of c2a-linker joint if Vesicle=1
	double rtmd[NSegMax][3];	// position of tmd on vesicle if Vesicle=1
	
	double rc1[NSegMax][3];	// position of centers of left face
	double rc2[NSegMax][3];	// position of centers of right face
	
	double dr1[NSegMax][3];	// vector stretching between rcl of i & rcr of i+1
	double dr2[NSegMax][3];	// vector stretching between rcr of i & rcl of i-1
	
	double rsm[NSegMax][3];	// membrane interaction sites
	double fsm[NSegMax][3];
	
	double nx[NSegMax][3];	// local orientation vectors from i to i+1
	double ny[NSegMax][3];	// nx: tangential, ny: radial, nz = nx x ny in z direction
	double nz[NSegMax][3];
	
	double nnxt[NSegMax][3];	// center-to-center vector from i to i+1
	double nxa[NSegMax][3];	// average nx: direction pointing to the next segment
	double nya[NSegMax][3];	// average ny, nz
	double nza[NSegMax][3];
	
	double f[NSegMax][3];	// force of each segment
	double fb[NSegMax][3];	// bending force for short chains
	double ftmd[NSegMax][3];	// force on tmd if Vesicle=1
	double nrmtmd[NSegMax][3];	// normal direction of tmd

	double taum[NSegMax][3];	// torque from syt-membrane attraction
	double taust[NSegMax][3];	// stretching torque
	double taurp[NSegMax][3];	// repulsion torque
	double taubd[NSegMax][3];	// bending torque for self-orientation (both radial and pitch-wise)
	double tauts[NSegMax][3];	// torsion torque for self-orientation
	double tautmd[NSegMax][3];	// torque from tmd if Vesicle=1
	double tauorient[NSegMax][3];	// torque for self-orientation = taubd+tauts+taum
	double tautot[3];	// total torque on this chain
	
	double frep[NSegMax][4][3], fatt[NSegMax][3];
	double fext[NSegMax][3], ftot[3];	// external (from tmd, mnode & ves) and total force
	
	double Eb[NSegMax], Em[NSegMax], Ev[NSegMax];	// energy for each segment: bending, syt-mem, syt-ves
	double E[NSegMax];	// total energy of each segment
	double Poff[NSegMax];	// Poff of each segment
	
	//CHAIN *prev;
	CHAIN *next;	// next chain
};

CHAIN *chain=NULL;

#if !LGV_MC
	CHAIN *chain_backup=NULL;
#elif BD_2ND
	CHAIN *chain_try=NULL;
#endif

//----- target membrane -----

typedef struct MNODE MNODE;
typedef struct MFACE MFACE;

struct MNODE {
	int id[2];
	int nnbr, nbr[NnbrMax][2];	// neighbor list
	int nface, iface[NnbrMax];
	int nfaceplus, ifaceplus[NnbrMaxPlus];	// larger neighborhood
	int circbnd;
	
	double r[3], rnbr[NnbrMax][3];	// r: current location in lab
	double lnbr[NnbrMax], drnbr[NnbrMax][3], dirnbr[NnbrMax][3];
	double area, areanbr[NnbrMax], nrmnbr[NnbrMax][3];	// areanbr is for area-weighted normal
	double Av;	// Voronoi area, for Vallet-Levy's Laplace algorithm
	double C1, C2;	// principal curvatures
	double H, K;	// mean and Gaussian curvatures
	double LaplaceH;	// Laplacian of H
	double e1[3], e2[3], e3[3];	// principal directions and normal
	double f[3];	// tension, pressure, bending, total
	double drb[3];	// displacement due to bending
	double Eb, Ea, Ec, Es, E;	// energy: bending, areal, carbon-adhesion, syt, total
};

struct MFACE {
	int id, inode[3][2];
	double nrmsign;	// nrmsign=+1 if r12 x r13 is along outward normal, -1 if not
	double nrm[3], area, vol;
	double rc[3], r1[3], r2[3], r3[3];	// coord. of center of weight & 3 vertices
};

MNODE **mnode=NULL;
MFACE *mface=NULL;

#if !LGV_MC
	MNODE **mnode_backup=NULL;
	MFACE *mface_backup=NULL;
#elif BD_2ND
	MNODE **mnode_try=NULL;
	MFACE *mface_try=NULL;
#endif

//----- membrane zones -----

typedef struct MZONE MZONE;
typedef struct VZONE VZONE;

struct MZONE {
	int nnode, nface;
	int inode[NodeZoneMax][2], iface[NodeZoneMax];
};

struct VZONE {
	int nnode, nface;
	int nnodedown, nfacedown;
	int inode[NodeZoneMax], iface[NodeZoneMax];
	int inodedown[NodeZoneMax], ifacedown[NodeZoneMax];
};

MZONE **mzone=NULL;
VZONE **vzone=NULL;


//----- vesicle mesh -----
// initial icosahedron
#define Xico 0.525731112119133606
#define Zico 0.850650808352039932
static double Vico[12][3] = {
{-Xico, 0.0, Zico}, {Xico, 0.0, Zico}, {-Xico, 0.0, -Zico}, {Xico, 0.0, -Zico},
{0.0, Zico, Xico}, {0.0, Zico, -Xico}, {0.0, -Zico, Xico}, {0.0, -Zico, -Xico},
{Zico, Xico, 0.0}, {-Zico, Xico, 0.0}, {Zico, -Xico, 0.0}, {-Zico, -Xico, 0.0}
};
static int Fico[20][3] = {
{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

typedef struct VNODE VNODE;
typedef struct VEDGE VEDGE;
typedef struct VFACE VFACE;

struct VNODE {
	int id, nnbr, nbr[NnbrMax];
	int nface, iface[NnbrMax];
	int nfaceplus, ifaceplus[NnbrMaxPlus];	// larger neighborhood
	int nedge, iedge[NnbrMax];
	double r0[3], r[3], rc[3];	// r0: initial location, r: current location in lab, rc: current location in RC's frame
	double lnbr[NnbrMax], drnbr[NnbrMax][3];
	double area;
	double Av;	// Voronoi area, for Vallet-Levy's Laplace algorithm
	double C1, C2;	// principal curvatures
	double H, K;	// mean and Gaussian curvatures
	double LaplaceH;	// Laplacian of H
	double e1[3], e2[3], e3[3];	// principal directions and normal
	double f[3];	// tension, pressure, bending -> total
	double fext[3];	// external force (excluding Helfrich forces)
	double Eb, Ea, Es, E;	// energy: bending, areal, syt, total
};

struct VEDGE {
	int id, inode[2];
	double l0, l, dl, dir[3];	// dir is a normalized vector from node_0 to node_1
};

struct VFACE {
	int id, inode[3];
	double nrmsign;	// nrmsign=+1 if r12 x r13 is along outward normal, -1 if not
	double nrm[3], area, vol;
	double rclab[3], rcrc[3], r1[3], r2[3], r3[3];	// coord. of center & 3 vertices
};

VNODE *vnode=NULL;
VEDGE *vedge=NULL;
VFACE *vface=NULL;

#if LGV_MC && BD_2ND
VNODE *vnode_try=NULL;
VEDGE *vedge_try=NULL;
VFACE *vface_try=NULL;
#endif


//------- MT19937p ----------------

typedef struct mt19937p mt19937p;

struct mt19937p {
	unsigned long mt[624];
	int mti;
	unsigned long mag01[2];
};

mt19937p *mtseeds=NULL;


// -------------

#define max2(x,y) (x>y?x:y)
#define min2(x,y) (x<y?x:y)
#define max3(x,y,z) (max2(x,max2(y,z)))
#define min3(x,y,z) (min2(x,min2(y,z)))
#define sign(x) (x>=0.0?1.0:-1.0)
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))


#if RNG_64
	#define RND() (rand64())	// fast
	#define RNDMP(i) (rand64mp(i))
#else
	//#define RND() (rand32())	// slow
	//#define RNDMP(i) (rand32mp(i))
	#define RND() (genrand(&(mtseeds[0])))	// faster than rand64
	#define RNDMP(i) (genrand(&(mtseeds[i])))
#endif

void initview(void);
void reset(void);
void DoSimu1(void);
void DoSimu2(void);

void SimTime(void);
void terminate(void);

double ran2(long *idum);
double ran3(long *idum);

void sgenrand(unsigned long seed, mt19937p *config);
double genrand(mt19937p *config);

double Ei(double x);
double expint(int n, double x);

void NewSegPos(CHAIN *p, int tipid, int growpos, double r[3], double nx[3], double ny[3], double nz[3]);
void NewSegPosRing(CHAIN *p, int tipid, int growpos, double r[3], double nx[3], double ny[3], double nz[3]);
void GetSytStretch(CHAIN *p);

double getEfieldNear(double r, double cs, double radius);
double getEfieldFar(double r, double cs, double radius);

double getAmix_k(double dr_ki[3], double dr_kk2[3], double r2_ki, double r2_kk2, double r2_k2i, double cotbk, double cotak2);

double getSMdistance(MNODE **MN, double rsm[3]);
void getSMdistanceplus(MNODE **MN, double rsm[3], double *dist, double nrm[3], int mi[3], int mj[3]);

double getEsm(double rsm[3]);
double getEsv(double rsm[3]);
double getEmc(double z, double area);
double getEsmPF(double r, double cs, double a);

void testfatt(void);

void UpdateVesFaceCoord(VNODE *VN, VFACE *VF);
void CalcVesAreaVolume(VNODE *VN, VFACE *VF, int status);

void UpdateSytSitesAll(CHAIN *CH);
void UpdateSytStretchAll(CHAIN *CH);
void UpdateChain(CHAIN *CH);

void GetMemArea(MNODE **MN, MFACE *MF, int status);

void readcontour(void);
double getcontour(double r);

int ReadData(void);

double LambertW(double x);
double zbrent(double (*func)(double), double x1, double x2, double tol);
double zbrent2(double (*func)(double, double), double val, double x1, double x2, double tol);
double zbrent3(double (*func)(double, double, double), double v1, double v2, double x1, double x2, double tol);









