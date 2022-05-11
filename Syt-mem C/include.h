
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if OPENMP
#include <omp.h>
#endif

#if OPENGL
//#include <GL/glut.h>
#include<GL/freeglut.h>
	#if SAVPIC
	#include <IL/il.h>
	#endif
#endif


#include "defs.h"
#include "readpara.c"
#include "func.c"

#include "rng.c"

#include "initmem.c"
#include "initves.c"

#include "initsyt1.c"
#include "initsyt2.c"

#include "dynsyt1.c"
#include "dynsyt2.c"

#include "memprop.c"
#include "vesprop.c"

#include "fmem.c"
#include "fves.c"
#include "fsyt.c"

#include "fsytmem.c"
#include "fsytves.c"
#include "fmemves.c"

#include "tsyt.c"

#include "flipmem.c"
#include "flipves.c"

#include "movmem.c"
#include "movves.c"
#include "movsyt.c"
#include "movtmd.c"

#include "energy.c"

#if !LGV_MC
	#include "mcaux.c"
	#include "mcenergy.c"
	#include "mc.c"
#endif

#if LGV_MC && BD_2ND
	#include "bd2nd.c"
#endif

#include "savedata.c"
#include "simu.c"

#if OPENGL
#include "glfunc.c"
#include "graphfunc.c"
#include "graphics.c"
#endif


