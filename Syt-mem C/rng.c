
void getmseed(void)
{
#if OPENGL && SAVPIC
	mseed=-2;
#else
	struct tm *newtime;
	time_t long_time;

	time(&long_time);	// RNG seed
	newtime=localtime(&long_time);
	mseed=-(int) fabs((newtime->tm_hour+1)*(2*newtime->tm_min+3)*(111*newtime->tm_sec+5)+
		191*(newtime->tm_min)+901*(newtime->tm_sec)-1);
		
	//############
	//mseed=-1;
	//############
#endif
}


void initseedarray(void)	// after obtaining Nthread
{
	int i;
	mseedarray=calloc(Nthread, sizeof(long));	// for parallel processing
	mtseeds=calloc(Nthread, sizeof(mt19937p));	// for MT19937p
	
#if RNG_64
	U64MP=calloc(Nthread, sizeof(Ullong));
	randSeed64();	// uses mseed
#else
	/*
	U32MP=calloc(Nthread, sizeof(Uint));
	V32MP=calloc(Nthread, sizeof(Uint));
	W132MP=calloc(Nthread, sizeof(Uint));
	W232MP=calloc(Nthread, sizeof(Uint));
	randSeed32();	// uses mseed
	*/
#endif
	
	for(i=0; i<Nthread; i++) {
		mseedarray[i] = (long)(1000000*ran3(&mseed));
		sgenrand(mseedarray[i], &(mtseeds[i]));	// MT19937p
	#if RNG_64
		U64MP[i] = 4101842887655102017LL;
		randSeed64mp(i);	// uses mseedarray[i]
	#else
		/*
		V32MP[i] = 2244614371U;
		W132MP[i] = 521288629U;
		W232MP[i] = 362436069U;
		randSeed32mp(i);	// uses mseedarray[i]
		*/
	#endif
	}
}



void testRNG(void)
{
	int i, j;
	double r;
#if OPENMP
	int tid;
#endif
	
	time(&start_time);
	
	for(i=0; i<10; i++) {
		printf("%d", i+1);
		
#if OPENMP
#pragma omp parallel for private(j,r)
#endif

		for(j=0; j<1000000000; j++) {
		#if OPENMP
			tid = omp_get_thread_num();
			r=RNDMP(tid);
		#else
			r=RND();
		#endif
			if(r<0 || r>1) {
				printf("\nrnd=%.4g\n", r);
				exit(0);
			}
		}
		printf("\n");
	}
	
	SimTime();
	exit(0);
}

