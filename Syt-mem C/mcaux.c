void InitChainBackup(void)
{
	CHAIN *p, *q;

	if(!chain_backup) {
		p = (CHAIN *)malloc(sizeof(CHAIN));
		p->next = NULL;
		chain_backup = p;
	}
	else {
		p = chain_backup;
		while(p->next) p = p->next;	// p is the end of chain
		q = (CHAIN *)malloc(sizeof(CHAIN));
		q->next = NULL;
		p->next = q;
	}
}



void InitMemBackup(void)
{
	int i;
	
	mnode_backup=(MNODE **)calloc(Nxp1, sizeof(MNODE *));
	for(i=0; i<=Nx; i++) {
		mnode_backup[i]=(MNODE *)calloc(Nyp1, sizeof(MNODE));
		if(!mnode_backup[i]) { printf("Cannot initialize mnode_backup\n"); exit(0); }
	}
	mface_backup=(MFACE *)calloc(NMface, sizeof(MFACE));
}




void PrintChain(CHAIN *p)
{
	int i;
	
	printf("Chains:\n");
	while(p) {
		printf("%d\n", p->nseg);
		printf("%d: (%.2f, %.2f) -- (%.2f, %.2f)\n", p->id, p->rhead[0], p->rhead[1], p->rtail[0], p->rtail[1]);
		for(i=0; i<(p->nseg); i++) {
			//printf("\t%d", i);
			printf("%d:\t%.2f, %.2f, %.2f\n", i, p->r[i][0], p->r[i][1], p->r[i][2]);
			//printf("%d:\t%.2f, %.2f, %.2f\n", i, p->rtmd[i][0], p->rtmd[i][1], p->rtmd[i][2]);
			//printf("%d:\t%.2f, %.2f, %.2f\n", i, p->nx[i][0], p->nx[i][1], p->nx[i][2]);
			//printf("%d:\t%.2f, %.2f, %.2f\n", i, p->ny[i][0], p->ny[i][1], p->ny[i][2]);
			//printf("%d:\t%.2f, %.2f, %.2f\n", i, p->nz[i][0], p->nz[i][1], p->nz[i][2]);
		}
		printf("--------------------------------\n");
		p = p->next;
	}
}



void CleanSytBackup(void)
{
	CHAIN *p, *q;
	
	if(chain_backup) {
		p = chain_backup;
		while(p) {
			q = p->next;
			free(p);
			chain_backup = p = q;
		}
		chain_backup=NULL;
	}
}



void CleanMemBackup(void)
{
	int i;
	
	if(mnode_backup) {
		for(i=0; i<=Nx; i++) free(mnode_backup[i]);
		free(mnode_backup);
	}
	mnode_backup=NULL;
	
	if(mface_backup) free(mface_backup);
	mface_backup=NULL;
}

