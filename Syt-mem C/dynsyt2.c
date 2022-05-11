// for Limit=1

int checkcontact(CHAIN *p, CHAIN *q, int i, int j)
{
	int z;
	double r[3], d2, dotp, dotp1, dotp2;
	
	z=0;
	
	vecsub(q->r[j], p->r[i], r);	// r is from p_i to q_j
	d2=norm2(r);
	if(d2>Distclose2) return 0;	// 1: should be close enough
	
	dotp = dotprod(p->nx[i], q->nx[j]);
	if(dotp<=0) return 0;	// 2: tangents are in the same direction
	
	normalize(r, r);	// r from p_i to q_j
	dotp1 = dotprod(p->nx[i], r);
	dotp2 = dotprod(q->nx[j], r);
	
	if(dotp1*dotp2<=0) return 0;	// 3: connection between p_i and q_j is not in line with the 2 nx's
	//else if(fabs(dotp1)<csdTHclose || fabs(dotp2)<csdTHclose) return 0;
	//else if(fabs(dotp1)<0.5 || fabs(dotp2)<0.5) return 0;
	
	z=CheckCloseOrientation(p,q,i,j);
	
	return z;
}


int CheckChainContact(CHAIN *p, CHAIN *q)
{
	int n1m1, n2m1;
	int z;
	
	n1m1 = p->nseg - 1;
	n2m1 = q->nseg - 1;
	
	z=0;
	if(n1m1+n2m1==0) z=checkcontact(p, q, 0, 0);	// z=1: head-head, only for 2 monomers
	if(z==0 && n2m1!=0) z=2*checkcontact(p, q, 0, n2m1);	// z=2: p_head-q_tail -> q will be new head
	if(z==0 && n1m1!=0) z=3*checkcontact(p, q, n1m1, 0);	// z=3: p_tail-q_head -> p will be new head
	if(z==0 && (n1m1+n2m1)!=0) z=4*checkcontact(p, q, n1m1, n2m1);	// z=4: tail-tail
	return z;
}



void Append2End(CHAIN *p, CHAIN *q)	// append q's segments to p's end
{
	int i, j, n1, n2, n1m1, n2m1;
	
	n1 = p->nseg;
	n2 = q->nseg;
	n1m1 = n1-1;
	n2m1 = n2-1;
	
	for(i=0; i<n2; i++) {	// copy q's segments to p's head
		j = n1+i;
		CopySegments(q, i, p, j);
	}
	
	p->nseg = n1+n2;
	veccopy(p->r[n1+n2-1], p->rtail);
}



void Insert2Front(CHAIN *p, CHAIN *q)	// insert q's segments to p's front
{
	int i, j, n1, n2, n1m1;
	
	n1 = p->nseg;
	n2 = q->nseg;
	n1m1 = n1-1;
	
	for(i=n1m1; i>=0; i--) {	// shift back p's segments by n2
		j = i+n2;
		CopySegments(p, i, p, j);	// copy p-i to p-j
	}
	for(i=0; i<n2; i++) {	// copy q's segments to p's head
		CopySegments(q, i, p, i);
	}
	
	p->nseg = n1+n2;
	veccopy(p->r[0], p->rhead);
}



void Join2Chains(CHAIN *p, CHAIN *q, int contact)
{
	int dir;
	double r[3], crs[3];
	
	//printf("%d + %d: %d\n", p->id, q->id, contact);
	if(contact==1) {	// head-head, only for 2 monomers
		vecsub(q->r[0], p->r[0], r);	// r is from p to q
		crossprod(r, p->ny[0], crs);	// crs = r x np
		if(crs[2]>=0) dir=1;	// p-q is cw, p will be the new head
		else dir=-1;	// p-q is ccw, q will be the new head
	}
	else if(contact==2) dir=-1;	// p_head-q_tail: q will be new head
	else if(contact==3) dir=1;	// p_tail-q_head: p will be new head
	
	if(dir==1) Append2End(p, q);	// p is the new head: append q segments to p's end, then remove q
	else Insert2Front(p, q);	// q is the new head: insert q segments to p's front, then remove q
}



void Annealing(void)
{
	int contact;
	CHAIN *p, *q, *s;
	
	p = chain;
	while(p) {
		if(p->closed == 0) {
			s = p;	// s is previous to q
			q = p->next;
			while(q) {
				if(q->closed == 0) {
					contact=CheckChainContact(p, q);
					if(contact!=0) {
						Join2Chains(p, q, contact);	// add q to p, then remove q in the following
						s->next = q->next;
						free(q);
						Nchain--;
						q = s->next;
						//CheckChainMono();
						//printf("----------\n");
					}
					else {
						s = q;
						q = q->next;
					}
				}	// end of q->closed == 0
				else q = q->next;
			}	// end of while(q)
		}	// end of p->closed == 0
		p = p->next;
	}
}




void DynamicsSyt2(void)
{
	Annealing();
	CloseChain();
	//ChainBreak();
	//RemoveOrphans();
	//RemoveCoils();
}
