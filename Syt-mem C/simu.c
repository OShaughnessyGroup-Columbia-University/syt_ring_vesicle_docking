//---------- Langevin ----------

void SimuCore1(void)	// Langevin
{
	TrackFsyt=TrackFves=TrackFmem=0;
	
	if(SYT) {
		if(Limit==0) DynamicsSyt1();
		else DynamicsSyt2();
	}

	if(Membrane) UpdateMem(mnode, mface);
	if(Vesicle) UpdateVes(vnode, vface);
	
	if(SYT) ResetSytForceTorque(chain);

	if(Membrane) ResetMemForce(mnode);
	if(Vesicle) ResetVesForce(vnode);
	
	if(Membrane) ForceMem(mnode, mface);
	if(Vesicle) ForceVes(vnode);
	if(SYT) ForceSyt(chain);
	
	if(Membrane && SYT) FTSytMem(mnode, chain);

	if(Vesicle) {
		if(SYT) FTSytVes(chain, vnode, vface);
		if(Membrane) FTMemVes(vnode, vface, mnode, mface);
		TotalVesForce();	// get Fves
	}
	if(SYT) {
		AllSytForce(chain);
		TorqueSyt(chain);
		//AllEnergy();
		MotionSyt(chain);
	}

	if(Membrane) MotionMem(mnode, mface);
	if(Vesicle) MotionVes(vnode, vface, mnode, chain);
	
	if(Vesicle && TMD) MotionTMD(chain);
	
	t+=dt;
}



#if LGV_MC && BD_2ND

void SimuCore1BD2nd(void)	// Langevin w/ 2nd order BD
{
	TrackFsyt=TrackFves=TrackFmem=0;
	
	if(SYT) {
		if(Limit==0) DynamicsSyt1();
		else DynamicsSyt2();
	}

	if(Membrane) UpdateMem(mnode, mface);
	if(Vesicle) UpdateVes(vnode, vface);
	
	// --- 1st run: get forces & torques at old position ---
	
	if(SYT) ResetSytForceTorque(chain);	// reset before calculating all forces
	if(Membrane) ResetMemForce(mnode);
	if(Vesicle) ResetVesForce(vnode);

	if(SYT) ForceSyt(chain);
	if(Membrane) ForceMem(mnode, mface);
	if(Vesicle) ForceVes(vnode);
	
	if(Membrane && SYT) FTSytMem(mnode, chain);

	if(Vesicle) {
		if(SYT) FTSytVes(chain, vnode, vface);
		if(Membrane) FTMemVes(vnode, vface, mnode, mface);
		TotalVesForce();	// get Fves
	}
	
	if(SYT) {
		AllSytForce(chain);
		TorqueSyt(chain);
	}
	
	// make a copy of everything
	if(SYT) CopyChain(chain, chain_try);
	if(Membrane) {
		CopyMemNode(mnode, mnode_try);
		//CopyMemFace(mface, mface_try);
	}
	if(Vesicle) {
		CopyVesNode(vnode, vnode_try);
		//CopyVesEdge(vedge, vedge_try);
		//CopyVesFace(vface, vface_try);
		veccopy(Fves, Fves_old);	// save current Fves & Tauves
		veccopy(Tauves, Tauves_old);
	}
	
	// trial move on the copy
	if(SYT) MotionSyt2ndBD(chain_try);
	if(Membrane) MotionMem2ndBD(mnode_try, mface);	// use old mface
	if(Vesicle) MotionVes2ndBD(vnode_try, vface, mnode_try, chain_try);	// use old vface
	if(Vesicle && TMD) MotionTMD2ndBD(chain_try);

	// --- 2nd run: get forces & torques at new position ---

	if(SYT) ResetSytForceTorque(chain_try);	// reset before calculating all forces
	if(Membrane) ResetMemForce(mnode_try);
	if(Vesicle) ResetVesForce(vnode_try);

	if(SYT) ForceSyt(chain_try);
	if(Membrane) ForceMem(mnode_try, mface);	// use old mface
	if(Vesicle) ForceVes(vnode_try);
	
	if(Membrane && SYT) FTSytMem(mnode_try, chain_try);

	if(Vesicle) {
		if(SYT) FTSytVes(chain_try, vnode_try, vface);	// use old vface
		if(Membrane) FTMemVes(vnode_try, vface, mnode_try, mface);	// use old vface & mface
		TotalVesForce();	// get Fves
	}
	
	if(SYT) {
		AllSytForce(chain_try);
		TorqueSyt(chain_try);
	}
	
	// get average forces & torques from the 2 runs
	if(SYT) AveSytFT2bd(chain, chain_try);	// update chain with averaged forces & torques
	if(Membrane) AveMemFT2bd(mnode, mnode_try);
	if(Vesicle) AveVesFT2bd(vnode, vnode_try);
	
	// motion
	if(SYT) MotionSyt(chain);
	if(Membrane) MotionMem(mnode, mface);
	if(Vesicle) MotionVes(vnode, vface, mnode, chain);
	
	if(Vesicle && TMD) MotionTMD(chain);

	// ------

	t+=dt;
}

#endif



#if LGV_MC

void DoSimu1(void)	// Langevin
{
#if SAVDAT
	if(Nfile==0) {
		if(Membrane) UpdateMem(mnode, mface);
		if(Vesicle) UpdateVes(vnode, vface);
		SaveData();	// initial state
		if(Membrane) SaveContour();
		AllEnergy();
		SaveEnergy();
		Nfile++;
	}
#endif


#if !BD_2ND
	SimuCore1();
#else
	SimuCore1BD2nd();
#endif

	
#if SAVDAT
	SimCnt2++;	// use n=n1*1e6+cnt2
	if(SimCnt2>=1e6) {
		SimCnt1++;
		SimCnt2=0;
	}
	if(SimCnt1 >= SaveCntMax1 && SimCnt2 >= SaveCntMax2) {
		if(Nfile <= Nsave) {
			SaveData();
			if(Membrane) SaveContour();
			AllEnergy();
			SaveEnergy();
			Nfile++;
		}
		SimCnt1=SimCnt2=0;
	}
#endif
	
	if(Membrane && !Vesicle) {
		if(fabs(dHmem-dHmemprev)/max2(dHmem,eps)<1e-6 && t>5e-7) terminate();
	}
}

#endif


//---------- Monte-Carlo ----------

#if !LGV_MC

void SimuCore2(void)	// Monte-Carlo
{
	int i;
	double rnd;

	if(SYT) {
		if(Limit==0) DynamicsSyt1();
		else DynamicsSyt2();
	}
	if(Membrane) UpdateMem(mnode, mface);
	if(Vesicle) UpdateVes(vnode, vface);

	printf("Loop = %d\n", MCloopcnt);
	
	for(i=0; i<MCperloop; i++) {
		
		rnd = RND();

		if(rnd < MCp[0]) {
			MCmove1();
		}
		else if(rnd < MCp[1]) {
			MCmove2();
		}
		else {
			MCmove3();
		}
	}
}

#endif



void DoSimu2(void)	// Monte-Carlo
{
#if SAVDAT
	if(Nfile==0) {
		if(Membrane) UpdateMem(mnode, mface);
		if(Vesicle) UpdateVes(vnode, vface);
		SaveData();	// initial state
		if(SYT && Membrane) SaveContour();
		AllEnergy();
		SaveEnergy();
		Nfile++;
	}
#endif

#if !LGV_MC
	SimuCore2();
	MCloopcnt++;
#endif

#if SAVDAT
	SaveData();
	if(Membrane) SaveContour();
	AllEnergy();
	SaveEnergy();
	Nfile++;
#endif
	
	temperature *= 0.85;

	if(MCloopcnt>=MCloop) {
		MCloopcnt=0;
		MCcyclecnt++;
		temperature=0.01;
	}
}
