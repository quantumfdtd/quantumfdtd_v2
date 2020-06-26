/*

   mpisolve.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __mpisolve_h__
#define __mpisolve_h__

typedef complex<double> dcomp; 

extern int NUMX,NUM,STEPS,UPDATE,SNAPUPDATE,POTENTIAL,INITCONDTYPE,INITSYMMETRY,NF,SAVEWAVEFNCS,DUMPSNAPS;
extern double A,EPS,SIG,MASS,T,TC,SIGMA,XI,TOLERANCE,eB,Kx,SPINEXP;

extern int nodeID, numNodes, debug;

extern dcomp energyCollect;
extern dcomp normalizationCollect;
extern dcomp vInfinityCollect;
extern dcomp rRMS2Collect;    
extern dcomp xAvgCollect;           
extern dcomp yAvgCollect;          
extern dcomp zAvgCollect;        

/* debug values */
#define DEBUG_OFF		0
#define DEBUG_ON		1
#define DEBUG_FULL		2

/* message tags */
#define HELLO			1
#define DONE			2
#define SYNC_LEFT		3
#define SYNC_LEFT_MESSAGE	4
#define SYNC_RIGHT		5
#define SYNC_RIGHT_MESSAGE	6

// the main solve routines
void solveInitialize();
void solve();
void solveFinalize();
void evolve(int);
void findExcitedStates();

// boundary sync stuff
void syncBoundaries(dcomp*** wfnc);
void sendLeftBoundary(dcomp*** wfnc);
void sendRightBoundary(dcomp*** wfnc);
void receiveLeftBoundary();
void receiveRightBoundary();
void loadRightBoundaryFromBuffer(dcomp*** wfnc);
void loadLeftBoundaryFromBuffer(dcomp*** wfnc);

#endif /* __mpisolve_h__ */
