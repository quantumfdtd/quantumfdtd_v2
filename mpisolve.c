/*
 
 Parallelized Finite Difference Time Domain Schrodinger Eq Solver
 
 mpisolve.c
 
 Copyright (c) Michael Strickland
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <complex>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

#include "mpi.h" 
#include "mpisolve.h"
#include "grid.h"
#include "initialconditions.h"
#include "potential.h"
#include "outputroutines.h"
#include "paramreader.h"

// these global vars are initialized from parameters file
// defaults set here are overridden by that file
int    NUMX=20,NUM=20,STEPS=40000,UPDATE=100,SNAPUPDATE=1000;
int    POTENTIAL=0,INITCONDTYPE=0,INITSYMMETRY=0,NF=2,SAVEWAVEFNCS=0,DUMPSNAPS=0;
double  A=0.05,EPS=0.001,SIG=0.06,MASS=1.0,T=1.0,TC=0.192,SIGMA=0.223,XI=0.0,TOLERANCE=1.e-10,eB=0,Kx=0,SPINEXP=0;

// mpi vars
int nodeID,numNodes;

// files
fstream debug_out;

// debug flag; options are DEBUG_{OFF,ON,FULL}
int debug = DEBUG_OFF;

// used for MPI non-blocking sends
double *leftSendBuffer,*rightSendBuffer;
MPI_Status leftSendStatus,rightSendStatus;
MPI_Status leftSendMessageStatus,rightSendMessageStatus;
MPI_Request leftSend,rightSend;
MPI_Request leftMessageSend,rightMessageSend;

// used for MPI non-blocking receives
double *leftReceiveBuffer,*rightReceiveBuffer;
MPI_Status leftReceiveStatus,rightReceiveStatus;
MPI_Status leftReceiveMessageStatus,rightReceiveMessageStatus;
MPI_Request leftReceive,rightReceive;
MPI_Request leftMessageReceive,rightMessageReceive;

// used for worker only intercommunication
MPI_Comm workers_comm;
MPI_Group workers_group;

// variables which will be loaded and reduced across nodes
dcomp energy=0;			// the local node energy
dcomp energyCollect=0;		// the total energy
dcomp normalization=0;		// the local node normalization squared
dcomp normalizationCollect=0;  	// the total normalization squared
dcomp vInfinity=0;		// the local node expectation value of v_infty
dcomp vInfinityCollect=0;      	// the total expectation value of v_infty
dcomp rRMS2=0;                 	// the local node <r^2>  
dcomp rRMS2Collect=0;		// the total <r^2>  
dcomp xAvg=0;			// the local <x>
dcomp xAvgCollect=0;		// the total <x>
dcomp yAvg=0;			// the local <y>
dcomp yAvgCollect=0;		// the total <y>
dcomp zAvg=0;			// the local <z>
dcomp zAvgCollect=0;		// the total <z>

// ground state energy and final time saved in global var after convergence
dcomp	EGrnd, timef;

// counter used for recording snapshots
int snapcnt = -1;

int main( int argc, char *argv[] ) 
{ 
    int done=0,checksum=0;
    char message[64]; 
    char fname[32];
    int exclude[1];
    struct tms starttime,endtime;

    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&numNodes); 
    MPI_Comm_rank(MPI_COMM_WORLD,&nodeID); 
    MPI_Status status;
    MPI_Group all_group;
	
    // setup group consisting of computational nodes only for internal communications
    MPI_Comm_group(MPI_COMM_WORLD, &all_group);
    exclude[0]=0;
    MPI_Group_excl(all_group, 1, exclude, &workers_group);
    MPI_Comm_create(MPI_COMM_WORLD, workers_group, &workers_comm);
	
    if (debug) {
		sprintf(fname,"debug/debug_%d.txt",nodeID);
		debug_out.open(fname, ios::out);
		debug_out << "==> Node " << nodeID << " is ready" << endl; 
    }
	
    // node 0 is the master
    if (nodeID == 0) {
		times(&starttime); // load start time into starttime structure
		print_line();
		cout << "Parameters from file" << endl;
		print_line();
		readParametersFromFile((char *)"params.txt",1);
		if (argc>1) {
			print_line();
			cout << "Parameters from commandline" << endl;
			print_line();
			readParametersFromCommandLine(argc,argv,1);
		}
    }
    else {
		readParametersFromFile((char *)"params.txt",0);
		readParametersFromCommandLine(argc,argv,0);
    }
	
    if (NUM%(numNodes-1)!=0) {
		if (nodeID==0)
			print_line();
        	cout << "ERROR: Unable to partition lattice ... exiting" << endl;
        	print_line();
		if (debug) {
			debug_out << "==> Goodbye from node " << nodeID << endl; 
			debug_out.close();
		}
        	MPI_Finalize(); 
		exit(0);
    } else {
		NUMX = NUM/(numNodes-1);	
    }
	
    if (nodeID == 0) {
		
		//	
		// master node
		//	
		
		for (int node=1;node<numNodes;node++) {
			sprintf(message,"Hello to node %d",node); 
			MPI_Send(message, strlen(message), MPI_CHAR, node, HELLO, MPI_COMM_WORLD); 
		}
		
		// master loops and waits for children to report that they are ready to start
		checksum=0; 
		do {
			MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, HELLO, MPI_COMM_WORLD, &status); 
			checksum += done;
			if (debug) debug_out << "Received: hello from computational node" << endl;
		} while( checksum < numNodes-1 ); 
		
		// cluster is ready
		if (debug) debug_out << "==> Cluster ready" << endl;
		
		// Currently the master process does nothing.  
		// It simply starts and waits for the others
		
		// master loops and waits for children to report that they are done
		checksum=0;
		do {
			MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status); 
			checksum += done;
			if (debug) debug_out << "Received: checkout from computational node" << endl;
			sleep(1.0); // sleep 0.1 seconds between checks in order to reduce CPU usage of master
		} while( checksum < numNodes-1 ); 
		
		times(&endtime); // load end time into endtime structure
		cout << "==> User time: " << (endtime.tms_utime - starttime.tms_utime)/((double)sysconf(_SC_CLK_TCK)) << " seconds"<< endl;
		cout << "==> System time: " << (endtime.tms_stime - starttime.tms_stime)/((double)sysconf(_SC_CLK_TCK)) << " seconds" << endl;
		print_line();
		cout << "Done." << endl;
		print_line();
		
    } else {
		
		//	
		// computational node
		//	
		
		// set initial conditions and get ready for computation
		solveInitialize();
		
		// blocking receive to wait for master to fire up and say hello
		MPI_Recv(message, 20, MPI_CHAR, 0, HELLO, MPI_COMM_WORLD, &status); 
		if (debug == DEBUG_FULL) debug_out << "==> Received : "<< message << endl;
		
		// the master said hello so now let's get to work
		solve();
		
		// done with main computation, now do any analysis required
		solveFinalize();
		
		// send message to master that we are done
		done = 1;
		MPI_Send( &done, 1, MPI_INT, 0, DONE, MPI_COMM_WORLD );
		
		cout.flush();
		
    }
	
    if (debug) { 
		debug_out << "==> Goodbye from node " << nodeID << endl; 
		debug_out.close();
    }
	
    MPI_Group_free(&workers_group);
    MPI_Finalize(); 
    return 0; 
} 

// solve initialize
void solveInitialize() {
	
	// allocate memory
	allocateMemory();
	
	// load the potential
	if (nodeID==1) {
		print_line();
		cout << "==> Loading Potential Arrays" << endl;
		flush(cout);
	}
	loadPotentialArrays();
	
	if (nodeID==1) print_line();
	
	// set initial conditions
	setInitialConditions(nodeID+1);
	
	// output some summary information
	if (nodeID==1) { 
		print_line();
      		cout << "==> Number of computational nodes : " << numNodes-1 << endl; 
      		cout << "==> NUMX : " << NUMX << endl; 
      		print_line();
      		cout << "Spatial Step Size (A): " << A << endl;
      		cout << "Temporal Step Size (EPS): " << EPS << endl;
      		cout << "Standard Deviation of initial wavefunction noise (SIG): " << SIG << endl;
      		cout << "Box Size (A*" << NUM << "): " << A*NUM << endl;
      		print_line();
      		cout.width(dwidth); cout << "Time";
      		cout.width(dwidth); cout << "Energy";
      		cout.width(dwidth); cout << "Binding Energy";
      		cout.width(dwidth); cout << "r_RMS";   
      		cout << endl;
      		print_line();
	}
}

// reduce observables across nodes to first worker node
void computeObservables(dcomp*** wfnc) {
	
	// sum energy across nodes
	double energy_re=0.,energy_im=0.;
	double energy_re_collect=0.,energy_im_collect=0.;
	energy = wfncEnergy(wfnc);
	energy_re = real(energy);
	energy_im = imag(energy);
	MPI_Reduce(&energy_re,&energy_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&energy_im,&energy_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	energyCollect = dcomp(energy_re_collect,energy_im_collect);
	
	// sum normalization squared across nodes
	double normalization_re=0.,normalization_im=0.;
	double normalization_re_collect=0.,normalization_im_collect=0.;
	normalization = wfncNorm2(wfnc);
	normalization_re = real(normalization);
	normalization_im = imag(normalization);
	MPI_Reduce(&normalization_re,&normalization_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&normalization_im,&normalization_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	normalizationCollect = dcomp(normalization_re_collect,normalization_im_collect);
	
	// sum expectation across nodes
	double vInfinity_re=0.,vInfinity_im=0.;
	double vInfinity_re_collect=0.,vInfinity_im_collect=0.;
	vInfinity = vInfinityExpectationValue(wfnc);
	vInfinity_re = real(vInfinity);
	vInfinity_im = imag(vInfinity);	
	MPI_Reduce(&vInfinity_re,&vInfinity_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&vInfinity_im,&vInfinity_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	vInfinityCollect = dcomp(vInfinity_re_collect,vInfinity_im_collect);
	
	// sum r-squared across nodes
	double rRMS2_re=0.,rRMS2_im=0.;
	double rRMS2_re_collect=0.,rRMS2_im_collect=0.;
	rRMS2 = r2ExpectationValue(wfnc);
	rRMS2_re = real(rRMS2);
	rRMS2_im = imag(rRMS2);
	MPI_Reduce(&rRMS2_re,&rRMS2_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&rRMS2_im,&rRMS2_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	rRMS2Collect = dcomp(rRMS2_re_collect,rRMS2_im_collect);

	// sum x across nodes
	double x_re=0.,x_im=0.;
	double x_re_collect=0.,x_im_collect=0.;
	xAvg = xExpectationValue(wfnc);
	x_re = real(xAvg);
	x_im = imag(xAvg);
	MPI_Reduce(&x_re,&x_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&x_im,&x_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	xAvgCollect = A*dcomp(x_re_collect,x_im_collect);

	// sum y across nodes
	double y_re=0.,y_im=0.;
	double y_re_collect=0.,y_im_collect=0.;
	yAvg = yExpectationValue(wfnc);
	y_re = real(yAvg);
	y_im = imag(yAvg);
	MPI_Reduce(&y_re,&y_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&y_im,&y_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	yAvgCollect = A*dcomp(y_re_collect,y_im_collect);

	// sum z across nodes
	double z_re=0.,z_im=0.;
	double z_re_collect=0.,z_im_collect=0.;
	zAvg = zExpectationValue(wfnc);
	z_re = real(zAvg);
	z_im = imag(zAvg);
	MPI_Reduce(&z_re,&z_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	MPI_Reduce(&z_im,&z_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
	zAvgCollect = A*dcomp(z_re_collect,z_im_collect);
	
}

// main computational solve routine
void solve() {
	
	dcomp energytot,lastenergy = 1.0e10;
	int step=0,done=1;
	char label[64]; 
	
	leftSendBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
	rightSendBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
	leftReceiveBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
	rightReceiveBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
	
	// send message to master that we are starting
	MPI_Send( &done, 1, MPI_INT, 0, HELLO, MPI_COMM_WORLD );

	// evolve lattice in steps of UPDATE and output data along the way
	do {
		
		// sync boundaries
		syncBoundaries(w);
		
		// reduce observables across nodes
		computeObservables(w);
		
		// output 2d snapshots of the wavefunction for inspection
		// and check convergence of ground state energy
		if (step%SNAPUPDATE==0) {
			// broadcast observables
			MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			MPI_Bcast(&energyCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
			// force symmetry
			symmetrizeWavefunction();
            		// normalize wavefunction
            		normalizeWavefunction(w);
            		// check convergence and break if tolerance is achieved
			// otherwise, record snapshot for use in excited state 
			// computation and keep going
			energytot =  energyCollect/normalizationCollect;
			if (abs(energytot-lastenergy)<TOLERANCE) {
				if (nodeID==1) outputMeasurements(step*EPS);
				break;
			} else {
				lastenergy = energytot;
				// record and output snapshot
				if (step!=STEPS) { 
					snapcnt = (snapcnt+1)%2; // assume only two snapshots for now, so cycle
					recordSnapshot(w,snapcnt);
					if (DUMPSNAPS)  {
						sprintf(label,"%d_%d",nodeID,step); 
						outputSnapshot(w,label);
					}
				}
			}
		}
		
		if (nodeID==1) outputMeasurements(step*EPS);
		
		if (step<STEPS) evolve(UPDATE);
		step += UPDATE;
		
	} while (step<=STEPS);
	
	// save ground state energy and tau_f in global variables
	EGrnd = energytot;
	timef = step*EPS;
	
	if (nodeID==1) outputSummaryData("Ground State");
	
	if (debug) {
		debug_out << "==> Unnormalized Energy : " << energy << endl;
		debug_out << "==> Normalization2 : " << normalization << endl;
	}
	
	free(rightSendBuffer);
	free(leftSendBuffer);
	free(rightReceiveBuffer);
	free(leftReceiveBuffer);
	
}

// solve finalize
void solveFinalize() {
	
	// this routine currently computes the first excited state energy and wavefunction
	findExcitedStates();
	
}

// evolves solution nsteps
void evolve(int nsteps) {
	
	int leftTest,rightTest;
	
	for (int i=1;i<=nsteps;i++) {
		
		// receive boundary sync
		leftTest=1; 
		rightTest=1;
		if (i>1) {
			if (nodeID+1 < numNodes) { 
				receiveRightBoundary(); 
				rightTest = 0;
            		}
			if (nodeID-1 >= 1) {
				receiveLeftBoundary();
				leftTest = 0;
            		}
			while (!leftTest || !rightTest) {
				if (!rightTest) {
					MPI_Test(&rightReceive,&rightTest,&rightReceiveStatus);
					if (rightTest) loadRightBoundaryFromBuffer(w);
				}
				if (!leftTest) {
					MPI_Test(&leftReceive,&leftTest,&leftReceiveStatus);
					if (leftTest) loadLeftBoundaryFromBuffer(w);
				}
			}
		}
		
		// first update boundary so that the send can be happening while we update the interior
		updateBoundaries(EPS);
		
		// wait and make sure send buffers are ready
		if (i>1) {
			if (nodeID+1 < numNodes) MPI_Wait(&rightSend,&rightSendStatus);
			if (nodeID-1 >= 1) MPI_Wait(&leftSend,&leftSendStatus);
		}
		
		// send boundary sync 
		leftTest=1; 
		rightTest=1;
		if (i==1) {
			if (nodeID+1 < numNodes) sendRightBoundary(W);
			if (nodeID-1 >= 1) sendLeftBoundary(W);
		}
		else if (i!=nsteps && i>1) {
			if (nodeID+1 < numNodes) rightTest=0;
			if (nodeID-1 >= 1) leftTest=0;
			while (!leftTest || !rightTest) {
				if (!rightTest) {
					MPI_Test(&rightSend,&rightTest,&rightSendStatus);
					if (rightTest) sendRightBoundary(W);
				}
				if (!leftTest) {
					MPI_Test(&leftSend,&leftTest,&leftSendStatus);
					if (leftTest) sendLeftBoundary(W);
				}
			}
		}
		
		// evolve interior of the grid forward in time by EPS storing updated fields in capital vars
		updateInterior(EPS);
		
		// copy fields from capital vars (updated) down to lowercase vars (current)
		copyDown();
		
	}
	
}

/*---------------------------------------------------------------------------*/
/* Compute Wavefunction Overlap                                              */
/*---------------------------------------------------------------------------*/

dcomp computeOverlap(dcomp*** wfnc1, dcomp*** wfnc2){
        dcomp overlap=0,overlapCollect=0;
        for (int sx=1;sx<=NUMX;sx++)
          for (int sy=1;sy<=NUM;sy++)
            for (int sz=1; sz<=NUM;sz++)
              overlap += conj(wfnc1[sx][sy][sz])*wfnc2[sx][sy][sz];

        MPI_Reduce(&overlap,&overlapCollect,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,workers_comm);
        MPI_Bcast(&overlapCollect,1,MPI_DOUBLE_COMPLEX,0,workers_comm);
        return overlapCollect;
}

/*---------------------------------------------------------------------------*/
/* Find excited states                                                       */
/*---------------------------------------------------------------------------*/

void findExcitedStates() {
	
	char label[64]; 
	dcomp ener,dEtau;

	/*---------------------------------------------------------------------------*/
	/* Find first excited state                                                  */
	/*---------------------------------------------------------------------------*/
	
	int snap = snapcnt;
	
	// compute overlap
	dcomp overlap = computeOverlap(w,wstore[snap]);

	// subtract overlap
	for (int sx=0;sx<NUMX+2;sx++) 
		for (int sy=0;sy<NUM+2;sy++)
       	    		for (int sz=0; sz<NUM+2;sz++) 
				W[sx][sy][sz] = wstore[snap][sx][sy][sz] - overlap*w[sx][sy][sz];
	
	// compute observables
	computeObservables(W);
	
	// normalize
	MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
	normalizeWavefunction(W);

	if (DUMPSNAPS) {	
		// output snapshot of excited states
		sprintf(label,"first_excited_%d",nodeID); 
		outputSnapshot(W,label);
	}
	
	if (nodeID==1) {

		outputSummaryData("First Excited State");

		// consistency check and warning if fail
		ener = energyCollect/normalizationCollect;
		dEtau = (ener-EGrnd)*timef;     
		if (real(dEtau)<4.) {
			print_line();
			cout << "==> WARNING: states nearly degenerate, tau_f too small!" << endl;
			print_line();
		}
		
	}
	EGrnd = ener; // redfine EGrnd to hold energy of first excited state for comparison below

	/*---------------------------------------------------------------------------*/
	/* Find second excited state                                                 */
	/*---------------------------------------------------------------------------*/
	
	int save = snap;
	snap = (snapcnt+1)%2;
	
	// compute overlap
	overlap = computeOverlap(w,wstore[snap]);
	dcomp overlap2 = computeOverlap(W,wstore[snap]);

	// subtract overlap
	for (int sx=0;sx<NUMX+2;sx++) 
		for (int sy=0;sy<NUM+2;sy++)
       	    		for (int sz=0; sz<NUM+2;sz++) 
				wstore[save][sx][sy][sz] = wstore[snap][sx][sy][sz] - overlap*w[sx][sy][sz] - overlap2*W[sx][sy][sz];
	
	// compute observables
	computeObservables(wstore[save]);
	
	// normalize
	MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
	normalizeWavefunction(wstore[save]);

	if (DUMPSNAPS) {	
		// output snapshot of excited states
		sprintf(label,"second_excited_%d",nodeID); 
		outputSnapshot(wstore[save],label);
	}
	
	if (nodeID==1) {

		outputSummaryData("Second Excited State");

		// consistency check and warning if fail
		ener = energyCollect/normalizationCollect;
		dEtau = (ener-EGrnd)*timef;     
		if (real(dEtau)<4.) {
			print_line();
			cout << "==> WARNING: states nearly degenerate, tau_f too small!" << endl;
			print_line();
		}
		
	}

	/*---------------------------------------------------------------------------*/
	/* Save extracted wavefunctions                                              */
	/*---------------------------------------------------------------------------*/
	
	if (SAVEWAVEFNCS) {
		// save 3d wavefunction for extracted states
		sprintf(label,"0_%d",nodeID); 
		outputWavefunction(w,label);
		// For now only output ground state, NFS causes delays for writes and this eats time
		/*
		sprintf(label,"1_%d",nodeID); 
		outputWavefunction(W,label);
		sprintf(label,"2_%d",nodeID); 
		outputWavefunction(wstore[save],label);
		*/
	}
	
	return;
}

/*---------------------------------------------------------------------------*/
/* Boundary sync routines                                                    */
/*---------------------------------------------------------------------------*/

void syncBoundaries(dcomp ***wfnc) {
	
	// initiate sends and receives
	if (nodeID+1 < numNodes) { 
		sendRightBoundary(wfnc);
		receiveRightBoundary(); 
	}
	if (nodeID-1 >= 1) {
		sendLeftBoundary(wfnc);
		receiveLeftBoundary(); 
	}
	
	// now wait for communications to complete and sync wfnc when they do
	if (nodeID+1 < numNodes) { 
		MPI_Wait(&rightReceive,&rightReceiveStatus);
		loadRightBoundaryFromBuffer(wfnc);
		MPI_Wait(&rightSend,&rightSendStatus);
	}
	if (nodeID-1 >= 1) {
		MPI_Wait(&leftReceive,&leftReceiveStatus);
		loadLeftBoundaryFromBuffer(wfnc);
		MPI_Wait(&leftSend,&leftSendStatus);
	}
}

void sendRightBoundary(dcomp*** wfnc) {
	char message[64]; 
	if (debug == DEBUG_FULL) {
		sprintf(message,"%d -> %d",nodeID,nodeID+1); 
		debug_out << "==> Sending : " << message << endl;
		MPI_Isend(message, strlen(message), MPI_CHAR, nodeID+1, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &rightMessageSend); 
	}
	for (int sy=0;sy<NUM+2;sy++)
		for (int sz=0;sz<NUM+2;sz++) {
			rightSendBuffer[sy*(NUM+2)+sz%(NUM+2)] = real(wfnc[NUMX][sy][sz]);
			rightSendBuffer[sy*(NUM+2)+sz%(NUM+2) + (NUM+2)*(NUM+2)] = imag(wfnc[NUMX][sy][sz]);
		}
	MPI_Isend(rightSendBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID+1, SYNC_RIGHT, MPI_COMM_WORLD, &rightSend); 
}

void sendLeftBoundary(dcomp*** wfnc) {
	char message[64]; 
	if (debug == DEBUG_FULL) {
		sprintf(message,"%d -> %d",nodeID,nodeID-1); 
		debug_out << "==> Sending : " << message << endl;
		MPI_Isend(message, strlen(message), MPI_CHAR, nodeID-1, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &leftMessageSend); 
	}
	for (int sy=0;sy<NUM+2;sy++)
		for (int sz=0;sz<NUM+2;sz++) { 
			leftSendBuffer[sy*(NUM+2)+sz%(NUM+2)] = real(wfnc[1][sy][sz]);
			leftSendBuffer[sy*(NUM+2)+sz%(NUM+2) + (NUM+2)*(NUM+2)] = imag(wfnc[1][sy][sz]);
		}
	MPI_Isend(leftSendBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID-1, SYNC_LEFT, MPI_COMM_WORLD, &leftSend);
}

void receiveRightBoundary() {
	char message[64]; 
	if (debug == DEBUG_FULL) {
		MPI_Irecv(message, 255, MPI_CHAR, nodeID+1, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &rightMessageReceive); 
		debug_out << "==> Received : " << message << endl;
	}
	MPI_Irecv(rightReceiveBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID+1, SYNC_LEFT, MPI_COMM_WORLD, &rightReceive); 
}

inline void loadRightBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array right boundary
	for (int sy=0;sy<NUM+2;sy++)
		for (int sz=0;sz<NUM+2;sz++)
			wfnc[NUMX+1][sy][sz] = dcomp(rightReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)],rightReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)+(NUM+2)*(NUM+2)]);
}

void receiveLeftBoundary() {
	char message[64];
	if (debug == DEBUG_FULL) {
		MPI_Irecv(message, 255, MPI_CHAR, nodeID-1, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &leftMessageReceive); 
		debug_out << "==> Received : "<< message << endl;
	}
	MPI_Irecv(leftReceiveBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID-1, SYNC_RIGHT, MPI_COMM_WORLD, &leftReceive); 
}

inline void loadLeftBoundaryFromBuffer(dcomp ***wfnc) {
	// update w array left boundary
	for (int sy=0;sy<NUM+2;sy++)
		for (int sz=0;sz<NUM+2;sz++)
			wfnc[0][sy][sz] = dcomp(leftReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)],leftReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)+(NUM+2)*(NUM+2)]);
}
