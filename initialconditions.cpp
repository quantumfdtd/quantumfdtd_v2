/*

   initialconditions.cpp

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
#include <sys/stat.h>
#include <sys/time.h>
#include <string>
#include <vector>
#include <sstream>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "initialconditions.h"
#include "random.h"

// initializes the variables 
void setInitialConditions(int seedMult)
{
  	double sig=SIG; // standard deviation
	int sx,sy,sz,tx,ty,tz;
	double dx,dy,dz,r,costheta,cosphi,md,temp,temp2;
	fstream input;
	char fname[32];
	string line;
	vector<string> lines;
	int inputLatticeSize,inputLatticeSizeX,stridein=1,strideout=1,linenumber;
	
	
	//cout << "==> Initializing variables\n";

	srand(((unsigned int)time(0))*seedMult);
	//srand(19691123);

	switch (INITCONDTYPE) { 
	  case 0: 
		// read from file
	    sprintf(fname,"data/wavefunction_0_%d.dat",nodeID);
	    input.open(fname, ios::in);
		if (nodeID==1) cout << "==> Initial wavefunction : From file" << endl;
		if (!input) {
			cout << "==> Error : Unable to open wavefunction file " << fname << ". Using random Gaussian instead." << endl;
			for (sx=0;sx<NUMX+2;sx++)
				for (sy=0;sy<NUM+2;sy++)
					for (sz=0; sz<NUM+2;sz++)
						w[sx][sy][sz] = randGauss(sig);
		}
		while( getline( input, line ) ) lines.push_back( line ) ;
		inputLatticeSize = round(pow((numNodes-1)*lines.size(),1/3.));
		inputLatticeSizeX = inputLatticeSize/(numNodes-1);
		if (NUMX > inputLatticeSizeX) strideout = NUMX/inputLatticeSizeX;
		if (NUMX < inputLatticeSizeX) stridein = inputLatticeSizeX/NUMX;
		for (sx=1;sx<=NUMX;sx++)
			for (sy=1;sy<=NUM;sy++)
				for (sz=1; sz<=NUM;sz++) {
					if (debug && nodeID==1) cout << "Mark : " << sx << ", " << sy << ", " << sz << endl;
			        if (strideout==1 && strideout==1) {
						linenumber  = (sx-1)*inputLatticeSize*inputLatticeSize + (sy-1)*inputLatticeSize + (sz-1);
					}					
			        if (strideout>1) { 
						// If input wavefunction has lower resolution, spread it out
						tx = ceil(sx/((double)strideout));
						ty = ceil(sy/((double)strideout));
						tz = ceil(sz/((double)strideout));
						linenumber  = (tx-1)*inputLatticeSize*inputLatticeSize + (ty-1)*inputLatticeSize + (tz-1);
						if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}
			        if (stridein>1) {
						// If input wavefunction has higher resolution, sample it						
						tx = sx*stridein;
						ty = sy*stridein;
						tz = sz*stridein;
						linenumber  = (tx-1)*inputLatticeSize*inputLatticeSize + (ty-1)*inputLatticeSize + (tz-1);
						if (debug && nodeID==1) cout << "Respond : " << tx << ", " << ty << ", " << tz << endl;
					}					
					line = lines.at(linenumber);
					int space_index = line.find_last_of("\t");
					string linemod = line.substr(0,space_index-1);	
					int space_index2 = linemod.find_last_of("\t");
					std::istringstream stream;
					stream.str(line.substr(space_index,line.length()-space_index));
					stream >> temp; // imaginary part
					std::istringstream stream2;
					stream2.str(line.substr(space_index2,space_index-space_index2));
					stream2 >> temp2; // real part
					if (debug && nodeID==1) cout << "line #" << linenumber << " : " << line << endl;
					if (debug && nodeID==1) cout << "real : " << temp2 << endl;
					if (debug && nodeID==1) cout << "imag : " << temp << endl;
					w[sx][sy][sz] = dcomp(temp2,temp);
				}
		input.close();
		break;
	  case 1:
		// random
		if (nodeID==1) cout << "==> Initial wavefunction : Random" << endl;
        	for (sx=0;sx<NUMX+2;sx++)
          	  for (sy=0;sy<NUM+2;sy++)
            	    for (sz=0; sz<NUM+2;sz++)
                      w[sx][sy][sz] = dcomp(randGauss(sig),0.);
		break;
	  case 2:
		// coulomb like
		if (nodeID==1) cout << "==> Initial wavefunction : Coulomb" << endl;
        	for (sx=0;sx<=NUM+1;sx++)
          	  for (sy=0;sy<=NUM+1;sy++)
            	    for (sz=0; sz<=NUM+1;sz++) {
		      // coordinate system is centered in simulation volume
		      dx = sx - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
		      dy = sy - ((double)NUM+1.)/2.;
		      dz = sz - ((double)NUM+1.)/2.;
		      r = A*sqrt(dx*dx+dy*dy+dz*dz);
		      costheta = A*dz/r;
		      cosphi = A*dx/r;
                      w[sx][sy][sz]  = exp(-MASS*r);				// n=1
                      w[sx][sy][sz] += (2-MASS*r)*exp(-MASS*r/2);		// n=2,l=0
                      w[sx][sy][sz] += MASS*r*exp(-MASS*r/2)*costheta;		// n=2,l=1,m=0
                      w[sx][sy][sz] += MASS*r*exp(-MASS*r/2)*sqrt(1-costheta*costheta)*cosphi; //n=2,l=1,m=+-1
	    	    }
		break;
	  case 3:
		// constant
		if (nodeID==1) cout << "==> Initial wavefunction : Constant" << endl;
        	for (sx=0;sx<=NUMX+1;sx++)
          	  for (sy=0;sy<=NUM+1;sy++)
            	    for (sz=0; sz<=NUM+1;sz++)
                      w[sx][sy][sz] = 0.1;
		break;
	  case 4:
		// test grid
		if (nodeID==1) cout << "==> Initial wavefunction : Boolean Test" << endl;
        	for (sx=0;sx<=NUMX+1;sx++)
          	  for (sy=0;sy<=NUM+1;sy++)
            	    for (sz=0; sz<=NUM+1;sz++)
                      w[sx][sy][sz] = (sx%2)*(sy%2)*(sz%2);
		break;
	  default:
		cout << "Invalid initial condition type" << endl;
		exit(0);
		break;
	}

	// enforce BCs
        for (sx=0;sx<=NUMX+1;sx++)
          for (sy=0;sy<=NUM+1;sy++) {
		w[sx][sy][0] = 0;
		w[sx][sy][NUM+1] = 0;
	  }

        for (sz=0;sz<=NUM+1;sz++)
          for (sy=0;sy<=NUM+1;sy++) {
		w[0][sy][sz] = 0;
		w[NUMX+1][sy][sz] = 0;
	  }

        for (sx=0;sx<=NUMX+1;sx++)
          for (sz=0;sz<=NUM+1;sz++) {
		w[sx][0][sz] = 0;
		w[sx][NUM+1][sz] = 0;
	  }

	// zero out updated wavefnc for safety's sake
        for (sx=0;sx<NUMX+2;sx++)
          for (sy=0;sy<NUM+2;sy++)
            for (sz=0;sz<NUM+2;sz++)
		W[sx][sy][sz] = 0;

	// symmetrize the intial condition
	symmetrizeWavefunction();
}

void symmetrizeWavefunction()
{
	int sx,sy,sz,x,y,z,sign=1;
	double check=0;

	// set sign for inversions; convention is that odd numbers 
	// are symmetric and even numbers antisymmetric 
	sign = 2*(INITSYMMETRY%2)-1;

	switch (INITSYMMETRY) { 
	  case 0:
		// no symmetry
		break;
	  case 1:
		// symmetric in z
	  case 2:
		// antisymmetric in z
        	for (sx=0;sx<NUMX+2;sx++)
		{
		  x=sx; 	
          	  for (sy=1;sy<=NUM;sy++)
		  {
		    y=sy; 	
            	    for (sz=1;sz<=NUM;sz++)
		    {
			z=sz;
		    	if (z>NUM/2)
		      	  z = NUM + 1 - z;
			if (sz>NUM/2)
			  w[sx][sy][sz] = ((dcomp)sign)*w[x][y][z];
		    }
		  }
		}
		break;
	  case 3:
		// symmetric in y
	  case 4:
		// antisymmetric in y
        	for (sx=0;sx<NUMX+2;sx++)
		{
		  x=sx; 	
          	  for (sy=1;sy<=NUM;sy++)
		  {
		    y=sy; 	
		    if (y>NUM/2)
		      y = NUM + 1 - y;
            	    for (sz=1;sz<=NUM;sz++)
		    {
			z=sz;
			if (sy>NUM/2)
			  w[sx][sy][sz] = ((dcomp)sign)*w[x][y][z];
		    }
		  }
		}
		break;
	  default:
		cout << "Invalid symmetry type" << endl;
		exit(0);
		break;
	}
}
