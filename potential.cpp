/*

   potential.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <cmath>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "specialfunctions.h"
#include "potential.h"
#include "intde2.h"

#define Power(a,b)	pow((double)a,(double)b)
#define Sqrt(a)		sqrt((double)a)

// used for numerical integration when necessary
int     lenaw=LENAW;
double  aw[LENAW];

// global variable useful for subroutines
double  dx,dy,dz;
double  r;
double  md;
double  mdp;

// determines x distance to center of simulation volume in lattice units
double distx(int sx) 
{ 
	return sx - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
}

// determines y distance to center of simulation volume in lattice units
double disty(int sy) 
{ 
	return sy - ((double)NUM+1.)/2.;
}

// determines z distance to center of simulation volume in lattice units
double distz(int sz) 
{ 
	return sz - ((double)NUM+1.)/2.;
}

// determines square of distance to center of simulation volume in lattice units
double distsq(int sx,int sy, int sz) 
{
	double dx,dy,dz,r2;

	// coordinate system is centered in simulation volume 
	dx = sx - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
	dy = sy - ((double)NUM+1.)/2.;
	dz = sz - ((double)NUM+1.)/2.;
	r2 = (dx*dx+dy*dy+dz*dz);
	return r2;
}



dcomp potential(int sx,int sy, int sz) 
{
	double temp,iV,rV;
	double res,err;
	double m12,B,wc,e,rho;
	
	// coordinate system is centered in simulation volume 
	dx = ((double) sx) - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
	dy = ((double) sy) - ((double)NUM+1.)/2.;
	dz = ((double) sz) - ((double)NUM+1.)/2.;
	r = A*sqrt(dx*dx+dy*dy+dz*dz);
	rho = A*sqrt(dx*dx+dy*dy);
	
	switch(POTENTIAL) {
	  case 0:
		// none
          	return 0.;
		break;
	  case 1:
		// cubic well
		if ( (sx>NUM/4 && sx<=3*NUM/4) && (sy>NUM/4 && sy<=3*NUM/4) && (sz>NUM/4 && sz<=3*NUM/4) )
			return -10.0;
		else
			return 0.0;
		break;
	  case 2:
		// quadrilateral-well in center of cube with short side in z direction
		if ( (sx>NUM/4 && sx<=3*NUM/4) && (sy>NUM/4 && sy<=3*NUM/4) && (sz>3*NUM/8 && sz<=5*NUM/8) )
			return -10.0;
		else
			return 0.0;
		break;
	  case 3:
		// 3d periodic
		temp  = sin(2*M_PI*(sx-1)/(NUM-1))*sin(2*M_PI*(sx-1)/(NUM-1));
		temp *= sin(2*M_PI*(sy-1)/(NUM-1))*sin(2*M_PI*(sy-1)/(NUM-1));
		temp *= sin(2*M_PI*(sz-1)/(NUM-1))*sin(2*M_PI*(sz-1)/(NUM-1));
		return -temp+1;
		break;
	  case 4:
		// coulomb
		if (r < A)
		  return 0.0;
		else
		  return -1./r + 1./A;
		break;
	  case 5:
		// elliptical coulomb
		dz *= 2;
		r = A*sqrt(dx*dx+dy*dy+dz*dz);
		if (r < A)
		  return 0.0;
		else
		  return -1./r + 1./A;
		break;
	  case 6:
		// cornell plus spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4);
		if (r < A)
		  return 4*MASS;
		if (r>5.5745)
		  r = 5.5745;
		//return -0.385/r + SIGMA*r + 4*MASS;
		return -0.385/r + SIGMA*r - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 7:
		// screened cornell
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md + 4*MASS;
		break;
	  case 8:
		// screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 9:
		// anisotropically screened short distance piece + isotropic cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-mu(T,1.4)*r))/mu(T,1.4)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 10:
		// anisotropically screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/mu(T,1.4)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 11:
		// fully anisotropic screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 12:
		// fully anisotropic screened cornell using small xi expression for mu + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1));
		if (r < A)
		  return 4*MASS;
		else
		  return -alphas(2*M_PI*T)*(4./3.)*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 13:
		// modified fully anisotropic screened cornell + spin correction
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(T,1.4)*(pow(1+1.85*pow(XI,1.27),-0.20)+(pow(1+0.74*pow(XI,1.20),-0.23)-pow(1+1.85*pow(XI,1.27),-0.20))*(1-A*A*dz*dz/(r*r)));
		if (r < A)
			return 4*MASS;
		else
			return -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 14:
		// newpotential add entropy contribution
		// units here are GeV for energy/momentum and GeV^(-1) for distance
		md = mu(3,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		if (r < A)
			return 4*MASS;
		else
			return -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 15:
		// 3d harmonic oscillator
          	return r*r/2;
		break;			
	  case 16:
		// Mickey Mouse's Head
		double Dx, Dy, Dz, R;
		if (r/A <= NUM/4) return -100.0; // head
		Dx = ((double) sx) - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
		Dy = sy - ((double)NUM+1.)/2. - (1/sqrt(2.)+0.25)*NUM/4;
		Dz = sz - ((double)NUM+1.)/2. - (1/sqrt(2.)+0.25)*NUM/4;
		R = sqrt(4*Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUM/8) return -105.0; // ear
		Dy = sy - ((double)NUM+1.)/2. + (1/sqrt(2.)+0.25)*NUM/4;
		Dz = sz - ((double)NUM+1.)/2. - (1/sqrt(2.)+0.25)*NUM/4;
		R = sqrt(4*Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUM/8) return -105.0; // ear
		Dx = ((double) sx) - ((double)NUMX+1.)/2. - ((double)NUM/8.)+ ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
		Dy = sy - ((double)NUM+1.)/2.;
		Dz = sz - ((double)NUM+1.)/2.;
		R = sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
		if (R < NUM/6) return -100.0; // nose
		return 0;
		break;
	  case 17:
		// Dodecahedron
		double x, y, z;
		x = dx/((NUM-1)/2);
		y = dy/((NUM-1)/2);
		z = dz/((NUM-1)/2);
		if (12.70820393249937 + 11.210068307552588*x >= 14.674169922690343*z && 11.210068307552588*x <= 12.70820393249937 + 14.674169922690343*z && 5.605034153776295*(3.23606797749979*x - 1.2360679774997896*z) <= 6.*(4.23606797749979 + 5.23606797749979*y) && 18.1382715378281*x + 3.464101615137755*z <= 12.70820393249937 && 9.06913576891405*x + 15.70820393249937*y <= 12.70820393249937 + 3.464101615137755*z && 9.70820393249937*y <= 12.70820393249937 + 5.605034153776294*x + 14.674169922690343*z && 12.70820393249937 + 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z >= 0. && 15.70820393249937*y + 3.464101615137755*z <= 12.70820393249937 + 9.06913576891405*x && 5.605034153776295*(-6.47213595499958*x - 1.2360679774997896*z) <= 25.41640786499874 && 3.464101615137755*z <= 9.06913576891405*x + 3.*(4.23606797749979 + 5.23606797749979*y) && 1.7320508075688772*(3.23606797749979*x + 8.47213595499958*z) <= 3.*(4.23606797749979 + 3.23606797749979*y) && 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z <= 12.70820393249937) return -100.0;
		else
		  return 0.;
	  case 18:
		// Complex 3d harmonic oscillator
		return dcomp(1.,1.)*dcomp(r*r/2,0.);
		break;
	  case 19:			
			// complex coulomb
			if (r < A)
				return dcomp(0.,1.)*dcomp(-1./A,0.);
			else
				return dcomp(0.,1.)*dcomp(-1./r,0.);
			break;			
	  case 20:
		// Complex Dodecahedron
		// double x, y, z;
		x = dx/((NUM-1)/2);
		y = dy/((NUM-1)/2);
		z = dz/((NUM-1)/2);
		if (12.70820393249937 + 11.210068307552588*x >= 14.674169922690343*z && 11.210068307552588*x <= 12.70820393249937 + 14.674169922690343*z && 5.605034153776295*(3.23606797749979*x - 1.2360679774997896*z) <= 6.*(4.23606797749979 + 5.23606797749979*y) && 18.1382715378281*x + 3.464101615137755*z <= 12.70820393249937 && 9.06913576891405*x + 15.70820393249937*y <= 12.70820393249937 + 3.464101615137755*z && 9.70820393249937*y <= 12.70820393249937 + 5.605034153776294*x + 14.674169922690343*z && 12.70820393249937 + 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z >= 0. && 15.70820393249937*y + 3.464101615137755*z <= 12.70820393249937 + 9.06913576891405*x && 5.605034153776295*(-6.47213595499958*x - 1.2360679774997896*z) <= 25.41640786499874 && 3.464101615137755*z <= 9.06913576891405*x + 3.*(4.23606797749979 + 5.23606797749979*y) && 1.7320508075688772*(3.23606797749979*x + 8.47213595499958*z) <= 3.*(4.23606797749979 + 3.23606797749979*y) && 5.605034153776294*x + 9.70820393249937*y + 14.674169922690343*z <= 12.70820393249937) return dcomp(-100.,-100.);
		else
		  return 0.;
	  case 21:
		// Anisotropically Screened Quarkonium Potential with Imaginary Part
		if (r < A)
		  return dcomp(4*MASS,0.);
		// real part
		md = mu(3,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		rV = -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.* SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		// imaginary part
		md = mu(3,1.4); // do not include angular modification of md for imaginary part since this is already acctd for!
		iV = ImV(r*md,acos(A*dz/r),XI); // ImV defined in specialfunctions.cpp
		iV *= 0.385*T*TC;
		return dcomp(rV,iV);
		break;
	  case 22: 
		// Include Large-XI Behavior of md
		x = A*dz/r; // Cos(theta) 
		double b;
		b = 9./16.; // phenomenological value
		md = mu(3,1.4)*Power(1 + XI*(1 + Power(2,1 + b)*(-1 + Power(x,2))*Power(1 + XI,2)*Power(2 + XI,-2 - b)),-0.25)*T/3;
		if (r < A)
			return 4*MASS;
		else 
			return -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.*SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS; 
		break;
	  case 23:
		// Real Part
		// Include Large-XI Behavior of md - PRL
		x = A*dz/r; // Cos(theta) 
		b = 9./16.; // phenomenological value
		md = mu(3,1.4)*Power(1 + XI*(1 + Power(2,1 + b)*(-1 + Power(x,2))*Power(1 + XI,2)*Power(2 + XI,-2 - b)),-0.25)*T/3;
		if (r < A)
			rV = 4*MASS;
		else 
			rV =  -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.*SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS; 

		// imaginary part
		md = mu(3,1.4); // do not include angular modification of md for imaginary part since this is already acctd for! // this choice also fixes md to md @ 3 Tc. :|
		iV = ImV(r*md,acos(A*dz/r),XI); // ImV defined in specialfunctions.cpp
		iV *= 0.385*T*TC;
		return dcomp(rV,iV);
		break;
	  case 24:
		// Real Part
		// Include new model of Large-XI Behavior of md -- NPA B
		// Running Coupling 
		x = A*dz/r; // Cos(theta) 
		md = mu(T,1.4)*Power(1. + XI*(1. - (0.0944049*(2.16919 - 29.6088*Power(x,2))*Power(1. + XI,1.5))/(3. + Power(XI,2)))*(1.62114 - (1.*(0.878423 + Power(1. + XI,0.125)))/Sqrt(3. + XI)),-0.25);
		if (r < A)
			rV = 4*MASS;
		else 
			rV =  -0.385*exp(-md*r)/r*(1.0 + md*r) + 2.*SIGMA*(1. - exp(-md*r))/md - SIGMA*r*exp(-md*r)  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS; 
		// imaginary part
		md = mu(T,1.4); // do not include angular modification of md for imaginary part since this is already acctd for!
		iV = ImV(r*md,acos(A*dz/r),XI); // ImV defined in specialfunctions.cpp
		iV *= 4*alphas(2*M_PI*T)*T*TC/3;
		return dcomp(rV,iV);
		break;
	  case 25:
		// Real Part
		// Include new model of Large-XI Behavior of md -- NPA A
		// Running Coupling 
		x = A*dz/r; // Cos(theta) 
		md = mu(T,1.4)*Power(1. + XI*(1. - (0.0944049*(2.16919 - 29.6088*Power(x,2))*Power(1. + XI,1.5))/(3. + Power(XI,2)))*(1.62114 - (1.*(0.878423 + Power(1. + XI,0.125)))/Sqrt(3. + XI)),-0.25);
		if (r < A)
			rV = 4*MASS;
		else 
			rV =  -0.385*exp(-md*r)/r + SIGMA*(1. - exp(-md*r))/md  - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS; 
		// imaginary part
		md = mu(T,1.4); // do not include angular modification of md for imaginary part since this is already acctd for!
		iV = ImV(r*md,acos(A*dz/r),XI); // ImV defined in specialfunctions.cpp
		iV *= 4*alphas(2*M_PI*T)*T*TC/3;
		return dcomp(rV,iV);
		break;
          case 26:
                // cornell + magnetic field
                // units here are GeV for energy/momentum and GeV^(-1) for distance
                m12 = 4.2048; // m1 + m2 for B0 meson
		wc = eB/MASS/3; // factor of 1/3 due to charge of bottom quark
                if (r < A)
                  return m12;
                  //return b;
                else
                return -0.385/r + SIGMA*r + MASS*wc*wc*rho*rho/8 + m12;
                break;
          case 27:
		// 3d harmonic oscillator + magnetic field - w0 = 1 - shifted
		wc  = eB/MASS;
		double a,bb,c,d,kx,ky; 
		kx = Kx;
		ky = 0;
		a = MASS*(1+wc*wc/4);	
		bb = wc*ky/4;
		c = wc*kx/4;
		d = MASS;
		return a*rho*rho/2 + d*A*A*dz*dz/2;
          case 28:
		// 3d harmonic oscillator + magnetic field - w0 = 1 - not shifted
		wc  = eB/MASS;
		kx = Kx;
		ky = 0;
		a = MASS*(1+wc*wc/4);	
		bb = wc*ky/4;
		c = wc*kx/4;
		d = MASS;
		return a*rho*rho/2 -bb*dx*A + c*dy*A + d*A*A*dz*dz/2;
          case 29:
                // coulomb + magnetic field
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = eB/MASS;
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
                return -bb*dx*A + c*dy*A - 1./r + MASS*wc*wc*rho*rho/8.;
                break;
          case 30:
                // cornell + magnetic field (no tune)
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = eB/MASS/3;
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
                return -bb*dx*A + c*dy*A - 0.385/r + SIGMA*r + MASS*wc*wc*rho*rho/8.;
                break;
          case 31:
                // cornell + magnetic field + spin-spin interaction : J/Psi (J/Psi tuned)
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = 2*eB/MASS/3; // J/Psi
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
		a = 2.06; // J/Psi
		d = 1.982; // J/Psi
                return -bb*dx*A + c*dy*A - 0.312/r + SIGMA*r + MASS*wc*wc*rho*rho/8. + SPINEXP*a*exp(-d*r);
                break;
          case 32:
                // cornell + magnetic field + spin-spin interaction : J/Psi (Upsilon tuned)
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = 2*eB/MASS/3; // J/Psi
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
		a = 0.825; // J/Psi
		d = 1.982; // J/Psi
                return -bb*dx*A + c*dy*A - 0.42059/r + SIGMA*r + MASS*wc*wc*rho*rho/8. + SPINEXP*a*exp(-d*r);
                break;
          case 33:
                // cornell + magnetic field + spin-spin interaction : Upsilon (Upsilon tuned)
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = -eB/MASS/3; // Upsilon
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
		a = 0.318; // Upsilon
		d = 1.982; // Upsilon
                return -bb*dx*A + c*dy*A - 0.42059/r + SIGMA*r + MASS*wc*wc*rho*rho/8. + SPINEXP*a*exp(-d*r);
                break;
          case 34:
                // cornell + magnetic field + spin-spin interaction : J/Psi -- SHIFTED ALONG Y AXIS (Upsilon tuned)
                // units here are GeV for energy/momentum and GeV^(-1) for distance
		wc = 2*eB/MASS/3; // J/Psi
		kx = Kx;
		ky = 0;
		bb = wc*ky/4;
		c = wc*kx/4;
		a = 0.825; // J/Psi
		d = 1.982; // J/Psi	
		double y0;
		y0 = round((4*SIGMA*MASS - 2*kx*eB/3)/(4*eB*eB/9)/A);
		dy += y0;
		r = A*sqrt(dx*dx+dy*dy+dz*dz);
		rho = A*sqrt(dx*dx+dy*dy);
                return -bb*dx*A + c*dy*A - 0.42059/r + SIGMA*r + MASS*wc*wc*rho*rho/8. + SPINEXP*a*exp(-d*r);
                break;
	  default:
		return 0.;
		break;
	}
}

// returns value of potential which should be subtracted when computing binding energies
dcomp potentialSub(int sx, int sy, int sz) 
{
	double iV,rV;

	// coordinate system is centered in simulation volume 
	dx = ((double) sx) - ((double)NUMX+1.)/2. + ( ((double)nodeID) - ((double)numNodes)/2. )*NUMX;
	dy = ((double) sy) - ((double)NUM+1.)/2.;
	dz = ((double) sz) - ((double)NUM+1.)/2.;
	r = A*sqrt(dx*dx+dy*dy+dz*dz);

	switch(POTENTIAL) {
	  case 0:
	  case 1:
	  case 2:
	  case 3:
		return 0.;
		break;
	  case 4:
	  case 5:
		return 1./A;
		break;
	  case 6:
		r = 5.5745;
		//return -0.385/r + SIGMA*r + 4*MASS;
		return -0.385/r + SIGMA*r - 0.8*SIGMA/(4*MASS*MASS*r) + 4*MASS;
		break;
	  case 7:
	  case 8:
	  case 9:
	  case 10:
		return SIGMA/mu(T,1.4) + 4*MASS;
		break;
	  case 11:
		md = mu(T,1.4)*(1 + 0.07*pow(XI,0.2)*(1-A*A*dz*dz/(r*r)))*pow(1+XI,-0.29);
		return SIGMA/md + 4*MASS;
		break;
	  case 12:
		md = mu(T,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1));
		//cout << SIGMA <<  ", " << md << ", " << MASS << endl;
		return SIGMA/md + 4*MASS;
		break;
	  case 13:
		md = mu(T,1.4)*(pow(1+1.85*pow(XI,1.27),-0.20)+(pow(1+0.74*pow(XI,1.20),-0.23)-pow(1+1.85*pow(XI,1.27),-0.20))*(1-A*A*dz*dz/(r*r)));
		return SIGMA/md + 4*MASS;
		break;
	  case 14:
		md = mu(3,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		break;
	  case 15:
		return 0.;
		break;			
	  case 16:
		return 0.;
		break;			
	  case 17:
		return 0.;
		break;			
	  case 18:
		return 0.;
		break;	
	  case 19:
		return 0.;
		break;	
	  case 20:
		return 0.;
		break;		
	  case 21:
		md = mu(3,1.4)*(1 - 0.125*XI*(A*A*dz*dz/(r*r)+1))*T/3;
		rV = 2*SIGMA/md + 4*MASS;
		//iV = -1 + XI/6;
		//iV *= 0.385*T*TC;
		iV = 0;
		return dcomp(rV,iV);
		break;				
	  case 22:
		double x,b;
		x = A*dz/r; // Cos(theta) 
		b = 9./16.; // phenomenological value
		md = mu(3,1.4)*Power(1 + XI*(1 + Power(2,1 + b)*(-1 + Power(x,2))*Power(1 + XI,2)*Power(2 + XI,-2 - b)),-0.25)*T/3;
		return 2*SIGMA/md + 4*MASS;
		break;
	  case 23:
		x = A*dz/r; // Cos(theta) 
		b = 9./16.; // phenomenological value
		md = mu(3,1.4)*Power(1 + XI*(1 + Power(2,1 + b)*(-1 + Power(x,2))*Power(1 + XI,2)*Power(2 + XI,-2 - b)),-0.25)*T/3;
		rV = 2*SIGMA/md + 4*MASS;
		iV = 0;
		return dcomp(rV,iV);
		break;
	  case 24:
		x = A*dz/r; // Cos(theta) 
		md = mu(T,1.4)*Power(1. + XI*(1. - (0.0944049*(2.16919 - 29.6088*Power(x,2))*Power(1. + XI,1.5))/(3. + Power(XI,2)))*(1.62114 - (1.*(0.878423 + Power(1. + XI,0.125)))/Sqrt(3. + XI)),-0.25);
		rV = 2*SIGMA/md + 4*MASS;
		iV = 0;
		return dcomp(rV,iV);
		break;
	  case 25:
		x = A*dz/r; // Cos(theta) 
		md = mu(T,1.4)*Power(1. + XI*(1. - (0.0944049*(2.16919 - 29.6088*Power(x,2))*Power(1. + XI,1.5))/(3. + Power(XI,2)))*(1.62114 - (1.*(0.878423 + Power(1. + XI,0.125)))/Sqrt(3. + XI)),-0.25);
		rV = SIGMA/md + 4*MASS;
		iV = 0;
		return dcomp(rV,iV);
		break;
	  case 26:
		return 0.;
		break;
	  case 27:
		return 0.;
		break;
	  case 28:
		return 0.;
		break;
	  case 29:
		return 0.;
		break;
	  case 30:
		return 0.;
		break;
	  case 31:
		return 0.;
		break;
	  case 32:
		return 0.;
		break;
	  case 33:
		return 0.;
		break;
	  case 34:
		return 0.;
		break;
	  default:
		return 0.;
		break;
	}
	return 0.;
}

// phir integrand
double phir(double z)
{
	return 2*z*(1-sin(z*r*md)/(z*r*md))/(z*z+1)/(z*z+1);
}

// utility function for ps1 and ps2
double psig(double r, double z)
{	// note r here is rhat
	return (r*z*cos(r*z)-sin(r*z))/(r*r*r*z*z*z);
}

// psi1 integrand
double psi1(double z)
{
	double cos2theta = A*A*dz*dz/(r*r);
	double sin2theta = 1. - cos2theta;
	return z*(1 - 1.5*(sin2theta*sin(z*md*r)/(z*md*r) + (1-3*cos2theta)*psig(r*md, z)))/(z*z+1)/(z*z+1);
}

// psi2 integrand
double psi2(double z)
{
	double cos2theta = A*A*dz*dz/(r*r);
	double sin2theta = 1. - cos2theta;
	return -4*z*(1 - 3*((2./3. - cos2theta)*sin(z*md*r)/(z*md*r) + (1-3*cos2theta)*psig(r*md, z)))/(z*z+1)/(z*z+1)/(z*z+1)/3;

}


// three loop running coupling
// mu is assumed to be in units of TC
double alphas(double mu)
{
        double b0,b1,b2,L,R,nc,nf,lambda_ms,t; 

	nc = 3;
	nf = (double) NF;

	b0 = (11*nc-2*nf)/(12*M_PI);
	b1 = (17*nc*nc-nf*(10*nc+6*(nc*nc-1)/(2*nc))/2)/(24*M_PI*M_PI);
	b2 = (2857 - 5033*nf/9 + 325*nf*nf/27)/(128*M_PI*M_PI*M_PI);

	lambda_ms = 0.344;
	t = 2*log(mu*TC/lambda_ms);

	return (1 - (b1*log(t))/(Power(b0,2)*t) + (b0*b2 + Power(b1,2)*(-1 - log(t) + Power(log(t),2)))/(Power(b0,4)*Power(t,2)) - (3*b0*b1*b2*log(t) + Power(b1,3)*(0.5 - 2*log(t) - (5*Power(log(t),2))/2. + Power(log(t),3)))/(Power(b0,6)*Power(t,3)))/(b0*t);
}

// derivative of three loop running coupling with respect to mu
// mu is assumed to be in units of TC
double alphasp(double mu)
{
        double b0,b1,b2,L,R,nc,nf,lambda_ms,t; 

	nc = 3;
	nf = (double) NF;

	b0 = (11*nc-2*nf)/(12*M_PI);
	b1 = (17*nc*nc-nf*(10*nc+6*(nc*nc-1)/(2*nc))/2)/(24*M_PI*M_PI);
	b2 = (2857 - 5033*nf/9 + 325*nf*nf/27)/(128*M_PI*M_PI*M_PI);

	lambda_ms = 0.344;
	t = 2*log(mu*TC/lambda_ms);

	return ((-2*(-4*Power(b1,3) + 3*b0*b1*b2 + b1*log(t)*(3*(Power(b1,2) - 4*b0*b2) + Power(b1,2)*(13 - 4*log(t))*log(t)) + Power(b0,2)*t*(-2*Power(b1,2) + 3*b0*b2 + Power(b1,2)*log(t)*(-5 + 3*log(t)) + Power(b0,2)*t*(b1 + Power(b0,2)*t - 2*b1*log(t)))))/(Power(b0,7)*mu*Power(t,5)));
}

// debye screening mass
double mu(double t, double fac)
{
	return fac*sqrt( (1+((double)NF)/6)*4*M_PI*alphas(2*M_PI*t) )*t*TC;
}

// derivative of debye screening mass wrt T
double mup(double t, double fac)
{
	return (fac*(6 + NF)*Sqrt((2*M_PI)/3.)*(alphas(2*M_PI*t) + M_PI*t*alphasp(2*M_PI*t)))/Sqrt((6 + NF)*alphas(2*M_PI*t));
}
