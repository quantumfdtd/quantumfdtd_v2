/*

   potential.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __potential_h__
#define __potential_h__

dcomp potential(int sx, int sy, int sz);
dcomp potentialSub(int sx, int sy, int sz);
double alphas(double mu);
double alphasp(double mu);
double mu(double t, double fac);
double mup(double t, double fac);

double distx(int sx);
double disty(int sy);
double distz(int sz);
double distsq(int sx, int sy, int sz);

double phir(double z);
double psi1(double z);
double psi2(double z);

#endif /* __potential_h__ */
