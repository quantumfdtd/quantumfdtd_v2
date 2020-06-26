/*

   ouputroutines.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "outputroutines.h"
#include "potential.h"

void outputMeasurements(const double time) {

	dcomp ener = energyCollect/normalizationCollect;
	dcomp binding = ener - vInfinityCollect/normalizationCollect;
	dcomp rRMS2 = rRMS2Collect/normalizationCollect;  
        //dcomp yAvg = yAvgCollect/normalizationCollect;  

	// output to screen

	cout.precision(8);
	cout.width(dwidth); cout << time;
	cout.width(dwidth); cout << setprecision (7) << ener;
	cout.width(dwidth); cout << setprecision (7) << binding;
	cout.width(dwidth); cout << setprecision (7) << sqrt(real(rRMS2));  
	//cout.width(dwidth); cout << setprecision (7) << real(yAvg);  
	cout << endl;
}

void outputSummaryData(string label) {


      dcomp ener = energyCollect/normalizationCollect;
      dcomp binding = ener - vInfinityCollect/normalizationCollect;
      dcomp rRMS2 = rRMS2Collect/normalizationCollect;  
      dcomp xAvg = xAvgCollect/normalizationCollect;  
      dcomp yAvg = yAvgCollect/normalizationCollect;  
      dcomp zAvg = zAvgCollect/normalizationCollect;  

      print_line();
      cout << "==> " << label << " Energy : " << ener << endl;
      cout << "==> " << label << " Binding Energy : " << binding << endl;
      cout << "==> " << label << " r_RMS : " << sqrt(real(rRMS2)) << endl;  
      cout << "==> " << label << " <x> : " << xAvg << endl;  
      cout << "==> " << label << " <y> : " << yAvg << endl;  
      cout << "==> " << label << " <z> : " << zAvg << endl;  
      cout << "==> " << label << " L/r_RMS : " << float(NUM)/sqrt(real(rRMS2)) << endl;  
      print_line();

      // convert label to lower case and convert spaces to underscores for file name
      const int length = label.length();
      for(int i=0; i < length; ++i)
      {
        label[i] = std::tolower(label[i]);
	if (label[i]==' ') label[i]='_';
      }

      // output summary data to output file for later use
      fstream out;
      char fname[255];
      sprintf(fname,"data/%s.out",label.c_str());
      out.open(fname, ios::out);
      out.precision(8);
      //out << T << "\t";
      //out << XI << "\t";
      out << eB << "\t";
      out << Kx << "\t";
      out << real(binding) << "\t";
      //out << imag(binding) << "\t";
      out << real(xAvg) << "\t";
      out << real(yAvg) << "\t";
      out << real(zAvg) << endl;
      out.close();
}

void outputSnapshot(dcomp ***wfnc, char* label) {

  int x;
  static int h=NUM/2;
  static int hx=NUMX/2;
  dcomp data;

  fstream out;
  char fname[255];

  // dump wavefunction

  // output slices suitable for 2d viewing
  sprintf(fname,"data/snapshot/wavefunction_%s.dat",label);
  out.open(fname, ios::out);
  out.precision(10);
  for (int s=0;s<=NUMX+1;s++) {
    x=(nodeID-1)*NUMX + s;	  
    out << x << "\t";
    data = 0.5*(wfnc[s][h][h]+wfnc[s][h+1][h+1]);
    out << scientific << real(data) << "\t";
    out << scientific << imag(data);
    out << endl;
  }
  out << "&&" << endl;
  for (int s=0;s<=NUM+1;s++) {
    out << s << "\t";
    data = 0.5*(wfnc[hx][s][h]+wfnc[hx+1][s][h+1]);
    out << scientific << real(data) << "\t";
    out << scientific << imag(data);
    out << endl;
  }
  out << "&&" << endl;
  for (int s=0;s<=NUM+1;s++) {
    out << s << "\t";
    data = 0.5*(wfnc[hx][h][s]+wfnc[hx+1][h+1][s]);
    out << scientific << real(data) << "\t";
    out << scientific << imag(data);
    out << endl;
  }
  out.close();
  
  return;
}

void outputWavefunction(dcomp ***wfnc, char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"data/wavefunction_%s.dat",label);

  cout << "==> Dumping wave function to " << fname << endl;

  out.open(fname, ios::out);
  out.precision(8);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
                out << real(wfnc[sx][sy][sz]) << "\t";
                out << imag(wfnc[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v 3d
void outputPotential(char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"data/potential_%s.dat",label);

  cout << "==> Dumping potential to " << fname << endl;

  out.open(fname, ios::out);
  out.precision(8);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
                out << real(v[sx][sy][sz]) << "\t";
                out << imag(v[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v along principal axes
void dumpPotential() {

  int h=NUM/2;
  fstream out;

  out.open("data/potential.dat", ios::out);
  for (int s=0;s<=NUM+1;s++) {
                out << s << "\t";
                out << v[s][h][h] << "\t";
                out << v[h][s][h] << "\t";
                out << v[h][h][s] << "\t";
                out << endl;
  }
  out.close();
  return;

}

void print_line() {
        for (int i=0;i<4*dwidth;i++) cout << "-"; cout << endl;
        return;
}

