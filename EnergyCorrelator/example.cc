// Example showing usage of energy correlator classes.
//
// Compile it with "make example" and run it with
//
//   ./example < ../data/single-event.dat
//
// Copyright (c) 2013
// Andrew Larkoski, Gavin Salam, and Jesse Thaler
//
// $Id: example.cc 758 2014-11-13 15:45:06Z larkoski $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <ctime>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <sstream>
#include "EnergyCorrelator.hh" // In external code, this should be fastjet/contrib/EnergyCorrelator.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how this EnergyCorrelator contrib works

  analyze(event);

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}

////////
//
//  Main Routine for Analysis 
//
///////

void analyze(const vector<PseudoJet> & input_particles) {

   /////// EnergyCorrelator /////////////////////////////
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm; 
   double jet_rad = 1.00; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   ClusterSequence clust_seq(input_particles,jetDef);
   vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());
   
   for (int j = 0; j < 2; j++) { // Two hardest jets per event
      if (antikt_jets[j].perp() > 200) {
         
         PseudoJet myJet = antikt_jets[j];
         
         // various values of beta
         vector<double> betalist;
         betalist.push_back(0.1);
         betalist.push_back(0.2);
         betalist.push_back(0.5);
         betalist.push_back(1.0);
         betalist.push_back(1.5);
         betalist.push_back(2.0);

         // checking the two energy/angle modes
         vector<EnergyCorrelator::Measure> measurelist;
         measurelist.push_back(EnergyCorrelator::pt_R);
         measurelist.push_back(EnergyCorrelator::E_theta); 

         vector<string> modename;
         modename.push_back("pt_R");
         modename.push_back("E_theta");

         for (unsigned int M = 0; M < measurelist.size(); M++) {
            
            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelator:  ECF(N,beta) with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s %15s\n","beta", "N=1 (GeV)", "N=2 (GeV^2)", "N=3 (GeV^3)", "N=4 (GeV^4)", "N=5 (GeV^5)");
            
            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];
               
               EnergyCorrelator ECF0(0,beta,measurelist[M]);
               EnergyCorrelator ECF1(1,beta,measurelist[M]);
               EnergyCorrelator ECF2(2,beta,measurelist[M]);
               EnergyCorrelator ECF3(3,beta,measurelist[M]);
               EnergyCorrelator ECF4(4,beta,measurelist[M]);
               EnergyCorrelator ECF5(5,beta,measurelist[M]);

               printf("%7.3f %14.2f %14.2f %14.2f %14.2f %15.2f \n",beta,ECF1(myJet),ECF2(myJet),ECF3(myJet),ECF4(myJet),ECF5(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorRatio:  r_N^(beta) = ECF(N+1,beta)/ECF(N,beta) with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s %15s \n","beta", "N=0 (GeV)", "N=1 (GeV)", "N=2 (GeV)", "N=3 (GeV)","N=4 (GeV)");
            
            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];
               
               EnergyCorrelatorRatio r0(0,beta,measurelist[M]);
               EnergyCorrelatorRatio r1(1,beta,measurelist[M]);
               EnergyCorrelatorRatio r2(2,beta,measurelist[M]);
               EnergyCorrelatorRatio r3(3,beta,measurelist[M]);
               EnergyCorrelatorRatio r4(4,beta,measurelist[M]);

               printf("%7.3f %14.4f %14.4f %14.4f %14.4f %15.4f \n",beta,r0(myJet),r1(myJet),r2(myJet),r3(myJet),r4(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorDoubleRatio:  C_N^(beta) = r_N^(beta)/r_{N-1}^(beta) with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3", "N=4");
            
            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];
               
               EnergyCorrelatorDoubleRatio C1(1,beta,measurelist[M]);
               EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M]);
               EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M]);
               EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M]);

               printf("%7.3f %14.6f %14.6f %14.6f %14.6f \n",beta,C1(myJet),C2(myJet),C3(myJet),C4(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl;
 
            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC1:  C_1^(beta) = ECF(2,beta)/ECF(1,beta)^2 with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C1 obs");

            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];

               EnergyCorrelatorC1 c1(beta,measurelist[M]);

               printf("%7.3f %14.6f \n",beta,c1(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorC2:  C_2^(beta) = ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2 with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","C2 obs");

            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];

               EnergyCorrelatorC2 c2(beta,measurelist[M]);

               printf("%7.3f %14.6f \n",beta,c2(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl;

           
            cout << "-------------------------------------------------------------------------------------" << endl;
            cout << "EnergyCorrelatorD2:  D_2^(beta) = ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3 with " << modename[M] << endl;
            cout << "-------------------------------------------------------------------------------------" << endl;
            printf("%7s %14s \n","beta","D2 obs");

            for (unsigned int B = 0; B < betalist.size(); B++) {
               double beta = betalist[B];

               EnergyCorrelatorD2 d2(beta,measurelist[M]);

               printf("%7.3f %14.6f \n",beta,d2(myJet));
            }
            cout << "-------------------------------------------------------------------------------------" << endl << endl; 

                       
            // timing tests for the developers
            double do_timing_test = false;
            if (do_timing_test) {
            
               cout << "jet with pt = " << myJet.pt() << " and " << myJet.constituents().size() << " constituents" << endl;

               clock_t clock_begin, clock_end;
               double num_iter;
               double beta = 0.5;

               cout << setprecision(6);

               // test C1
               num_iter = 20000;
               clock_begin = clock();
               EnergyCorrelatorDoubleRatio C1s(1,beta,measurelist[M],EnergyCorrelator::slow);
               EnergyCorrelatorDoubleRatio C1f(1,beta,measurelist[M],EnergyCorrelator::storage_array);
               cout << "timing " << C1s.description() << endl;
               cout << "timing " << C1f.description() << endl;
               for (int t = 0; t < num_iter; t++) {
                  C1s(myJet);
               }
               clock_end = clock();
               cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C1"<< endl;

               num_iter = 20000;
               clock_begin = clock();
               for (int t = 0; t < num_iter; t++) {
                  C1f(myJet);
               }
               clock_end = clock();
               cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C1"<< endl;


               // test C2
               num_iter = 1000;
               clock_begin = clock();
               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M],EnergyCorrelator::slow);
                  C2(myJet);
               }
               clock_end = clock();
               cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C2"<< endl;

               num_iter = 10000;
               clock_begin = clock();
               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M],EnergyCorrelator::storage_array);
                  C2(myJet);
               }
               clock_end = clock();
               cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C2"<< endl;

               // test C3
               num_iter = 100;
               clock_begin = clock();

               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M],EnergyCorrelator::slow);
                  C3(myJet);
               }
               clock_end = clock();
               cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C3"<< endl;

               num_iter = 3000;
               clock_begin = clock();
               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M],EnergyCorrelator::storage_array);
                  C3(myJet);
               }
               clock_end = clock();
               cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C3"<< endl;

              // test C4
               num_iter = 10;
               clock_begin = clock();

               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M],EnergyCorrelator::slow);
                  C4(myJet);
               }
               clock_end = clock();
               cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C4"<< endl;

               num_iter = 300;
               clock_begin = clock();
               for (int t = 0; t < num_iter; t++) {
                  EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M],EnergyCorrelator::storage_array);
                  C4(myJet);
               }
               clock_end = clock();
               cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C4"<< endl;

            }
         }
      }
   }
}



