//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2013
//  David Krohn, Jesse Thaler, and Lian-Tao Wang
//
//  $Id: example.cc 599 2014-04-18 14:29:54Z jthaler $
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

#include <iostream>
#include <sstream>
#include <stdio.h>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "VariableRPlugin.hh" // In external code, this should be fastjet/contrib/VariableRPlugin.hh
//#include "VariableR.hh" // This header is still available for backwards compatibility.

using namespace std;
using namespace fastjet;
using namespace contrib;

void print_jets (const fastjet::ClusterSequence &,
                 const vector<fastjet::PseudoJet> &);

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){
   
   //----------------------------------------------------------
   // read in input particles
   vector<PseudoJet> event;
   read_event(event);
   cout << "# read an event with " << event.size() << " particles" << endl;
   
   //----------------------------------------------------------
   // illustrate how VariableR contrib works
   // anti-kT variable R
   //----------------------------------------------------------
   
   // defining parameters
   double rho(2000.), min_r(0.0), max_r(2.0),ptmin(5.0);
   
   VariableRPlugin lvjet_pluginAKT(rho, min_r, max_r, VariableRPlugin::AKTLIKE);
   fastjet::JetDefinition jet_defAKT(&lvjet_pluginAKT);
   fastjet::ClusterSequence clust_seqAKT(event, jet_defAKT);
   
   // tell the user what was done
   cout << "Ran " << jet_defAKT.description() << endl;
   
   // extract the inclusive jets with pt > 5 GeV
   vector<fastjet::PseudoJet> inclusive_jetsAKT = clust_seqAKT.inclusive_jets(ptmin);
   
   // print them out
   cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
   cout << "---------------------------------------\n";
   print_jets(clust_seqAKT, inclusive_jetsAKT);
   cout << endl;
   
   //----------------------------------------------------------
   // Same for Cambridge-Aachen variable R
   //----------------------------------------------------------
   
   VariableRPlugin lvjet_pluginCA(rho, min_r, max_r, VariableRPlugin::CALIKE);
   fastjet::JetDefinition jet_defCA(&lvjet_pluginCA);
   fastjet::ClusterSequence clust_seqCA(event, jet_defCA);
   
   // tell the user what was done
   cout << "Ran " << jet_defCA.description() << endl;
   
   // extract the inclusive jets with pt > 5 GeV
   vector<fastjet::PseudoJet> inclusive_jetsCA = clust_seqCA.inclusive_jets(ptmin);
   
   // print them out
   cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
   cout << "---------------------------------------\n";
   print_jets(clust_seqCA, inclusive_jetsCA);
   cout << endl;
   
   //----------------------------------------------------------
   // Illustrating preclustering feature new in v1.1
   // Need to have minimum jet radius for preclustering to make sense
   //----------------------------------------------------------
   
   min_r = 0.4;  // change small radius to allow for preclustering
   VariableRPlugin lvjet_pluginAKT_precluster(rho, min_r, max_r, VariableRPlugin::AKTLIKE,true);
   fastjet::JetDefinition jet_defAKT_precluster(&lvjet_pluginAKT_precluster);
   fastjet::ClusterSequence clust_seqAKT_precluster(event, jet_defAKT_precluster);

   // tell the user what was done
   cout << "Ran " << jet_defAKT_precluster.description() << endl;
   
   // extract the inclusive jets with pt > 5 GeV
   vector<fastjet::PseudoJet> inclusive_jetsAKT_precluster = clust_seqAKT_precluster.inclusive_jets(ptmin);
   
   // print them out
   cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
   cout << "---------------------------------------\n";
   print_jets(clust_seqAKT_precluster, inclusive_jetsAKT_precluster);
   cout << endl;

   //----------------------------------------------------------
   // As a cross check, same as above, but with no preclustering
   // (Results should be nearly identical)
   //----------------------------------------------------------
   
   VariableRPlugin lvjet_pluginAKT_noprecluster(rho, min_r, max_r, VariableRPlugin::AKTLIKE,false);
   fastjet::JetDefinition jet_defAKT_noprecluster(&lvjet_pluginAKT_noprecluster);
   fastjet::ClusterSequence clust_seqAKT_noprecluster(event, jet_defAKT_noprecluster);
   
   // tell the user what was done
   cout << "Ran " << jet_defAKT_noprecluster.description() << endl;
   
   // extract the inclusive jets with pt > 5 GeV
   vector<fastjet::PseudoJet> inclusive_jetsAKT_noprecluster = clust_seqAKT_noprecluster.inclusive_jets(ptmin);
   
   // print them out
   cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
   cout << "---------------------------------------\n";
   print_jets(clust_seqAKT_noprecluster, inclusive_jetsAKT_noprecluster);
   cout << endl;

   
   
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

//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const fastjet::ClusterSequence & clust_seq,
                 const vector<fastjet::PseudoJet> & jets) {
   
   // sort jets into increasing pt
   vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);
   
   // label the columns
   printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity",
          "phi", "pt","m","e", "n constituents");
   
   // print out the details for each jet
   for (unsigned int i = 0; i < sorted_jets.size(); i++) {
      int n_constituents = clust_seq.constituents(sorted_jets[i]).size();
      printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %8u\n",
             i, sorted_jets[i].rap(), sorted_jets[i].phi(),
             sorted_jets[i].perp(),sorted_jets[i].m(),sorted_jets[i].e(), n_constituents);
   }
}
