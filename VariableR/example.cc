//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2013
//  David Krohn, Jesse Thaler, and Lian-Tao Wang
//
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

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "VariableR.hh" // In external code, this should be fastjet/contrib/VariableR.hh

using namespace std;
using namespace fastjet;
using namespace contrib;

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

  AKTVR lvjet_pluginAKT(2000,2.0);
  fastjet::JetDefinition jet_defAKT(&lvjet_pluginAKT);
  fastjet::ClusterSequence clust_seqAKT(event, jet_defAKT);
  
  // tell the user what was done
  cout << "Ran " << jet_defAKT.description() << endl;

  // extract the inclusive jets with pt > 5 GeV
  double ptminAKT = 5.0;
  vector<fastjet::PseudoJet> inclusive_jetsAKT = clust_seqAKT.inclusive_jets(ptminAKT);
  
  // print them out
  cout << "Printing inclusive jets with pt > "<< ptminAKT <<" GeV\n";
  cout << "---------------------------------------\n";
  fastjet::contrib::print_jets(clust_seqAKT, inclusive_jetsAKT);
  cout << endl;

  //----------------------------------------------------------
  // illustrate how VariableR contrib works
  // Cambridge-Aachen variable R
  //----------------------------------------------------------

  CAVR lvjet_pluginCA(2000,2.0);
  fastjet::JetDefinition jet_defCA(&lvjet_pluginCA);
  fastjet::ClusterSequence clust_seqCA(event, jet_defCA);
  
  // tell the user what was done
  cout << "Ran " << jet_defCA.description() << endl;

  // extract the inclusive jets with pt > 5 GeV
  double ptminCA = 5.0;
  vector<fastjet::PseudoJet> inclusive_jetsCA = clust_seqCA.inclusive_jets(ptminCA);
  
  // print them out
  cout << "Printing inclusive jets with pt > "<< ptminCA <<" GeV\n";
  cout << "---------------------------------------\n";
  fastjet::contrib::print_jets(clust_seqCA, inclusive_jetsCA);
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
