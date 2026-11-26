// To run this example, use the following command:
//
//  ./example < ../data/single-epDIS-event.dat
//
//----------------------------------------------------------------------
// $Id: example.cc 1556 2026-06-03 09:43:28Z silviaferrarioravasio $
//
// Copyright (c) 2026-, Melissa van Beekveld, Silvia Ferrario Ravasio, 
// Alexander Karlberg and Darcy Peake.
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
#include "fastjet/contrib/DISGenkt.hh" 

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void output_clust(const vector<PseudoJet> & jets, contrib::DISGenktPlugin * plugin);

//----------------------------------------------------------------------
int main(){
  // read in input particles -- they must be in the Breit Frame
  vector<PseudoJet> event;
  read_event(event);
  cout << "# Read an event with " << event.size() << " particles" << endl;
  
  // set the sign of the (proton) beam
  int beam_sign = 1;

  // set the power p
  double p = 0.; //< CA variant
  // choose the jet radius parameter
  double R = M_PI/2.;

  // run jet clustering
  contrib::DISGenktPlugin * disjets_plugin = new contrib::DISGenktPlugin(p, beam_sign, R);
  JetDefinition jet_def_p0(disjets_plugin);
  ClusterSequence clust(event, jet_def_p0);

  // output
  cout << "# Output of clustering with p = 0 and R = pi/2" << endl;
  vector<PseudoJet> jets = clust.inclusive_jets();
  output_clust(jets, disjets_plugin);

  // redo clustering with p = 1 (kt variant)
  p = 1;
  disjets_plugin = new contrib::DISGenktPlugin(p, beam_sign, R);
  JetDefinition jet_def_p1(disjets_plugin);
  clust = ClusterSequence(event, jet_def_p1);
  cout << "# Output of clustering with p = 1 and R = pi/2" << endl;
  jets = clust.inclusive_jets();
  output_clust(jets, disjets_plugin);

  // redo clustering with p = -1 (anti-kt variant), R = Pi/2
  p = -1;
  disjets_plugin = new contrib::DISGenktPlugin(p, beam_sign, R);
  JetDefinition jet_def_pm1(disjets_plugin);
  clust = ClusterSequence(event, jet_def_pm1);
  cout << "# Output of clustering with p = -1 and R = pi/2" << endl;
  jets = clust.inclusive_jets();
  output_clust(jets, disjets_plugin);

  // A value for R <= Pi/3 is recommended to 
  // avoid clustering too many beam remnant
  // particles into the macrojet.
  // Redo clustering with p = -1 (anti-kt variant), R = Pi/3 
  R = M_PI/3.;
  disjets_plugin = new contrib::DISGenktPlugin(p, beam_sign, R);
  jet_def_pm1 = JetDefinition(disjets_plugin);
  clust       = ClusterSequence(event, jet_def_pm1);
  cout << "# Output of clustering with p = -1 and R = pi/3" << endl;
  jets = clust.inclusive_jets();
  output_clust(jets, disjets_plugin);

  // sorting the jets by their light-cone projection along the original struct quark ("zjet")
  vector<PseudoJet> jets_sorted = disjets_plugin->sorted_by_zjet(jets);
  cout << "# Output of clustering with p = -1 and R = pi/3, sorted by largest z-projection" << endl;
  output_clust(jets_sorted, disjets_plugin);


  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there is extra "pollution"
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

void output_clust(const vector<PseudoJet> & jets, contrib::DISGenktPlugin * plugin){
  // find the macrojet
  int idx_macro = plugin->find_idx_macrojet(jets);
  for (unsigned int i=0; i<jets.size(); i++){
    const PseudoJet &jet = jets[i];
    vector<fastjet::PseudoJet> constituents = jet.constituents();
    if(i == idx_macro) cout << "Macrojet with ";
    else cout << "jet with ";
    cout << "E = " << jet.e() << ", kT = " << sqrt(jet.px()*jet.px() + jet.py()*jet.py()) << ", rap = " << jet.rap() << ", n constituents " << constituents.size() << endl;
  }
  cout << endl;
}
