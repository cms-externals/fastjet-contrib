// To run this example, use the following command:
//
//   ./example-GHS [nevmax [njetmax]] < ../data/pythia8_Zq_vshort.dat
//
//----------------------------------------------------------------------
// $Id: example-GHS.cc 1507 2025-01-30 10:33:36Z gstagn $
//
// Copyright (c) 2025, Rhorry Gauld, Alexander Huss, Giovanni Stagnitto
//
// based on initial version by Fabrizio Caola, Radoslaw Grabarczyk,
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
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
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/GHSAlgo.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(int iargc, char **argv){

  // give user control over printout (mainly relevant for make check)
  // usage: "./example [nevmax [njetmax]] < data/pythia8_Zq_vshort.dat"
  unsigned int nevmax = 2;
  unsigned int njetmax = 1;
  if (iargc > 1) nevmax  = stoi(argv[1]);
  if (iargc > 2) njetmax = stoi(argv[2]);

  // print banner for FastJet at the start, so it doesn't mix
  // into the other output
  ClusterSequence::print_banner();

  // we start with a base jet definition (here antikt)
  double R = 0.4;
  JetDefinition base_jet_def(antikt_algorithm, R);
  // enable it to track flavours (default is net flavour)
  FlavRecombiner flav_recombiner;
  base_jet_def.set_recombiner(&flav_recombiner);

  // GHS parameters:
  double GHS_alpha = 1.0; // < flav-kt distance parameter alpha
  double GHS_omega = 2.0; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)

  double ptcut = 15.0; // < overall ptcut

  cout << "! this analysis considers only b-flavour from input events !" << endl;
  cout << "base jet definition: " << base_jet_def.description() << endl;
  cout << "GHS jet definition:  " << endl
       << " alpha = "  << GHS_alpha << endl
       << " omega = "  << GHS_omega << endl
       << " (ptcut = " << ptcut << " GeV)" << endl;

  // loop over some number of events
  unsigned int n_events = 10;
  for (unsigned int iev = 0; iev < n_events && iev < nevmax; iev++) {

    // read in input particles: see that routine for info
    // on how to set up the PseudoJets with flavour information
    vector<PseudoJet> event;
    read_event(event);
    cout << "\n#---------------------------------------------------------------\n";
    cout << "# read event " << iev << " with " << event.size() << " particles" << endl;

    // run the jet clustering with the base jet definition
    vector<PseudoJet> base_jets = base_jet_def(event);

    // run the GHS algorithm: require base jets & a hardness cut that should be chosen to match the fiducial selection
    vector<PseudoJet> GHS_jets  = run_GHS(base_jets, ptcut, GHS_alpha, GHS_omega, flav_recombiner);

    // make sure the sizes are the same (after the ptcut)
    Selector select_pt = SelectorPtMin(ptcut);
    vector<PseudoJet> selected_base_jets = select_pt(base_jets);

    assert(selected_base_jets.size() == GHS_jets.size());

    // ----------------------------------------------------
    // loop over the two leading jets and print out their properties
    for (unsigned int ijet = 0; ijet < selected_base_jets.size() && ijet < njetmax; ijet++) {
      // first print out the original anti-kt jets and the GHS jets
      const auto & base_jet = selected_base_jets[ijet];
      const auto & GHS_jet  = GHS_jets [ijet];
      cout << endl;
      cout << "base jet " << ijet << ": ";
      cout << "pt=" << base_jet.pt() << " rap=" << base_jet.rap() << " phi=" << base_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(base_jet).description() << endl;
      cout << "GHS jet  " << ijet << ": ";
      cout << "pt=" << GHS_jet.pt() << " rap=" << GHS_jet.rap() << " phi=" << GHS_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(GHS_jet).description() << endl;

      // for the first event, print out the jet constituents' pt and initial and final flavours
      cout << "constituents:" << endl;
      for (const auto & c: sorted_by_pt(base_jet.constituents())) {
        cout << "  pt = " << setw(10) << c.pt();
        cout << ", orig. flav = " << setw(8) << FlavHistory::initial_flavour_of(c).description();
        cout << ", final flav = " << setw(8) << FlavHistory::current_flavour_of(c).description();
        cout << endl;
      }
    }
  }

  return 0;
}

// read in input particles and set up PseudoJets with flavour information
void read_event(vector<PseudoJet> &event){
    // read in the input particles and their PDG IDs
    string line;
    double px, py, pz, E;
    int    pdg_id;
    event.resize(0);
    while(getline(cin,line)) {
      if(line[0] == '#') continue;

      istringstream iss(line);
      iss >> px >> py >> pz >> E >> pdg_id;
      // create a fastjet::PseudoJet with these components and put it onto
      // back of the input_particles vector
      PseudoJet p(px,py,pz,E);

      // assign information about flavour (will be deleted automatically)
      FlavInfo flav_info_init(pdg_id);
      // we take only b-quarks
      flav_info_init.reset_all_but_flav(5);
      // need to cast to const for constructor
      p.set_user_info(new FlavHistory(const_cast<const FlavInfo&>(flav_info_init)));
      event.push_back(p);

      if (cin.peek() == '\n' || cin.peek() == EOF) {
        getline(cin,line);
        break;
      }
    }
}
