// To run this example, use the following command:
//
//   ./example < data/pythia8_Zq_vshort.dat
//
// NB: the example file reads in a file with 6 light flavours, and
//     an extra high density of quarks, to help with "make check"
//     test things more thoroughly
//----------------------------------------------------------------------
// $Id: example-IFN.cc 1473 2024-12-12 14:36:33Z salam $
//
// Copyright (c) 2023, Fabrizio Caola, Radoslaw Grabarczyk, 
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
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/IFNPlugin.hh" // In external code, this may become fastjet/contrib/IFNPlugin.hh

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

  // we start with a base jet definition (should be either
  // antikt_algorithm or cambridge_algorithm, or their e+e- variants)
  JetDefinition base_jet_def(antikt_algorithm, 0.4);
  // enable it to track flavours (default is net flavour)
  FlavRecombiner flav_recombiner;
  base_jet_def.set_recombiner(&flav_recombiner);

  // And then we set up the IFNPlugin that builds on the base_jet_def
  // The main free parameter, alpha, in the uij distance, 
  //   uij = max(pt_i, pt_j)^alpha min(pt_i, pt_j)^(2-alpha) Omega_ij
  double alpha = 2.0;

  // The parameter that sets the nature of the Omega rapidity term;
  // only change the default of 3-alpha if you are sure you know what you are doing
  double omega = 3.0 - alpha;

  // The flavour summation scheme; should be one of 
  //   - FlavRecombiner::net
  //   - FlavRecombiner::modulo_2
  FlavRecombiner::FlavSummation flav_summation = FlavRecombiner::net;

  // then construct the IFNPlugin jet definition
  auto ifn_plugin = new IFNPlugin(base_jet_def, alpha, omega, flav_summation);
  JetDefinition IFN_jet_def(ifn_plugin);
  IFN_jet_def.delete_plugin_when_unused();

  cout << "base jet definition: " << base_jet_def.description() << endl;
  cout << "IFN jet definition:  " << IFN_jet_def.description() << endl;

  // loop over some number of events
  int n_events = 10;
  for (int iev = 0; iev < n_events && iev < nevmax; iev++) {

    // read in input particles: see that routine for info 
    // on how to set up the PseudoJets with flavour information
    vector<PseudoJet> event;
    read_event(event);
    cout << "\n#---------------------------------------------------------------\n";
    cout << "# read event " << iev << " with " << event.size() << " particles" << endl;

    // run the jet clustering with the base jet definition and the
    // IFNPlugin-based jet definition
    vector<PseudoJet> base_jets = base_jet_def(event);
    vector<PseudoJet> IFN_jets  = IFN_jet_def(event);

    // make sure the sizes are the same
    assert(base_jets.size() == IFN_jets.size());

    // ----------------------------------------------------
    // loop over the two leading jets and print out their properties
    for (unsigned int ijet = 0; ijet < base_jets.size() && ijet < njetmax; ijet++) {
      // first print out the original anti-kt jets and the IFN jets
      const auto & base_jet = base_jets[ijet];
      const auto & IFN_jet  = IFN_jets [ijet];
      cout << endl;
      cout << "base jet " << ijet << ": ";
      cout << "pt=" << base_jet.pt() << " rap=" << base_jet.rap() << " phi=" << base_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(base_jet).description() << endl;
      cout << "IFN jet  " << ijet << ": ";
      cout << "pt=" << IFN_jet.pt() << " rap=" << IFN_jet.rap() << " phi=" << IFN_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(IFN_jet).description() << endl;
      
      // for the first event, print out the jet constituents' pt and initial and final flavours
      cout << "constituents:" << endl;
      for (const auto & c: sorted_by_pt(IFN_jet.constituents())) {
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
      p.set_user_info(new FlavHistory(pdg_id));
      event.push_back(p);

      if (cin.peek() == '\n' || cin.peek() == EOF) {
        getline(cin,line);
        break;
      }
    }
}
