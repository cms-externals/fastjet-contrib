// To run this example, use the following command:
//
//   ./example-CMP < ../data/pythia8_Zq_vshort.dat
//
// NB: the example file reads in a file with 6 light flavours
//     (though for this algorithm we need to strip off all but one), and
//     an extra high density of quarks, to help with "make check"
//     test things more thoroughly
//----------------------------------------------------------------------
// $Id$
//
// Copyright (c) 2024, Michal Czakon, Alexander Mitov, and Rene Poncelet
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
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/CMPPlugin.hh"

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

  // CMP parameters:
  // CMP 'a' parameter in
  //   kappa_ij = 1/a * (kT_i^2 + kT_j^2) / (2*kT_max^2)
  double CMP_a = 0.1;
  // correction to original CMP algo: do not change this if you want IRC safety!
  CMPPlugin::CorrectionType CMP_corr =
    CMPPlugin::CorrectionType::SqrtCoshyCosPhiArgument_a2;
  // Dynamic definition of ktmax
  CMPPlugin::ClusteringType CMP_clust = CMPPlugin::ClusteringType::DynamicKtMax;
  // CMP plugin
  JetDefinition CMP_jet_def(new CMPPlugin(R, CMP_a, CMP_corr, CMP_clust));

  // enable it to track flavours (default is net flavour)
  CMP_jet_def.set_recombiner(&flav_recombiner);
  CMP_jet_def.delete_plugin_when_unused();

  cout << "! this analysis considers only b-flavour from input events !" << endl;
  cout << "base jet definition: " << base_jet_def.description() << endl;
  cout << "CMP jet definition:  " << CMP_jet_def.description() << endl;

  // loop over some number of events
  unsigned int n_events = 10;
  for (unsigned int iev = 0; iev < n_events && iev < nevmax; iev++) {

    // read in input particles: see that routine for info 
    // on how to set up the PseudoJets with flavour information
    vector<PseudoJet> event;
    read_event(event);
    cout << "\n#---------------------------------------------------------------\n";
    cout << "# read event " << iev << " with " << event.size() << " particles" << endl;

    // run the jet clustering with the base jet definition and the
    // CMPPlugin-based jet definition
    vector<PseudoJet> base_jets = base_jet_def(event);
    vector<PseudoJet> CMP_jets  = CMP_jet_def(event);

    // ----------------------------------------------------
    // loop over the two leading jets and print out their properties
    for (unsigned int ijet = 0; ijet < base_jets.size() && ijet < njetmax; ijet++) {
      // first print out the original anti-kt jets
      const auto & base_jet = base_jets[ijet];
      cout << endl;
      cout << "base jet " << ijet << ": ";
      cout << "pt=" << base_jet.pt() << " rap=" << base_jet.rap() << " phi=" << base_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(base_jet).description() << endl;

      // for the first event, print out the jet constituents' pt and initial and final flavours
      cout << "constituents:" << endl;
      for (const auto & c: sorted_by_pt(base_jet.constituents())) {
        cout << "  pt = " << setw(10) << c.pt();
        cout << ", orig. flav = " << setw(8) << FlavHistory::initial_flavour_of(c).description();
        cout << ", final flav = " << setw(8) << FlavHistory::current_flavour_of(c).description();          
        cout << endl;
      }
    }
    // ----------------------------------------------------
    // loop over the two leading jets and print out their properties
    for (unsigned int ijet = 0; ijet < CMP_jets.size() && ijet < njetmax; ijet++) {
      // and the CMP jets
      const auto & CMP_jet  = CMP_jets [ijet];
      cout << endl;
      cout << "CMP jet  " << ijet << ": ";
      cout << "pt=" << CMP_jet.pt() << " rap=" << CMP_jet.rap() << " phi=" << CMP_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(CMP_jet).description() << endl;
      
      // for the first event, print out the jet constituents' pt and initial and final flavours
      cout << "constituents:" << endl;
      for (const auto & c: sorted_by_pt(CMP_jet.constituents())) {
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
