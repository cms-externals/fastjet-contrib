//--------------------------------------------------------------------
//
// file: example.cc
// Example program demonstrating how to run the DynamicR FastJet
// contrib plugin on input events and print the resulting jets.
// To run this example, use the following command:
//   make example
//   ./example < ../data/pythia8_Zq_vshort.dat
//
//--------------------------------------------------------------------
// $Id: example.cc 1567 2026-06-24 12:19:44Z tousik $
//
// Copyright (c) 2023-, Tousik Samui 
//
//--------------------------------------------------------------------
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
//--------------------------------------------------------------------
#include <iomanip>
#include <iostream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/DynamicR.hh"

// Header and Namespace.
using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);


int main() {

  // No. of events.
  int nEvent = 10 ;             // Number of events to process

  // FastJet parameters
  double radius = 0.5;          // Initial jet radius parameter
  double pTjetMin = 10.0;       // Minimum jet transverse momentum for output
  vector<PseudoJet> particles;  // Container holding input particles for one event

  // some useful variables
  double Rd;

  ClusterSequence::set_fastjet_banner_stream(NULL);
  // Loop over events 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    particles.resize(0);        // Clear the particle container for the next event

    //----------------------------------------------------------
    // read in input particles
    read_event(particles);      // Fill the particle vector with one event from standard input

    cout << "############################################################" << endl;
    cout << "Event # " << iEvent+1 << "." << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "# read an event with " << particles.size() << " particles" << endl;
    cout << "------------------------------------------------------------" << endl;

    // continue to next event if no final state particle.
    if (int(particles.size())==0) continue; 

    // DR-AK algorithm
    DRAK DRAKJP(radius);                                  // Construct the Dynamic-R anti-kt-like plugin
    JetDefinition jet_def_DRAK(&DRAKJP);                  // Build the FastJet jet definition from the plugin
    ClusterSequence cs_DRAK(particles, jet_def_DRAK);    // Run clustering on the event particles
    vector<PseudoJet> DRAKjets = sorted_by_pt(cs_DRAK.inclusive_jets(pTjetMin));  // Keep jets above pT threshold and sort by pT

    cout << "Output of " << jet_def_DRAK.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    // Loop over the clustered DR-AK jets and print their kinematics and dynamic radius
    for (unsigned int iD = 0; iD < DRAKjets.size(); iD++) {
      Rd = radius;
      if (DRAKjets[iD].has_user_info<DynamicRJetInfo>()) {
        Rd = DRAKjets[iD].user_info<DynamicRJetInfo>().Rd();   // Retrieve the jet's stored dynamic radius
      }

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRAKjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRAKjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-CA algorithm
    DRCA DRCAJP(radius);                                  // Construct the Dynamic-R Cambridge/Aachen-like plugin
    JetDefinition jet_def_DRCA(&DRCAJP);                 // Build the corresponding FastJet jet definition
    ClusterSequence cs_DRCA(particles, jet_def_DRCA);    // Run clustering on the same event particles
    vector<PseudoJet> DRCAjets = sorted_by_pt(cs_DRCA.inclusive_jets(pTjetMin));  // Keep jets above pT threshold and sort by pT

    cout << "Output of " << jet_def_DRCA.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    // Loop over the clustered DR-CA jets and print their kinematics and dynamic radius
    for (unsigned int iD = 0; iD < DRCAjets.size(); iD++) {
      Rd = radius;
      if (DRCAjets[iD].has_user_info<DynamicRJetInfo>()) {
        Rd = DRCAjets[iD].user_info<DynamicRJetInfo>().Rd();   // Retrieve the jet's stored dynamic radius
      }

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRCAjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRCAjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

    // DR-KT algorithm
    DRKT DRKTJP(radius);                                  // Construct the Dynamic-R kt-like plugin
    JetDefinition jet_def_DRKT(&DRKTJP);                 // Build the corresponding FastJet jet definition
    ClusterSequence cs_DRKT(particles, jet_def_DRKT);    // Run clustering on the same event particles
    vector<PseudoJet> DRKTjets = sorted_by_pt(cs_DRKT.inclusive_jets(pTjetMin));  // Keep jets above pT threshold and sort by pT

    cout << "Output of " << jet_def_DRKT.description() << "." << endl;
    cout << setw(5) << right << "jet #" << setw(10) << right << "rapidity" << setw(10) << right << "phi" << setw(15) << right << "pt (GeV)" << setw(10) << right << "Rd" << endl; 

    // Loop over the clustered DR-KT jets and print their kinematics and dynamic radius
    for (unsigned int iD = 0; iD < DRKTjets.size(); iD++) {
      Rd = radius;
      if (DRKTjets[iD].has_user_info<DynamicRJetInfo>()) {
        Rd = DRKTjets[iD].user_info<DynamicRJetInfo>().Rd();   // Retrieve the jet's stored dynamic radius
      }

      cout << setw(5) << right << iD+1;
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].rap();
      cout << setw(10) << right << fixed << setprecision(2) << DRKTjets[iD].phi();
      cout << setw(15) << right << fixed << setprecision(2) << DRKTjets[iD].pt();
      cout << setw(10) << right << fixed << setprecision(2) << Rd << endl;
    }
    cout << "------------------------------------------------------------" << endl;
    cout << endl;

  }

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){
  string line;   // Holds one line of event input at a time
  while (getline(cin, line)) {
    istringstream linestream(line);   // Stream used to parse momentum components from the line

    // take substrings to avoid problems when there is extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}   // End of current event marker
    if (line.substr(0,1) == "\n") {return;}     // Stop if an empty-line-like delimiter is found
    if (line.substr(0,1) == "") {return;}       // Stop if the line is empty
    if (line.substr(0,1) == "#") {continue;}    // Skip comment lines

    double px,py,pz,E;   // Particle four-momentum
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);   // Add the particle to the current event
  }
}

