// $Id: example_whole_event_using_charged_info.cc 995 2017-01-19 21:59:55Z berta $
//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_whole_event < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//----------------------------------------------------------------------
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


#include "functions.hh"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  // set up before event loop:
  contrib::ConstituentSubtractor subtractor; // no need to provide background estimator in this case
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There  are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta  
  subtractor.set_max_distance(1); // free parameter for the maximal allowed distance between particle i and ghost k
  subtractor.set_alpha(2);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
  subtractor.set_ghost_area(0.01); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_do_mass_subtraction(true);   // specify if also the mass term sqrt(pT^2+m^2)-pT should be corrected or not 
  double CBS=1.0;  // choose the scale for scaling the background charged particles
  double CSS=1.0;  // choose the scale for scaling the signal charged particles


  // event loop
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  vector<PseudoJet> *hard_event_charged=new vector<PseudoJet>;
  vector<PseudoJet> *background_event_charged=new vector<PseudoJet>;

  read_event(hard_event, full_event, hard_event_charged, background_event_charged);

  double maxEta=4;  // specify the maximal rapidity for the particles used in the subtraction
  double maxEta_jet=3; // the maximal rapidity for selected jets

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsEtaMax(maxEta)(hard_event);
  full_event = SelectorAbsEtaMax(maxEta)(full_event);
  *hard_event_charged = SelectorAbsEtaMax(maxEta)(*hard_event_charged);
  *background_event_charged = SelectorAbsEtaMax(maxEta)(*background_event_charged);

  cout << "# read an event with " << hard_event.size() << " signal particles, " << full_event.size() - hard_event.size() << " background particles, " << hard_event_charged->size() << " signal charged particles, and " << background_event_charged->size() << " background charged particles" << " with pseudo-rapidity |eta|<4" << endl;

 
  // do the clustering with ghosts and get the jets
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.7);
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(maxEta,1)); // the area definiton is used only for the jet backgroud estimator. It is not important for the ConstituentSubtractor when subtracting the whole event - this is not true when subtracting the individual jets

  ClusterSequenceArea clust_seq_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea clust_seq_full(full_event, jet_def, area_def);

  Selector sel_jets = SelectorNHardest(3) * SelectorAbsRapMax(maxEta_jet);

  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

 // create what we need for the background estimation
  //----------------------------------------------------------
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsEtaMax(maxEta-0.4);
  //  ClusterSequenceArea clust_seq_rho(full_event, jet_def, area_def);  

  cout << subtractor.description() << endl;

  vector<PseudoJet> corrected_event=subtractor.subtract_event_using_charged_info(full_event,CBS,*background_event_charged,CSS,*hard_event_charged,maxEta);
  ios::fmtflags f( cout.flags() );
  cout << setprecision(4) << fixed;
  cout << endl << "Corrected particles in the whole event:" << endl;
  for (unsigned int i=0; i<corrected_event.size(); i++){
    const PseudoJet &particle = corrected_event[i];
    cout << "pt = " << particle.pt()
	 << ", phi = " << particle.phi()
	 << ", rap = " << particle.rap()
	 << ", |mass| = " << fabs(particle.m()) << endl;
  }
  cout << endl;


  ClusterSequenceArea clust_seq_corr(corrected_event, jet_def, area_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
  cout << setprecision(4);
  cout << "Original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "Unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "Subtracted full jets" << endl;
  for (unsigned int i=0; i<corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];

    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
	 }
  cout << endl;

  return 0;
}



