// $Id: example_whole_event.cc 995 2017-01-19 21:59:55Z berta $
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
  contrib::ConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  subtractor.set_max_distance(1); // free parameter for the maximal allowed distance between particle i and ghost k
  subtractor.set_alpha(2);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
  subtractor.set_ghost_area(0.01); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.

  // event loop
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  double maxEta=4;  // specify the maximal rapidity for the particles used in the subtraction
  double maxEta_jet=3; // the maximal rapidity for selected jets

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsEtaMax(maxEta)(hard_event);
  full_event = SelectorAbsEtaMax(maxEta)(full_event);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with pseudo-rapidity |eta|<4" << endl;

 
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

  JetMedianBackgroundEstimator bge_rho(rho_range, jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity=new BackgroundJetScalarPtDensity();
  bge_rho.set_jet_density_class(scalarPtDensity); // this changes the computation of pt of patches from vector sum to scalar sum. The scalar sum seems more reasonable.
  bge_rho.set_particles(full_event);
  subtractor.set_background_estimator(&bge_rho);

  // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
  subtractor.set_common_bge_for_rho_and_rhom(true); // for massless input particles it does not make any difference (rho_m is always zero)
  cout << subtractor.description() << endl;


  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event,maxEta);
  ClusterSequenceArea clust_seq_corr(corrected_event, jet_def, area_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());


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

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
  cout << setprecision(4);
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
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



