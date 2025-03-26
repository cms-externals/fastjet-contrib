//
//----------------------------------------------------------------------
// Example on how to estimate background (rho) using the SignalFreeBackgroundEstimator. This background estimate is then used within the Iterative Constituent Subtraction.
//
// run it with
//  ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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


#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "SignalFreeBackgroundEstimator.hh"
//#include "fastjet/contrib/IterativeConstituentSubtractor.hh"

#include "functions.hh"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  // This should be set up before event loop:
  double max_eta=4;   // specify the maximal pseudorapidity for the input particles. It is used for the subtraction. Particles with eta>|max_eta| are removed and not  used during the subtraction (they are not returned). The same parameter should be used for the SignalFreeBackgroundEstimator as it is demonstrated in this example.
  double max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not important for the subtraction.

  // Signal-free background estimator
  fastjet::contrib::SignalFreeBackgroundEstimator bge_rho(max_eta,0.55); // maximal pseudo-rapidity cut is used inside ConstituentSubtraction, but in SignalFreeBackgroundEstimator, the range is specified by maximal rapidity cut. Therefore, it is important to apply the same pseudo-rapidity cut also for particles used for background estimation (using function "set_particles") and also derive the rho dependence on rapidity using this max pseudo-rapidity cut to get the correct rescaling function!
  bge_rho.set_signal_seed_parameters(0.3,0.4);  // Set the parameters for signal seeds (the first is the distance parameter for anti-kt clustering of signal particles into jets, and the second is the deltaR distance around signal jets to exclude areas with signal).
  bge_rho.set_jet_rho_min(5,8);    // Set the parameters for the minimal jet rho when selecting jets from the signal jets (intercept and slope of the linear function as a function of some measure of pileup specified in the set_particles function)
  bge_rho.set_jet_rho_min_charged(20); // If tracking info available, set the minimal jet rho when selecting jets obtained from charged signal particles
  bge_rho.set_window_parameters(0.5,0.4,0.1);   // Set parameters for the window to be used to compute the weighted median 

  /*  Commenting event-wide pileup corrector (ICS) in this example, because currently no other contribs can be used.
  contrib::IterativeConstituentSubtractor subtractor;
  subtractor.set_max_eta(max_eta); // parameter for maximal |eta| cut. It is specified above.
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  vector<double> max_distances;
  max_distances.push_back(0.15);
  max_distances.push_back(0.2);
  vector<double> alphas;
  alphas.push_back(1);
  alphas.push_back(1);
  subtractor.set_parameters(max_distances,alphas); // in this example, 2 CS corrections will be performed: 1. correction with max_distance of 0.1 and alpha of 0, 2. correction with max_distance of 0.15 and alpha of 0. After the first correction, the scalar sum of pt from remaining ghosts is evaluated and redistributed uniformly accross the event.
  subtractor.set_ghost_removal(true);  // set to true if the ghosts (proxies) which were not used in the previous CS procedure should be removed for the next CS procedure
  subtractor.set_ghost_area(0.004); // parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_background_estimator(&bge_rho); // specify the background estimator to estimate rho. You can use also specify background estimator for the mass term as an additional parameter
  subtractor.initialize();  // this is new compared to previous usages of ConstituentSubtractor! It should be used after specifying all the parameters and before event loop.
  // print info (optional)
  cout << subtractor.description() << endl;
  */

  // Grid-median background estimator and alternative ICS method to obtain input to the signal-free background estimator
  GridMedianBackgroundEstimator bge_rho_grid(max_eta,0.5);
  /*  contrib::IterativeConstituentSubtractor subtractor_to_get_signal_seeds;
  subtractor_to_get_signal_seeds.set_parameters(max_distances,alphas);  // can be arbitrary
  subtractor_to_get_signal_seeds.set_ghost_removal(true);
  subtractor_to_get_signal_seeds.set_max_eta(max_eta);
  subtractor_to_get_signal_seeds.set_background_estimator(&bge_rho_grid);*/


  // in event loop
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  vector<PseudoJet> *hard_event_charged=new vector<PseudoJet>;
  vector<PseudoJet> *background_event_charged=new vector<PseudoJet>;
  read_event(hard_event, full_event, hard_event_charged, background_event_charged);

  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);
  *hard_event_charged = SelectorAbsEtaMax(max_eta)(*hard_event_charged);
  *background_event_charged = SelectorAbsEtaMax(max_eta)(*background_event_charged);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with pseudo-rapidity |eta|<4" << endl;

  // jet definition and selector for jets
  JetDefinition jet_def(antikt_algorithm, 0.7);
  Selector sel_jets = SelectorNHardest(3) * SelectorAbsEtaMax(max_eta_jet);

  // clustering of the hard-scatter event and uncorrected event
  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());


  // There are several options how to define signal seeds for the signal-free background estimator:
  // 1. Perform an event-wide pileup correction with your preferred method to get estimated signal particles. No tracking information is needed. For example, you can use ICS with grid-median background estimator:
  //bge_rho_grid.set_particles(full_event);
  //vector<PseudoJet> estimated_signal=subtractor_to_get_signal_seeds.subtract_event(full_event);
  //fastjet::PseudoJet jet_zeroEta(5,0,0,10); // px,py,pz,E
  //double estimated_rho=bge_rho_grid.rho(jet_zeroEta);

  //bge_rho.set_particles(full_event,estimated_signal,estimated_rho);
  //vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);

  // 2. An additional improvement can be done by doing an additional iteration: The corrected event from option 1 can be used as input to the signal-free background estimator:
  //  bge_rho_grid.set_particles(full_event);
  //  vector<PseudoJet> estimated_signal=subtractor_to_get_signal_seeds.subtract_event(full_event);
  //  fastjet::PseudoJet jet_zeroEta(5,0,0,10); // px,py,pz,E
  //  double estimated_rho=bge_rho_grid.rho(jet_zeroEta);
  //  bge_rho.set_particles(full_event,estimated_signal,estimated_rho); // first iteration
  //  vector<PseudoJet> estimated_signal2=subtractor.subtract_event(full_event);
  //  double estimated_rho2=bge_rho.rho(jet_zeroEta);

  //  bge_rho.set_particles(full_event,estimated_signal2,estimated_rho2);  // second iteration
  //  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);


  // 3. If tracking information is available, then the charged signal particles can be used to find the signal seeds, without the need to do any event-wide pileup correction in advance.
  vector<PseudoJet> estimated_signal;   // empty vector should be provided in this case
  bge_rho.set_particles(full_event,estimated_signal,-1,*hard_event_charged);  // hard_event_charged can be obtained from Charged Hadron Subtraction (not done in this example, instead assiming prefect knowledge about chargen particles)
  cout << "obtained rho with SignalFreeBackgroundEstimator: " << bge_rho.rho()  << endl;
  // vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);  // use this to get the corrected event using ICS

  // 4. The strategy to find signal seeds described in the above steps 1, 2, and 3 can be arbitrarily combined. For example, the combination of steps 1 and 3 can be done as follows:
  // bge_rho_grid.set_particles(full_event);   // Alternatively, one can use only the neutral component, i.e. perform Charged Hadron Subtraction in advance, and then use the output here and on the next line.
  // vector<PseudoJet> estimated_signal=subtractor_to_get_signal_seeds.subtract_event(full_event);
  // fastjet::PseudoJet jet_zeroEta(5,0,0,10); // px,py,pz,E
  // double estimated_rho=bge_rho_grid.rho(jet_zeroEta);

  //bge_rho.set_particles(full_event,estimated_signal,estimated_rho,*hard_event_charged);
  //vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);

  // 5. The user can provide his/her own set of signal seeds:
  // vector<PseudoJet> additional_seeds=....; // the user obtaines the signal seeds based on his/her preferred method
  // bge_rho.add_seeds_from_user(additional_seeds);  // must be done before calling "set_particles"
  // bge_rho.set_particles(full_event);   // no pileup-corrected particles or charged signal particles are provided here. Alternatively, the user can also provide them, and then a combination of signal seeds will be used.

  // 6. It is also possible to not provide any input for signal seeds. Then no areas will be excluded, and the algortithm is similar to GridMedianBackgroundEstimator, except the fact that the size and the position of the window to compute rho can be specified using the function SignalFreeBackgroundEstimator::set_window_parameters. In this case use:
  //  bge_rho.set_particles(full_event);
  //  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);


  /*   Clustering and jet variable printing (in case pileup correction is applied)
  // clustering of the corrected event
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
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
  //  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
  cout << setprecision(4);
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
      //	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
      //	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
  for (unsigned int i=0; i<corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];

    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
      //	 << ", width = " << width(jet) << endl;
	 }
  cout << endl;
  */
  return 0;
}



