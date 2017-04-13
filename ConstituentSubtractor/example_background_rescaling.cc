//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_rescaling_TH1 < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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
#include "RescalingClasses.hh"
//#include "TH1D.h"
//#include "TF1.h"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  ///**** rescaling in heavy ion events: ****
  // Set the seven parameters for rescaling using function
  //  BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2)
  /// which is parametrized as
  ///  f(y,phi) = phi_term(phi) * rap_term(y)                                                                                                           
  /// where                                                                                                                                             
  ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))                                     
  ///  with four parameters v2, v3, v4, and psi.                                                                                                        
  /// rap_term(y) = a1*exp(-pow(y,2)/(2*sigma1^2)) + a2*exp(-pow(y,2)/(2*sigma2^2))         
  ///  with four parameters sigma1, sigma2, 1a, and a2.                                     
  ///
  /// You need to set the parameters event-by-event. Example:
  contrib::BackgroundRescalingYPhi rescaling(0.1,0.1,0.001,0,1,5,0,10);
  rescaling.use_rap_term(false);    // set to true, if you have derived also the rapidity terms for the rescaling



  ///**** rescaling using rapidity dependence stored in root TH1 object: ****
  // find the rapidity distribution of pileup particles from minimum bias events in a separate run. Fill a root TH1 histogram with this distribution.
  // Here as an example where the TH1 histogram is filled with a random distribution. Do not use this!!! Fill it with all particles (topoclusters) weighted by their pt using minimum bias events.
  /* TH1D* hist=new TH1D("hist","",100,-5,5);
  TF1 polynom("polynom","1+0.03*x*x-0.003*x*x*x*x",-5,5);
  for (int i=1;i<hist->GetNbinsX()+1;i++){
    hist->SetBinContent(i,polynom(hist->GetBinCenter(i)));
  }

  contrib::BackgroundRescalingYFromRoot<TH1D> rescaling(hist);
*/




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
 
  // setting the TH1 rescaling:
  bge_rho.set_rescaling_class(&rescaling);

  bge_rho.set_particles(full_event);
  subtractor.set_background_estimator(&bge_rho);

  // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
  subtractor.set_common_bge_for_rho_and_rhom(true); // for massless input particles it does not make any difference (rho_m is always zero)
  cout << subtractor.description() << endl;


  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event,maxEta);
  ClusterSequenceArea clust_seq_corr(corrected_event, jet_def, area_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
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



