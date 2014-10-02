// $Id: ConstituentSubtractor.cc 587 2014-04-06 13:44:24Z berta $
//
// ConstituentSubtractor package
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz
//
// Copyright (c) 2014-, Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner
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

#include "ConstituentSubtractor.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{


  /// Constructor that takes a pointer to a background estimator for
  /// rho and optionally a pointer to a background estimator for
  /// rho_m.
  ConstituentSubtractor::ConstituentSubtractor(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom, double alpha, double maxDeltaR) :
    _bge_rho(bge_rho), _bge_rhom(bge_rhom), _common_bge(false), 
    _externally_supplied_rho_rhom(false), _alpha(alpha), _max_deltaR(maxDeltaR){
    if (_max_deltaR>0) _use_max_deltaR=true;
    else _use_max_deltaR=false;
  }


  // Constructor that takes an externally supplied value for rho and, optionally, for rho_m.
  ConstituentSubtractor::ConstituentSubtractor(double rho, double rhom, double alpha, double maxDeltaR) :      
    _bge_rho(0), _bge_rhom(0), _common_bge(false), _rho(rho), _rhom(rhom), _externally_supplied_rho_rhom(true), _alpha(alpha), _max_deltaR(maxDeltaR) {
    if (_max_deltaR>0) _use_max_deltaR=true;
    else _use_max_deltaR=false;
    assert(_rho  >= 0);
    assert(_rhom >= 0);
  }




  // the action on a given jet
  fastjet::PseudoJet ConstituentSubtractor::result(const PseudoJet &jet) const{
    // make sure we have a BGE or a rho value
    if (!_bge_rho && !_externally_supplied_rho_rhom){
      throw Error("ConstituentSubtractor::result() constituent subtraction needs a BackgroundEstimator or a value for rho.");
    }


    //----------------------------------------------------------------------
    // sift ghosts and particles in the input jet
    std::vector<PseudoJet> particles, ghosts;
    SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    unsigned long nGhosts=ghosts.size();
    unsigned long nParticles=particles.size();


    //----------------------------------------------------------------------
    // determinating what to do with rho_m and storing the rho and rho_m values for each ghost. Since the estimated background densities, rho and rho_m, may be position-dependent, separate value of background densities is stored for each ghost. 
    bool use_rhom=false;
    std::vector<double> rho;
    std::vector<double> rhom;
    if (_externally_supplied_rho_rhom){
      use_rhom=true;
      for (unsigned int j=0;j<nGhosts; j++){
	rho.push_back(_rho);
	rhom.push_back(_rhom);
      }
    }
    else{
      for (unsigned int j=0;j<nGhosts; j++){
	rho.push_back(_bge_rho->rho(ghosts[j]));
      }
      if (_bge_rhom){
	use_rhom=true;
	for (unsigned int j=0;j<nGhosts; j++){
	  rhom.push_back( _bge_rhom->rho(ghosts[j]));
	}
      } else if (_common_bge){
	use_rhom=true;
	BackgroundJetPtMDensity deltaMass_density;
	JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
	const FunctionOfPseudoJet<double> * orig_density = jmbge->jet_density_class();
      
	jmbge->set_jet_density_class(&deltaMass_density);
	for (unsigned int j=0;j<nGhosts; j++){
	  rhom.push_back(jmbge->rho(ghosts[j]));
	}
	jmbge->set_jet_density_class(orig_density);
      }
    }

    //----------------------------------------------------------------------
    // computing and sorting the distances, deltaR
    double maxDeltaR_squared=pow(_max_deltaR,2); 
    double alpha_times_two=_alpha*2.;
    std::vector<std::pair<double,int> > deltaRs;  // the first element is deltaR, the second element is only the index in the vector used for sorting
    std::vector<int> particle_indices_unsorted;
    std::vector<int> ghost_indices_unsorted;
    for (unsigned int i=0;i<nParticles; i++){
      double pt_factor=1.;
      if (fabs(alpha_times_two)>1e-5) pt_factor=pow(particles[i].pt(),alpha_times_two);
      for (unsigned int j=0;j<nGhosts; j++){
	double deltaR_squared = ghosts[j].squared_distance(particles[i])*pt_factor;
	if (!_use_max_deltaR || deltaR_squared<=maxDeltaR_squared){
	  particle_indices_unsorted.push_back(i);
	  ghost_indices_unsorted.push_back(j);
	  int deltaRs_size=deltaRs.size();  // current position
	  deltaRs.push_back(std::make_pair(deltaR_squared,deltaRs_size));
	}
      }
    }
    std::sort(deltaRs.begin(),deltaRs.end(),ConstituentSubtractor::_function_used_for_sorting);
    unsigned long nStoredPairs=deltaRs.size();

    //----------------------------------------------------------------------
    // the iterative process. Here, only finding the fractions of pt or deltaM to be corrected. The actual correction of particles is done later.
    std::vector<double> ghosts_fraction_of_pt(nGhosts,1.);
    std::vector<double> particles_fraction_of_pt(nParticles,1.);
    std::vector<double> ghosts_fraction_of_mtMinusPt(nGhosts,1.);
    std::vector<double> particles_fraction_of_mtMinusPt(nParticles,1.);
    for (unsigned long iindices=0;iindices<nStoredPairs;++iindices){
      int particle_index=particle_indices_unsorted[deltaRs[iindices].second];
      int ghost_index=ghost_indices_unsorted[deltaRs[iindices].second];

      if (ghosts_fraction_of_pt[ghost_index]>0 && particles_fraction_of_pt[particle_index]>0){
	double ratio_pt=particles[particle_index].pt()*particles_fraction_of_pt[particle_index]/rho[ghost_index]/ghosts[ghost_index].area()/ghosts_fraction_of_pt[ghost_index];
	if (ratio_pt>1){
	  particles_fraction_of_pt[particle_index]*=1-1./ratio_pt;
	  ghosts_fraction_of_pt[ghost_index]=-1;
	}
	else {
	  ghosts_fraction_of_pt[ghost_index]*=1-ratio_pt;
	  particles_fraction_of_pt[particle_index]=-1;
	}
      }
      if (use_rhom && ghosts_fraction_of_mtMinusPt[ghost_index]>0 && particles_fraction_of_mtMinusPt[particle_index]>0){
	double ratio_mtMinusPt=(particles[particle_index].mt()-particles[particle_index].pt())*particles_fraction_of_mtMinusPt[particle_index]/rhom[ghost_index]/ghosts[ghost_index].area()/ghosts_fraction_of_mtMinusPt[ghost_index];
	if (ratio_mtMinusPt>1){
	  particles_fraction_of_mtMinusPt[particle_index]*=1-1./ratio_mtMinusPt;
	  ghosts_fraction_of_mtMinusPt[ghost_index]=-1;
	}
	else{
	  ghosts_fraction_of_mtMinusPt[ghost_index]*=1-ratio_mtMinusPt;
	  particles_fraction_of_mtMinusPt[particle_index]=-1;
	}
      }
    }

    //----------------------------------------------------------------------
    // do the actual correction for particles:
    std::vector<PseudoJet> subtracted_particles;
    for (unsigned int i=0;i<particles_fraction_of_pt.size(); i++){
      if (particles_fraction_of_pt[i]<=0) continue;  // particles with zero pt are not used (but particles with zero mtMinusPt are used)
      double rapidity=particles[i].rap();
      double azimuth=particles[i].phi();
      double subtracted_pt=0;
      if (particles_fraction_of_pt[i]>0) subtracted_pt=particles[i].pt()*particles_fraction_of_pt[i];
      double subtracted_mtMinusPt=0;
      if (particles_fraction_of_mtMinusPt[i]>0) subtracted_mtMinusPt=(particles[i].mt()-particles[i].pt())*particles_fraction_of_mtMinusPt[i];
      PseudoJet subtracted_const(subtracted_pt*cos(azimuth),subtracted_pt*sin(azimuth),(subtracted_pt+subtracted_mtMinusPt)*sinh(rapidity),(subtracted_pt+subtracted_mtMinusPt)*cosh(rapidity));
      subtracted_particles.push_back(subtracted_const);
    }
    fastjet::PseudoJet subtracted_jet=join(subtracted_particles);

    return subtracted_jet;
  }





  void ConstituentSubtractor::use_common_bge_for_rho_and_rhom(bool value){ 
    if (value){
      if (_bge_rhom)  // checking if a BackGroundEstimator for rho_m is already provided
	throw Error("ConstituentSubtractor::use_common_bge_for_rho_and_rhom() is not allowed in the presence of an existing background estimator for rho_m.");
      if (_externally_supplied_rho_rhom)   // if the rho is externally supplied, the rho_m cannot be estimated (it should be supplied externally as well)
	throw Error("ConstituentSubtractor::use_common_bge_for_rho_and_rhom() is not allowed when supplying externally the values for rho and rho_m.");
      JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
      if (!jmbge)    //  currently, only for JetMedianBackgroundEstimator is allowed to change the density class with set_jet_density_class (for GridMedianBackgroundEstimator not).
	throw Error("ConstituentSubtractor::use_common_bge_for_rho_and_rhom() is currently only allowed for background estimators of JetMedianBackgroundEstimator type.");
    }
    _common_bge=value;
  }



  std::string ConstituentSubtractor::description() const{
    std::ostringstream descr;
    if ( _externally_supplied_rho_rhom){
      descr << "ConstituentSubtractor using externally supplied rho = " << _rho << " and rho_m = " << _rhom << " to describe the background";
    } else {
      if (_bge_rhom) {
	descr << "ConstituentSubtractor using [" << _bge_rho->description() << "] and [" << _bge_rhom->description() << "] to estimate the background";
      } else {
	descr << "ConstituentSubtractor using [" << _bge_rho->description() << "] to estimate the background";
      }
    }  
    return descr.str();
  }


  void ConstituentSubtractor::set_max_deltaR(double max_deltaR){
    if (max_deltaR>0){
      _use_max_deltaR=true;
      _max_deltaR=max_deltaR;
    }
    else _use_max_deltaR=false; 
  }


  void ConstituentSubtractor::set_alpha(double alpha){
    _alpha=alpha;
  }

 
  bool ConstituentSubtractor::_function_used_for_sorting(std::pair<double,int> i,std::pair<double, int> j){
    return (i.first < j.first);
  }




} // namespace contrib

FASTJET_END_NAMESPACE
