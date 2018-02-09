// $Id: ConstituentSubtractor.cc 1094 2017-12-18 13:18:20Z berta $
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

LimitedWarning ConstituentSubtractor::_warning_unused_rhom;


  /// Constructor that takes a pointer to a background estimator for
  /// rho and optionally a pointer to a background estimator for
  /// rho_m.
  ConstituentSubtractor::ConstituentSubtractor(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom, double alpha, double max_distance, Distance distance) :
    _bge_rho(bge_rho), _bge_rhom(bge_rhom), _common_bge(false), _rhom_from_bge_rhom(false), _externally_supplied_rho_rhom(false), _distance(distance), _alpha(alpha), _max_distance(max_distance){
    if (_max_distance>0) _use_max_distance=true;
    else _use_max_distance=false;
    _polarAngleExp=0;
    _ghost_area=0.01;
    _remove_zero_pt_and_mtMinusPt_particles=true;
    _remove_all_zero_pt_particles=false;
    _max_eta=0;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _do_mass_subtraction=false;
    if (_common_bge || bge_rhom) _do_mass_subtraction=true;
  }

  ///
  /// Constructor that takes an externally supplied value for rho and, optionally, for rho_m.
  ConstituentSubtractor::ConstituentSubtractor(double rho, double rhom, double alpha, double max_distance, Distance distance) :      
    _bge_rho(0), _bge_rhom(0), _common_bge(false), _rhom_from_bge_rhom(false), _rho(rho), _rhom(rhom), _externally_supplied_rho_rhom(true), _distance(distance), _alpha(alpha), _max_distance(max_distance) {
    if (_max_distance>0) _use_max_distance=true;
    else _use_max_distance=false;
    assert(_rho  >= 0);
    assert(_rhom >= 0);
    _polarAngleExp=0;
    _ghost_area=0.01;
    _remove_zero_pt_and_mtMinusPt_particles=true;
    _remove_all_zero_pt_particles=false;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _do_mass_subtraction=true;
    _max_eta=0;
  }


  void ConstituentSubtractor::set_background_estimator(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom){
    _bge_rho=bge_rho;
    _bge_rhom=bge_rhom;
    if (_common_bge || bge_rhom) _do_mass_subtraction=true;
  }


  void ConstituentSubtractor::set_scalar_background_density(double rho, double rhom){
    _rho=rho;
    _rhom=rhom;
    assert(_rho  >= 0);
    assert(_rhom >= 0);
    _externally_supplied_rho_rhom=true;
    _common_bge=false;
    _do_mass_subtraction=true;
  }

  ///----------------------------------------------------------------------
  /// the action on a given jet
  fastjet::PseudoJet ConstituentSubtractor::result(const PseudoJet &jet) const{
    // make sure we have a BGE or a rho value
    if (!_bge_rho && !_externally_supplied_rho_rhom){
      throw Error("ConstituentSubtractor::result() constituent subtraction needs a BackgroundEstimator or a value for rho.");
    }
    if (_ghosts_constructed) throw Error("ConstituentSubtractor::result() The ghosts are constructed, but they are not needed when using this function. When you want to perform jet-by-jet correction, initialize a new ConstituentSubtractor without construction of ghosts.");

    ///----------------------------------------------------------------------
    /// sift ghosts and particles in the input jet
    std::vector<fastjet::PseudoJet> particles, ghosts;
    SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    std::vector<double> ghosts_area;
    unsigned long nGhosts=ghosts.size();
    for (unsigned int j=0;j<nGhosts; ++j){
      ghosts_area.push_back(ghosts[j].area());
    }

    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(ghosts,ghosts_area);
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(particles,backgroundProxies);
    fastjet::PseudoJet subtracted_jet=join(subtracted_particles);
    subtracted_jet.set_user_index(jet.user_index());

    return subtracted_jet;    
  }


  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::subtract_event(std::vector<fastjet::PseudoJet> const &particles, double max_eta){
    if (fabs(_max_eta/max_eta-1)>1e-5) _ghosts_constructed=false;
    if (!_ghosts_constructed) this->construct_ghosts_uniformly(max_eta);
    _ghosts_rapidity_sorted=true; // the ghosts are now sorted according to rapidity. This variable needs to be set to true to be able to use faster algorithm in "do_subtraction".
    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    std::vector<fastjet::PseudoJet> selected_particles;
    for (unsigned int i=0;i<particles.size();++i){
      if (fabs(particles[i].eta())<max_eta) selected_particles.push_back(particles[i]);
    }
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(selected_particles,backgroundProxies);
    _ghosts_rapidity_sorted=false;
    return subtracted_particles;
  }




  std::vector<fastjet::PseudoJet> ConstituentSubtractor::get_background_proxies_from_ghosts(std::vector<fastjet::PseudoJet> const &ghosts,std::vector<double> const &ghosts_area) const{
    std::vector<fastjet::PseudoJet> proxies;
    unsigned long nGhosts=ghosts.size();
    std::vector<double> pt;
    std::vector<double> mtMinusPt;
    if (_externally_supplied_rho_rhom){
      for (unsigned int j=0;j<nGhosts; ++j){
	pt.push_back(_rho*ghosts_area[j]);
	mtMinusPt.push_back(_rhom*ghosts_area[j]);
      }
    }
    else{
      for (unsigned int j=0;j<nGhosts; ++j) pt.push_back(_bge_rho->rho(ghosts[j])*ghosts_area[j]);
      if (_bge_rhom){
	if (_rhom_from_bge_rhom){
#if FASTJET_VERSION_NUMBER >= 30100
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rhom->rho_m(ghosts[j])*ghosts_area[j]);
#else
	  throw(Error("ConstituentSubtractor:: _rhom_from_bge_rhom not allowed for FJ<3.1"));
#endif  // end of code specific to FJ >= 3.1
	} else {
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rhom->rho(ghosts[j])*ghosts_area[j]);
	}
      }
      else if (_common_bge){
	// since FJ 3.1.0, some background estimators have an automatic internal calculation of rho_m
#if FASTJET_VERSION_NUMBER >= 30100
	// check if the BGE has internal support for rho_m
	if (_bge_rho->has_rho_m()){
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rho->rho_m(ghosts[j])*ghosts_area[j]);
	} else {
#endif  // end of code specific to FJ >= 3.1
	  BackgroundJetPtMDensity m_density;
	  JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
	  const FunctionOfPseudoJet<double> * orig_density = jmbge->jet_density_class();
	  jmbge->set_jet_density_class(&m_density);
	  for (unsigned int j=0;j<nGhosts; ++j)  mtMinusPt.push_back(jmbge->rho(ghosts[j])*ghosts_area[j]);
	  jmbge->set_jet_density_class(orig_density);
#if FASTJET_VERSION_NUMBER >= 30100
	}
#endif
      }
      else { // a single bge, only rho requested
	for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(1e-200);
#if FASTJET_VERSION_NUMBER >= 30100
	// In FJ3.1 and BGE with rho_m support, add a warning, similar to that in Subtractor
	double const rho_m_warning_threshold = 1e-5;
	if (_bge_rho->has_rho_m() && _bge_rho->rho_m()>rho_m_warning_threshold*_bge_rho->rho()){
	  _warning_unused_rhom.warn("ConstituentSubtractor:: Background estimator indicates non-zero rho_m, but the constituent subtractor does not use rho_m information; consider calling set_common_bge_for_rho_and_rhom(true) to include the rho_m information");
	}
#endif     
      }
    }


    fastjet::PseudoJet proxy(0,0,0,1);
    for (unsigned int j=0;j<nGhosts; ++j){
      double mass_squared=pow(mtMinusPt[j]+pt[j],2)-pow(pt[j],2);
      double mass=0;
      if (mass_squared>0) mass=sqrt(mass_squared);
      proxy.reset_momentum_PtYPhiM(pt[j],ghosts[j].rap(),ghosts[j].phi(),mass);
      proxies.push_back(proxy);
    }
    return proxies;
  }
  
  



  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::do_subtraction(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const &backgroundProxies,std::vector<fastjet::PseudoJet> *remaining_backgroundProxies) const{
    unsigned int nBackgroundProxies=backgroundProxies.size();
    unsigned int nParticles=particles.size();
    double max_distance_transformed=-1;
    if (_distance==ConstituentSubtractor::deltaR) max_distance_transformed=pow(_max_distance,2); 
    if (_distance==ConstituentSubtractor::angle) max_distance_transformed=-cos(_max_distance);

    ///
    /// sort particles according to rapidity to achieve faster performance for the whole event subtraction
    std::vector<fastjet::PseudoJet>  particles_sorted=particles;
    std::sort(particles_sorted.begin(),particles_sorted.end(),ConstituentSubtractor::_rapidity_sorting);

    ///
    /// get the kinematic variables for particles and background proxies in advance to achieve faster performance
    std::vector<double> particles_phi,particles_rap,particles_pt,particles_mt,pt_factors,polarAngle_factors,particle_factors;
    for (unsigned int i=0;i<nParticles; ++i){
      particles_phi.push_back(particles_sorted[i].phi());
      particles_rap.push_back(particles_sorted[i].rap());
      particles_pt.push_back(particles_sorted[i].pt());
      particles_mt.push_back(particles_sorted[i].mt());
      if (fabs(_alpha)>1e-5) pt_factors.push_back(pow(particles_pt[i],_alpha));
      else  pt_factors.push_back(1);
      if (fabs(_polarAngleExp)>1e-5) polarAngle_factors.push_back(pow(particles_pt[i]/sqrt(particles_sorted[i].pt2()+particles_sorted[i].pz()*particles_sorted[i].pz()),_polarAngleExp));
      else  polarAngle_factors.push_back(1);
      particle_factors.push_back(pt_factors[i]*polarAngle_factors[i]);
    }
    std::vector<double> backgroundProxies_phi,backgroundProxies_rap,backgroundProxies_pt,backgroundProxies_mt;
    for (unsigned int j=0;j<nBackgroundProxies;++j){
      backgroundProxies_phi.push_back(backgroundProxies[j].phi());
      backgroundProxies_rap.push_back(backgroundProxies[j].rap());
      backgroundProxies_pt.push_back(backgroundProxies[j].pt());
      backgroundProxies_mt.push_back(backgroundProxies[j].mt());
    }

    ///
    /// finding the rapidity range of particles for each ghost to achieve faster performance in the double loop over particles and proxies below
    double max_number_pairs=0;
    std::vector<unsigned int> backgroundProxies_minParticleIndex,backgroundProxies_maxParticleIndex;
    if (_use_max_distance && _distance==ConstituentSubtractor::deltaR && _ghosts_rapidity_sorted){
      for (unsigned int j=0;j<_ghosts_rapidities.size();++j){
	unsigned int min=this->_find_index_after(_ghosts_rapidities[j]-_max_distance,particles_rap);
	unsigned int max=this->_find_index_before(_ghosts_rapidities[j]+_max_distance,particles_rap);
	for (int k=0;k<_n_ghosts_phi;++k){
 	  backgroundProxies_minParticleIndex.push_back(min);
	  backgroundProxies_maxParticleIndex.push_back(max);
	}
	max_number_pairs+=max-min;
      }
      max_number_pairs=max_number_pairs*_n_ghosts_phi*_max_distance/3.1415;
    }
    else{
      for (unsigned int j=0;j<nBackgroundProxies; ++j){
	backgroundProxies_minParticleIndex.push_back(0);
	backgroundProxies_maxParticleIndex.push_back(nParticles);      
      }
      if (_use_max_distance && _ghosts_constructed){
	if (_distance==ConstituentSubtractor::deltaR) max_number_pairs=nParticles*nBackgroundProxies*_max_distance*_max_distance/4./_max_eta;
	if (_distance==ConstituentSubtractor::angle) max_number_pairs=nParticles*nBackgroundProxies*_max_distance/3.141593;
      }
      else max_number_pairs=nParticles*nBackgroundProxies;
    }
    
    ///
    /// computation of the CS distances
    std::vector<std::pair<double,std::pair<int,int> > > CS_distances;  // storing three elements: the CS distance, and corresponding particle and proxy indexes
    CS_distances.reserve(max_number_pairs);

    bool skip_particles_outside_phi_range=false; // used for speed optimization
    if (_distance==ConstituentSubtractor::deltaR && _use_max_distance && _max_distance<twopi/2.*0.9999 && _ghosts_constructed) skip_particles_outside_phi_range=true;
    double distance_transformed = 0;
    bool switched=false;
    double particle_phi_max=0;
    double particle_phi_min=0;
    for (unsigned int j=0;j<nBackgroundProxies; ++j){
      if (skip_particles_outside_phi_range){
	switched=false;
	particle_phi_max=backgroundProxies_phi[j]+_max_distance;
	particle_phi_min=backgroundProxies_phi[j]-_max_distance;
	if (particle_phi_max>twopi){
	  particle_phi_min=particle_phi_max-twopi;
	  particle_phi_max=backgroundProxies_phi[j]-_max_distance;
	  switched=true;
	}
	if (particle_phi_min<0){
	  particle_phi_max=particle_phi_min+twopi;
	  particle_phi_min=backgroundProxies_phi[j]+_max_distance;
	  switched=true;
	}
      }
      for (unsigned int i=backgroundProxies_minParticleIndex[j];i<backgroundProxies_maxParticleIndex[j];++i){
	if (_distance==ConstituentSubtractor::deltaR){
	  if (skip_particles_outside_phi_range) if ((switched && particles_phi[i]>particle_phi_min && particles_phi[i]<particle_phi_max) || (!switched && (particles_phi[i]<particle_phi_min || particles_phi[i]>particle_phi_max))) continue;  // speed optimization only, this line has no effect on the subtraction performance
	  double deltaPhi=fabs(backgroundProxies_phi[j]-particles_phi[i]);
	  if (deltaPhi>pi) deltaPhi=twopi-deltaPhi;
	  double deltaRap=backgroundProxies_rap[j]-particles_rap[i];
	  distance_transformed = deltaPhi*deltaPhi+deltaRap*deltaRap;
	}
	if (_distance==ConstituentSubtractor::angle){
	  double particle_magnitude=sqrt(particles_sorted[i].pt2()+particles_sorted[i].pz()*particles_sorted[i].pz());
	  double backgroundProxy_magnitude=sqrt(backgroundProxies[j].pt2()+backgroundProxies[j].pz()*backgroundProxies[j].pz());
	  double scalar_product = particles_sorted[i].px()*backgroundProxies[j].px() + particles_sorted[i].py()*backgroundProxies[j].py() + particles_sorted[i].pz()*backgroundProxies[j].pz();
	  distance_transformed = -scalar_product/particle_magnitude/backgroundProxy_magnitude;
	}
	if (!_use_max_distance || distance_transformed<=max_distance_transformed){
	  double CS_distance = distance_transformed*particle_factors[i];
	  CS_distances.push_back(std::make_pair(CS_distance,std::make_pair(i,j))); // have tried to use emplace_back and tuple here - did not lead to any speed improvement
	}
      }
    }

    ///
    /// sorting of the CS distances
    std::sort(CS_distances.begin(),CS_distances.end(),ConstituentSubtractor::_function_used_for_sorting);

    ///----------------------------------------------------------------------
    /// the iterative process. Here, only finding the fractions of pt or delta_m=mt-pt to be corrected. The actual correction of particles is done later.
    unsigned long nStoredPairs=CS_distances.size();

    std::vector<double> backgroundProxies_fraction_of_pt(nBackgroundProxies,1.);
    std::vector<double> particles_fraction_of_pt(nParticles,1.);
    std::vector<double> backgroundProxies_fraction_of_mtMinusPt(nBackgroundProxies,1.);
    std::vector<double> particles_fraction_of_mtMinusPt(nParticles,1.);
    for (unsigned long iindices=0;iindices<nStoredPairs;++iindices){
      int particle_index=CS_distances[iindices].second.first;
      int proxy_index=CS_distances[iindices].second.second;

      if (backgroundProxies_fraction_of_pt[proxy_index]>0 && particles_fraction_of_pt[particle_index]>0 && particles_pt[particle_index]>0 && backgroundProxies[proxy_index].pt()>0){
	double ratio_pt=particles_pt[particle_index]*particles_fraction_of_pt[particle_index]/backgroundProxies_pt[proxy_index]/backgroundProxies_fraction_of_pt[proxy_index];
	if (ratio_pt>1){
	  particles_fraction_of_pt[particle_index]*=1-1./ratio_pt;
	  backgroundProxies_fraction_of_pt[proxy_index]=-1;
	}
	else {
	  backgroundProxies_fraction_of_pt[proxy_index]*=1-ratio_pt;
	  particles_fraction_of_pt[particle_index]=-1;
	}
      }
      if (_do_mass_subtraction && backgroundProxies_fraction_of_mtMinusPt[proxy_index]>0 && particles_fraction_of_mtMinusPt[particle_index]>0 && particles_mt[particle_index]>particles_pt[particle_index] && backgroundProxies_mt[proxy_index]>backgroundProxies_pt[proxy_index]){
	double ratio_mtMinusPt=(particles_mt[particle_index]-particles_pt[particle_index])*particles_fraction_of_mtMinusPt[particle_index]/(backgroundProxies_mt[proxy_index]-backgroundProxies_pt[proxy_index])/backgroundProxies_fraction_of_mtMinusPt[proxy_index];
	if (ratio_mtMinusPt>1){
	  particles_fraction_of_mtMinusPt[particle_index]*=1-1./ratio_mtMinusPt;
	  backgroundProxies_fraction_of_mtMinusPt[proxy_index]=-1;
	}
	else{
	  backgroundProxies_fraction_of_mtMinusPt[proxy_index]*=1-ratio_mtMinusPt;
	  particles_fraction_of_mtMinusPt[particle_index]=-1;
	}
      }
    }

    ///----------------------------------------------------------------------
    /// do the actual correction for particles:
    std::vector<fastjet::PseudoJet> subtracted_particles;
    for (unsigned int i=0;i<nParticles; ++i){
      ///
      /// particles with zero pt and zero mtMinusPt are removed. Particles with zero pt and non-zero mtMinusPt (i.e. massive particles in rest) are kept with very small pt (1e-200 GeV) - the user can decide in his/her code if to remove or not to remove them.

      /*      if (particles_phi[i]>5.9556 && particles_phi[i]<5.9558){
	std::cout << "fraction of pt: " << particles_fraction_of_pt[i] << "  pt: " << particles_pt[i] << "  fraction of deltam: " << particles_fraction_of_mtMinusPt[i] << "  deltam: " << particles_mt[i] << "  rap: " << particles_rap[i] << std::endl;
	}*/
      bool particle_pt_larger_than_zero=(particles_fraction_of_pt[i]>0 && particles_pt[i]>0);
      bool particle_mtMinusPt_larger_than_zero=(_do_mass_subtraction && particles_fraction_of_mtMinusPt[i]>0 && particles_mt[i]>particles_pt[i]);
      if (!particle_pt_larger_than_zero && !particle_mtMinusPt_larger_than_zero && _remove_zero_pt_and_mtMinusPt_particles) continue;
      if (!particle_pt_larger_than_zero && _remove_all_zero_pt_particles) continue;
      double subtracted_pt=1e-50;
      if (particle_pt_larger_than_zero) subtracted_pt=particles_pt[i]*particles_fraction_of_pt[i];
      double  subtracted_mtMinusPt=0;
      if (particle_mtMinusPt_larger_than_zero) subtracted_mtMinusPt=(particles_mt[i]-particles_pt[i])*particles_fraction_of_mtMinusPt[i];
      PseudoJet subtracted_const(0,0,0,1);
      double mass_squared=pow(subtracted_pt+subtracted_mtMinusPt,2)-pow(subtracted_pt,2);
      if (mass_squared<0) mass_squared=0;
      subtracted_const.reset_momentum_PtYPhiM(subtracted_pt,particles_rap[i],particles_phi[i],sqrt(mass_squared));
      subtracted_const.set_user_index(particles_sorted[i].user_index());
      subtracted_particles.push_back(subtracted_const);
      /*if (particles_phi[i]>5.9556 && particles_phi[i]<5.9558){
	std::cout << "pt: " << subtracted_const.pt() << "  rap: " << subtracted_const.rap() << "  mass: " << subtracted_const.m() << std::endl;
	}*/
    }
    ///
    /// get the remaining background proxies if requested:
    if (remaining_backgroundProxies){
      for (unsigned int i=0;i<nBackgroundProxies; ++i){
      /// keeping all background proxies
	bool proxy_pt_larger_than_zero=(backgroundProxies_fraction_of_pt[i]>0 && backgroundProxies_pt[i]>0);
	bool proxy_mtMinusPt_larger_than_zero=(_do_mass_subtraction && backgroundProxies_fraction_of_mtMinusPt[i]>0 && backgroundProxies_mt[i]>backgroundProxies_pt[i]);
	double subtracted_pt=1e-200;
	if (proxy_pt_larger_than_zero) subtracted_pt=backgroundProxies_pt[i]*backgroundProxies_fraction_of_pt[i];
	double subtracted_mtMinusPt=0;
	if (proxy_mtMinusPt_larger_than_zero) subtracted_mtMinusPt=(backgroundProxies_mt[i]-backgroundProxies_pt[i])*backgroundProxies_fraction_of_mtMinusPt[i];
	PseudoJet subtracted_const(0,0,0,1);
	double mass_squared=pow(subtracted_pt+subtracted_mtMinusPt,2)-pow(subtracted_pt,2);
	if (mass_squared<0) mass_squared=0;
	subtracted_const.reset_momentum_PtYPhiM(subtracted_pt,backgroundProxies_rap[i],backgroundProxies_phi[i],sqrt(mass_squared));
	remaining_backgroundProxies->push_back(subtracted_const);
      }
    }

    return subtracted_particles;
  }




  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::sequential_subtraction(std::vector<fastjet::PseudoJet> const &particles, double max_eta){
    if (fabs(_max_eta/max_eta-1)>1e-5) _ghosts_constructed=false;
    if (!_ghosts_constructed) this->construct_ghosts_uniformly(max_eta);
    _ghosts_rapidity_sorted=true; // the ghosts are now sorted according to rapidity. This variable needs to be set to true to be able to use faster algorithm in "do_subtraction". 
    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    std::vector<fastjet::PseudoJet> selected_particles;
    for (unsigned int i=0;i<particles.size();++i){
      if (fabs(particles[i].eta())<max_eta) selected_particles.push_back(particles[i]);
    }
    double total_background_pt=0,total_background_mt=0;
    for (unsigned int i=0;i<backgroundProxies.size();++i){
      total_background_pt+=backgroundProxies[i].pt();
      total_background_mt+=backgroundProxies[i].mt();
    }
    this->set_max_distance(_max_distance_sequential);
    std::vector<fastjet::PseudoJet> *remaining_backgroundProxies=new std::vector<fastjet::PseudoJet>;
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(selected_particles,backgroundProxies,remaining_backgroundProxies);
    double remaining_background_pt=0,remaining_background_mt=0;
    for (unsigned int i=0;i<remaining_backgroundProxies->size();++i){
      remaining_background_pt+=(remaining_backgroundProxies->at(i)).pt();
      remaining_background_mt+=(remaining_backgroundProxies->at(i)).mt();
    }
    std::vector<fastjet::PseudoJet> backgroundProxies2;
    for (unsigned int i=0;i<backgroundProxies.size();++i){
      double rapidity=backgroundProxies[i].rap();
      double azimuth=backgroundProxies[i].phi();
      double pt= backgroundProxies[i].pt()*remaining_background_pt/total_background_pt;
      double mtMinusPt= (backgroundProxies[i].mt()-backgroundProxies[i].pt())*(remaining_background_mt-remaining_background_pt)/(total_background_mt-total_background_pt);
      PseudoJet proxy(pt*cos(azimuth),pt*sin(azimuth),(pt+mtMinusPt)*sinh(rapidity),(pt+mtMinusPt)*cosh(rapidity));
      backgroundProxies2.push_back(proxy);
    }
    delete remaining_backgroundProxies;
    this->set_max_distance(_max_distance);
    std::vector<fastjet::PseudoJet> subtracted_particles2=this->do_subtraction(subtracted_particles,backgroundProxies2);
    _ghosts_rapidity_sorted=false;
    return subtracted_particles2;
  }



  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::subtract_event_using_charged_info(std::vector<fastjet::PseudoJet> const &particles, double charged_background_scale, std::vector<fastjet::PseudoJet> const &charged_background, double charged_signal_scale, std::vector<fastjet::PseudoJet> const &charged_signal, double max_eta, fastjet::FunctionOfPseudoJet<double> *rescaling){
    if (fabs(_max_eta/max_eta-1)>1e-5) _ghosts_constructed=false;
    if (!_ghosts_constructed) this->construct_ghosts_uniformly(max_eta);
    _ghosts_rapidity_sorted=false;  // no speed optimization implemented for this function yet

    std::vector<fastjet::PseudoJet>  scaled_charged_all;
    std::vector<fastjet::PseudoJet>  scaled_charged_signal;
    std::vector<fastjet::PseudoJet>  scaled_charged_background;
    for (unsigned int i=0;i<charged_background.size();++i){
      if (fabs(charged_background[i].eta())>max_eta) continue;
      scaled_charged_background.push_back(charged_background_scale*charged_background[i]);
      scaled_charged_all.push_back(scaled_charged_background[scaled_charged_background.size()-1]);
    }
    for (unsigned int i=0;i<charged_signal.size();++i){
      if (fabs(charged_signal[i].eta())>max_eta) continue;
      scaled_charged_signal.push_back(charged_signal_scale*charged_signal[i]);
      scaled_charged_all.push_back(scaled_charged_signal[scaled_charged_signal.size()-1]);
    }
    std::vector<fastjet::PseudoJet> selected_particles;
    for (unsigned int i=0;i<particles.size();++i){
      if (fabs(particles[i].eta())<max_eta) selected_particles.push_back(particles[i]);
    }
    std::vector<fastjet::PseudoJet> *remaining_charged_background= new std::vector<fastjet::PseudoJet>;
    double maxDeltaR=this->get_max_distance();
    if (maxDeltaR<=0) maxDeltaR=0.5;
    this->set_max_distance(0.2);
    std::vector<fastjet::PseudoJet> subtracted_particles_using_scaled_charged_signal=this->do_subtraction(selected_particles,scaled_charged_signal); 
    std::vector<fastjet::PseudoJet> subtracted_particles_using_scaled_charged_all=this->do_subtraction(subtracted_particles_using_scaled_charged_signal,scaled_charged_background,remaining_charged_background);  // remaining neutral background particles
    std::vector<fastjet::PseudoJet> scaled_charged_background_used_for_subtraction=this->do_subtraction(scaled_charged_background,*remaining_charged_background); 
    _bge_rho= new fastjet::GridMedianBackgroundEstimator(max_eta, 0.6);
    if (_do_mass_subtraction) this->set_common_bge_for_rho_and_rhom(true);
    _bge_rho->set_rescaling_class(rescaling);
    _bge_rho->set_particles(subtracted_particles_using_scaled_charged_all);

    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    backgroundProxies.insert(backgroundProxies.end(), scaled_charged_background_used_for_subtraction.begin(), scaled_charged_background_used_for_subtraction.end());

    this->set_max_distance(maxDeltaR);
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(selected_particles,backgroundProxies);
    delete remaining_charged_background; 
    delete _bge_rho;
    return subtracted_particles;
  }



  void ConstituentSubtractor::set_common_bge_for_rho_and_rhom(bool value){ 
    if (!_bge_rho)  throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed when _bge_rho is not set!");
    if (value){
      if (_bge_rhom) throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed in the presence of an existing background estimator for rho_m.");
      if (_externally_supplied_rho_rhom) throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed when supplying externally the values for rho and rho_m.");
    }
#if FASTJET_VERSION_NUMBER >= 30100
    if (!_bge_rho->has_rho_m()){
#endif
      JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
      if (!jmbge)	throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is currently only allowed for background estimators of JetMedianBackgroundEstimator type.");
#if FASTJET_VERSION_NUMBER >= 30100
    }
#endif
    _common_bge=value;
    _do_mass_subtraction=value;
  }


  // setting this to true will result in rho_m being estimated using bge_rhom->rho_m() instead of bge_rhom->rho()
  void ConstituentSubtractor::set_use_bge_rhom_rhom(bool value){
    if (!value){
      _rhom_from_bge_rhom=false;
      return;
    }
 
#if FASTJET_VERSION_NUMBER < 30100
    throw Error("ConnstituentSubtractor::use_rhom_from_bge_rhom() can only be used with FastJet >=3.1.");
#else	
    if (!_bge_rhom) throw Error("ConstituentSubtractor::use_rhom_from_bge_rhom() requires a background estimator for rho_m.");
    
    if (!(_bge_rhom->has_rho_m())) throw Error("ConstituentSubtractor::use_rhom_from_bge_rhom() requires rho_m support for the background estimator for rho_m.");
#endif	
    _rhom_from_bge_rhom=true; 
    _do_mass_subtraction=true;
  }


  void ConstituentSubtractor::set_do_mass_subtraction(bool do_mass_subtraction){
    _do_mass_subtraction=do_mass_subtraction;
  }



  std::string ConstituentSubtractor::description() const{
    std::ostringstream descr;
    if ( _externally_supplied_rho_rhom){
      descr << "ConstituentSubtractor using externally supplied rho = " << _rho << " and rho_m = " << _rhom << " to describe the background";
    } else {
      if (_bge_rhom && _bge_rho) {
	descr << "ConstituentSubtractor using [" << _bge_rho->description() << "] and [" << _bge_rhom->description() << "] to estimate the background";
      } else {
	if (_bge_rho) descr << "ConstituentSubtractor using [" << _bge_rho->description() << "] to estimate the background";
	else descr << "ConstituentSubtractor: no externally supplied rho, nor background estimator";
      }
    }  
    descr << std::endl << "perform mass subtraction: " << _do_mass_subtraction << std::endl;
    return descr.str();
  }



  void ConstituentSubtractor::set_distance_type(Distance distance){
    _distance=distance;
  }


  void ConstituentSubtractor::set_max_distance(double max_distance){
    if (max_distance>0){
      _use_max_distance=true;
      _max_distance=max_distance;
    }
    else _use_max_distance=false; 
  }


  void ConstituentSubtractor::set_max_standardDeltaR(double max_distance){
    this->set_max_distance(max_distance);
  }


  void ConstituentSubtractor::set_max_distance_sequential(double max_distance_sequential){
    _max_distance_sequential=max_distance_sequential;
  }

  double ConstituentSubtractor::get_max_distance(){
    return _max_distance;
  }


  void ConstituentSubtractor::set_alpha(double alpha){
    _alpha=alpha;
  }

  void ConstituentSubtractor::set_polarAngleExp(double polarAngleExp){
    _polarAngleExp=polarAngleExp;
  }

 
  void ConstituentSubtractor::set_ghost_area(double ghost_area){
    _ghost_area=ghost_area;
    _ghosts_constructed=false;
  }


  void ConstituentSubtractor::set_remove_zero_pt_and_mtMinusPt_particles(bool value){
    _remove_zero_pt_and_mtMinusPt_particles=value;
  }


  void ConstituentSubtractor::set_remove_all_zero_pt_particles(bool value){
    _remove_all_zero_pt_particles=value;
  }


  bool ConstituentSubtractor::_function_used_for_sorting(std::pair<double,std::pair<int,int> >  const &i,std::pair<double,std::pair<int,int> >  const &j){
    return (i.first < j.first);
  }
  
  bool ConstituentSubtractor::_rapidity_sorting(fastjet::PseudoJet const &i,fastjet::PseudoJet  const &j){
    return (i.rap() < j.rap());
  }

  unsigned int ConstituentSubtractor::_find_index_after(double const &value, std::vector<double> const &vec) const{
    int size=vec.size();
    if (size==0) return -1;
    int nIterations=log(size)/log(2)+2;
    unsigned int lowerBound=0;
    unsigned int upperBound=size-1;
    if (value<vec[0]) return 0;
    if (value>vec[size-1]) return size-1;
    for (int i=0;i<nIterations;++i){
      //      std::cout << i << "  lowerBound: " << lowerBound << std::endl;
      unsigned int test=(upperBound+lowerBound)/2;
      if (value>vec[test]){
	if (value<vec[test+1]) return test+1;
	lowerBound=test;
      }
      else{
	if (value>vec[test-1]) return test;
	upperBound=test;
      }
    }
    return lowerBound;
  }

  unsigned int ConstituentSubtractor::_find_index_before(double const &value, std::vector<double> const &vec) const{
    int size=vec.size();
    if (size==0) return -1;
    int nIterations=log(size)/log(2)+1;
    unsigned int lowerBound=0;
    unsigned int upperBound=size-1;
    if (value<vec[0]) return 1;  // it is higher by one to account for the "<" comparison in the for loop
    if (value>vec[size-1]) return size;  // it is higher by one to account for the "<" comparison in the for loop
    for (int i=0;i<nIterations;++i){
      unsigned int test=(upperBound+lowerBound)/2;
      if (value>vec[test]){
	if (value<vec[test+1]) return test+1;  // it is higher by one to account for the "<" comparison in the for loop
	lowerBound=test;
      }
      else{
	if (value>vec[test-1]) return test;  // it is higher by one to account for the "<" comparison in the for loop
	upperBound=test;
      }
    }
    return upperBound+1;
  }




  void ConstituentSubtractor::construct_ghosts_uniformly(double max_eta){
    _ghosts.clear();
    _ghosts_rapidities.clear();
    _ghosts_area.clear();
    _max_eta=max_eta;
    double a=sqrt(_ghost_area);
    _n_ghosts_phi=(2*3.14159265/a+0.5); // rounding 
    int nRap=(2*max_eta/a+0.5); // rounding 
    _grid_size_phi=2*3.14159265/(double)_n_ghosts_phi;
    _grid_size_rap=2*max_eta/(double)nRap;
    double used_ghost_area=_grid_size_phi*_grid_size_rap;
    fastjet::PseudoJet ghost(0,0,0,1);
    for (int iRap=0;iRap<nRap;++iRap){
      double rapidity=_grid_size_rap*(iRap+0.5)-max_eta;
      _ghosts_rapidities.push_back(rapidity);
      for (int iPhi=0;iPhi<_n_ghosts_phi;++iPhi){
	ghost.reset_momentum_PtYPhiM(1,rapidity,_grid_size_phi*(iPhi+0.5),1e-200);
	_ghosts.push_back(ghost);
	_ghosts_area.push_back(used_ghost_area);
      }
    }
    _ghosts_constructed=true;
  }


  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::get_ghosts(){
    return _ghosts;
  }


  std::vector<double>  ConstituentSubtractor::get_ghosts_area(){
    return _ghosts_area;
  }

} // namespace contrib


FASTJET_END_NAMESPACE
