//
// SignalFreeBackgroundEstimator package
// Questions/comments: peter.berta@cern.ch
//
// Copyright (c) 2024-, Peter Berta
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

#include "SignalFreeBackgroundEstimator.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  void SignalFreeBackgroundEstimator::set_particles(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& signalParticles, const double& measureOfPileup, const std::vector<fastjet::PseudoJet>& chargedSignalParticles){

    assert(all_tiles_equal_area());

    _cache_available = false;
    _cached_estimate.reset();
    _cached_estimate.set_has_sigma(false);
    _cached_estimate.set_mean_area(mean_tile_area());

    // Cluster signal particles
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(_rapidity_max));
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, _signal_jets_R);
    fastjet::ClusterSequenceArea clusSeq(signalParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> allSignalJets = sorted_by_pt(clusSeq.inclusive_jets());

    std::vector<fastjet::PseudoJet> signal_seeds;
    double jet_rho_min=_jet_rho_min_intercept;
    if (measureOfPileup>0) jet_rho_min=_jet_rho_min_intercept+sqrt(measureOfPileup)*_jet_rho_min_slope;
    for (size_t i = 0; i < allSignalJets.size(); ++i) {
      if (allSignalJets.at(i).pt()/allSignalJets.at(i).area() > jet_rho_min) signal_seeds.emplace_back(allSignalJets.at(i));
    }
    if (chargedSignalParticles.size()>0){
      fastjet::ClusterSequenceArea clusSeq_charged(chargedSignalParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> chargedJets = sorted_by_pt(clusSeq_charged.inclusive_jets());
      for (size_t i = 0; i < chargedJets.size(); ++i) {
        if (chargedJets.at(i).pt()/chargedJets.at(i).area() > _jet_rho_min_charged) signal_seeds.emplace_back(chargedJets.at(i));
      }
    }
    // Add signal seeds provided by the user:
    signal_seeds.insert(signal_seeds.end(),_signal_seeds_from_user.begin(),_signal_seeds_from_user.end());

    double scalarPtSumSignal=0;
    for (size_t i = 0; i < signalParticles.size(); ++i) scalarPtSumSignal+=signalParticles[i].pt();
    _signal_fraction=scalarPtSumSignal/14000.;


    // Determine tile area outside the exclusion area
    std::vector<bool> tileStates(n_tiles(),false);
    std::vector<double> tileAreas(n_tiles());
    int nGhosts_phi=dphi()/_ghost_size;
    int nGhosts_rap=drap()/_ghost_size;
    fastjet::PseudoJet ghost(0,0,0,1);
    int nPhi = std::round(fastjet::twopi / dphi());
    for (int i = 0; i < n_tiles(); ++i) {
      double areaFraction=1;
      if (_use_weighted_tiles) {
        int tileIndexPhi = i % nPhi;
        int tileIndexRap = i / nPhi;
        int nGhostsCloseToHardScatter=0;
        for (int iphi = 0; iphi < nGhosts_phi; ++iphi) {
          double phi=(tileIndexPhi + (iphi+0.5)/(double)nGhosts_phi)*dphi();
          for (int irap = 0; irap < nGhosts_rap; ++irap) {
            double rap=(tileIndexRap + (irap+0.5)/(double)nGhosts_rap)*drap() + rapmin();
            ghost.reset_momentum_PtYPhiM(1,rap,phi,1e-200);
            for (size_t isignal = 0; isignal < signal_seeds.size(); ++isignal) {
              if (signal_seeds.at(isignal).delta_R(ghost) < _exclusion_deltaR) {
                nGhostsCloseToHardScatter++;
                break;
              }
            }
          }
        }
        areaFraction=1-nGhostsCloseToHardScatter/(double)(nGhosts_phi*nGhosts_rap);
      }

      tileAreas[i]=mean_tile_area()*areaFraction;
      if (tileAreas[i]>_tile_area_min) tileStates[i]=true;
    }

    std::vector<double> tilePt(n_tiles(), 0.);
    std::vector<double> tileMtMinusPt(n_tiles(), 0.);
    for (size_t i = 0; i < particles.size(); ++i) {
      // Find index of the tile the particle belongs to
      int tileIndex = tile_index(particles[i]);
      if (tileIndex < 0) {
        std::cout << "WARNING SignalFreeBackgroundEstimator::set_particles: Tile index is negative. Particle (pt,rap,phi,m): " <<  particles[i].pt() << " " << particles[i].rap() << " " <<  particles[i].phi() << " " << particles[i].m() << std::endl;
        continue;
      }

      bool closeToHardScatter = false;
      for (size_t j = 0; j < signal_seeds.size(); ++j) {
        if (signal_seeds.at(j).delta_R(particles[i]) < _exclusion_deltaR) {
          closeToHardScatter = true;
          break;
        }
      }
      if (closeToHardScatter) continue;

      if (_rescaling_class) {
        double r = (*_rescaling_class)(particles[i]);
        tilePt[tileIndex] += particles[i].pt() / r;
        if (_enable_rho_m) tileMtMinusPt[tileIndex] += (particles[i].mt()-particles[i].pt()) / r;
      }
      else{
        tilePt[tileIndex] += particles[i].pt();
        if (_enable_rho_m) tileMtMinusPt[tileIndex] += particles[i].mt()-particles[i].pt();
      }
    }

    /// Compute rho and rho_m for each tile, and assign weights based on tile area
    std::vector<std::pair<double,double> > tileRho, tileRhom;
    _weights_sum=0;
    for (int i = 0; i < n_tiles(); ++i) {
      if (!tileStates[i]) continue;
      tileRho.emplace_back(std::make_pair(tilePt[i]/tileAreas[i],tileAreas[i]));
      _weights_sum+=tileAreas[i];
      if (_enable_rho_m) tileRhom.emplace_back(std::make_pair(tileMtMinusPt[i]/tileAreas[i],tileAreas[i]));
    }
    if (tileRho.size()==0){
      std::cout << "WARNING SignalFreeBackgroundEstimator::set_particles: All tiles are excluded in this event, and it is not possible to estimate rho with SignalFreeBackgroundEstimator! Estimating it with GridMedianBackgroundEstimator instead. Crosscheck whether it is expected that all tiles are excluded in this event and try to optimize the parameters of SignalFreeBackgroundEstimator, if needed." << std::endl;
      std::unique_ptr<fastjet::GridMedianBackgroundEstimator> gridMedian(new fastjet::GridMedianBackgroundEstimator(_rapidity_max,0.55));
      if (_rescaling_class) gridMedian->set_rescaling_class(_rescaling_class);
      gridMedian->set_particles(particles);
      _cached_estimate=gridMedian->estimate();
      _cache_available = true;
      return;
    }

    /// Compute and set the final rho and rho_m
    double rho=_compute_weighted_median(tileRho);
    _cached_estimate.set_rho(rho);

    if (_enable_rho_m){
      double rhom=_compute_weighted_median(tileRhom);;
      _cached_estimate.set_has_rho_m(true);
      _cached_estimate.set_rho_m(rhom);
    }

    _cache_available = true;
  }


  /// Compute weighted median based on the weights corresponding to tile areas.
  double SignalFreeBackgroundEstimator::_compute_weighted_median(std::vector<std::pair<double,double> > &tileRho) const{
    sort(tileRho.begin(),tileRho.end());

    // Get the actual center in case the user requested to use floating center:
    double shift=0;
    if (_regulator_for_floating_center>=0) shift=(_signal_fraction-_signal_fraction_min)*_regulator_for_floating_center;
    if (shift<0) shift=0;
    if (shift>_shift_max) shift=_shift_max;
    if (shift>_center-_half_window) shift=_center-_half_window;
    double center_shifted=_center-shift;

    // Compute rho from the window specified by the user:
    double rho_nominator=0;
    double rho_denominator=0;
    bool insideWindow=false;
    double fractions_partialSum=0;
    for (size_t i = 0; i < tileRho.size(); ++i) {
      double fraction_i=tileRho[i].second/_weights_sum;
      double remaining_to_start=center_shifted-_half_window-fractions_partialSum;
      double remaining_to_stop=center_shifted+_half_window-fractions_partialSum;
      if (remaining_to_start<fraction_i && !insideWindow){
        rho_nominator=tileRho[i].first*(fraction_i-remaining_to_start);
        rho_denominator=fraction_i-remaining_to_start;
        insideWindow=true;
        if (remaining_to_stop<fraction_i) break;
      }
      else if (remaining_to_stop<fraction_i){
        rho_nominator+=tileRho[i].first*remaining_to_stop;
        rho_denominator+=remaining_to_stop;
        break;
      }
      else if (insideWindow){
        rho_nominator+=tileRho[i].first*fraction_i;
        rho_denominator+=fraction_i;
      }
      fractions_partialSum+=fraction_i;
    }

    return rho_nominator/rho_denominator;
  }



  /// Get the full set of background properties.
  fastjet::BackgroundEstimate SignalFreeBackgroundEstimator::estimate() const {
    verify_particles_set();

    return _cached_estimate;
  }


  /// Get the full set of background properties for a given reference jet.
  fastjet::BackgroundEstimate SignalFreeBackgroundEstimator::estimate(const fastjet::PseudoJet& jet) const {
    verify_particles_set();

    if (_rescaling_class) {
      fastjet::BackgroundEstimate local_estimate = _cached_estimate;
      local_estimate.apply_rescaling_factor((*_rescaling_class)(jet));

      return local_estimate;
    }

    return _cached_estimate;
  }


  /// Get rho, the median background density per unit area
  double SignalFreeBackgroundEstimator::rho() const {
    verify_particles_set();
    return _cached_estimate.rho();
  }


  /// Returns rho, the average background density per unit area, locally at the position of a given jet
  double SignalFreeBackgroundEstimator::rho(const fastjet::PseudoJet & jet) {
    if (_rescaling_class) return (*_rescaling_class)(jet)*_cached_estimate.rho();
    return _cached_estimate.rho();
  }

  // returns rho_m (particle-masses contribution to the 4-vector density)
  double SignalFreeBackgroundEstimator::rho_m() const {
    if (!_enable_rho_m){
      throw fastjet::Error("SignalFreeBackgroundEstimator: rho_m requested but rho_m calculation has been disabled.");
    }
    verify_particles_set();
    return _cached_estimate.rho_m();
  }

  // returns rho_m locally at the position of a given jet.
  double SignalFreeBackgroundEstimator::rho_m(const PseudoJet & jet)  {
    double rescaling = (_rescaling_class == 0) ? 1.0 : (*_rescaling_class)(jet);
    return rescaling*rho_m();
  }



  /// Get estimator description
  std::string SignalFreeBackgroundEstimator::description() const {
    std::ostringstream desc;
    desc << "SignalFreeBackgroundEstimator, with " << RectangularGrid::description();
    return desc.str();
  }


  // verify that particles have been set and throw an error if not
  void SignalFreeBackgroundEstimator::verify_particles_set() const {
    if (!_cache_available) throw fastjet::Error("SignalFreeBackgroundEstimator::verify_particles_set: rho() or sigma() called without particles having been set");
  }


  /// Set a pointer to a class that calculates the rescaling factor as a function of the jet position.
  void SignalFreeBackgroundEstimator::set_rescaling_class(const fastjet::FunctionOfPseudoJet<double>* rescaling_class_in) {
    if (!_cache_available) {
      _warning_rescaling.warn("SignalFreeBackgroundEstimator::set_rescaling_class: Found cached result. Set particles again to obtain correct calculation!");
    }

    BackgroundEstimatorBase::set_rescaling_class(rescaling_class_in);
  }

  /// User can add his/her own signal seeds (either additional seeds to those obtained from jet clustering or only these user-defined seeds can be used). This function must be called before "set_particles"
  void SignalFreeBackgroundEstimator::add_seeds_from_user(const std::vector<fastjet::PseudoJet> &seeds_from_user){
    _signal_seeds_from_user=seeds_from_user;
  }

}  // contrib namespace


FASTJET_END_NAMESPACE
