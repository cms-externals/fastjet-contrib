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


#ifndef __SIGNAL_FREE_BACKGROUND_ESTIMATOR_H__
#define __SIGNAL_FREE_BACKGROUND_ESTIMATOR_H__

#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/RectangularGrid.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib{

  /// \class SignalFreeBackgroundEstimator
  /// Background estimator class based on the method described in arXiv:2304.08383. The design of this class is based on the GridMedianBackgroundEstimator class located in the official FastJet software: ./tools/fastjet/tools/GridMedianBackgroundEstimator.hh


  class SignalFreeBackgroundEstimator : public fastjet::BackgroundEstimatorBase,
                                        public fastjet::RectangularGrid {
  public:
    SignalFreeBackgroundEstimator(double rapidity_max,
                                  double requested_grid_spacing): fastjet::RectangularGrid(rapidity_max, requested_grid_spacing), _rapidity_max(rapidity_max) {}


    /// Return a pointer to a copy of this estimator
    SignalFreeBackgroundEstimator * copy() const override {
      return new SignalFreeBackgroundEstimator(*this);
    };

    ///
    /// default dtor
    virtual ~SignalFreeBackgroundEstimator(){}

    /// Tell the background estimator that it has a new event, composed of the specified particles.
    void set_particles(const std::vector<fastjet::PseudoJet>& particles) override {
      std::vector<fastjet::PseudoJet> empty;
      this->set_particles(particles,empty);
    }

    /// Tell the background estimator that it has a new event, composed of the specified particles and also include estimated signal particles. If the variable measureOfPileup is set to -1, then the minimal jet rho threshold is constant and set to the parameter _jet_rho_min_intercept. If it is >0, then it is not constant and varies linearly with measureOfPileup. If signalChargedParticles are also provided, then the signal jets are obtained also from clustering of these particles (separate jet clustering is used).
    void set_particles(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& signalParticles, const double& measureOfPileup=-1 /* can be nPU, nPV, mu, rho*/, const std::vector<fastjet::PseudoJet>& signalChargedParticles=std::vector<fastjet::PseudoJet>());

    /// Define signal seeds from the user. This function must be called before "set_particles"
    void add_seeds_from_user(const std::vector<fastjet::PseudoJet> &additional_seeds);

    /// Get the full set of background properties
    fastjet::BackgroundEstimate estimate() const override;

    /// Get the full set of background properties for a given reference jet
    fastjet::BackgroundEstimate estimate(const fastjet::PseudoJet& jet) const override;

    /// Returns rho, the average background density per unit area
    double rho() const override;

    /// Returns rho, the average background density per unit area, locally at the position of a given jet
    double rho(const fastjet::PseudoJet & jet) override;

    /// Returns rho_m, the purely longitudinal, particle-mass-induced component of the background density per unit area
    double rho_m() const FASTJET_OVERRIDE;

    /// Returns rho_m locally at the jet position. As for rho(jet), it is non-const.
    double rho_m(const PseudoJet & jet) FASTJET_OVERRIDE;

    /// determine whether the automatic calculation of rho_m (by default true)
    void set_compute_rho_m(bool enable){ _enable_rho_m = enable; }

    /// Returns true if this background estimator has support for determination of rho_m.
    bool has_rho_m() const FASTJET_OVERRIDE {return _enable_rho_m;}

    /// Returns a textual description of the background estimator.
    std::string description() const override;

    /// Set a pointer to a rescaling function
    virtual void set_rescaling_class(const fastjet::FunctionOfPseudoJet<double>* rescaling_class) override;

    /// Set the parameters for signal seeds (distance parameter for clustering of signal particles into jets, and deltaR distance around signal jets to exclude areas with signal)
    void set_signal_seed_parameters(const double &signal_jets_R, const double &exclusion_deltaR){ _signal_jets_R = signal_jets_R; _exclusion_deltaR = exclusion_deltaR; }

    /// Set the parameters for the minimal jet rho when selecting jets from the signal jets
    void set_jet_rho_min(const double &jet_rho_min_intercept, const double &jet_rho_min_slope){ _jet_rho_min_intercept = jet_rho_min_intercept; _jet_rho_min_slope=jet_rho_min_slope;}

    /// Set the parameter for the minimal jet rho when selecting jets obtained from charged signal particles
    void set_jet_rho_min_charged(const double &jet_rho_min_charged){ _jet_rho_min_charged = jet_rho_min_charged;}

    /// Set the parameter for the spacing of ghosts
    void set_ghost_size(const double &ghost_size){ _ghost_size = ghost_size; }

    /// Set the parameter for the minimal tile area
    void set_tile_area_min(const double &tile_area_min){ _tile_area_min = tile_area_min; }

    /// Set parameters for the window to be used to compute the weighted median
    void set_window_parameters(const double &center, const double &half_window, const double &regulator_for_floating_center=-1, const double &signal_fraction_min=0.005, const double &shift_max=0.15){ _center = center; _half_window = half_window; _regulator_for_floating_center=regulator_for_floating_center; _signal_fraction_min=signal_fraction_min; _shift_max=shift_max; }

    /// Set to use weighted tiles for rho computation
    void set_use_weighted_tiles(const bool &use_weighted_tiles){ _use_weighted_tiles = use_weighted_tiles; }

  private:
    // function to compute the weighted median
    double _compute_weighted_median(std::vector<std::pair<double,double> > &tileRho) const;

    /// Store eta extend.
    double _rapidity_max=-1;

    /// Distance parameter for jet clustering with anti-kt algorithm to obtain the signal seeds
    double _signal_jets_R=0.3;

    /// Parameter for the exclusion delta R around the signal seeds
    double _exclusion_deltaR=0.4;

    /// Parameters to compute the minimal jet rho threshold: intercept and slope of the linear function as a function of some measure of pileup specified in the set_particles function:
    double _jet_rho_min_intercept=5;
    double _jet_rho_min_slope=8;

    /// Parameter for the minimal jet rho threshold for the jets clustered from charged particles. The default value is 20 GeV/(unit area), if the input charged particles are provided in GeV.
    double _jet_rho_min_charged=20;

    /// Parameter for the minimal tile area after excluding signal-contaminated areas. Tiles with smaller areas are excluded from the weighted median computation. Can be modified with function set_tile_area_min
    double _tile_area_min=0.00001;

    /// Parameter for the spacing of ghosts to evaluate tile areas after excluding signal-contaminated areas. The smaller ghostSize, the more precise is the computation, but also slower. The default is 0.04. Can be modified with function set_ghost_size
    double _ghost_size=0.04;

    /// Window parameters:
    /// Parameters for the window center and the size of half window used to compute rho
    double _center=0.5;
    double _half_window=0.1;
    /// Parameters for the floating center of the window to compute rho
    double _regulator_for_floating_center=-1;
    double _signal_fraction_min=0.005;
    double _shift_max=0.15;
    /// Fraction of scalar pT from signal particles used to move the window to get the median. It is computed based on the provided signal particles.
    double _signal_fraction;

    /// Sum of all tile areas which are used as weights to compute the weighted median
    double _weights_sum;

    /// bool to determine whether the rho should be computed using weigted tiles
    bool _use_weighted_tiles=true;

    /// verify that particles have been set and throw an error if not
    void verify_particles_set() const;

    /// Rescaling warning
    fastjet::LimitedWarning _warning_rescaling;

    bool _enable_rho_m=true;
    std::vector<fastjet::PseudoJet> _signal_seeds_from_user;

  };

}

FASTJET_END_NAMESPACE

#endif // __SIGNAL_FREE_BACKGROUND_ESTIMATOR_H__

