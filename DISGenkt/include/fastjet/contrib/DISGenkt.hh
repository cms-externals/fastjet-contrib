// $Id: DISGenkt.hh 1554 2026-05-18 15:33:57Z silviaferrarioravasio $
//
// Copyright (c) 2026-, Melissa van Beekveld, Silvia Ferrario Ravasio, 
// Alexander Karlberg and Darcy Peake.
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

#ifndef __FASTJET_CONTRIB_DISGenktPlugin_HH__
#define __FASTJET_CONTRIB_DISGenktPlugin_HH__

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{
  /// \class DISjetInfo
  /// \brief Store the fixed parameters needed during clustering.
  ///
  /// This helper object passes the distance-measure power and the beam
  /// orientation to the nearest-neighbour machinery.
  class DISjetInfo {
  public:
    /// \brief Build the clustering parameter container.
    /// \param p Power entering the DIS generalised-\f$k_t\f$ distance.
    /// \param beam_sign Sign of the proton-beam direction along the
    /// \param R jet radius parameter
    /// \f$z\f$ axis.
    DISjetInfo(double p, int beam_sign, double R) {
      _p = p;
      _beam_sign = beam_sign;
      _R = R;
    }

    double p() { return _p; }
    double R() { return _R; }
    int beam_sign() { return _beam_sign; }

  private:
    double _p, _R;
    int _beam_sign;
  };

  /// \class DISBriefJet
  /// \brief Lightweight jet representation used by the nearest-neighbour
  /// search.
  ///
  /// The pairwise distance is
  /// \f$d_{ij} = \min(E_i^{2p}, E_j^{2p})(1-\cos\theta_{ij})\f$, while
  /// the beam distance is
  /// \f$d_{iB} = (1-\mathrm{beam\_sign}\, n_{iz}) E_i^{2p}\f$, with
  /// \f$n_{iz}\f$ the normalised longitudinal direction of the jet
  /// momentum.
  class DISBriefJet {
  public:
    /// \brief Initialise the cached quantities used in distance
    /// calculations.
    /// \param jet Input pseudojet.
    /// \param info Clustering parameters shared by all jets.
    void init(const PseudoJet & jet, DISjetInfo * info);

    /// \brief Compute the pairwise distance to another cached jet.
    /// \param jet Jet to compare against.
    /// \return The DIS generalised-\f$k_t\f$ pairwise distance.
    double distance(const DISBriefJet * jet) const;

    /// \brief Compute the distance to the proton beam.
    /// \return The beam distance for this jet.
    double beam_distance() const;

  private:
    // normalised vector of three-momentum and E^{2p}
    double _nx, _ny, _nz, _E2p;
    double _dij_norm; 
    int _beam_sign;
  };

//------------------------------------------------------------------------
/// \class DISGenktPlugin
/// \brief FastJet plugin for the inclusive DIS generalised-kt
/// algorithm in the Breit frame.
///
/// This plugin implements the spherically invariant generalised-kt
/// clustering for Deep Inelastic Scattering in the Breit frame. The 
/// distance is controlled by the power parameter p>=0, while
/// beam_sign selects the proton-beam direction along the z axis.
class DISGenktPlugin : public JetDefinition::Plugin{
public:
  /// \brief Construct a DIS generalised-\f$k_t\f$ plugin instance.
  /// \param p Power entering the pairwise distance
  ///   \f$d_{ij} = 2\min(E_i^{2p}, E_j^{2p})(1-\cos\theta_{ij})\f$
  /// and the beam distance \f$d_{iB} = 2E_i^{2p}(1-\cos\theta_{iB})\f$
  /// \param beam_sign Sign of the proton-beam direction along the \f$z\f$
  ///   axis; must be either \c +1 or \c -1.
  /// \param R optional radius parameter, used to normalise the pairwise
  ///   distance. 
  DISGenktPlugin(double p, int beam_sign, double R = M_PI/2.)
      : _p(p), _R(R), _beam_sign(beam_sign) {
    _Q2 = -1;
    if(_p < -1)
      throw std::runtime_error("p should be larger or equal than 0");
    if(abs(_beam_sign) != 1)
      throw std::runtime_error(
          "Expected a beam sign of +/- 1, got" +
          std::to_string(_beam_sign));
  }

  /// \brief Destructor.
  ~DISGenktPlugin(){}

  /// \brief Return a human-readable description of the jet definition.
  virtual std::string description() const override;

  /// \brief Run the clustering for the supplied cluster sequence.
  /// \param cs Cluster sequence to be filled by the plugin.
  virtual void run_clustering(ClusterSequence &) const override;

  /// \brief Return a dummy jet radius.
  /// \return _R
  virtual double R() const override { return _R; }

  /// \brief Report that exclusive jets are meaningful for this sequence.
  /// \return Always returns \c true.
  virtual bool exclusive_sequence_meaningful() const override { return p() >= 0;}

  /// \brief Report whether the algorithm is spherical.
  /// Intended for use on a spherical geometry - DIS detectors
  /// do use a spherical geometry. 
  /// \return Always returns \c true
  virtual bool is_spherical() const override { return true; }

  /// \brief Return the power parameter of the clustering measure.
  /// \return The current value of \f$p\f$.
  double p() const { return _p; }

  /// \brief Return the sign of the proton-beam direction.
  /// \return Either \c +1 or \c -1.
  int beam_sign() const { return _beam_sign; }

  /// \brief Update the beam direction after construction.
  /// \param beam_sign New beam direction sign along the \f$z\f$ axis.
  void reset_beam_sign(const int beam_sign) { _beam_sign = beam_sign; }

  /// \brief reset Q2 (= DIS invariant)
  /// \param Q2 new value of Q2
  void reset_Q2(const double Q2) { _Q2 = Q2; }

  /// \brief Find the macrojet with the largest light-cone component
  /// along the direction of the original struct quark in the Breit
  /// frame (i.e. anti-parallel to the proton).
  /// \param jets Jets to inspect.
  /// \return Index of the jet maximising
  ///   \f$E - \mathrm{beam\_sign}\, p_z \f$, or \c -1 if no such jet is found.
  int find_idx_macrojet(const std::vector<PseudoJet>& jets) const;

  /// \brief Sort jets by decreasing light-cone component along the
  /// direction of the original struct quark in the Breit frame
  /// (i.e. anti-parallel to the proton).
  /// \param jets Jets to sort.
  /// \return A reordered copy of \p jets, sorted by decreasing
  ///   \f$E - \mathrm{beam\_sign}\, p_z\f$.
  std::vector<PseudoJet> sorted_by_zjet(
      const std::vector<PseudoJet> & jets) const;


  private:
    double _p, _R, _Q2;     
    int    _beam_sign;
};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_DISGenktPlugin_HH__
