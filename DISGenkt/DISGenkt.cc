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

#include "fastjet/contrib/DISGenkt.hh"
#include "fastjet/NNH.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{
  /// \brief Cache the normalised momentum components and energy weight
  /// for a pseudojet.
  /// \param jet Input pseudojet.
  /// \param info Clustering parameters shared by all jets.
  void DISBriefJet::init(const PseudoJet & jet, DISjetInfo * info) {
    double norm = 1.0/sqrt(jet.modp2());
    _nx  = jet.px() * norm;
    _ny  = jet.py() * norm;
    _nz  = jet.pz() * norm;
    _E2p = (info->p() == 0) ? 1 : std::pow(jet.E(), 2*info->p());
    _beam_sign = info->beam_sign();
    // compute norm for distances between final-state particles
    if(info->R() < M_PI) 
      _dij_norm = 1 - cos(info->R());
    else 
      _dij_norm = 3 + cos(info->R()); 
  }

  /// \brief Evaluate the pairwise DIS generalised-\f$k_t\f$ distance.
  /// \param jet Jet to compare against.
  /// \return The pairwise distance used by the nearest-neighbour search.
  double DISBriefJet::distance(const DISBriefJet * jet) const {
    double dij = 1 - _nx*jet->_nx
                   - _ny*jet->_ny
                   - _nz*jet->_nz;
    dij = dij * fmin(_E2p, jet->_E2p) / _dij_norm;
    return dij;
  }

  /// \brief Evaluate the distance between the jet and the proton beam.
  /// \return The beam distance used by the clustering sequence.
  double DISBriefJet::beam_distance() const {
    return (1 - _beam_sign * _nz)*_E2p;
  }

  /// \brief Build the text description reported by the plugin.
  /// \return A summary string containing the clustering parameters.
  std::string DISGenktPlugin::description () const {
    std::ostringstream desc;
    desc << "DISCambridge with"
         << " p = " << p() 
         << ", R = " << R()
         << ", beam_sign = " << beam_sign(); 
    return desc.str();
  }

  /// \brief Perform the DIS generalised-\f$k_t\f$ clustering.
  /// \param cs Cluster sequence to populate with the recombination
  ///   history.
  void DISGenktPlugin::run_clustering(fastjet::ClusterSequence & cs) const{
    int njets = cs.jets().size();
    DISjetInfo info(p(), beam_sign(), R());
    NNH<DISBriefJet,DISjetInfo> nnh(cs.jets(), &info);

    while (njets > 0) {
      int i, j, k;
      double dij = nnh.dij_min(i, j);
      // now merge or record jet as final 
      if (j >= 0) {
        cs.plugin_record_ij_recombination(i, j, dij, k);
        nnh.merge_jets(i, j, cs.jets()[k], k);
      } else {
        cs.plugin_record_iB_recombination(i, dij);
        nnh.remove_jet(i);
      }
      njets--;
    }
  }

  /// \brief Find the jet with the largest DIS longitudinal projection.
  /// \param jets Jets to inspect.
  /// \return Index of the jet maximising
  ///   \f$E - \mathrm{beam\_sign}\, p_z\f$.
  int DISGenktPlugin::find_idx_macrojet(
      const std::vector<PseudoJet>& jets) const {
    double best = 0.0;
    int idx = -1;
    for (size_t i = 0; i < jets.size(); ++i) {
      double proj = jets[i].E() - beam_sign() * jets[i].pz();
      if (proj > best) {
        best = proj;
        idx  = static_cast<int>(i);
      }
    }
    return idx;
  }

  /// \brief Sort jets by decreasing DIS longitudinal projection.
  /// \param jets Jets to sort.
  /// \return A reordered copy of \p jets from largest to smallest in
  ///   \f$E - \mathrm{beam\_sign}\, p_z\f$.
  std::vector<PseudoJet> DISGenktPlugin::sorted_by_zjet(
      const std::vector<PseudoJet> & jets) const {
    std::vector<double> zproj(jets.size());
    for (size_t i = 0; i < jets.size(); i++) {
      zproj[i] = -(jets[i].E() - beam_sign() * jets[i].pz());
    }
    return objects_sorted_by_values(jets, zproj);
  }
} // namespace contrib

FASTJET_END_NAMESPACE
