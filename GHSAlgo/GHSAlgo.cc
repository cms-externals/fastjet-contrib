// $Id: GHSAlgo.cc 1507 2025-01-30 10:33:36Z gstagn $
//
// Copyright (c) 2025, Rhorry Gauld, Alexander Huss, Giovanni Stagnitto
//
// based on initial version by Fabrizio Caola, Radoslaw Grabarczyk,
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
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

#include "fastjet/contrib/GHSAlgo.hh"

#include <limits>
// to facilitate use with fjcore
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/ClusterSequence.hh"
#include "fastjet/NNH.hh"
#endif
#include <iomanip>
#include <iostream>

#include "fastjet/contrib/FlavInfo.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  using namespace std;

  const double _deltaR2_handover =
      pow(std::numeric_limits<double>::epsilon(), 0.5);

#ifdef VERBOSEDIJ
  void print_PJ(ostream * ostr, const PseudoJet &p, unsigned precision,
                bool short_version, bool final_flav) {
    (*ostr).precision(precision);
    (*ostr) << setw(precision + 8) << p.px() << setw(precision + 8) << p.py()
            << setw(precision + 8) << p.pz() << setw(precision + 8) << p.E()
            << setw(precision + 8) << p.pt() << setw(precision + 8) << p.rap()
            << setw(precision + 8) << p.phi() << setw(precision + 8) << p.m()
            << setw(6) << p.user_index();
#ifdef PRINTMASSFLAVINFO
    (*ostr) << setw(6) << p.has_user_info<MassFlavHistory>();
#else
    (*ostr) << setw(6) << "";
#endif
    if (p.has_user_info<FlavHistory>()) {
      if (final_flav)
        (*ostr) << setw(14) << FlavHistory::current_flavour_of(p).description();
      else
        (*ostr) << setw(14) << FlavHistory::initial_flavour_of(p).description();
      // cout << json(FlavHistory::current_flavour_of(p)) << endl;
    }
    if (!short_version) {
      (*ostr) << setw(6) << p.cluster_hist_index();
    }
  }
#endif

  struct GHSInfo {
    double alpha;
    double omega;
    unsigned njets;
    vector<PseudoJet> jets;
    FlavRecombiner flav_recombiner;
  };

  class GHSBriefJet {
   public:
    void init(const PseudoJet &jet, GHSInfo *info) {
      _jet = jet;
      _info = info;
      _pt2 = jet.pt2();
      double pt = sqrt(_pt2);
      _nx = jet.px() / pt;
      _ny = jet.py() / pt;
      _rap = jet.rap();
      _phi = jet.phi();
      constexpr double rap_transition = 0.1;
      if (fabs(_rap) < rap_transition) {
        // for moderately small rapidities switch to special rapidity formula
        //
        // rap = 1/2 log[(E+pz)/(E-pz)]
        //     = 1/2 log[1 + 2pz/(E-pz)]
        //
        // and use log1p for the latter
        _rap = 0.5 * log1p(2 * jet.pz() / (jet.E() - jet.pz()));
      }

      // now, evalute the beam distance for the PseudoJet
      if (is_jet()) {
        // no beam distance
        _diB = numeric_limits<double>::max();
      } else if (is_particle()) {
        // beam distance is defined only if the particle has no assigned jet
        if (associated_jet() == -1) {
          double ptBp = 0, ptBm = 0;
          for (const auto &j : _info->jets) {
            // _info->jets are PJs not BriefJets (as they must be to have an
            // associated cs) so we must recalculate the precise rapidity here.
            constexpr double rap_transition = 0.1;
            double jrap = j.rap();
            if (fabs(jrap) < rap_transition) {
              jrap = 0.5 * log1p(2 * j.pz() / (j.E() - j.pz()));
            }
            double jetrap = jet.rap();
            if (fabs(jetrap) < rap_transition) {
              jetrap = 0.5 * log1p(2 * jet.pz() / (jet.E() - jet.pz()));
            }
            double dyj = jrap - jetrap;
            ptBp += j.pt() * exp(min(0.0, -dyj));
            ptBm += j.pt() * exp(min(0.0, dyj));
          }
          double alpha = _info->alpha;
          double diBp = max(pow(jet.pt(), alpha), pow(ptBp, alpha)) *
                        min(pow(jet.pt(), 2 - alpha), pow(ptBp, 2 - alpha));
          double diBm = max(pow(jet.pt(), alpha), pow(ptBm, alpha)) *
                        min(pow(jet.pt(), 2 - alpha), pow(ptBm, 2 - alpha));
          _diB = min(diBp, diBm);
        } else {
          _diB = numeric_limits<double>::max();
        }

      } else {
        assert(false &&
               "the PseudoJet should either be a particle or a jet, "
               "but was neither");
      }
    }
    //> idx in jet if associated [0,...,jets.size()];  -1 if unassociated; 0 if
    // jet
    int associated_jet() const { return abs(_jet.user_index()) - 1; }
    bool is_jet() const { return _jet.user_index() == 1; }
    bool is_particle() const { return _jet.user_index() <= 0; }

    static double geometrical_distance_squared(const GHSBriefJet *first,
                                               const GHSBriefJet *other) {
      // do straight rapidity difference, because we already took
      // care of making the rapidity as accurate as possible
      double delta_y = first->_rap - other->_rap;
      double delta_phi = std::fabs(first->_phi - other->_phi);
      if (delta_phi > pi) delta_phi = twopi - delta_phi;

      // transition is somewhat arbitrary, but should be such that
      // we are in a region where arcsin() is unambiguous; can be
      // O(1), but must be < pi/2
      constexpr double phi_transition = 0.1;

      if (delta_phi < phi_transition) {
        // take a cross product of the n's (normalised), which
        // is simply equal to sin(delta_phi)
        double cross = first->_nx * other->_ny - other->_nx * first->_ny;
        // the sign can come out negative, but this isn't an issue
        // because we will use it in a context where the sign
        // disappears
        delta_phi = asin(cross);
      }

      // omega = 0: we use deltaR^2 explicitly
      double omega = first->_info->omega;
      if (omega == 0.0) {
        return delta_y * delta_y + delta_phi * delta_phi;
      } else {
        double deltaR2 = delta_y * delta_y + delta_phi * delta_phi;
        if (deltaR2 > _deltaR2_handover)
          return 2 * ((cosh(omega * delta_y) - 1) / (omega * omega) -
                      (cos(delta_phi) - 1));
        else
          return deltaR2;
      }
    }

    double distance(const GHSBriefJet *other_in) {
      // make sure that if at least one of them is a particle, then
      // the first is always a particle
      const GHSBriefJet *first = this;
      const GHSBriefJet *other = other_in;
      if (other->is_particle()) std::swap(first, other);

      // this can only happen if both are jets
      if (first->is_jet()) return numeric_limits<double>::max();

      // return the particle-particle distance (REVISIT FLAVOUR SUM?)
      if (other->is_particle()) {
        // give a flav dependence (only opposite flavours can annihilate)
        FlavInfo flavA = this->_jet.user_info<FlavHistory>().current_flavour();
        FlavInfo flavB =
            other_in->_jet.user_info<FlavHistory>().current_flavour();
        int assA = this->associated_jet();
        int assB = other_in->associated_jet();
        if ((flavA.is_flavourless() || flavB.is_flavourless()) &&
            (assA != assB)) {
          return numeric_limits<double>::max();
        }
        FlavInfo flavsum = flavA + flavB;
        _info->flav_recombiner.apply_summation_choice(flavsum);
        if (!flavA.is_flavourless() && !flavB.is_flavourless() &&
            !flavsum.is_flavourless()) {
          return numeric_limits<double>::max();
        } else {
          return dij(first, other);
        }
      } else {
        // ---- other must be a jet -----
        // first check to see if this particle is associated with any jet at all
        // (if not there is no distance)
        if (first->associated_jet() < 0) return numeric_limits<double>::max();
        // then see if it's associated with specifically with other
        if (_info->jets[first->associated_jet()].cluster_hist_index() ==
            other->_jet.cluster_hist_index()) {
          return dij(first, other);
        } else {
          return numeric_limits<double>::max();
        }
      }
    }
    double beam_distance() const { return _diB; }

    static double dij(const GHSBriefJet *first, const GHSBriefJet *other) {
      double alpha = first->_info->alpha;
      double ptf = sqrt(first->_pt2);
      double pto = sqrt(other->_pt2);
#ifdef VVERBOSE
      cout << " i (pt rap phi ui) = " << sqrt(first->_pt2) << " " << first->_rap
           << " " << first->_phi << " " << first->is_jet() << endl
           << " j (pt rap phi ui) = " << sqrt(other->_pt2) << " " << other->_rap
           << " " << other->_phi << " " << other->is_jet() << endl
           << " with distance = "
           << geometrical_distance_squared(first, other) *
                  max(pow(ptf, alpha), pow(pto, alpha)) *
                  min(pow(ptf, 2 - alpha), pow(pto, 2 - alpha))
           << endl;
#endif
      return geometrical_distance_squared(first, other) *
             max(pow(ptf, alpha), pow(pto, alpha)) *
             min(pow(ptf, 2 - alpha), pow(pto, 2 - alpha));
    }

   private:
    PseudoJet _jet;
    GHSInfo *_info;
    double _diB, _pt2, _rap, _phi, _nx, _ny;
  };

  LimitedWarning ghs_fiducial_warning(10);

  std::vector<PseudoJet> run_GHS(const std::vector<PseudoJet> &jets_base,
                                 double ptcut, double alpha, double omega,
                                 const FlavRecombiner &flav_recombiner) {
    if (jets_base.size() == 0) return jets_base;

#ifdef VERBOSE
    cout << " -- apply pt threshold & define particles -- " << endl;
#endif

    const ClusterSequence &cs = *(jets_base[0].associated_cs());
    vector<PseudoJet> inputs(cs.jets().begin(),
                              cs.jets().begin() + cs.n_particles());

    // ghs_fiducial_warning.warn("Watch out: in run_GHS a fiducial selector is
    // being applied to jets_base (this should be removed; but then run_GHS needs
    // to get a list of all particles)");

    Selector select_pt = SelectorPtMin(ptcut);
    vector<PseudoJet> jets_hard = select_pt(jets_base);

    if (jets_hard.size() == 0) {
      return jets_hard;
    }

#ifdef VERBOSE
    cout << " -- dressing jets -- " << endl;
#endif

    vector<PseudoJet> all;
    vector<PseudoJet> final_jets = jets_hard;
    vector<FlavInfo> jet_flavs(final_jets.size(), 0);
    int njets = final_jets.size();

    GHSInfo ghs_info;
    ghs_info.jets = final_jets;
    ghs_info.njets = final_jets.size();
    ghs_info.alpha = alpha;
    ghs_info.omega = omega;
    ghs_info.flav_recombiner = flav_recombiner;

    // make sure that jets have no flavour, and give them an index meaning that
    // they are a jet
    for (auto &j : final_jets) {
      j.set_user_info(new FlavHistory(FlavInfo(0)));
      // FlavInfo flavj = j.user_info<FlavHistory>().current_flavour();
      // flavj.update_flavourless_attribute();
      j.set_user_index(1);  //< the signal that it is a jet
      all.push_back(j);
    }
    for (auto &c : inputs) {
      //> value <= 0 indicates that we're dealing with a flavour cluster
      //>   0:  flavour cluster not associated with a jet
      //>  -i:  flavour cluster associated with jet `i`
      c.set_user_index(0);  //> start with an unassociated cluster
      for (unsigned i = 0; i < final_jets.size(); i++) {
        if (c.is_inside(final_jets[i])) {
          c.set_user_index(-1 -
                           i);  //> need -1 offset since counting starts at 0
          break;
        }
      }
      all.push_back(c);
    }

    if (all.size() == 0) {
      return all;
    }

#ifdef VERBOSEDIJ
    cout << "initial jets + inputs: (0 = cluster, 1 = jet)" << endl;
    for (unsigned int i = 0; i < all.size(); i++) {
      cout << i << " ";
      print_PJ(&cout, all[i], 5, true, true);
      cout << endl;
    }
#endif
    // all.insert(all.end(), inputs.begin(), inputs.end());
    //  set up nnh
    NNH<GHSBriefJet, GHSInfo> nnh(all, &ghs_info);
    int iA, iB;
    while (njets >
           0) {  // the loop does not change njets, but njets > 0 is necessary
                 // given that the selector could cut off all jets

      double dij = nnh.dij_min(iA, iB);
#ifdef VERBOSEDIJ
      cout << setprecision(12);
      cout << "dij = " << dij << " between " << iA << " and " << iB << endl;
#endif

      // LS-2023-02-10: not sure this is very safe...
      // if (dij > 0.9*numeric_limits<double>::max()) {
      if (dij == numeric_limits<double>::max()) {
        break;
      }
      if (iB >= 0) {
        if (iA > iB) std::swap(iA, iB);
        // we must never have two jets
        assert(iB >= njets && "second entry must be a particle");
        if (iA < njets) {
          // if the first is a jet, assign B's flavour to A and then remove B
          // (note that through the shared pointer, this also affects the
          // flavour of the objects in the NNH object -- which is dangerous --
          // one should really remove the jet and add it back in)
          // FlavInfo * flav_info =
          // dynamic_cast<FlavInfo*>(final_jets[iA].user_info_shared_ptr().get());
          //*flav_info = *flav_info + all[iB].user_info<FlavInfo>();
          // FlavInfo flavA =
          // final_jets[iA].user_info<FlavHistory>().current_flavour();
          FlavInfo flavB = all[iB].user_info<FlavHistory>().current_flavour();
          jet_flavs[iA] = jet_flavs[iA] + flavB;
          flav_recombiner.apply_summation_choice(jet_flavs[iA]);
          // final_jets[iA].set_user_info(new FlavHistory(flavA + flavB));
          nnh.remove_jet(iB);
#ifdef VERBOSE
          cout << "adding flavour of " << iB << " to " << iA << ", removing "
               << iB << endl;
#endif
        } else {
          //> iA & iB are both flavour inputs
          // merge
          PseudoJet new_pseudojet = all[iA];
          new_pseudojet.reset_momentum(
              all[iA] + all[iB]);  //<- resetting only the momentum keeps the
          //> determine the jet association for the merged cluster:
          //> can only be associated with a jet if *both* inputs were
          // associated with the *same* jet
          if (all[iA].user_index() == all[iB].user_index()) {
            new_pseudojet.set_user_index(all[iA].user_index());
          } else {
            new_pseudojet.set_user_index(0);
          }
          FlavInfo flav = FlavHistory::current_flavour_of(all[iA]) +
                          FlavHistory::current_flavour_of(all[iB]);
          flav_recombiner.apply_summation_choice(flav);
          /// set FlavInfo attribute
          new_pseudojet.set_user_info(new FlavHistory(flav));
#ifdef VERBOSE
          cout << "two flavoured inputs combine." << endl;
          cout << "iA = " << iA;
          print_PJ(&cout, all[iA], 12, true, true);
          cout << endl;
          cout << "iB = " << iB;
          print_PJ(&cout, all[iB], 12, true, true);
          cout << endl;
          cout << "into ";
          print_PJ(&cout, new_pseudojet, 12, true, true);
          cout << endl;
#endif
          all.push_back(new_pseudojet);
          nnh.merge_jets(iA, iB, new_pseudojet, all.size() - 1);
        }
      } else {
        //   assert(iA >= njets && "for beam clustering, iA must be a
        //   particle");
#ifdef VERBOSE
        cout << "removing " << iA << endl;
#endif
        nnh.remove_jet(iA);
      }
    }
    for (unsigned i = 0; i < final_jets.size(); i++) {
      jet_flavs[i].update_flavourless_attribute();
      final_jets[i].set_user_info(new FlavHistory(FlavInfo(jet_flavs[i])));
      // restore user index to what it was
      final_jets[i].set_user_index(jets_hard[i].user_index());
    }

#ifdef VERBOSE
    cout << "final dressed jets = " << endl;
    for (auto &a : final_jets) {
      print_PJ(&cout, a, 5, true, true);
      cout << endl;
    }
#endif

    return final_jets;
  }

} // namespace contrib

FASTJET_END_NAMESPACE
