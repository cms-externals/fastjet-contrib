#include "fastjet/contrib/FlavInfo.hh"
#include "fastjet/contrib/MassFlav.hh"
#include "fastjet/contrib/FlavNeutraliser.hh"

#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SharedPtr.hh"
#endif

#include <limits>

#include <numeric>
#include <tuple>

using namespace std;

//#define DEBUG

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

LimitedWarning _warn_old_measure(100);
LimitedWarning _warn_aktlike_measure(100);

namespace contrib{
  
// Set the scale at which to switch from standard measures (e.g. cosphi_coshy)
// to delta_R measure for small delta_R
const double FlavNeutraliser::_deltaR2_handover =
    pow(std::numeric_limits<double>::epsilon(), 0.5);

/// returns true if there is flavour to neutralise betwen jets j & k
bool have_flavour_to_neutralise(const PseudoJet & j, const PseudoJet & k, bool modulo_2) {
  const FlavInfo & flav1 = FlavHistory::current_flavour_of(j);
  const FlavInfo & flav2 = FlavHistory::current_flavour_of(k);

  if (!modulo_2) {
    for (unsigned i = 1; i <= 6; i++) {
      if (flav2.operator[](i) * flav1.operator[](i) < -0) return true;
    }
  } else {
    for (unsigned i = 1; i <= 6; i++) {
      if (flav2.operator[](i) == 1 && flav1.operator[](i) == 1) return true;
    }
  }
  return false;
}


/// The distance measure
double FlavNeutraliser::neutralisation_distance(
       const PseudoJet &p1, const PseudoJet &p2, double ref_scale) 
                                    const {
  if (_spherical_algo) {
    // Tried to do with single normalisation but sqrt() precision gives
    // failures, ie sqrt(p1.modp2()*p2.modp2()) !=
    // sqrt(p1.modp2())*sqrt(p2.modp2()). The below still uses 2 sqrt()
    // functions in modp()

    double norm = 1.0 / (p1.modp() * p2.modp());
    double one_minus_cos_theta =
        1.0 - (p1.px()*p2.px() + p1.py()*p2.py() + p1.pz()*p2.pz()) * norm;
    // safety checks when the angle gets small
    if (one_minus_cos_theta*one_minus_cos_theta < std::numeric_limits<double>::epsilon()) {
      double cx = p1.py() * p2.pz() - p2.py() * p1.pz();
      double cy = p1.pz() * p2.px() - p2.pz() * p1.px();
      double cz = p1.px() * p2.py() - p2.px() * p1.py();
      double sin2theta = (cx*cx + cy*cy + cz*cz) * norm*norm;
      one_minus_cos_theta = sin2theta/2;
    }

    double u;
    if (_measure == aktlike_pair_refratio) {
      _warn_aktlike_measure.warn("FlavNeutraliser::neutralisation_distance: using aktlike_pair_refratio, which is not validated");
      // default to the anti-kt like measure
      u = 2*one_minus_cos_theta/pow(std::max(p1.E(), p2.E()),2);
      // if both particles are flavoured, multiply by Emax^4/Q^4
      if ( !FlavHistory::current_flavour_of(p1).is_flavourless()
        && !FlavHistory::current_flavour_of(p2).is_flavourless()) {
          u *= pow(std::max(p1.E(),p2.E())/ref_scale, 4);
      }
    } else if (_measure == jade || _measure == jadea2) {
      u = 2 * p1.E() * p2.E() * one_minus_cos_theta;
    } else if (_measure == maxscale) {
      u = 2 * pow(std::max(p1.E(), p2.E()), 2) * one_minus_cos_theta;
    } else if (_measure == general) {
      u = 2 * pow(std::max(p1.E(), p2.E()), 2*_p) * pow(std::min(p1.E(), p2.E()), 2*_q) * one_minus_cos_theta;
    } else {
      // all other options 
      _warn_old_measure.warn("FlavNeutraliser::neutralisation_distance: using deprecated old ratio measure");
      u = 2 * pow(std::max(p1.E(), p2.E()) / std::min(p1.E(), p2.E()), 2*_pp) 
          * one_minus_cos_theta;
    }
    return u;

  } else {
    assert(_pp==1);
    // suggestion; first calculate the ratio of squared pts
    // (this avoids square roots -- hopefully should not
    // trigger under/overflow)
    double p1t2 = p1.pt2();
    double p2t2 = p2.pt2();
    double pt2ratio = p1t2 > p2t2 ? p1t2 / p2t2 : p2t2 / p1t2;
    double u, maxpt2, maxpt2p, minpt2q, p1tp2t;
    double drap = fabs(accurate_rap(p1) - accurate_rap(p2));
    double dphi = accurate_absdphi(p1,p2);

    //double old_deltaR2 = p1.squared_distance(p2);
    double deltaR2 = drap*drap + dphi*dphi;

    switch (_measure) {  
      case general:
        maxpt2p = pow(std::max(p1.pt2(), p2.pt2()),_p);
        minpt2q = pow(std::min(p1.pt2(), p2.pt2()),_q);
        if (deltaR2 > _deltaR2_handover) return maxpt2p*minpt2q * 2*((cosh(_a*drap) - 1)/(_a*_a) 
                                                                    - (cos(dphi) - 1));
        else                             return maxpt2p*minpt2q * deltaR2;  

      case sinh_delta_R:
        _warn_old_measure.warn("FlavNeutraliser::neutralisation_distance: using deprecated old ratio measure");
        if (deltaR2 > _deltaR2_handover) return pt2ratio * pow(sinh(sqrt(deltaR2)), 2);
        else                             return pt2ratio * deltaR2;
      case delta_R:
        _warn_old_measure.warn("FlavNeutraliser::neutralisation_distance: using deprecated old ratio measure");
        return pt2ratio * deltaR2;
      case jade_delta_R:
        p1tp2t = sqrt(p1t2*p2t2);
        return p1tp2t * deltaR2;
      case maxscale_delta_R:
        maxpt2 = std::max(p1.pt2(), p2.pt2());
        return maxpt2 * deltaR2;
      case phi2_coshy:
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        if (deltaR2 > _deltaR2_handover) return pt2ratio * (pow(dphi, 2) + 2*(cosh(drap) - 1));
        else                             return pt2ratio * deltaR2;
      case cosphi_coshy:
        _warn_old_measure.warn("FlavNeutraliser::neutralisation_distance: using deprecated old ratio measure");
        //cout << "accr: " << drap << " " << dphi << " " << new_deltaR2 << "\n";
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        //cout << "orig: " << drap << " " << dphi << " " << deltaR2 << "\n\n";
        // _deltaR2_handover is defined to be sqrt(epsilon), where epsilon is 
        // machine precision
        if (deltaR2 > _deltaR2_handover) return pt2ratio * 2*(cosh(drap)-cos(dphi));
        else                             return pt2ratio * deltaR2;
      case aktlike_pair_refratio:
      case aktlike_pair_dynrefratio:
        _warn_aktlike_measure.warn("FlavNeutraliser::neutralisation_distance: using aktlike_pair_refratio, which is not validated");
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        maxpt2 = std::max(p1.pt2(), p2.pt2());
        u = 1.0/maxpt2;
        if (deltaR2 > _deltaR2_handover) u *= 2*(cosh(drap)-cos(dphi));
        else                             u *= deltaR2;
        if ( !FlavHistory::current_flavour_of(p1).is_flavourless()
          && !FlavHistory::current_flavour_of(p2).is_flavourless()) {
            u *= pow(maxpt2/pow(ref_scale,2), 2);
        }
        return u;
      case jade:
        // a hadron-collider version of the JADE (without the correction factor a)
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        p1tp2t = sqrt(p1t2*p2t2);
        if (deltaR2 > _deltaR2_handover) return p1tp2t * 2*(cosh(drap)-cos(dphi));
        else                             return p1tp2t * deltaR2;

      case jadea2:
        // a hadron-collider version of the JADE (with a = 2)
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        p1tp2t = sqrt(p1t2*p2t2);
        //cout << "p1tp2t: " << p1tp2t << ", deltaR2 = " << deltaR2 << ", u=" << p1tp2t * deltaR2 << "\n";
        //cout << "drap = " << drap << ", dphi = " << dphi << "\n";
        if (deltaR2 > _deltaR2_handover) return p1tp2t * 2*( (cosh(2*drap)-1)/4 + (1-cos(dphi)) );
        else                             return p1tp2t * deltaR2;
        
      case maxscale:
        // use the max(pt^2) option
        // dphi = p1.delta_phi_to(p2);
        // drap = std::abs(p1.rap() - p2.rap());
        maxpt2 = std::max(p1.pt2(), p2.pt2());
        if (deltaR2 > _deltaR2_handover) return maxpt2 * 2*(cosh(drap)-cos(dphi));
        else                             return maxpt2 * deltaR2;

      default:
        throw Error("Unrecognised neutralisation measure");
    }
  }
}


//----------------------------------------------------------------------
void neutralise_flavour(PseudoJet & j, PseudoJet & k, int hist_step,
                        bool modulo_2) {
  FlavInfo flav1 = FlavHistory::current_flavour_of(j);
  FlavInfo flav2 = FlavHistory::current_flavour_of(k);

  if (!modulo_2) {
    // if using normal flavour (non-mod2) then check that flavours are
    // opposite, i.e. flav1*flav2 < 0 in order to neutralise
    for (unsigned i = 1; i <= 6; i++) {
      // Change this to account for both cases at once
      if (flav2.operator[](i) * flav1.operator[](i) < -0 &&
          (abs(flav2.operator[](i)) <= abs(flav1.operator[](i)))) {
        flav1._flav_content[i] =
            flav2._flav_content[i] + flav1._flav_content[i];
        flav2._flav_content[i] = 0;
      }
      if (flav2.operator[](i) * flav1.operator[](i) < -0 &&
          (abs(flav2.operator[](i)) > abs(flav1.operator[](i)))) {
        flav2._flav_content[i] =
            flav2._flav_content[i] + flav1._flav_content[i];
        flav1._flav_content[i] = 0;
      }
    }
  } else {
    // if using mod2 flavour then each flavour is either 1 or 0 so just
    // require flavours to both be 1 in order to neutralise
    for (unsigned i = 1; i <= 6; i++) {
      if (flav2.operator[](i) == 1 && flav1.operator[](i) == 1) {
        flav1._flav_content[i] = 0;
        flav2._flav_content[i] = 0;
      }
    }
  }
  flav1.update_flavourless_attribute();
  flav2.update_flavourless_attribute();

  dynamic_cast<FlavHistory *>(j.user_info_shared_ptr().get())
      ->update_flavour_history(flav1, hist_step);
  dynamic_cast<FlavHistory *>(k.user_info_shared_ptr().get())
      ->update_flavour_history(flav2, hist_step);
}

//----------------------------------------------------------------------
// Compare flavour of each jet from a pair of vector<fastjet::PseudoJet>
// objects. If there are >max_jets jets in one vector, only the first max_jets
// terms are compared (in order of decreasing hardness).

vector<PseudoJet> sorted_by_px(const vector<PseudoJet> & jets) {
   vector<double> px(jets.size());
   for (size_t i = 0; i < jets.size(); i++) {px[i] = jets[i].px();}
   return objects_sorted_by_values(jets, px);
 }

bool jet_flavour_compare(const vector<fastjet::PseudoJet> &j,
                         const vector<fastjet::PseudoJet> &k,
                         const int max_jets, bool sort_by_px) {
    vector<fastjet::PseudoJet> jets_j = sorted_by_E(j);
    vector<fastjet::PseudoJet> jets_k = sorted_by_E(k);
    if(sort_by_px){
      jets_j = sorted_by_px(j);
      jets_k = sorted_by_px(k);
    }
  unsigned max_count =
      max_jets < 0
          ? std::min({int(jets_j.size()), int(jets_k.size())})
          : std::min({int(jets_j.size()), int(jets_k.size()), max_jets});
  for (unsigned i = 0; i < max_count; i++) {
    assert(jets_j[i].has_user_info() && jets_k[i].has_user_info());
    FlavInfo flav1 = jets_j[i].user_info<FlavHistory>().current_flavour();
    FlavInfo flav2 = jets_k[i].user_info<FlavHistory>().current_flavour();
    flav1.update_flavourless_attribute();
    flav2.update_flavourless_attribute();
    if (flav1 != flav2) {
      return false;
    }
  }
  return true;
}

//----------------------------------------------------------------------
// Compare net flavour of each jet from a pair of vector<fastjet::PseudoJet>
// objects.
bool jet_net_flavour_compare(vector<fastjet::PseudoJet> &j,
                             vector<fastjet::PseudoJet> &k) {
  FlavInfo flav_j, flav_k; //net flavours

  for (unsigned i = 0; i < j.size(); i++) {
      assert(j[i].has_user_info());
      FlavInfo fl = FlavHistory::current_flavour_of(j[i]);
      for (unsigned l=1; l<=6; l++){
        flav_j._flav_content[l] += fl._flav_content[l];
      }
  }

  for (unsigned i = 0; i < k.size(); i++) {
      assert(k[i].has_user_info());
      FlavInfo fl = FlavHistory::current_flavour_of(k[i]);
      for (unsigned l=1; l<=6; l++){
        flav_k._flav_content[l] += fl._flav_content[l];
      }
  }

  if (flav_j!=flav_k) {
    return 0;
  }
  return 1;
}

//----------------------------------------------------------------------
// The method takes a ClusterSequence as argument and goes through
// each step of the declustering to neutralise flavour
std::vector<PseudoJet> FlavNeutraliser::neutralise(ClusterSequence & cs) const {

  // Setup the FlavRecombiner
  unique_ptr<JetDefinition::Recombiner> flav_recombiner;
  flav_recombiner.reset(_use_mass_flav ? new MassFlavRecombiner()
                                       : new FlavRecombiner());
  //flav_recombiner.reset(new FlavRecombiner());


  // Copy across a the list of all pseudojets
  vector<PseudoJet> jets = cs.jets();
  // The clustering history
  const vector<ClusterSequence::history_element> & hist = cs.history();

  // shorthand for the hardness (as of 2022-07-26, not yet being used everywhere)
  auto hardness = _spherical_algo  
      ? [](const PseudoJet & j) {return  j.E();}
      : [](const PseudoJet & j) {return  j.pt();};


  // NOTE -- pp ref scale should evolve to become something
  // less sensitive to UE & pileup, but for now take something simple
  // It is used only for the aktlike_pair_refratio option
  double ref_scale = 0;
  if (_measure == aktlike_pair_refratio) {
    for (unsigned i = 0; i < cs.n_particles(); i++) {
      ref_scale += hardness(jets[i]);
    }
  } else if (_measure == aktlike_pair_dynrefratio) {
    for (unsigned i = 0; i < cs.n_particles(); i++) {
      ref_scale = max(ref_scale, hardness(jets[i]));
    }
  }

  // Loop over each step of the clustering
  for (unsigned ih_step = 0; ih_step < hist.size(); ++ih_step) {

    const ClusterSequence::history_element & hist_step = hist[ih_step];

    // Get the two parent jets' history indices
    int index1 = hist_step.parent1;
    int index2 = hist_step.parent2;

    // Check that it's a valid recombined jet
    if (index1 < 0 || index2 < 0) continue;


    // arranges the indices such that index1 corresponds to the softer jet
    // or correspondingly jet_i below is the softer one
    if (hardness(jets[hist[index1].jetp_index]) 
        > hardness(jets[hist[index2].jetp_index])) std::swap(index1, index2);
    //if (!_spherical_algo) {
    //  if (jets[hist[index1].jetp_index].pt() >
    //      jets[hist[index2].jetp_index].pt())
    //    std::swap(index1, index2);
    //} else {
    //  if (jets[hist[index1].jetp_index].E() >
    //      jets[hist[index2].jetp_index].E())
    //    std::swap(index1, index2);
    //}
    
    PseudoJet & jet_i = jets[hist[index1].jetp_index];
    PseudoJet & jet_j = jets[hist[index2].jetp_index];

    // If i is flavourless, get the child, do the standard recombination and
    // move on to the next history step
    const FlavInfo & flav_i = FlavHistory::current_flavour_of(jet_i);
    if (flav_i.is_flavourless()) {
      int cluster_hist_index_temp = jets[hist_step.jetp_index].cluster_hist_index();
      flav_recombiner->recombine(jet_i, jet_j, jets[hist_step.jetp_index]);
      jets[hist_step.jetp_index].set_cluster_hist_index(cluster_hist_index_temp);
      dynamic_cast<FlavHistory *>(
          jets[hist_step.jetp_index].user_info_shared_ptr().get())
          ->amend_last_history_index(ih_step);
      if (_modulo_2) {
        dynamic_cast<FlavHistory *>(
            jets[hist_step.jetp_index].user_info_shared_ptr().get())
            ->apply_modulo_2();
      }
      // if relevant, update the dynamic reference scale before continuing with the loop
      if (_measure == aktlike_pair_dynrefratio) 
           ref_scale = max(ref_scale, hardness(jets[hist_step.jetp_index]));
      continue;
    }

    const double uij = neutralisation_distance(jet_i, jet_j, ref_scale);
    if (_writeout_uijs) cout << "uij for pair about to cluster is " 
                             << uij << ", i,j=" << jet_i.cluster_hist_index() 
                             << ", " << jet_j.cluster_hist_index() << endl;

    // If not, identify all flavoured PseudoJets and their distance to
    // the current jet i getting recombined, and put them in a list K
    std::vector<std::pair<PseudoJet*, double>> flavour_candidates;

    // note: this loop is not especially efficient -- it may be
    // worth exploring some kind of dynamically maintained list of 
    // flavoured objects (on the other hand the loop probbaly isn't
    // used all that often)
    for (unsigned k = 0; k < jets.size(); ++k) {

      PseudoJet & jet_k = jets[k];

      // Check k != i (i.e. not itself)
      if (int(k) == hist[index1].jetp_index) continue;     

      // Check k != j to stop needless neutralisation of i and j before 
      // recombination
      if (int(k) == hist[index2].jetp_index) continue;

      // If the jet is not currently around, skip
      // First check it has already been created
      if (jets[k].cluster_hist_index() >= int(ih_step)) continue;

      // then check that it wasn't recombined meanwhile
      if ( (hist[hist[jets[k].cluster_hist_index()].child].parent2 != -1)
           && (hist[jets[k].cluster_hist_index()].child < int(ih_step)) ) continue;

      // Get the flavour info
      const FlavInfo & flav = FlavHistory::current_flavour_of(jet_k);

      // If it carries flavour, compute the distance u_ik and add it to the list
      if (!flav.is_flavourless()) {
        //double uik = neutralisation_distance(jet_i, *j);
        //flavour_candidates.push_back(std::make_pair(j, uik));
        flavour_candidates.push_back(std::make_pair(&jet_k, 0.0));
      }
    }

    if (recursive()) {
      // it will be useful to have jet j in the list as part of the recursion
      flavour_candidates.push_back(std::make_pair(&jet_j, 0.0));
      use_neutralisation_candidates_recursive(jet_i, uij, ih_step, flavour_candidates, ref_scale, &jet_j);
    } else {
      use_neutralisation_candidates(jet_i, uij, ih_step, flavour_candidates, ref_scale);
    }

//    // Sort the list in decreasing order of u_ik
//    sort(flavour_candidates.begin(), flavour_candidates.end(),
//        [](std::pair<PseudoJet*,double> & a, std::pair<PseudoJet*,double> & b) {return a.second > b.second;});
//
//    // Now, loop over all jets in the list K
//    while (!flavour_candidates.empty()) {
//
//      // Get the last element (the one with the smaller u_ik)
//      const std::pair<PseudoJet*,double> current_neutraliser = flavour_candidates.back();
//
//      // If u_ik >= u_ij, we're done
//      if (current_neutraliser.second >= uij) break;
//
//      // Otherwise, take the corresponding pseudojet k, and neutralise as much flavour from jet_i
//      // as one can with k
//      neutralise_flavour(jet_i, *current_neutraliser.first, ih_step, _modulo_2);
//
//      FlavInfo flav_i = FlavHistory::current_flavour_of(jet_i);
//      if (!flav_i.is_flavourless()) {
//        flavour_candidates.pop_back();
//      } else {
//        // i is flavourless, so go to recombiner step below while loop
//        break;
//      }
//
//    }

    // If i is flavourless, the set is empty, or all remaining k have u_ik > u_ij, recombine
    int cluster_hist_index_temp = jets[hist_step.jetp_index].cluster_hist_index();
    assert(cluster_hist_index_temp == int(ih_step));
    flav_recombiner->recombine(jet_i, jet_j, jets[hist_step.jetp_index]);
    jets[hist_step.jetp_index].set_cluster_hist_index(cluster_hist_index_temp);
    dynamic_cast<FlavHistory *>(
        jets[hist_step.jetp_index].user_info_shared_ptr().get())
        ->amend_last_history_index(ih_step);
    if (_modulo_2) {
      dynamic_cast<FlavHistory *>(
          jets[hist_step.jetp_index].user_info_shared_ptr().get())
          ->apply_modulo_2();
    }

    // if relevant, update the dynamic reference scale before continuing with the loop
    if (_measure == aktlike_pair_dynrefratio) 
          ref_scale = max(ref_scale, hardness(jets[hist_step.jetp_index]));
  }

  return jets;

}

void FlavNeutraliser::use_neutralisation_candidates(
    PseudoJet & jet_i,
    double uij,
    int ih_step,
    std::vector<std::pair<PseudoJet*, double>> & flavour_candidates, double ref_scale) const {

  // set up the uik neutralisation distances
  for (auto & fc: flavour_candidates) {
    fc.second = neutralisation_distance(jet_i, *(fc.first), ref_scale);
  }

  // Sort the list in decreasing order of u_ik
  sort(flavour_candidates.begin(), flavour_candidates.end(),
      [](std::pair<PseudoJet*,double> & a, std::pair<PseudoJet*,double> & b) {return a.second > b.second;});

  // Now, loop over all jets in the list K
  while (!flavour_candidates.empty()) {

    // Get the last element (the one with the smaller u_ik)
    const std::pair<PseudoJet*,double> current_neutraliser = flavour_candidates.back();

    // If u_ik >= u_ij, we're done
    if (current_neutraliser.second >= uij) break;

    // Otherwise, take the corresponding pseudojet k, and neutralise as much flavour from jet_i
    // as one can with k
    //if (have_flavour_to_neutralise(jet_i, *current_neutraliser.first, _modulo_2)) {
    neutralise_flavour(jet_i, *current_neutraliser.first, ih_step, _modulo_2);
    //}

    FlavInfo flav_i = FlavHistory::current_flavour_of(jet_i);
    if (!flav_i.is_flavourless()) {
      flavour_candidates.pop_back();
    } else {
      // i is flavourless, so go to recombiner step below while loop
      break;
    }

  }

}

/// This is an attempt to think about a potential more global
/// neutralisation procedure.
///
/// Suppose we have a situation  
///
///  E| 0                               0 = u
///   |                                 1 = ubar
///   |                                 2 = ubar
///   |    1     2 3                    3 = u
///   +---------------------------------
///                                  rap
///
///  - jet_i = 1, jet_j that it will cluster with is 0 and u10 is quite large
///  - 1 would find it favourable to neutralise with 3, u13 < u10
///    (u12 < u13, but neutralisation not possible)
///
///  - before performing the 13 neutralisation, we check if 3 has
///    a preferred neutralisation candidate, one with u3x < u13, and
///    in this case u23 < u13, so 2 & 3 neutralise each other, and then
///    we return to the search for a neutralisation partner for 1
///
///  A slightly non-trivial question is the set of particles over which
///  we do the search for 3's potential neutralisation partners. Should
///  it be all flavoured particles, including 0? How do we organise the
///  corresponding information? The conclusion on 2022-05-29
///  is that 
///
///  - when we examine neutralisation candidates for 1, we should
///    not consider 0 (because neutralisation will happen when they
///    recombine); 
///
///  - when we examine neutralisation candidates for 3, we add 0
///    back into the mix, but we do not consider 1 (and if we find neutralisation
///    candidates for 3 and recurse again, we do not consider 1 or 3)
///
/// The exclude argument serves solely to exclude the original "0"
void FlavNeutraliser::use_neutralisation_candidates_recursive(
    PseudoJet & jet_i,
    double uij,
    int ih_step,
    std::vector<std::pair<PseudoJet*, double>> & flavour_candidates,
    double ref_scale, 
    const PseudoJet * jet_j_exclude) const {

  // set up the uik neutralisation distances
  for (auto & fc: flavour_candidates) {
    fc.second = neutralisation_distance(jet_i, *(fc.first), ref_scale);
    if (_writeout_uijs) {
      cout << "in recursive step: u" << jet_i.cluster_hist_index() << "," 
           << fc.first->cluster_hist_index() << " = " << fc.second << std::endl;
    }
  }

  // Sort the list in order of increasing u_ik
  sort(flavour_candidates.begin(), flavour_candidates.end(),
      [](std::pair<PseudoJet*,double> & a, std::pair<PseudoJet*,double> & b) {return a.second < b.second;});

  // Now, loop over all jets in the list K
  for (auto current_neutraliser: flavour_candidates) {

    //
    if (current_neutraliser.first == jet_j_exclude) continue;

    // If u_ik >= u_ij, we're done
    if (current_neutraliser.second >= uij) break;

    // Otherwise, take the corresponding pseudojet k, and neutralise as
    // much flavour from jet_i as one can with k
    if (have_flavour_to_neutralise(jet_i, *current_neutraliser.first, _modulo_2)) {
      // but before doing that, we will allow k to explore its potential neutralisers
      // so that we avoid stealing the flavour of k from something that should have
      // preferentially have gone to k
      //
      // start by copying the candidates and removing the current one
      std::vector<std::pair<PseudoJet*,double>> flavour_candidates_copy;
      flavour_candidates_copy.reserve(flavour_candidates.size() - 1);
      for (auto & fc: flavour_candidates) {
        if (fc.first != current_neutraliser.first) flavour_candidates_copy.push_back(fc);
      }
      // then run the neutralisation procedure on anything that remains
      use_neutralisation_candidates_recursive(*current_neutraliser.first, 
        current_neutraliser.second, ih_step, flavour_candidates_copy, ref_scale);

      // once we are done with the recursive call neutralise anything that 
      // is still there to be neutralised
      neutralise_flavour(jet_i, *current_neutraliser.first, ih_step, _modulo_2);
    }

    const FlavInfo & flav_i = FlavHistory::current_flavour_of(jet_i);
    // if i is flavourless, there's nothing left to do, so we exit
    if (flav_i.is_flavourless()) break;
  }

}


} // namespace contrib
FASTJET_END_NAMESPACE
