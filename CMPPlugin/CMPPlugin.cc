#include "fastjet/contrib/CMPPlugin.hh"

#include <sstream>
#include "fastjet/contrib/FlavInfo.hh"
// to facilitate use with fjcore
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/NNH.hh"
#endif

FASTJET_BEGIN_NAMESPACE  // defined in fastjet/internal/base.hh

using namespace std;

// Transition points to use small rap/phi formulae
constexpr double rap_transition = 0.1;
constexpr double phi_transition = 0.1;

// Transition point to switch from Omega_ik^2 -> DeltaR_ik^2
const double CMPPlugin::_deltaR2_handover =
    pow(std::numeric_limits<double>::epsilon(), 0.5);

/// Info class for particle-independent information
class CMPNNInfo {
public:
  CMPNNInfo(double a, double R, double deltaR2_handover,
            CMPPlugin::CorrectionType correction_type,
            CMPPlugin::ClusteringType clustering_type,
            bool spherical)
     : _a(a),
       _R(R),
       _deltaR2_handover(deltaR2_handover),
       _correction_type(correction_type),
       _clustering_type(clustering_type),
       _spherical(spherical) {}
 double R() const { return _R; }
 double a() const { return _a; }
 double deltaR2_handover() const { return _deltaR2_handover; }
  CMPPlugin::CorrectionType correction_type() const { return _correction_type; }
  CMPPlugin::ClusteringType clustering_type() const { return _clustering_type; }
 bool is_spherical() const { return _spherical; }
private:
  double _a, _R, _deltaR2_handover;
  CMPPlugin::CorrectionType _correction_type;
  CMPPlugin::ClusteringType _clustering_type;
  bool _spherical;
};

/// Convenience class to run CMP plugin algorithm
class CMPBriefJet {
public:
  void init(const PseudoJet & jet, CMPNNInfo * info) {
    // relevant jet kinematics
    _E = jet.E();
    _px = jet.px();
    _py = jet.py();
    _pz = jet.pz();
    _modp2 = jet.modp2();
    _m2 = jet.m2();

    _kt = jet.pt();
    _phi = jet.phi();
    _rap = jet.rap();
    _pt2 = jet.pt2();
    double pt = sqrt(_pt2);
    _nx = jet.px() / pt;
    _ny = jet.py() / pt;
    if (fabs(_rap) < rap_transition) {
      // for moderately small rapidities switch to special rapidity formula
      //
      // rap = 1/2 log[(E+pz)/(E-pz)]
      //     = 1/2 log[1 + 2pz/(E-pz)]
      // 
      // and use log1p for the latter
        _rap = 0.5 * log1p(2*_pz/(_E - _pz));
  
      }

    // algorithm parameters
    _R = info->R();
    _a = info->a();
    _deltaR2_handover = info->deltaR2_handover();
    _correction_type = info->correction_type();
    _clustering_type = info->clustering_type();
    _spherical = info->is_spherical();

    // jet flavour info
    _flavoured = (!jet.user_info<contrib::FlavHistory>()
                       .current_flavour()
                       .is_flavourless());
    _flavour = (jet.user_info<contrib::FlavHistory>()
                       .current_flavour());
  }

  //----------------------------------------------------------------------
  /// Returns the 3-vector cross-product of p1 and p2. If lightlike is false
  /// then the energy component is zero; if it's true the the energy component
  /// is arranged so that the vector is lighlike
  void cross_product(const CMPBriefJet * p1, const CMPBriefJet * p2,
                     PseudoJet & retp, bool lightlike=false) const {
    double px = p1->_py * p2->_pz - p2->_py * p1->_pz;
    double py = p1->_pz * p2->_px - p2->_pz * p1->_px;
    double pz = p1->_px * p2->_py - p2->_px * p1->_py;

    double E;
    if (lightlike) {
      E = sqrt(px*px + py*py + pz*pz);
    } else {
      E = 0.0;
    }

    // is this exactly what we want? The masses seem to change, need to check!
    retp.reset(px, py, pz, E);

    return;
  }

  double dot_product_3d(const CMPBriefJet * a, const CMPBriefJet * b) const {
    return a->_px*b->_px + a->_py*b->_py + a->_pz*b->_pz;
  }
  double dot_product(const CMPBriefJet * a, const CMPBriefJet * b) const {

    return a->_E*b->_E - dot_product_3d(a,b);
  }
  /// Returns (1-cos theta) where theta is the angle between p1 and p2
  double one_minus_costheta(const CMPBriefJet * p1, const CMPBriefJet * p2) const {

    if (p1->_m2 == 0 && p2->_m2 == 0) {
      // use the 4-vector dot product.
      // For massless particles it gives us E1*E2*(1-cos theta)
      double res = dot_product(p1,p2) / (p1->_E * p2->_E);
      return res;
    } else {
      double p1mod = sqrt(p1->_modp2);
      double p2mod = sqrt(p2->_modp2);
      double p1p2mod = p1mod*p2mod;
      double dot = dot_product_3d(p1,p2);

      if (dot > (1-std::numeric_limits<double>::epsilon()) * p1p2mod) {
        PseudoJet cross_result;
        cross_product(p1, p2, cross_result, false);
        // the mass^2 of cross_result is equal to
        // -(px^2 + py^2 + pz^2) = (p1mod*p2mod*sintheta_ab)^2
        // so we can get
        double res = -cross_result.m2()/(p1p2mod * (p1p2mod+dot));

        return res;
      }
      return 1.0 - dot/p1p2mod;

    }
  }

  /// Returns delta_R^2 between jets
  double geometrical_distance(const CMPBriefJet* jet) const {
    double dphi = std::abs(_phi - jet->_phi);
    double deta = (_rap - jet->_rap);
    if (dphi > pi) {
      dphi = twopi - dphi;
    }
    
    if (dphi < phi_transition) {
      // take a cross product of the n's (normalised), which 
      // is simply equal to sin(delta_phi)
      double cross = _nx * jet->_ny - jet->_nx * _ny;
      assert(cross >= -1.0 && cross <= 1.0);
      // the sign can come out negative, but this isn't an issue
      // because we will use it in a context where the sign 
      // disappears
      dphi = asin(cross);
    }
    return dphi * dphi + deta * deta;
  }

  // for now just make it like antikt to ensure everything runs right then
  // change to CMP distance
  double distance(const CMPBriefJet * jet) const {

    // With NNH, we have to recompute the correction factor for all interjet
    // distances, for oppositely flavoured pairs. We do this manually in
    // _NN_clustering(), and set it to infinity for now.
    if ((_flavoured && jet->_flavoured) && (_flavour + jet->_flavour).is_flavourless()){
      return std::numeric_limits<double>::max();
    }
    if (!_spherical) {
      double deltaR2 = geometrical_distance(jet);
      // dij is normal antikt distance measure:
      return (1 / (_kt > jet->_kt ? _kt*_kt : jet->_kt*jet->_kt)) *
                  deltaR2 / pow(_R, 2);
    }
    else {
      // e+e- generalised anti-kt 
      return one_minus_costheta(this, jet) / (1 - cos(_R))
          / (_E > jet->_E ? _E*_E : jet->_E*jet->_E);
    }

  }

  // standard antikt beam distance
  double beam_distance() const {
    if (!_spherical)
      return 1/pow(_kt, 2);
    else
      return std::numeric_limits<double>::max();
  }

private:
  double _E, _px, _py, _pz, _modp2, _m2, _kt, _rap,
         _phi, _R, _a, _deltaR2_handover, _nx, _ny, _pt2;
  bool _flavoured;
  contrib::FlavInfo _flavour;
  CMPPlugin::CorrectionType _correction_type;
  CMPPlugin::ClusteringType _clustering_type;
  bool _spherical;
};


string CMPPlugin::description() const {
  ostringstream desc;
  desc << "CMP plugin with R = " << R() << " and a = " << a();
  if        (_clustering_type == DynamicKtMax) {
    desc << ", reference scale is dynamic";
  } else if (_clustering_type == FixedKtMax) {
    desc << ", reference scale is largest anti-kt jet pt (for spherical, Eref about in flux)";
  } else {
    throw Error("Unrecognised value of _clustering_type");
  }
  if        (_correction_type == NoCorrection) {
    desc << ", original algorithm";
  } else if (_correction_type == CoshyCosPhi) {
    desc << ", with coshy-cosphi correction to the cos term";
  } else if (_correction_type == OverAllCoshyCosPhi) {
    desc << ", with overall coshy-cosphi correction";
  } else if (_correction_type == OverAllCoshyCosPhi_a2) {
    desc << ", with overall coshy-cosphi correction (a = 2)";
  } else if (_correction_type == SqrtCoshyCosPhiArgument) { 
    desc << ", with a sqrt coshy-cosphi correction to the cos argument";
  } else if (_correction_type == SqrtCoshyCosPhiArgument_a2) {
    desc << ", with a sqrt coshy-cosphi (a = 2) correction to the cos argument";
  } else {
    throw Error("Unrecognised value of _correction_type");
  }
  return desc.str();
}

/// calculates the precise squared distance; I think we want to have this for PseudoJets
/// as we must recalculate the distance_opposite_flavour for jets taken from the cs,
/// which are not CMPBriefJets
double CMPPlugin::precise_squared_distance(const PseudoJet & j1, const PseudoJet & j2) const {
  
  double rap1 = j1.rap();
  double rap2 = j2.rap();

  if (fabs(rap1) < rap_transition) {
    rap1 = 0.5*log1p(2*j1.pz()/(j1.E() - j1.pz()));
  }

  if (fabs(rap2) < rap_transition) {
    rap2 = 0.5*log1p(2*j2.pz()/(j2.E() - j2.pz()));
  }

  double drap = rap1 - rap2;
  double dphi = std::fabs(j1.phi() - j2.phi());
  if (dphi > pi) {
    dphi = twopi - dphi;
  }

  if (dphi < phi_transition) {
    double invj1pt = 1.0/j1.pt();
    double invj2pt = 1.0/j2.pt();
    double cross = (j1.px()*invj1pt) * (j2.py()*invj2pt) - (j2.px()*invj2pt) * (j1.py()*invj1pt);
    //double cross = (p1.px() * p2.py() - p2.px() * p1.py())/sqrt(p1.pt2() * p2.pt2());
    assert(cross <= 1.0 && cross >= -1.0);
    //double cross = (j1.px() * j2.py() - j2.px() * j1.py())/(j1.pt()*j2.pt());
    dphi = asin(cross);
  }

  return drap*drap + dphi*dphi;

}

/// returns the distance with CMP correction factor between
/// two pseudojets of opposite flavour
double CMPPlugin::distance_opposite_flavour(const PseudoJet & j1, const PseudoJet & j2,
                                            const double ktmax) const {
  
  FlavInfo flav1 = j1.user_info<FlavHistory>().current_flavour();
  FlavInfo flav2 = j2.user_info<FlavHistory>().current_flavour();
  assert((!flav1.is_flavourless() && !flav2.is_flavourless()) && (flav1+flav2).is_flavourless());

  // RG: numerical instability generated here
  
  double deltaR2 = precise_squared_distance(j1, j2);

  double dFij;
  if (!_spherical) {
    // dij is normal antikt distance measure:
    dFij = (1 / (j1.pt() > j2.pt() ? j1.pt()*j1.pt() : j2.pt()*j2.pt())) *
                deltaR2 / pow(_R, 2);
  }
  else {
    // e+e- generalised anti-kt 
    dFij = one_minus_costheta(j1, j2) / (1 - cos(_R))
        / (j1.E() > j2.E() ? j1.E()*j1.E() : j2.E()*j2.E());
  }

  // if i, j is an oppositely flavoured pair then distance gets multiplied by 
  // Sij factor
  //
  // re-structuring for dynamic ktmax definition
  // With NNH, we have to recompute the correction factor for all interjet
  // distances, for oppositely flavoured pairs. We do this manually in
  // _NN_clustering(), and set it to infinity for now.
  double x;
  if (!_spherical) {
    x = (j1.pt()*j1.pt() + j2.pt()*j2.pt()) / (2 * _a * ktmax*ktmax);
  }
  else {
    x = (j1.E()*j1.E() + j2.E()*j2.E()) / (2 * _a * ktmax*ktmax);
  }

  double dphi = std::fabs(j1.phi() - j2.phi());
  if (dphi > pi) {
    dphi = twopi - dphi;
  }
  double drap = std::fabs(j1.rap() - j2.rap());

  double Sij;
  if (isinf(x)) {
    // otherwise get Sij = nan
    Sij = 1;
  } else if (x*(pi/2) > 1e-4) {
    double correction_factor;
    Sij = 1 - (1 - x > 0 ? 1 : 0) * cos((pi / 2) * x);

    if (_correction_type == CoshyCosPhi) {
      // MLH CORRECTION TO CMP - can add exponential factor to Sij
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;

      Sij = 1 - (1 - x > 0 ? 1 : 0) * cos((pi / 2) * x) * correction_factor;
    } else if (_correction_type == OverAllCoshyCosPhi) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;
      
      Sij *= correction_factor;
    } else if (_correction_type == OverAllCoshyCosPhi_a2) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2*( (cosh(2*drap)-1)/4 + (1-cos(dphi)) ) / deltaR2;
      else
        correction_factor = 1;

      Sij *= correction_factor;
    } else if (_correction_type == SqrtCoshyCosPhiArgument) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;

      double y = sqrt(correction_factor) * x;
      Sij = 1 - (1 - y > 0 ? 1 : 0) * cos((pi / 2) * y);
    } else if (_correction_type == SqrtCoshyCosPhiArgument_a2) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2*( (cosh(2*drap)-1)/4 + (1-cos(dphi)) ) / deltaR2;
      else
        correction_factor = 1;

      double y = sqrt(correction_factor) * x;
      Sij = 1 - (1 - y > 0 ? 1 : 0) * cos((pi / 2) * y);
    }

  } else {
    // small angle expansion (here 1-x > 0 for certain so no need for Heaviside function)
    double correction_factor;
    Sij = (pow(pi * x / 2, 2) / 2);

    if (_correction_type == CoshyCosPhi) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;

      Sij = 1 - (1 - (pi*x/2)*(pi*x/2)/2) * correction_factor;
    } else if (_correction_type == OverAllCoshyCosPhi) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;

      Sij *= correction_factor;
    } else if (_correction_type == OverAllCoshyCosPhi_a2) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2*( (cosh(2*drap)-1)/4 + (1-cos(dphi)) ) / deltaR2;
      else
        correction_factor = 1;

      Sij *= correction_factor;
    } else if (_correction_type == SqrtCoshyCosPhiArgument) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2 * (cosh(drap) - cos(dphi)) / deltaR2;
      else
        correction_factor = 1;

      double y = sqrt(correction_factor) * x;
      // here we know x is small but y can be quite large so only want to
      // use small angle approx for the cosine if y is small also
      if (y*pi/2 > 1e-4) {
        Sij = 1 - (1 - y > 0 ? 1 : 0) * cos((pi / 2) * y);
      } else {
        Sij = (pi*y/2)*(pi*y/2) / 2;
      }
    } else if (_correction_type == SqrtCoshyCosPhiArgument_a2) {
      if (deltaR2 > _deltaR2_handover)
        correction_factor = 2*( (cosh(2*drap)-1)/4 + (1-cos(dphi)) ) / deltaR2;
      else
        correction_factor = 1;

      double y = sqrt(correction_factor) * x;
      // here we know x is small but y can be quite large so only want to
      // use small angle approx for the cosine if y is small also
      if (y*pi/2 > 1e-4) {
        Sij = 1 - (1 - y > 0 ? 1 : 0) * cos((pi / 2) * y);
      } else {
        Sij = (pi*y/2)*(pi*y/2) / 2;
      }
    }
  }
  dFij *= Sij;
  return dFij;

}


void CMPPlugin::run_clustering(ClusterSequence & cs) const {
  // make sure jets are set up with FlavHistory correctly
  for (unsigned i = 0; i < cs.jets().size(); i++) {
    const PseudoJet &jet = cs.jets()[i];
    int hist_index = jet.cluster_hist_index();
    if (jet.has_user_info<contrib::FlavInfo>()) {
      /// it can be useful to be able to start from a FlavInfo
      cs.plugin_non_const_jet(i).set_user_info(
          new contrib::FlavHistory(jet.user_info<contrib::FlavInfo>(), hist_index));
    } else if (jet.has_user_info<contrib::FlavHistory>()) {
      // if we start from a FlavHistory, make sure that we copy the
      // object, copy its current_flavour and assign this CS's
      // hist_index.
      cs.plugin_non_const_jet(i).set_user_info(new contrib::FlavHistory(
          jet.user_info<contrib::FlavHistory>().current_flavour(), hist_index));
    } else {
      throw fastjet::Error(
          "A PseudoJet being clustered with CMPPlugin had neither "
          "FlavInfo nor FlavHistory user_info.");
    }
  }
  CMPNNInfo info(this->a(), this->R(), this->_deltaR2_handover,
                 this->_correction_type, this->_clustering_type,
                 this->_spherical);
  NNH<CMPBriefJet, CMPNNInfo> nnh(cs.jets(), &info);
  _NN_clustering(cs, nnh);
}

template<typename NN>
void CMPPlugin::_NN_clustering(ClusterSequence & cs, NN & nnh) const {
  const vector<PseudoJet> & jets = cs.jets();
  int njets = jets.size();

  double current_ktmax;
  // run normal antikt here to find ktmax (don't need to worry about flavour)
  if (_clustering_type == FixedKtMax) { // fixed ktmax
    if (!_spherical) {
      JetDefinition antikt_jet_def(antikt_algorithm, _R);
      ClusterSequence antikt_CS(jets, antikt_jet_def);
      current_ktmax = sorted_by_pt(antikt_CS.inclusive_jets())[0].pt();
    }
    else {
      // we set ktmax to the energy of the hardest particle, for now
      current_ktmax = sorted_by_E(jets)[0].E();
    }
  }
  else if (_clustering_type == DynamicKtMax){ // dynamic ktmax
    if (!_spherical) current_ktmax = sorted_by_pt(cs.jets())[0].pt();
    else current_ktmax = sorted_by_E(jets)[0].E();
  }
  else {
    cerr << "ClusteringType not recognized. Exiting." << endl;
    exit(-1);
  }

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j);

    // The above are correct for all pairs of non-oppositely flavoured pseudojets.
    // We recompute the distances by hand for such pairs.
    unsigned iA = 0;
    const std::vector<ClusterSequence::history_element> & h = cs.history();
    while (iA < cs.jets().size()) {
      // if jet is not flavoured, continue
      FlavInfo flavA = cs.jets()[iA].user_info<FlavHistory>().current_flavour();
      if (flavA.is_flavourless()) { ++iA; continue;}
      // if jet has recombined already, continue
      if (h[cs.jets()[iA].cluster_hist_index()].child > 0) { ++iA; continue;}

      unsigned iB = iA + 1;
      while (iB < cs.jets().size()) {
        // if the two jets don't have opposite flavour, their
        // distance was computed correctly
        FlavInfo flavB = cs.jets()[iB].user_info<FlavHistory>().current_flavour();
        if (flavB.is_flavourless() || !(flavA + flavB).is_flavourless()) { ++iB; continue;}
        // if jet has recombined already, continue
        if (h[cs.jets()[iB].cluster_hist_index()].child > 0) { ++iB; continue;}

        // then the distance needs to be recomputed. If it's smaller than
        // dij, replace
        double flav_dij = distance_opposite_flavour(cs.jets()[iA], cs.jets()[iB], current_ktmax);
        if (flav_dij < dij) {
          dij = flav_dij;
          i = iA;
          j = iB;
        }
        ++iB;
      }
      ++iA;
    }

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);

      // Put this in its own function?
      // update ktmax: either pt(k) > ktmax,
      // or one of either i or j had pt = ktmax
      if (_clustering_type == DynamicKtMax) {
        if (!_spherical) {    // non-spherical
          if (cs.jets()[k].pt() > current_ktmax) // pt(k) > current_ktmax -> update
            current_ktmax = cs.jets()[k].pt();
          // LS: could also just keep an index...
          else if (fabs(cs.jets()[i].pt() - current_ktmax) < 1e-4
                || fabs(cs.jets()[j].pt() - current_ktmax) < 1e-4) {
            // Then we loop through current jets to find the one with highest ktmax
            current_ktmax = 0.0;
            const std::vector<ClusterSequence::history_element> & h = cs.history();
            for (unsigned iA = 0; iA < cs.jets().size(); iA++) {
              if (h[cs.jets()[iA].cluster_hist_index()].child > 0) continue;
              if (cs.jets()[iA].pt() > current_ktmax) {
                current_ktmax = cs.jets()[iA].pt();
              }
            }
          }

        }
        else {    // spherical
          if (cs.jets()[k].E() > current_ktmax) // pt(k) > current_ktmax -> update
            current_ktmax = cs.jets()[k].E();
          // LS: could also just keep an index...
          else if (fabs(cs.jets()[i].E() - current_ktmax) < 1e-4
                || fabs(cs.jets()[j].E() - current_ktmax) < 1e-4) {
            // Then we loop through current jets to find the one with highest ktmax
            current_ktmax = 0.0;
            const std::vector<ClusterSequence::history_element> & h = cs.history();
            for (unsigned iA = 0; iA < cs.jets().size(); iA++) {
              if (h[cs.jets()[iA].cluster_hist_index()].child > 0) continue;
              if (cs.jets()[iA].E() > current_ktmax) {
                current_ktmax = cs.jets()[iA].E();
              }
            }
          }
        }
      }
    } else {
      double diB = 1 / pow(cs.jets()[i].pt(), 2);
      cs.plugin_record_iB_recombination(i, diB);
      nnh.remove_jet(i);
    }
    njets--;
  }
}

FASTJET_END_NAMESPACE
