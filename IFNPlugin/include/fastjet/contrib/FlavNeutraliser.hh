#ifndef __FJCONTRIB_FLAVNEUTRALISER_HH__
#define __FJCONTRIB_FLAVNEUTRALISER_HH__

#include "fastjet/contrib/FlavInfo.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

// ----------------------------------------------------------------------
/// Compare the net flavour of two vector<fastjet::PseudoJet> objects.
bool jet_net_flavour_compare(std::vector<fastjet::PseudoJet> & j, std::vector<fastjet::PseudoJet> & k);

/// Compare the flavour of each jet of two vector<fastjet::PseudoJet> objects.
/// If either vector has >3 jets, only the 3 hardest (i.e. largest energy) jets
/// are compared.
bool jet_flavour_compare(const std::vector<fastjet::PseudoJet> & j, const std::vector<fastjet::PseudoJet> & k, 
                        const int max_jets=-1, bool sort_by_px = false);


/// ----------------------------------------------------------------
/// Class for the new flavour neutraliser
class FlavNeutraliser {
public:
  enum measure { sinh_delta_R, delta_R, jade_delta_R, maxscale_delta_R, phi2_coshy, cosphi_coshy,
                 aktlike_pair_refratio, aktlike_pair_dynrefratio, 
                 jade, jadea2, maxscale, general };

  /// Main constructor for the FlavNeutraliser class.
  FlavNeutraliser(bool modulo_2,
                  measure measure_in,
                  bool use_mass_flav = false,
                  bool spherical_algo = false,
                  double pp = 1.0,
                  bool recursive_neutralisation = true) :
                  _modulo_2(modulo_2), _measure(measure_in), _use_mass_flav(use_mass_flav),
                  _spherical_algo(spherical_algo), _pp(pp), _recursive(recursive_neutralisation)
                  {}

  /// A constructor for the FlavNeutraliser class for 
  /// the general measure of the form:
  ///
  ///   u_ij = max( pti^2 , ptj^2 )^p min( pti^2 , ptj^2 )^q 
  ///          x 2[ 1/a^2 (cosh(a*Δy_ij) - 1) + (cos(Δφ_ij) - 1) ]
  FlavNeutraliser(double p, double q, double a,
                  bool modulo_2,
                  measure measure_in = general,
                  bool use_mass_flav = false,
                  bool spherical_algo = false,
                  double pp = 1.0,
                  bool recursive_neutralisation = true) :
                  _p(p), _q(q), _a(a), _modulo_2(modulo_2), _measure(measure_in), 
                  _use_mass_flav(use_mass_flav), _spherical_algo(spherical_algo), 
                  _pp(pp), _recursive(recursive_neutralisation)
                  {}

  void set_recursive(bool val) {_recursive = val;}
  bool recursive() const {return _recursive;}

  /// @brief Returns the neutralisation distance (u12) between the
  ///        two input jets
  /// @param p1 input pseudojet 1
  /// @param p2 input pseudojet 1
  /// @param ref_scale reference scale for distance measures that need such a scale
  /// @return u12
  double neutralisation_distance(const PseudoJet &p1, 
                                 const PseudoJet &p2, double ref_scale) const;

  /// Uses the ClusterSequence to build up neutralisation and returns a
  /// flavour-neutralised vector<PseudoJet> that is intended to replace
  /// cs.jets() (which contains the whole clustering history).
  std::vector<PseudoJet> neutralise(ClusterSequence &cs) const;

  /// returns the inclusive jets as obtained after flavour
  /// neutralisation (i.e. with the neutralise function).
  std::vector<PseudoJet> inclusive_jets(ClusterSequence &cs,
                                   const double ptmin = 0.0) const {
    const std::vector<ClusterSequence::history_element> & hist = cs.history();
    std::vector<PseudoJet> jets = neutralise(cs);
    std::vector<PseudoJet> jets_local;
    int u = int(hist.size()) - 1;
    while (u >= 0) {
      if (hist[u].parent2 == ClusterSequence::BeamJet) {
        int parent1 = hist[u].parent1;
        const PseudoJet & jet = jets[hist[parent1].jetp_index];
        if (jet.pt() >= ptmin) {jets_local.push_back(jet);}
      }
      u--;
    }
    return jets_local;
  }

  void use_neutralisation_candidates(
    PseudoJet & jet_i,
    double uij,
    int ih_step,
    std::vector<std::pair<PseudoJet*, double>> & flavour_candidates,
    double ref_scale) const;

  void use_neutralisation_candidates_recursive(
    PseudoJet & jet_i,
    double uij,
    int ih_step,
    std::vector<std::pair<PseudoJet*, double>> & flavour_candidates,
    double ref_scale,
    const PseudoJet * exclude = nullptr
    ) const;

  /// @brief  return a value for the jet's rapidity
  ///         that is FJ's default except at small rapidity
  ///         where a more accurate formula is used, allowing
  ///         representations of rapidities down to the 
  ///         numeric_limits<double>::min()
  ///
  /// @param p the input jet
  /// @return the input jet's rapidity
  static double accurate_rap(const PseudoJet & p) {
    double rap = p.rap();

    // switch to accurate rap if the momentum has moderately small
    // rapidity
    constexpr double rap_transition = 0.1;
    if (std::fabs(rap) < rap_transition) rap = 0.5 * log1p(2*p.pz()/(p.E() - p.pz()));
    return rap;
  }

  /// @brief  return a value for the absdphi between p1 and p2
  /// @param p1 input jet 1
  /// @param p2 input jet 2
  /// @return abs(dphi_{12})
  ///
  /// the result is FJ's default, except at small dphi values, where a 
  /// more accurate calculation is used based on a transverse cross
  /// product, which is more accurate in particular when both momenta
  /// are close to either the x or y axis
  static double accurate_absdphi(const PseudoJet & p1, const PseudoJet & p2) {
    double dphi = std::fabs(p1.delta_phi_to(p2));
    constexpr double phi_transition = 0.1;
    if (dphi < phi_transition) {
      // use a more accurate calculation based on a transverse cross product;
      // organise the cross product so as to avoid excessively large or small
      // numbers appearing in intermediate stages
      double invp1pt = 1.0/p1.pt();
      double invp2pt = 1.0/p2.pt();
      double cross = (p1.px()*invp1pt) * (p2.py()*invp2pt) - (p2.px()*invp2pt) * (p1.py()*invp1pt);
      //double cross = (p1.px() * p2.py() - p2.px() * p1.py())/sqrt(p1.pt2() * p2.pt2());
      assert(cross <= 1.0 && cross >= -1.0);
      dphi = std::fabs(asin(cross));
    }
    return dphi;
  }

private:
  /// the value of deltaR2 below which we replace 2*(cosh-cos) with
  /// (addition as of 2021-09-26, but not yet being used)
  double _p, _q, _a;
  static const double _deltaR2_handover;
  bool    _modulo_2;
  measure _measure;
  bool    _use_mass_flav;
  bool    _spherical_algo;
  double  _pp;
  bool    _recursive;

  bool    _writeout_uijs = false;


};

} // namespace contrib

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
#endif // __FJCONTRIB_FLAVNEUTRALISER_HH__
