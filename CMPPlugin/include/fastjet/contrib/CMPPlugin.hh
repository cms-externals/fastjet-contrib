/// Implementation of the flavoured anti-kt algorithm
/// by Michal Czakon, Alexander Mitov and Rene Poncelet,
/// as described in https://arxiv.org/pdf/2205.11879 (v2)
///
/// Authors: Michal Czakon, Alexander Mitov, and Rene Poncelet
///
/// based on initial version by
///
/// Authors: Fabrizio Caola, Radoslaw Grabarczyk, Maxwell Hutt,
///          Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
///

#ifndef __CMPPLUGIN_HH__
#define __CMPPLUGIN_HH__

#include "fastjet/contrib/FlavInfo.hh"
// to facilitate use with fjcore
#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/NNH.hh"
#endif

FASTJET_BEGIN_NAMESPACE  // defined in fastjet/internal/base.hh
using namespace std;
using namespace contrib;

/// Plugin for algorithm  proposed by Czakon, Mitov and Poncelet 
/// as described in arXiv:2205.11879 with the modification
class CMPPlugin : public JetDefinition::Plugin {
public:

 /// correction to original CMP
 enum CorrectionType {
  /// this uses the CMP distance measure factor as written in 2205.11879v1
  /// Eq. 2.9:
  ///
  ///     Sij = 1 - Theta(1-κ)cos(πκ/2)
  NoCorrection, 
  /// SqrtCoshyCosPhiArgument causes kappa in Eq.(2.9) to be
  /// mutiplied by
  ///   
  ///    sqrt(M), with M = 2*(cosh(drap) - cos(dphi))/deltaR2
  ///
  /// which tends to 1 in the small deltaR^2 limit, and fixes
  /// a joint infrared/collinear safety issue when deltaR is large
  /// 
  SqrtCoshyCosPhiArgument,
  /// Default: same as above, but with damping a = 2
  ///
  ///  sqrt(M), with M = 2*(1/a^2*[cosh(a*drap)-1] - [cos(dphi)-1])/deltaR2
  /// 
  SqrtCoshyCosPhiArgument_a2,
  // NOT RECOMMENDED: with M from above, uses
  //
  //     Sij = 1 - Theta(1-κ)cos(πκ/2) * M
  CoshyCosPhi, 
  // with M from above, uses
  //
  //     Sij = (1 - Theta(1-κ)cos(πκ/2)) * M
  OverAllCoshyCosPhi,
  // same as above, but with damping a = 2
  OverAllCoshyCosPhi_a2
 };

 /// definition of ktmax
 ///   0 : dynamic ktmax (pT of the hardest pseudojet lying around, default)
 ///   1 : fixed ktmax   (pT of the hardest pseudojet clustered with anti-kt)
 enum ClusteringType {DynamicKtMax, FixedKtMax};

  /// Main constructor for the class
  ///
  /// @param R: the jet radius parameter
  /// 
  /// @param a: CMP parameter a
  ///
  /// @param correction_type: the type of correction we use in the CMP distance
  ///
  /// @param clustering_type: the definition of ktmax in the CMP distance
  ///
  /// @param spherical: true if using the e+e- version of the algorithm
 CMPPlugin(double R, double a,
           CorrectionType correction_type = SqrtCoshyCosPhiArgument_a2,
           ClusteringType clustering_type = DynamicKtMax, bool spherical=false)
     : _R(R), _a(a), _correction_type(correction_type),
       _clustering_type(clustering_type), _spherical(spherical) {}

 /// copy constructor
 CMPPlugin(const CMPPlugin &plugin) { *this = plugin; }

 /// CMP parameter
 double a() const { return _a; }
 /// Jet radius of antikt distance
 double R() const { return _R; }

 /// whether to use the e+e- version of the algo
 bool is_spherical() const { return _spherical; }

 // Required by base class:
 virtual std::string description() const;

 virtual double precise_squared_distance(const PseudoJet & j1, const PseudoJet & j2) const; 

 virtual double distance_opposite_flavour(const PseudoJet & j1, const PseudoJet & j2,
                                          const double ktmax) const;

 virtual void run_clustering(ClusterSequence &) const;
 
 void cross_product(const PseudoJet & p1, const PseudoJet & p2,
                    PseudoJet & retp, bool lightlike=false) const {
   double px = p1.py() * p2.pz() - p2.py() * p1.pz();
   double py = p1.pz() * p2.px() - p2.pz() * p1.px();
   double pz = p1.px() * p2.py() - p2.px() * p1.py();
 
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
 double dot_product_3d(const PseudoJet & a, const PseudoJet & b) const {
   return a.px()*b.px() + a.py()*b.py() + a.pz()*b.pz();
 }
 /// Returns (1-cos theta) where theta is the angle between p1 and p2
 double one_minus_costheta(const PseudoJet & p1, const PseudoJet & p2) const {
 
   if (p1.m2() == 0 && p2.m2() == 0) {
     // use the 4-vector dot product.
     // For massless particles it gives us E1*E2*(1-cos theta)
     double res = dot_product(p1,p2) / (p1.E() * p2.E());
     return res;
   } else {
     double p1mod = sqrt(p1.modp2());
     double p2mod = sqrt(p2.modp2());
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

private:
  static const double _deltaR2_handover;
  double _R, _a;
  CorrectionType _correction_type;
  ClusteringType _clustering_type;
  bool _spherical;
  template <typename NN>
  void _NN_clustering(ClusterSequence &cs, NN &nn) const;
  //template <typename NN, typename CMPNNInfo>
  //void k_tmax_clustering(ClusterSequence &cs, NN &nn, CMPNNInfo info) const;
};

FASTJET_END_NAMESPACE

#endif  // __CMPPLUGIN_HH__
