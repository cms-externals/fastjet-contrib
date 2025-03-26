#ifndef __FJCONTRIB_IFNPLUGIN_HH__
#define __FJCONTRIB_IFNPLUGIN_HH__

// #ifdef __FJC_FLAVNEUT_USEFJCORE__
// #include "fjcore.hh"
// #define FASTJET_BEGIN_NAMESPACE "namespace fastjet {"
// #define FASTJET_END_NAMESPACE   "}"
// #else 
// #include "fastjet/JetDefinition.hh"
// #define precision_type double
// #endif 

#include "fastjet/contrib/FlavNeutraliser.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

/// Plugin that runs a specified jet algorithm in conjunction
/// with Flavour Neutralisation
class IFNPlugin : public JetDefinition::Plugin {
public:

 /// Main constructor for a IFNPlugin
 ///
 /// \param jet_def: the jet definition on which this plugin will be based
 /// \param alpha: the parameter alpha in the uij neutralisation measure
 /// \param omega: the parameter omega in the uij neutralisation measure (<=0 defaults to
 ///               3 - alpha as used in the paper); relevant only for pp algorithms
 /// \param flav_summation: the flavour summation method to be used (should be one of 
 ///        FlavRecombiner::net and FlavRecombiner::modulo_2)
 /// \param use_mass_flav: intended for IRC tests when dealing with large rapidities
 IFNPlugin(
    const JetDefinition & jet_def,
    double alpha,                  
    double omega = -1,
    FlavRecombiner::FlavSummation flav_summation = FlavRecombiner::net, 
    bool use_mass_flav = false
 ) :  _jet_def(jet_def), 
      _p(0.5*alpha), 
      _q(0.5*(2-alpha)), 
      _a(omega > 0 ? omega : 3 - alpha),
      _modulo_2(flav_summation == FlavRecombiner::modulo_2),
      _measure_in(FlavNeutraliser::general), 
      _use_mass_flav(use_mass_flav),
      _spherical_algo(_jet_def.is_spherical()),
       _pp(1), // dummy value
      _recursive(true) {
        if (flav_summation == FlavRecombiner::any_abs) {
          throw Error("IFNPlugin: FlavRecombiner::any_abs is not supported");
          }
          check_mod2_consistency();
      }

 /// Old main constructor for the class. Intended only for developers.
 ///
 /// \param jet_def: the jet definition on which this plugin will be based
 /// \param modulo_2: if true, flavour should be treated modulo 2
 /// \param measure_in: the distance measure to be used in the neutralisation
 /// \param use_mass_flav: only intended for tests
 /// \param pp: in [old, pre-release] e+e- alg we use (Emax^2/Emin^2)^pp
 /// \param recursive_neutralisation: if true (default), use recursive algo
 IFNPlugin(
     JetDefinition jet_def,
     bool modulo_2,
     FlavNeutraliser::measure measure_in,
     bool use_mass_flav = false,
     double pp = 1,
     bool recursive_neutralisation = true)
     : _jet_def(jet_def),
       _modulo_2(modulo_2),
       _measure_in(measure_in),
       _use_mass_flav(use_mass_flav),
       _spherical_algo(_jet_def.is_spherical()),
       _pp(pp),
       _recursive(recursive_neutralisation) {check_mod2_consistency();}
         


 /// Alternative constructor for the class that allows
 /// defining the neutraliser measure with general p, q and a
 /// according to
 ///   u_ij = max( pti^2 , ptj^2 )^p min( pti^2 , ptj^2 )^q 
 ///          x 2[ 1/a^2 (cosh(a*Δy_ij) - 1) + (cos(Δφ_ij) - 1) ]
 /// \param jet_def: the jet definition on which this plugin will be based
 /// \param p: p parameter in the definition above
 /// \param q: q parameter in the definition above
 /// \param a: a parameter in the definition above
 /// \param modulo_2: if true, flavour should be treated modulo 2
 /// \param measure_in: the distance measure to be used in the neutralisation
 /// \param use_mass_flav: only intended for tests
 /// \param pp: in e+e- alg we use (Emax^2/Emin^2)^pp (default = 1)
 /// \param recursive_neutralisation: if true (default), use recursive algo
 IFNPlugin(
     JetDefinition jet_def,
     double p, double q, double a,
     bool modulo_2 = false,
     FlavNeutraliser::measure measure_in = FlavNeutraliser::general,
     bool use_mass_flav = false,
     double pp = 1,
     bool recursive_neutralisation = true)
     : _jet_def(jet_def), 
       _p(p), _q(q), _a(a),
       _modulo_2(modulo_2),
       _measure_in(measure_in),
       _use_mass_flav(use_mass_flav),
       _spherical_algo(_jet_def.is_spherical()),
       _pp(pp),
       _recursive(recursive_neutralisation)
        {check_mod2_consistency();}

 typedef FlavNeutraliser::measure Measure;

 void set_recursive(bool val) {_recursive = val;}
 bool recursive() const {return _recursive;}

 void set_measure(FlavNeutraliser::measure measure) {_measure_in = measure;}
 FlavNeutraliser::measure measure() const {return _measure_in;}

 // Jet radius
 double R() const FASTJET_OVERRIDE { return _jet_def.R(); }

 // Required by base class:
 std::string description() const FASTJET_OVERRIDE;

 void run_clustering(ClusterSequence &) const FASTJET_OVERRIDE;

 bool is_spherical() const FASTJET_OVERRIDE {return _spherical_algo;}

private:

  /// checks that the mod2 setting is consistent with the base
  /// algorithm's flavour recombiner, if it has one
  void check_mod2_consistency() const;

  /// the base jet definition
  JetDefinition _jet_def;

  /// the general form of the uij distance measure is max(pti^2, ptj^2)^p min(pti^2, ptj^2)^q Omega_ij
  double _p;
  double _q;

  /// the internal variable that stores what is known as omega in the paper
  double _a;

  // whether modulo 2 should be used
  bool _modulo_2;
  FlavNeutraliser::measure _measure_in;
  bool _use_mass_flav;
  bool _spherical_algo;
  double _pp;
  bool _recursive;
};

// a typedef for backward compatibility
typedef IFNPlugin FlavNeutraliserPlugin;

} // namespace contrib


FASTJET_END_NAMESPACE

#endif // __FJCONTRIB_IFNPLUGIN_HH__
