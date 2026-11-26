///
/// \file MassFlav.hh
///
/// This file is mainly intended for IRC tests of the IFNPlugin
#ifndef __FJCONTRIB_MASSFLAV_HH__
#define __FJCONTRIB_MASSFLAV_HH__

#ifdef __FJC_FLAVINFO_USEFJCORE__
#include "FlavInfo.hh"
#else
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/FlavNeutraliser.hh"
#endif

#include <limits>
///
///

FASTJET_BEGIN_NAMESPACE  // defined in fastjet/internal/base.hh
namespace contrib{

  //class MassFlavHistory : public FlavInfo {
  class MassFlavHistory : public FlavHistory {
   public:
    // MassFlavHistory(FlavInfo flavinfo, double m_in = 0.0) : FlavInfo(flavinfo),
    // _m(m_in) {}
    MassFlavHistory(FlavHistory flavhist, double m_in = 0.0)
        : FlavHistory(flavhist), _m(m_in) {}
    void set_mass(double m_in) { _m = m_in; }
    double mass() const { return _m; }

   private:
    double _m;
  };

  class MassFlavRecombiner : public FlavRecombiner {
   public:
    MassFlavRecombiner() : FlavRecombiner() {}

    MassFlavRecombiner(const FlavRecombiner & flav_reco) : FlavRecombiner(flav_reco) {}


    /// Perform a recombination taking flavour into account
    /// as well as the known mass of the particles, so
    /// as to be able to correctly reconstruct the rapidity
    /// of the newly formed particle, even when pa and pb
    /// are at quite extreme rapidities.
    void recombine(const PseudoJet &pa, const PseudoJet &pb,
                   PseudoJet &pab) const {

      // old code for the case eta_a << 0, eta_b >> 0
      // not sure whether it's actually needed

      //PseudoJet pab_tilde;
      //FlavRecombiner::recombine(pa, pb, pab_tilde);

      //double pab_m2 = ma*ma + mb*mb + 2*dot_product(pa,pb);

      //// and set the correct rapidity
      //double rap = 0.5*log((max(0.0,pab_m2)+pab_tilde.kt2())/
      //                      pow(pab_tilde.E() + fabs(pab_tilde.pz()),2));
      //if (pab_tilde.pz() > 0) rap = -rap;

      //pab = PtYPhiM(pab_tilde.pt(), rap, pab_tilde.phi(), sqrt(pab_m2));
      //pab.set_user_info(new MassFlavHistory(pab_tilde.user_info<FlavHistory>(),
      //                  sqrt(pab_m2)));

      // special treatment if both rapidities are small, or if one
      // is small and it looks like the other will be small too
      constexpr double rap_lim = 0.1;

      // if both particles (or one particle and its result) are at
      // small rapidity, just use 4-vector sum, which preserves
      // accurate small x,y,z components more reliably
      bool pa_small_rap = std::fabs(pa.rap()) < rap_lim;
      bool pb_small_rap = std::fabs(pb.rap()) < rap_lim;
      double pab_approx_rap = (pa.pt()*pa.rap() + pb.pt()*pb.rap())/(pa.pt()+pb.pt());
      bool pab_small_rap = std::fabs(pab_approx_rap) < rap_lim;
      int n_small_rap = int(pa_small_rap ) + int(pb_small_rap) + int(pab_small_rap);
      // at least two of the above at small rapidity: revert to 4-vector sum
      if (n_small_rap >= 2) {
        // open question as to what we do about the mass here...
        // (i.e. we could try to make it more accurate, but will not for now)
        FlavRecombiner::recombine(pa, pb, pab);
        pab.set_user_info(new MassFlavHistory(pab.user_info<FlavHistory>(),
                                           pab.m()));
        return;
      }

      double ma = mass_of_particle(pa);
      double mb = mass_of_particle(pb);

      double avrap = (pa.rap() + pb.rap()) / 2;

      PseudoJet shifted_pa = PtYPhiM(pa.pt(), pa.rap() - avrap, pa.phi(), ma);
      PseudoJet shifted_pb = PtYPhiM(pb.pt(), pb.rap() - avrap, pb.phi(), mb);

      shifted_pa.set_user_info(
          new MassFlavHistory(pa.user_info<FlavHistory>(), ma));
      shifted_pb.set_user_info(
          new MassFlavHistory(pb.user_info<FlavHistory>(), mb));

      PseudoJet shifted_pab;

      // Recombine using the default recombiner
      FlavRecombiner::recombine(shifted_pa, shifted_pb, shifted_pab);

      pab = PtYPhiM(shifted_pab.pt(), shifted_pab.rap() + avrap,
                    shifted_pab.phi(), shifted_pab.m());
      pab.set_user_info(new MassFlavHistory(shifted_pab.user_info<FlavHistory>(),
                                        shifted_pab.m()));
    }

    double mass_of_particle(const PseudoJet &p) const {
      if (p.has_user_info<MassFlavHistory>())
        return p.user_info<MassFlavHistory>().mass();
      // otherwise we expect the mass to be zero to within rounding errors
      // so check this
      double m2 = p.m2();
      const double safety_factor = 10.0;
      if (std::fabs(m2) > safety_factor *
                              std::numeric_limits<double>::epsilon() *
                              pow(p.E(), 2)) {
        throw Error(
            "MassFlavRecombiner found particle without MassFlavHistory, whose "
            "mass is inconsistent with zero");
      }
      // if things are OK, then just return zero.
      return 0.0;
    }
  };

} // namespace contrib

FASTJET_END_NAMESPACE  // defined in fastjet/internal/base.hh

#endif // __FJCONTRIB_MASSFLAV_HH__
