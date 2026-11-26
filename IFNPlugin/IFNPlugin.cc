#include "fastjet/contrib/IFNPlugin.hh"
#include "fastjet/contrib/FlavNeutraliser.hh"

#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/ClusterSequence.hh"
#endif

// verification that we have a sufficiently recent version of FastJet
#ifdef FASTJET_VERSION_NUMBER
#if FASTJET_VERSION_NUMBER < 30401
#error "IFNPlugin requires FastJet version 3.4.1 or higher"
#endif
#endif

// for releases, this is commented out, though we might
// still need it for tests.
//#include "PseudoJetIO.hh"
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

using namespace std;

//---------------------------------------------------------------

string IFNPlugin::description () const {
  ostringstream desc;
  desc <<  "Interleaved Flavour Neutralisation (IFN) plugin based on " << _jet_def.description();
  if (_spherical_algo) {
    desc << ", using a spherical neutralisation measure of type ";
    switch (_measure_in) {
        case FlavNeutraliser::general:
          if (_p+_q == 1.0) desc << "standard uij, with alpha = " << 2*_p;
          else              desc << "general case with p = " << _p << " q = " << _q; 
          break;
      case FlavNeutraliser::jade:
          desc << "jade"; break;
      case FlavNeutraliser::maxscale:
          desc << "maxscale"; break;
      case FlavNeutraliser::aktlike_pair_refratio:
          desc << "aktlike_pair_refratio"; break;
      default:
          desc << "[deprecated, index=" << _measure_in<<"]"; 
          desc << ", with pp = " << _pp;
          break;
    }
  } else {
    desc << ", using a ";
    string nm = " neutralisation measure";
    switch (_measure_in) {
        case FlavNeutraliser::general:
          if (_p+_q == 1.0) desc << "standard uij " << nm << " with alpha = " << 2*_p << ", omega = " << _a;
          else              desc << "general case with p = " << _p << " q = " << _q << " omega = " << _a; 
          break;
        // older variants
        case FlavNeutraliser::sinh_delta_R:
          desc << "sinh_delta_R" << nm; break;
        case FlavNeutraliser::delta_R:
          desc << "delta_R" << nm; break;
        case FlavNeutraliser::jade_delta_R:
          desc << "jade_delta_R" << nm; break;
        case FlavNeutraliser::maxscale_delta_R:
          desc << "maxscale_delta_R" << nm; break;
        case FlavNeutraliser::phi2_coshy:
          desc << "phi2_coshy" << nm; break;
        case FlavNeutraliser::cosphi_coshy:
          desc << "cosphi_coshy" << nm; break;
        case FlavNeutraliser::aktlike_pair_refratio:
          desc << "aktlike_pair_refratio" << nm; break;
        case FlavNeutraliser::aktlike_pair_dynrefratio:
          desc << "aktlike_pair_dynrefratio" << nm; break;
        case FlavNeutraliser::jade:
          desc << "jade (without correction factor a)" << nm; break;
        case FlavNeutraliser::jadea2:
          desc << "jade (with a = 2)" << nm; break;
        case FlavNeutraliser::maxscale:
          desc << "maxscale"; break;
        default:
          desc << "UNRECOGNISED";
    }
  }
  desc << ", with modulo_2 = " << _modulo_2;
  desc << " and recursive = " << recursive();
  return desc.str();
}

void IFNPlugin::check_mod2_consistency() const {

  // a check that modulo 2 expectations are consistent with the flavour algorithm
  const auto * flav_recombiner = dynamic_cast<const FlavRecombiner*>(_jet_def.recombiner());
  if (flav_recombiner)  {
    if (_modulo_2) {
      if (flav_recombiner->flav_summation() != FlavRecombiner::modulo_2) throw Error(
                    "IFNPlugin modulo_2 is set to true, but base jet "
                    "definition (" + _jet_def.description() + 
                    ") has a FlavRecombiner with flav_summation != modulo_2");
    } else {
      if (flav_recombiner->flav_summation() != FlavRecombiner::net) throw Error(
                  "IFNPlugin modulo_2 is set to false, but base jet definition ()"
                   + _jet_def.description() + 
                  ") has a FlavRecombiner with flav_summation != net");
    }
  }

}

//---------------------------------------------------------------
void IFNPlugin::run_clustering(ClusterSequence & cs) const {

  // take the initial particles from the cs that gets passed to
  // (which is a cs that has not yet undergone any clustering)
  // and cluster them with the jet definition that was used
  // to construct the IFNPlugin.

  // need to make sure all jets have FlavHistory and that their indices are
  // correct
  for (unsigned i = 0; i < cs.jets().size(); i++) {

    const PseudoJet & jet = cs.jets()[i];
    int hist_index = jet.cluster_hist_index();
    if (jet.has_user_info<FlavInfo>()) {
      /// it can be useful to be able to start from a FlavInfo
      cs.plugin_non_const_jet(i).set_user_info(new FlavHistory(jet.user_info<FlavInfo>(), hist_index));
      if (_modulo_2) {
        dynamic_cast<FlavHistory *>(
            cs.plugin_non_const_jet(i).user_info_shared_ptr().get())
            ->apply_modulo_2();
      }
    } else if (jet.has_user_info<FlavHistory>()) {

      // if we start from a FlavHistory, make sure that we copy the
      // object, copy its current_flavour and assign this CS's
      // hist_index. 
      //
      // Taking a copy of the FlavHistory object ensures that if the
      // underlying PseudoJet is going to be used for other clusterings,
      // then we don't end up sharing a single FlavHistory object among
      // them (which would inevitably corrupt the FlavHistory of at
      // least one of the CS's).
      //
      // NB: if the user info is multiply derived (e.g. more than just
      // FlavHistory) then the other info will be lost; we still need to
      // think about how best to handle that scenario (possibly
      // introduce UserInfoCopyableBase in FJ with a generic copy() member?)
      cs.plugin_non_const_jet(i).set_user_info(new FlavHistory(
          jet.user_info<FlavHistory>().current_flavour(), hist_index));

      if (_modulo_2) {
        dynamic_cast<FlavHistory *>(
            cs.plugin_non_const_jet(i).user_info_shared_ptr().get())
            ->apply_modulo_2();
      }
    } else {
      throw fastjet::Error(
          "A PseudoJet being clustered with IFNPlugin had neither "
          "FlavInfo nor FlavHistory user_info.");
    }
  }

  ClusterSequence local_cs(cs.jets(), _jet_def);

  typedef ClusterSequence CS;

  /// get the neutralised jets
  FlavNeutraliser flav_neutraliser(_p, _q, _a, _modulo_2, _measure_in, _use_mass_flav, _spherical_algo, _pp);
  //FlavNeutraliser flav_neutraliser(_modulo_2, _measure_in, _use_mass_flav, _spherical_algo, _pp);
  flav_neutraliser.set_recursive(recursive());
  vector<PseudoJet> jets = flav_neutraliser.neutralise(local_cs);
  //const auto & jets = local_cs.jets();

  // transfer the neutralised version of the input particles to our
  // native CS.
  for (unsigned i = 0; i < cs.jets().size(); i++){
    cs.plugin_non_const_jet(i).set_user_info_shared_ptr(jets[i].user_info_shared_ptr());
  }

  const vector<ClusterSequence::history_element> & hist = local_cs.history();

  // Loop over each step of the clusterings to register the neutralised
  // PseudoJets and their clustering in our history.
  for (unsigned ih_step = 0; ih_step < hist.size(); ++ih_step) {
    const auto & h = hist[ih_step];
    if (h.parent1 == CS::InexistentParent && h.parent2 == CS::InexistentParent) {
      // this signifies an initial particle, so there is nothing to do.
      continue;
    } else if (h.parent1 >= 0 && h.parent2 == CS::BeamJet) {
      cs.plugin_record_iB_recombination(hist[h.parent1].jetp_index, h.dij);
    } else if (h.parent1 >= 0 && h.parent2 >= 0) {
      int newjet_k;
      cs.plugin_record_ij_recombination(
          hist[h.parent1].jetp_index,
          hist[h.parent2].jetp_index,
          h.dij,
          jets[h.jetp_index],
          newjet_k
          );
    } else {
      throw Error("Invalid h.parent1 and h.parent2 combination");
    }
  }
}

} // namespace contrib

FASTJET_END_NAMESPACE
