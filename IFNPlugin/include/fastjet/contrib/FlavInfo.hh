// Copyright (c) 2023, Fabrizio Caola, Radoslaw Grabarczyk,
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler

#ifndef __FJCONTRIB_FLAVINFO_HH__
#define __FJCONTRIB_FLAVINFO_HH__

#ifdef __FJC_FLAVINFO_USEFJCORE__
#define fastjet fjcore
#include "fjcore_local.hh"
#define FASTJET_BEGIN_NAMESPACE namespace fjcore {
#define FASTJET_OVERRIDE override
#define FASTJET_END_NAMESPACE }
#else
#include "fastjet/JetDefinition.hh"
//#include "fastjet/LimitedWarning.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/SharedPtr.hh"
#endif


#include <iostream>

//using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{
//----------------------------------------------------------------------
/// class to allow representation of flavour, including concepts
/// such as the fact that a particle is an incoming "beam", or
/// should be a "spectator" during the clustering.
///
/// The class also provides a facility to interpret PDG codes
/// and deduce the corresponding flavour content of the particle.
///
/// For jet clustering with IFN algorithms, it is often easier
/// to use the FlavHistory class, which is given below
class FlavInfo : public PseudoJet::UserInfoBase {

public:
  /// constructs a flavour info from a pdg code. This can handle most
  /// standard model codes without difficulty. A zero code produces a flavourless object.
  ///
  /// The flags argument is optional, and can be set either to
  /// FlavInfo::beam or FlavInfo::spectator
  ///
  FlavInfo (int pdg_code = 0, int flags = 0);
  /// constructs a flavour info object from the individual flavours
  FlavInfo (int n_d, int n_u, int n_s, int n_c, int n_b, int n_t, int flags = 0);

  /// apply modulo 2 to all flavours, for hadronic info
  void apply_modulo_2();

  /// apply any abs to all flavours (so any non-zero number becomes 1)
  void apply_any_abs();

  /// returns the net number of quarks of flavour iflv
  /// (iflv runs from 1=d to 6=top).
  int operator[](int iflv) const {return _flav_content[iflv];}

  /// sets the number of objects of flavour iflv (1=d to 6=top) to be n
  void set_flav(int iflv, int n) {_flav_content[iflv] = n; update_flavourless_attribute();}

  /// allows comparison of flavours
  bool operator==(const FlavInfo &) const;
  bool operator!=(const FlavInfo &) const;

  /// resets all flavours to zero except iflv; this should be called,
  /// for example, when considering b-flavour at hadron level, so that
  /// the algorithm doesn't get confused by the many other quark
  /// flavours that are around (and also because, experimentally,
  /// those other flavours may not be well known).
  void reset_all_but_flav(int iflv);

  /// returns the pdg_code, or 0 if it's unknown (e.g. due to result
  /// of recombination)
  int pdg_code() const {return _pdg_code;}

  /// label this particle as being an incoming beam particle
  void label_as_beam() {_flav_content[0] |= beam;}
  /// returns true if this particle is a beam particle.
  bool is_beam() const {return (_flav_content[0] & beam);}

  /// label this object as being a "spectator", such as a W, which is
  /// relevant for calculating the beam distance in flavour clusterings
  /// but does itself take part in the clustering
  void label_as_spectator() {_flav_content[0] |= spectator;}
  /// returns true if this particle is a spectator.
  bool is_spectator() const {return (_flav_content[0] & spectator);}

  /// returns true if the object has no net flavour
  bool is_flavourless() const {return (_flav_content[0] & _is_flavourless);}
  bool is_flavorless() const {return is_flavourless();}

  /// returns true if the object has more than one unit of flavour
  bool is_multiflavoured() const;
  bool is_multiflavored() const {return is_multiflavoured();}

  /// returns true if the particle has associated FlavInfo and if there
  /// exists a flavour for which (*this) and the particle have
  /// positive/negative (or vice versa) net amounts of that flavour. 
  ///
  /// BEWARE: do not use this to check if the flavours cancel. 
  /// Instead add the flavours and check that the result is flavourless
  bool has_opposite_flavour(const PseudoJet & particle) const;
  bool has_opposite_flavor(const PseudoJet & particle) const {return has_opposite_flavour(particle);}

  /// allows addition of flavour: note that beam, spectator and PDG status are lost
  FlavInfo operator+(const FlavInfo &) const;
  /// allows subtraction of flavour: note that beam, spectator and PDG status are lost
  FlavInfo operator-(const FlavInfo &) const;

  /// returns a string such as "u d", "cbar", etc.; "g" means gluon or anything
  /// else with no flavour
  std::string description() const;

  /// returns the flavour of a particle if that particle has flavour, otherwise
  /// just the default flavour
  static const FlavInfo & flavour_of(const PseudoJet & particle);
  // {
  //  if (particle.has_user_info<FlavInfo>()) return particle.user_info<FlavInfo>();
  //  else                                    return _no_flav;
  //}

  /// value of flag to indicate that the particle is an incoming beam particle
  static const int beam = 2;
  /// value of flag to indicate that the particle is a "spectator",
  /// such as a W, which is relevant for calculating the beam distance
  /// in flavour clusterings but does itself take part in the
  /// clustering
  static const int spectator = 4;
  int _flav_content[7];
  void update_flavourless_attribute();
  static const int _is_flavourless = 1;
private:
  int _pdg_code;
  static const FlavInfo _no_flav;

};

//----------------------------------------------------------------------
/// Class to keep track of flavour history of a given jet
///
class FlavHistory : public PseudoJet::UserInfoBase {
public:
  /// Constructor for generic initial history step from a PDG ID code
  FlavHistory(int initial_flavour_pdg_id) {
    _flavour_history.push_back(std::make_pair(-1, FlavInfo(initial_flavour_pdg_id)));
  }

  /// Constructor for generic initial history step from a FlavInfo object
  FlavHistory(const FlavInfo & initial_flavour) {
    _flavour_history.push_back(std::make_pair(-1, initial_flavour));
  }

  /// Constructor for known initial history step from a FlavInfo object
  FlavHistory(const FlavInfo & initial_flavour, 
              const int initial_hist_step) {
    _flavour_history.push_back(std::make_pair(initial_hist_step, initial_flavour));
  }

  /// returns the current flavour of the particle
  const FlavInfo & current_flavour() const {return _flavour_history.back().second;}
  /// returns the initial flavour of the particle
  const FlavInfo & initial_flavour() const {return _flavour_history.front().second;}

  /// Returns the index of the cluster_sequence history step at which this
  /// jet acquired its current flavour
  int current_hist_index() const {return _flavour_history.back().first;}

  /// Returns the index of the cluster_sequence history step at which this
  /// jet was created with its initial flavour
  int initial_hist_index() const {return _flavour_history.front().first;}





  /// Apply flavour modulo 2 to all elements of the history
  void apply_modulo_2() {
    for (unsigned i = 0; i < _flavour_history.size(); i++) {
      _flavour_history[i].second.apply_modulo_2();
    }
  }

  /// Label the most recent history element as a beam particle
  void label_as_beam() {
    _flavour_history.back().second.label_as_beam();
  }

  /// Add element to flavour history of a jet
  void update_flavour_history(FlavInfo new_flavour, int hist_step) {
    if (new_flavour != _flavour_history.back().second) {
      _flavour_history.push_back(std::make_pair(hist_step, new_flavour));
    }
  }

  /// Change the last history step of the flavour history. This could
  /// be used to, e.g., set the initial history step when the constructor
  /// with history step -1 is used.
  void amend_last_history_index(int new_hist_step) {
    _flavour_history.back().first = new_hist_step;
  }

  /// Return the flavour history vector
  const std::vector<std::pair<int, FlavInfo>> & history() const {return _flavour_history;}

  /// Return the current (final) flavour element of the history of a PseudoJet.
  /// For particles with just a FlavInfo, it returns that. 
  static const FlavInfo & current_flavour_of(const PseudoJet & particle);

  /// Return the final history step of the history of a PseudoJet
  static const int current_index_of(const PseudoJet &particle);

  /// Return the first flavour element of the history of a PseudoJet
  /// For particles with just a FlavInfo, it returns just that. 
  static const FlavInfo & initial_flavour_of(const PseudoJet &jet);

  /// Return the first history step of the history of a PseudoJet
  static const int initial_index_of(const PseudoJet &jet) {
    if (jet.has_user_info<FlavHistory>()) {
      return jet.user_info<FlavHistory>().history()[0].first;
    } else {
      throw fastjet::Error(
          "A particle without FlavHistory was searched for FlavHistory.");
    }
  }

  /// Return the flavour at a given history step
  const FlavInfo & flavour_at_step(int step) const {
    if (_flavour_history[0].first > step) {
      throw fastjet::Error(
          "A particle without FlavHistory was searched for FlavHistory.");
    }
    int index_needed = -1;
    for (unsigned i = 1; i < _flavour_history.size(); i++) {
      if ((_flavour_history[i].first > step) &&
          (_flavour_history[i - 1].first <= step))
        index_needed = i - 1;
    }
    if (index_needed == -1) {
      return _flavour_history.back().second;
    } else {
      return _flavour_history[index_needed].second;
    }
  }


  /// Reset the _flavour_history such that it has one element whose flavour is the initial
  /// flavour of the jet, and the hist_step is initiated to its initial value.
  void reset_flavour_history() {
    std::pair<int, FlavInfo> orig_flavour = _flavour_history[0];
    std::vector<std::pair<int, FlavInfo>> _new_history;
    _new_history.push_back(orig_flavour);
    _flavour_history = _new_history;
  }
  const std::vector<std::pair<int, FlavInfo>> & history() {return _flavour_history;}

  // US spelling
  static const FlavInfo & current_flavor_of(const PseudoJet & particle) {return current_flavour_of(particle);}
  static const FlavInfo & initial_flavor_of(const PseudoJet & particle) {return initial_flavour_of(particle);}  
  const FlavInfo & flavor_at_step(int step) const {return flavour_at_step(step);}
  const FlavInfo & current_flavor() const {return current_flavour();}
  const FlavInfo & initial_flavor() const {return initial_flavour();}


 private:
  std::vector<std::pair<int, FlavInfo>> _flavour_history;
};

/// ----------------------------------------------------------------
/// Class for a flavour recombiner
class FlavRecombiner : public JetDefinition::DefaultRecombiner {
public:

  enum FlavSummation {
    /// net flavour handling, 
    /// - b        = +1
    /// - bbar     = -1
    /// - b + bbar =  0
    /// - b + b    = +2
    net, // -> net_flav

    /// flavour is handled modulo 2
    /// - b        = 1
    /// - bbar     = 1
    /// - b + bbar = 0
    /// - b + b    = 0
    modulo_2,  // -> mod2_flav

    /// a given flavour is set to 1 if any non-zero number of quarks or anti-quarks is present. 
    /// So e.g. 
    /// - b        = 1
    /// - bbar     = 1
    /// - b + bbar = 1
    /// - b + b    = 1
    any_abs // -> any_flav
  };

  /// Constructor
  FlavRecombiner(FlavSummation flav_summation_in = net) : 
         DefaultRecombiner(), _flav_summation(flav_summation_in) {
    //assert(flav_summation == net && "handling of non-net flavour summation inside FlavRecombiner is still to come");
  }

  /// return the flavour summation choice
  FlavSummation flav_summation() const {return _flav_summation;}

  void preprocess(PseudoJet & p) const FASTJET_OVERRIDE {

    FlavInfo flav;
    if (p.has_user_info<FlavInfo>()) {
      flav = p.user_info<FlavInfo>();
    } else if (p.has_user_info<FlavHistory>()) {
      /// WARNING: not sure this is right -- depends very much on context
      /// (e.g. whether reclustering another algorithm's constituents, or clustering
      /// the flavoured jets that have come out of another algorithm)
      flav = p.user_info<FlavHistory>().initial_flavour();
    } else {
      throw Error("Could not identify FlavInfo or FlavHistory");
    }

    // make the initial flavour consistent with the summation choice
    apply_summation_choice(flav);
    p.set_user_info(new FlavHistory(flav));
  }

  /// Perform a recombination taking flavour into account
  void recombine(const PseudoJet &pa,
                 const PseudoJet &pb,
                 PseudoJet &pab) const FASTJET_OVERRIDE {

    // Recombine using the default recombiner
    DefaultRecombiner::recombine(pa, pb, pab);

    // put this condition early on
    assert(!pab.has_user_info<FlavHistory>());

    // Then, check if the resulting pseudojet is actually flavourless
    // If it isn't, add flavours. If it is, set it to a gluon
    FlavInfo flav = FlavHistory::current_flavour_of(pa) + FlavHistory::current_flavour_of(pb);

    // make the flavour consistent with the summation choice
    apply_summation_choice(flav);

    /// why don't we just do the following?
    pab.set_user_info(new FlavHistory(flav,pab.cluster_hist_index()));

  }

  // make the flavour consistent with the summation choice
  void apply_summation_choice(FlavInfo & flav) const {

    if      (_flav_summation == modulo_2) flav.apply_modulo_2();
    else if (_flav_summation == any_abs) flav.apply_any_abs();
    else if (_flav_summation != net) throw Error("FlavRecombiner: unknown FlavSummation");
  }

  /// returns the description of the recombiner
  std::string description() const FASTJET_OVERRIDE { return DefaultRecombiner::description() + " and " + to_string(_flav_summation) + " flavour recombination "; }

  std::string to_string(FlavSummation flav_summation) const {
    return (flav_summation == net ? "net_flav" : 
            flav_summation == modulo_2 ? "mod2_flav" : 
            flav_summation == any_abs ? "any_flav" : "unknown");
  }

private:
  FlavSummation _flav_summation;

};


} // namespace contrib
FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
#endif // __FJCONTRIB_FLAVINFO_HH__
