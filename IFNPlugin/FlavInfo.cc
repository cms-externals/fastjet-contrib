// Copyright (c) 2023, Fabrizio Caola, Radoslaw Grabarczyk,
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler

#include "fastjet/contrib/FlavInfo.hh"
#include <sstream>

#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SharedPtr.hh"
#endif

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

//----------------------------------------------------------------------
FlavInfo::FlavInfo(int n_d, int n_u, int n_s, int n_c, int n_b, int n_t, int flags) :
  _pdg_code(0) {
  _flav_content[0] = flags;
  _flav_content[1] = n_d;
  _flav_content[2] = n_u;
  _flav_content[3] = n_s;
  _flav_content[4] = n_c;
  _flav_content[5] = n_b;
  _flav_content[6] = n_t;
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
FlavInfo::FlavInfo(int pdg_code, int flags) : _pdg_code(pdg_code){
  _flav_content[0] = flags;
  for(unsigned i = 1; i <= 6; i++) _flav_content[i] = 0;

  // for particles with illicit (zero) pdg_code, no work to be done
  if (_pdg_code == 0) return;

  int netsign = (pdg_code >= 0 ? +1 : -1);
  pdg_code = abs(pdg_code);

  // extract digits of the pdg_code, since these contain information
  // on flavour of component quarks
  valarray<int> digit(4);
  int           ndigits = 0;
  for (int i = 0; i < 4; i++) {
    digit[i] = pdg_code % 10;
    if (digit[i] != 0) ndigits = i+1;
    pdg_code /= 10; // "shift" things along
  }

  // start this part with _flav_content already initialised to zero
  // in constructor
  if (ndigits == 1) { // a lone quark
    if (digit[0] > 6 || digit[0] == 0) {
      cerr << "FlavInfo failed to understand pdg_code = "<<_pdg_code<<endl; exit(-1);}
    _flav_content[digit[0]] = netsign;

  } else if (ndigits == 2) { // a lepton, photon or cluster [flav lost...]
    // do nothing...

  } else { // must be a meson, cluster or baryon
    // check sanity of codes
    for (int i=1; i < ndigits; i++) {
      if (digit[i] > 6) {cerr << "FlavInfo failed to understand pdg_code = "
			       <<_pdg_code<<endl; exit(-1);}}

    // now deal with different cases
    if (ndigits == 4) { // diquark [nm0x] or baryon [nmpx]
      for (int i=1; i < ndigits; i++) {
	if (digit[i] > 0) _flav_content[digit[i]] += netsign;}
    } else if (ndigits == 3) { // meson [nmx]
      // Beware of PDG convention that says that a K+ or B+ are a
      // particle and so have positive pdg_code (i.e. flavcodes > 1). So
      if (digit[2] == 3 || digit[2] == 5) netsign = -netsign;
      _flav_content[digit[2]] += netsign;
      _flav_content[digit[1]] -= netsign;
    } else {
      cerr << "FlavInfo failed to understand pdg_code = " <<_pdg_code<<endl; exit(-1);}
  }
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
void FlavInfo::apply_modulo_2() {
  for (unsigned iflv = 1; iflv <= 6; iflv++) {
    _flav_content[iflv] = abs(_flav_content[iflv] % 2);
  }
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
void FlavInfo::apply_any_abs() {
  for (unsigned iflv = 1; iflv <= 6; iflv++) {
    _flav_content[iflv] = _flav_content[iflv] == 0 ? 0 : 1;
  }
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
void FlavInfo::reset_all_but_flav(int iflv) {
  for (int i = 1; i <= 6; i++) {
    if (i != iflv) _flav_content[i] = 0;
  }
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
bool FlavInfo::operator==(const FlavInfo & other) const {
  for (unsigned i=0; i<=6; i++){
    if(operator[](i)!=other[i]) {return 0;};
  }
  return 1;
}

bool FlavInfo::operator!=(const FlavInfo & other) const {
  for (unsigned i=0; i<=6; i++){
    if(operator[](i)!=other[i]) {return 1;};
  }
  return 0;
}


//----------------------------------------------------------------------
FlavInfo FlavInfo::operator+(const FlavInfo & other) const {
  FlavInfo sum(operator[](1)+other[1],
               operator[](2)+other[2],
               operator[](3)+other[3],
               operator[](4)+other[4],
               operator[](5)+other[5],
               operator[](6)+other[6]
               );
  sum.update_flavourless_attribute();
  return sum;
}
//----------------------------------------------------------------------
FlavInfo FlavInfo::operator-(const FlavInfo & other) const {
  FlavInfo sum(operator[](1)-other[1],
               operator[](2)-other[2],
               operator[](3)-other[3],
               operator[](4)-other[4],
               operator[](5)-other[5],
               operator[](6)-other[6]
               );
  sum.update_flavourless_attribute();
  return sum;
}

//----------------------------------------------------------------------
void FlavInfo::update_flavourless_attribute() {
  for (unsigned i = 1; i <=6; i++) {
    if (_flav_content[i] != 0) {
      _flav_content[0] &= ~ _is_flavourless; // ~ is C++ bitwise not
      return;
    }
  }
  _flav_content[0] |= _is_flavourless;
}

//----------------------------------------------------------------------
bool FlavInfo::is_multiflavoured() const {
  int flavsum = 0;
  for (unsigned i = 1; i <=6; i++) flavsum += abs(_flav_content[i]);
  return flavsum > 1;
}

bool FlavInfo::has_opposite_flavour(const PseudoJet & particle) const {
 int t = 0;
 for(int i = 1; i<=6; i++){
  if((particle.has_user_info<FlavInfo>()) &&
    ((particle.user_info<FlavInfo>().operator[](i)*operator[](i) < -0))){
    t += 1;
  }
 }
  return t > 0;
}

//----------------------------------------------------------------------
// an object with no flavour, that we can conveniently point to
const FlavInfo FlavInfo::_no_flav;

//----------------------------------------------------------------------
string FlavInfo::description() const {
  const char * flavs = "duscbt";
  ostringstream result;

  result << "[";
  if (is_flavourless()) {
    result << "g";
  } else {
    for (int iflav = 1; iflav <= 6; iflav++) {
      int n = operator[](iflav);
      for (unsigned i = 0; i < abs(n); i++) {
        result << flavs[iflav-1];
        if (n<0) result << "bar";
        result << " ";
      }
    }
  }
  result << "]";
  if (is_beam()) result << "(beam) ";
  if (is_spectator()) result << "(spectator) ";
  return result.str();
}


/// returns the flavour of a particle if that particle has flavour, otherwise
/// just the default flavour
const FlavInfo & FlavInfo::flavour_of(const PseudoJet & particle) {
  if (particle.has_user_info<FlavInfo>()) {
    return particle.user_info<FlavInfo>();
  } else if (particle.has_user_info<FlavHistory>()) {
    throw Error("FlavInfo::flavour_of called on particle with FlavHistory. "
                "Use FlavHistory::current_flavour_of(...) or "
                "FlavHistory::initial_flavour_of(...) instead");
  } else {
    return _no_flav;
  }
}

/// Return the current (final) flavour element of the history of a PseudoJet
const FlavInfo & FlavHistory::current_flavour_of(const PseudoJet & particle) {
  if (particle.has_user_info<FlavHistory>()) {
    return particle.user_info<FlavHistory>().history().back().second;
  } else if (particle.has_user_info<FlavInfo>()) {
    return particle.user_info<FlavInfo>();
  } else {
    throw fastjet::Error("A particle without FlavHistory was searched for FlavHistory.");
  }
}

/// Return the first flavour element of the history of a PseudoJet
const FlavInfo & FlavHistory::initial_flavour_of(const PseudoJet &jet) {
  if (jet.has_user_info<FlavHistory>()) {
    return jet.user_info<FlavHistory>().history()[0].second;
  } else if (jet.has_user_info<FlavInfo>()) {
    return jet.user_info<FlavInfo>();
  } else {
    throw fastjet::Error(
        "A particle without FlavHistory was searched for FlavHistory.");
  }
}

//----------------------------------------------------------------------
} // namespace contrib
FASTJET_END_NAMESPACE
