// $Id: ConstituentSubtractor.hh 587 2014-04-06 13:44:24Z berta $
//
// ConstituentSubtractor package
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz
//
// Copyright (c) 2014-, Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner
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

#ifndef __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__
#define __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__



#include <fastjet/internal/base.hh>
#include <fastjet/ClusterSequenceAreaBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include "fastjet/tools/Transformer.hh" // to derive Subtractor from Transformer

#include <algorithm>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



//------------------------------------------------------------------------
/// \class ConstituentSubtractor
/// A class to perform subtraction of background, e.g. pileup, from a jet at particle-level. The output is a jet with corrected constituents.
///
/// This class corrects the input particles for background contamination with th algorithm described in:
/// Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner [arXiv:1403.3108]
/// Usage (see example.cc for more details):
/// 1. estimate the background with the help of object
///    fastjet::JetMedianBackgroundEstimator _bge_rho;
/// 2. create an object:
///    fastjet::PseudoJet _set_of_particles;
/// which contains real particles and ghosts among its constituents. There are two alternatives for this:
///    a) for the subtraction of a particular jet, define active jet area before clustering
///    b) for the subtraction of the whole event, add the ghosts manually among the jet constituents, uniformly distributed in y-phi plane. Alternatively, the  option a) can be used with large radius (>20) in the C/A clustering algorithm.
/// 3. use the ConstituentSubtractor class to subtract the background in the following way:
///   contrib::ConstituentSubtractor subtractor(&_bge_rho);
///   subtractor.use_common_bge_for_rho_and_rhom(true);  // recommended, useful only for massive input particles
///   subtractor.set_alpha(1);  // optional
///   subtractor.set_max_deltaR(2);  // optional
///   PseudoJet subtracted_jet = subtractor(_set_of_particles);
///
/// The class accounts for position-dependent (in rapidity-azimuth plane) background densities, rho and rho_m.
///
  class ConstituentSubtractor : public fastjet::Transformer{
public:
  /// default ctor
  ConstituentSubtractor() : 
    _bge_rho(0), _bge_rhom(0), _common_bge(false), _externally_supplied_rho_rhom(false), _alpha(0), _max_deltaR(-1), _use_max_deltaR(false){}

  /// Constructor that takes a pointer to a background estimator for
  /// rho and optionally a pointer to a background estimator for
  /// rho_m.  If the latter is not supplied, rho_m is assumed to
  /// always be zero (this behaviour can be changed by calling
  /// use_common_bge_for_rho_and_rhom).
  ConstituentSubtractor(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom=0, double alpha=0, double maxDeltaR=-1);

  /// Constructor that takes an externally supplied value for rho and rho_m.
  ConstituentSubtractor(double rho, double rhom=0, double alpha=0, double maxDeltaR=-1);
  
  /// default dtor
  virtual ~ConstituentSubtractor(){}

  /// a description of what this class does
  virtual std::string description() const;

  /// action of the correction on a given jet. The output is PseudoJet object with subtracted constituents
  virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const;



  /// when only one background estimator, bge_rho, is specified, calling this method with argument true, causes rho_m to be calculated from the same background estimator as rho, instead of being set to zero. Currently this only works if the estimator is a JetMedianBackgroundEstimator (or derived from it), and makes use of that class's set_jet_density_class(...) method.
  void use_common_bge_for_rho_and_rhom(bool value=true);


  /// method to change the alpha-parameter figuring in the distance measure deltaR. The larger the alpha, the more are preferred to be corrected the low pt particles. The default value is 0, i.e. by default the standard deltaR definition is used: deltaR=sqrt(deltay^2 + deltaphi^2)
  void set_alpha(double alpha);


  /// method to change the max_deltaR value. For particle-ghost pairs with deltaR>max_deltaR the iterative process stops. When max_deltaR<=0, the max_deltaR parameter is not used (no upper limit on deltaR). The default value is -1, i.e. by default there is no upper limit for possible deltaR values. 
  void set_max_deltaR(double max_deltaR);



/// This function returns true if the first argument is smaller than the second argument, otherwise returns false. The comparison is done only on the first element in the two pairs. This function is used to sort in ascending order the deltaR values for each pair particle-ghost while keeping track of particles and ghosts
  static bool _function_used_for_sorting(std::pair<double,int> i,std::pair<double, int> j);

protected:


  fastjet::BackgroundEstimatorBase *_bge_rho, *_bge_rhom;
  bool _common_bge;
  double _rho, _rhom;
  bool _externally_supplied_rho_rhom;
  double _alpha;
  double _max_deltaR;
  bool _use_max_deltaR;
};






} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__
