#ifndef __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
#define __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__

//  EnergyCorrelator Package
//  Questions/Comments?  larkoski@mit.edu gavin.salam@cern.ch jthaler@jthaler.net
//
//  Copyright (c) 2013
//  Andrew Larkoski, Gavin Salam, and Jesse Thaler
//
//  $Id: EnergyCorrelator.hh 759 2014-11-20 23:00:34Z larkoski $
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

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// \mainpage EnergyCorrelator contrib
/// 
/// The EnergyCorrelator contrib provides an implementation of energy
/// correlators and their ratios as described in arXiv:1305.0007 by
/// Larkoski, Salam and Thaler.  Additionally, the ratio observable
/// D2 described in arXiv:1409.6298 by Larkoski, Moult and Neill
/// is also included in this contrib.
///
/// <p>There are four main classes:
///
/// - EnergyCorrelator
/// - EnergyCorrelatorRatio
/// - EnergyCorrelatorDoubleRatio
/// - EnergyCorrelatorD2
///
/// each of which is a FastJet
/// FunctionOfPseudoJet. EnergyCorrelatorDoubleRatio in particular is
/// useful for quark/gluon discrimination and boosted object tagging.
/// EnergyCorrelationD2 has been shown to be the optimal discrimination
/// observable for boosted 2-prong jets.
///
/// See the file example.cc for an illustration of usage.

//------------------------------------------------------------------------
/// \class EnergyCorrelator
/// ECF(N,beta) is the N-point energy correlation function, with an angular exponent beta.
/// 
/// It is defined as follows 
///
///  - ECF(1,\f$ \beta)  = \sum_i E_i \f$
///  - ECF(2,\f$ \beta)  = \sum_{i<j} E_i E_j \theta_{ij}^\beta \f$
///  - ECF(3,\f$ \beta)  = \sum_{i<j<k} E_i E_j E_k (\theta_{ij} \theta_{ik} \theta_{jk})^\beta \f$
///  - ECF(4,\f$ \beta)  = \sum_{i<j<k<l} E_i E_j E_k E_l (\theta_{ij}  \theta_{ik} \theta_{il} \theta_{jk} \theta_{jl} \theta_{kl})^\beta \f$
///  - ...
///
/// The correlation can be determined with energies and angles (as
/// given above) or with transverse momenta and boost invariant angles
/// (the code's default). The choice is controlled by
/// EnergyCorrelator::Measure provided in the constructor.
///
/// The current implementation handles values of N up to and including 5.
/// Run times scale as n^N/N!, where n is the number of particles in a jet.
class EnergyCorrelator : public FunctionOfPseudoJet<double> {

public:

  enum Measure {
    pt_R,     ///< use transverse momenta and boost-invariant angles, 
              ///< eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} p_{ti} p_{tj} \Delta R_{ij}^{\beta} \f$
    E_theta   ///  use energies and angles, 
              ///  eg \f$\mathrm{ECF}(2,\beta) = \sum_{i<j} E_{i} E_{j}   \theta_{ij}^{\beta} \f$
  };

  enum Strategy {
    slow,          ///< interparticle angles are not cached. 
                   ///< For N>=3 this leads to many expensive recomputations, 
                   ///< but has only O(n) memory usage for n particles
    
    storage_array  /// the interparticle angles are cached. This gives a significant speed
                   /// improvement for N>=3, but has a memory requirement of (4n^2) bytes.
  };

public:
  
  /// constructs an N-point correlator with angular exponent beta,
  /// using the specified choice of energy and angular measure as well
  /// one of two possible underlying computational Strategy
  EnergyCorrelator(int N, 
                   double beta, 
                   Measure measure = pt_R, 
                   Strategy strategy = storage_array) : 
    _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

  /// destructor
  virtual ~EnergyCorrelator(){}
  
  /// returns the value of the energy correlator for a jet's
  /// constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;

  /// returns the the part of the description related to the parameters
  std::string description_parameters() const;
  std::string description_no_N() const;

private:

   int _N;
   double _beta;
   Measure _measure;
   Strategy _strategy;

   double energy(const PseudoJet& jet) const;
   double angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const;

};

// core EnergyCorrelator::result code in .cc file.



//------------------------------------------------------------------------
/// \class EnergyCorrelatorRatio
/// A class to calculate the ratio of (N+1)-point to N-point energy correlators, 
///     ECF(N+1,beta)/ECF(N,beta), 
/// called \f$ r_N^{(\beta)} \f$ in the publication. 
class EnergyCorrelatorRatio : public FunctionOfPseudoJet<double> {

public:

  /// constructs an (N+1)-point to N-point correlator ratio with
  /// angular exponent beta, using the specified choice of energy and
  /// angular measure as well one of two possible underlying
  /// computational strategies
  EnergyCorrelatorRatio(int N, 
                        double  beta, 
                        EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R, 
                        EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array) 
    : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

  virtual ~EnergyCorrelatorRatio() {}
   
  /// returns the value of the energy correlator ratio for a jet's
  /// constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;
  
private:

   int _N;
   double _beta;

   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;


};


inline double EnergyCorrelatorRatio::result(const PseudoJet& jet) const {

   double numerator = EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
   double denominator = EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet);

   return numerator/denominator;

}


//------------------------------------------------------------------------
/// \class EnergyCorrelatorDoubleRatio
/// Calculates the double ratio of energy correlators, ECF(N-1,beta)*ECF(N+1)/ECF(N,beta)^2.
///
/// A class to calculate a double ratio of energy correlators, 
///     ECF(N-1,beta)*ECF(N+1)/ECF(N,beta)^2, 
/// called \f$C_N^{(\beta)}\f$ in the publication, and equal to 
/// \f$ r_N^{(\beta)}/r_{N-1}^{(\beta)} \f$.  
///
/// Of the different energy correlator classes, this is the one
/// recommended for quark/gluon discrimination (N=1) and for boosted
/// N-prong object discrimination (N=2 for boosted W/Z/H, N=3 for
/// boosted top).
class EnergyCorrelatorDoubleRatio : public FunctionOfPseudoJet<double> {

public:

  EnergyCorrelatorDoubleRatio(int N, 
                              double beta, 
                              EnergyCorrelator::Measure measure = EnergyCorrelator::pt_R,  
                              EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array) 
    : _N(N), _beta(beta), _measure(measure), _strategy(strategy) {};

  virtual ~EnergyCorrelatorDoubleRatio() {}
   
   
  /// returns the value of the energy correlator double-ratio for a
  /// jet's constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;

private:

   int _N;
   double _beta;
   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;


};


inline double EnergyCorrelatorDoubleRatio::result(const PseudoJet& jet) const {

   double numerator = EnergyCorrelator(_N - 1, _beta, _measure, _strategy).result(jet) * EnergyCorrelator(_N + 1, _beta, _measure, _strategy).result(jet);
   double denominator = pow(EnergyCorrelator(_N, _beta, _measure, _strategy).result(jet), 2.0);

   return numerator/denominator;

}


//------------------------------------------------------------------------
/// \class EnergyCorrelatorD2
/// A class to calculate the observable formed from the ratio of the 
/// 3-point and 2-point energy correlators,
///     ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3,
/// called \f$ D_2^{(\beta)} \f$ in the publication.
class EnergyCorrelatorD2 : public FunctionOfPseudoJet<double> {

public:

  /// constructs an 3-point to 2-point correlator ratio with
  /// angular exponent beta, using the specified choice of energy and
  /// angular measure as well one of two possible underlying
  /// computational strategies
  EnergyCorrelatorD2(double  beta,
                        EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                        EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
    : _beta(beta), _measure(measure), _strategy(strategy) {};

  virtual ~EnergyCorrelatorD2() {}

  /// returns the value of the energy correlator ratio for a jet's
  /// constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;

private:

   double _beta;

   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;


};


inline double EnergyCorrelatorD2::result(const PseudoJet& jet) const {

   double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
   double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
   double denominator2 = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

   return numerator3*numerator1*numerator1*numerator1/denominator2/denominator2/denominator2;

}



//------------------------------------------------------------------------
/// \class EnergyCorrelatorC1
/// A class to calculate the normalized 2-point energy correlators,
///     ECF(2,beta)/ECF(1,beta)^2,
/// called \f$ C_1^{(\beta)} \f$ in the publication.
class EnergyCorrelatorC1 : public FunctionOfPseudoJet<double> {

public:

  /// constructs a 2-point correlator ratio with
  /// angular exponent beta, using the specified choice of energy and
  /// angular measure as well one of two possible underlying
  /// computational strategies
  EnergyCorrelatorC1(double  beta,
                        EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                        EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
    : _beta(beta), _measure(measure), _strategy(strategy) {};

  virtual ~EnergyCorrelatorC1() {}

  /// returns the value of the energy correlator ratio for a jet's
  /// constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;

private:

   double _beta;

   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;


};


inline double EnergyCorrelatorC1::result(const PseudoJet& jet) const {

   double numerator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);
   double denominator = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);

   return numerator/denominator/denominator;

}


//------------------------------------------------------------------------
/// \class EnergyCorrelatorC2
/// A class to calculate the double ratio of 3-point to 2-point 
/// energy correlators, 
///     ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2,
/// called \f$ C_2^{(\beta)} \f$ in the publication.
class EnergyCorrelatorC2 : public FunctionOfPseudoJet<double> {

public:

  /// constructs a 3-point to 2-point correlator double ratio with
  /// angular exponent beta, using the specified choice of energy and
  /// angular measure as well one of two possible underlying
  /// computational strategies
  EnergyCorrelatorC2(double  beta,
                        EnergyCorrelator::Measure  measure  = EnergyCorrelator::pt_R,
                        EnergyCorrelator::Strategy strategy = EnergyCorrelator::storage_array)
    : _beta(beta), _measure(measure), _strategy(strategy) {};

  virtual ~EnergyCorrelatorC2() {}

  /// returns the value of the energy correlator ratio for a jet's
  /// constituents. (Normally accessed by the parent class's
  /// operator()).
  double result(const PseudoJet& jet) const;

  std::string description() const;

private:

   double _beta;

   EnergyCorrelator::Measure _measure;
   EnergyCorrelator::Strategy _strategy;


};


inline double EnergyCorrelatorC2::result(const PseudoJet& jet) const {

   double numerator3 = EnergyCorrelator(3, _beta, _measure, _strategy).result(jet);
   double numerator1 = EnergyCorrelator(1, _beta, _measure, _strategy).result(jet);
   double denominator = EnergyCorrelator(2, _beta, _measure, _strategy).result(jet);

   return numerator3*numerator1/denominator/denominator;

}


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ENERGYCORRELATOR_HH__
