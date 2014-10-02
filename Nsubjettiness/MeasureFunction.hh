//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: MeasureFunction.hh 678 2014-06-12 20:43:03Z jthaler $
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

#ifndef __FASTJET_CONTRIB_MEASUREFUNCTION_HH__
#define __FASTJET_CONTRIB_MEASUREFUNCTION_HH__

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>
#include <limits>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// TauComponents
// (eventually we might want to put this in a separate header file)
//
///////
   
/// \class TauComponents
// This class creates a wrapper for the various tau/subtau values calculated in Njettiness. This class allows Njettiness access to these variables
// without ever having to do the calculation itself. It takes in subtau numerators and tau denominator from MeasureFunction
// and outputs tau numerator, and normalized tau and subtau.
class TauComponents {
   
public:
   // empty constructor necessary to initialize tau_components in Njettiness
   // later set correctly in Njettiness::getTau function
   TauComponents() {
      _jet_pieces_numerator.resize(1, 0.0);
      _beam_piece_numerator = 0.0;
      _denominator = 0;
      _numerator = 0;
      _jet_pieces.resize(1, 0.0);
      _beam_piece = 0.0;
      _tau = 0;
      _has_denominator = false;
      _has_beam = false;
   }
   
   // This constructor takes input vector and double and calculates all necessary tau components
   TauComponents(std::vector<double> jet_pieces_numerator, double beam_piece_numerator, double denominator, bool has_denominator, bool has_beam)
   : _jet_pieces_numerator(jet_pieces_numerator),
   _beam_piece_numerator(beam_piece_numerator),
   _denominator(denominator),
   _has_denominator(has_denominator),
   _has_beam(has_beam) {
      
      if (!_has_denominator) assert(_denominator == 1.0); //make sure no effect from _denominator if _has_denominator is false
      if (!_has_beam) assert (_beam_piece_numerator == 0.0); //make sure no effect from _beam_piece_numerator if _has_beam is false
      
      _numerator = _beam_piece_numerator;
      _jet_pieces.resize(_jet_pieces_numerator.size(),0.0);
      for (unsigned j = 0; j < _jet_pieces_numerator.size(); j++) {
         _jet_pieces[j] = _jet_pieces_numerator[j]/_denominator;
         _numerator += _jet_pieces_numerator[j];
      }
      
      _beam_piece = _beam_piece_numerator/_denominator;
      _tau = _numerator/_denominator;
   }
   
   
   // return values
   std::vector<double> jet_pieces_numerator() const { return _jet_pieces_numerator; }
   double beam_piece_numerator() const { return _beam_piece_numerator; }
   double denominator() const { return _denominator; }
   double numerator() const { return _numerator; }

   bool has_denominator() const { return _has_denominator; }
   bool has_beam() const { return _has_beam; }
   
   std::vector<double> jet_pieces() const { return _jet_pieces; }
   double beam_piece() const { return _beam_piece; }
   double tau() const { return _tau; }

private:
   
   // these values are input in the constructor
   std::vector<double> _jet_pieces_numerator;
   double _beam_piece_numerator;
   double _denominator;
   bool _has_denominator; //added so that TauComponents knows if denominator is used or not
   bool _has_beam; //added so that TauComponents knows if beam regions is used or not
   
   // these values are derived from above values
   std::vector<double> _jet_pieces;
   double _beam_piece;
   double _numerator;
   double _tau;
   
};

///////
//
// Measure Function
//
///////


//------------------------------------------------------------------------
/// \class MeasureFunction
// This class calculates the tau_N of a jet given a specific measure.
// It serves as a base class for calculating tau_N according to different measures that are implemented 
// in further derived classes, but does not define a particular measure itself.
class MeasureFunction {
   
public:
   //These functions define the measure by which tau_N is calculated,
   //and are overloaded by the various measures below
   
   // Distanes to axes.  These are called many times, so need to be as fast as possible
   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const = 0;
   virtual double beam_distance_squared(const fastjet::PseudoJet& particle) const = 0;
   
   // The actual measures used in N-(sub)jettiness
   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const  = 0;
   virtual double beam_numerator(const fastjet::PseudoJet& particle) const = 0;
   
   // a possible normalization factor
   virtual double denominator(const fastjet::PseudoJet& particle) const = 0;
   
   //------
   // The functions below call the above functions and are not virtual
   //------
   
   // Return all of the necessary TauComponents for specific input particles and axes
   // Also have optional pointers to get out information about partitioning
   TauComponents result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;

   // Just getting tau value if that is all that is needed
   double tau(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const {
      return result(particles,axes).tau();
   }
   
   // Create the partitioning and stores internally
   std::vector<fastjet::PseudoJet> get_partition(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes, PseudoJet * beamPartitionStorage = NULL) const;

   // Essentially same as get_partition, but in the form needed for the jet algorithm
   std::vector<std::list<int> > get_partition_list(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) const;

   // calculates the tau result using an existing partition
   TauComponents result_from_partition(const std::vector<fastjet::PseudoJet>& jet_partitioning, const std::vector<fastjet::PseudoJet>& axes, PseudoJet * beamPartitionStorage = NULL) const;

   // shorthand for squaring
   static inline double sq(double x) {return x*x;}

   //virtual destructor
   virtual ~MeasureFunction(){}
   
protected:
   //bool set by derived classes to choose whether or not to use the denominator
   bool _has_denominator;
   bool _has_beam;
   
   // This constructor allows _has_denominator to be set by derived classes
   MeasureFunction(bool has_denominator = true, bool has_beam = true) : _has_denominator(has_denominator), _has_beam(has_beam) {}
   
};


/// \class DefaultNormalizedMeasureFunction
// This class is the default measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined as the pT of the particle multiplied by deltaR 
// to the power of beta. This class includes the normalization factor determined by R0
class DefaultNormalizedMeasureFunction : public MeasureFunction {

public:

   DefaultNormalizedMeasureFunction(double beta, double R0, double Rcutoff, bool normalized = true)
   : MeasureFunction(normalized), _beta(beta), _R0(R0), _Rcutoff(Rcutoff) {}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      return particle.squared_distance(axis);
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet& /*particle*/) const {
      return sq(_Rcutoff);
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const{
      return particle.perp() * std::pow(jet_distance_squared(particle,axis),_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      return particle.perp() * std::pow(_Rcutoff,_beta);
   }

   virtual double denominator(const fastjet::PseudoJet& particle) const {
      return particle.perp() * std::pow(_R0,_beta);
   }
   
private:
   double _beta;
   double _R0;
   double _Rcutoff;

   
};

//------------------------------------------------------------------------
/// \class DefaultUnnormalizedMeasureFunction
// This class is the unnormalized default measure, inheriting from the class above. The only difference from above
// is that the denominator is defined to be 1.0 by setting _has_denominator to false.
class DefaultUnnormalizedMeasureFunction : public DefaultNormalizedMeasureFunction {

public:
   // Since all methods are identical, UnnormalizedMeasure inherits directly
   // from NormalizedMeasure. R0 is a dummy value since the value of R0 is unecessary for this class,
   // and the "false" flag sets _has_denominator in MeasureFunction to false so no denominator is used.
   DefaultUnnormalizedMeasureFunction(double beta, double Rcutoff)
   : DefaultNormalizedMeasureFunction(beta, std::numeric_limits<double>::quiet_NaN(), Rcutoff, false) {}
};

//------------------------------------------------------------------------
/// \class GeometricMeasureFunction
// This class is the geometic measure, inheriting from the class above. This class will calculate tau_N 
// of a jet according to this measure. This measure is defined by the Lorentz dot product between
// the particle and the axis. This class includes normalization of tau_N.
class GeometricMeasureFunction : public MeasureFunction {

public:
   // Right now, we are hard coded for beam_beta = 1.0, but that will need to change
   GeometricMeasureFunction(double jet_beta, double Rcutoff) : 
     MeasureFunction(false), // doesn't have denominator
     _jet_beta(jet_beta), _beam_beta(1.0), _Rcutoff(Rcutoff) {}

   virtual double jet_distance_squared(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double pseudoRsquared = 2.0*dot_product(lightFrom(axis),particle)/(lightAxis.pt()*particle.pt());
      return pseudoRsquared;
   }

   virtual double beam_distance_squared(const fastjet::PseudoJet&  /*particle*/) const {
      return sq(_Rcutoff);
   }

   virtual double jet_numerator(const fastjet::PseudoJet& particle, const fastjet::PseudoJet& axis) const {
      fastjet::PseudoJet lightAxis = lightFrom(axis);
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(lightAxis.pt(),_beam_beta - 1.0);
      return particle.pt() * weight * std::pow(jet_distance_squared(particle,axis),_jet_beta/2.0);
   }

   virtual double beam_numerator(const fastjet::PseudoJet& particle) const {
      double weight = (_beam_beta == 1.0) ? 1.0 : std::pow(particle.pt()/particle.e(),_beam_beta - 1.0);
      return particle.pt() * weight * std::pow(_Rcutoff,_jet_beta);
   }

   virtual double denominator(const fastjet::PseudoJet&  /*particle*/) const {
      return std::numeric_limits<double>::quiet_NaN();
   }
   
   
private:
   double _jet_beta;
   double _beam_beta;
   double _Rcutoff;
   
   // create light-like axis
   fastjet::PseudoJet lightFrom(const fastjet::PseudoJet& input) const {
      double length = sqrt(pow(input.px(),2) + pow(input.py(),2) + pow(input.pz(),2));
      return fastjet::PseudoJet(input.px()/length,input.py()/length,input.pz()/length,1.0);
   }

};
   
   
} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MEASUREFUNCTION_HH__
