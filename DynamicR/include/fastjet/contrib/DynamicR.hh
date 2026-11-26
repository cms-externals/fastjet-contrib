//--------------------------------------------------------------------
//
// file: DynamicR.hh
//
//----------------------------------------------------------------------
// $Id: DynamicR.hh 1567 2026-06-24 12:19:44Z tousik $
//
// Copyright (c) 2023-, Tousik Samui
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

#ifndef __FASTJET_CONTRIB_DYNAMICR_HH__
#define __FASTJET_CONTRIB_DYNAMICR_HH__

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"    // JetDefinition class

FASTJET_BEGIN_NAMESPACE      // FastJet namespace (macro defined in fastjet/internal/base.hh)

using namespace fastjet;

class PseudoJet;             // Forward declaration

namespace contrib{

//----------------------------------------------------------------------
/// \class DynamicR
/// Dynamic Radius jet algorithm based on 2301.13074
///
/// The Dynamic Radius (DR) jet algorithm starts from an initial radius
/// parameter, R0. It then adapts each evolving proto-jet's effective
/// radius based on the pt-weighted standard deviation of
/// inter-constituent distances in the rapidity-azimuth plane.
///
/// The main class is
///   DynamicR(R0, algorithm, recomb_scheme)
/// where
///     R0            - initial jet radius parameter, same for all jets
///                     in an event
///     algorithm     - one of DRAK_algorithm, DRCA_algorithm, or
///                     DRKT_algorithm corresponding to Dynamic-R
///                     anti-kt-like, C/A-like, or kt-like variatants,
///                     respectively.
///     recomb_scheme - recombination scheme
///                     (presently only E_scheme is implemented)
///
/// Supplementary classes:
///     DynamicRGlobalInfo - stores input R0 (globally) for further use
/// 
///     DynamicRJetInfo    - stores, for each jet, pt-weighted mean
///                          interconstituent distance, its mean-squared,
///                          and the associated total weight
///
/// The classes DRAK, DRCA, and DRKT are convenient wrapper classes for
/// Dynamic-R anti-kt, C/A, and kt type algorithms, respectively.
///
//----------------------------------------------------------------------

//----------------------------------------------------------------------

// Small helper class to store initial radius parameter R0
class DynamicRGlobalInfo {
  public:
  DynamicRGlobalInfo (double radin) {  // Constructor taking an input radius
    _r0   = radin;            // Store the input radius in the private member
  }
  double R0() const {return _r0;}  // Accessor returning the stored radius

  private:
  double _r0;                 // Stored radius value
};


// UserInfo class that can be attached to a PseudoJet to store additional radius statistics and total weight
class DynamicRJetInfo: public PseudoJet::UserInfoBase {
  public:
  DynamicRJetInfo(double mean_radin, double ms_radin, double wtin) {  // Constructor with mean radius, mean-squared radius, and total weight
  _mean_rad = mean_radin;    // Store mean radius
  _ms_rad = ms_radin;      // Store mean-squared radius
  _wt = wtin;                // Store weight
  }

  double mean_R() const {return _mean_rad;}  // Accessor for mean radius
  double ms_R() const {return _ms_rad;}    // Accessor for mean-squared radius
  double wt() const {return _wt;}            // Accessor for weight

  inline void set_Rd(const double rd_in) { _rd = rd_in;}  // set modified dynamic radius
  double Rfin() const {return _rd;}  // Final dynamic radius (same as Rd)
  double Rd() const {return _rd;}    // modified dynamic radius

  private:
  double _mean_rad, _ms_rad, _wt, _rd;  // Stored mean radius, ms radius, and weight
};


// Enumeration of the supported DynamicR jet algorithms
enum DynamicRJetAlgorithm {
  DRKT_algorithm = 0,        // Dynamic-R kt-like algorithm
  KTlike=DRKT_algorithm,     // Alias for DRKT_algorithm

  DRCA_algorithm = 1,        // Dynamic-R Cambridge/Aachen-like algorithm
  CAlike=DRCA_algorithm,     // Alias for DRCA_algorithm

  DRAK_algorithm = 2,        // Dynamic-R anti-kt-like algorithm
  AKlike=DRAK_algorithm      // Alias for DRAK_algorithm
};

//----------------------------------------------------------------------

/// @ingroup plugins
/// \class DynamicR
/// Implementation of the Dynamic-R jet algorithm (plugin for fastjet v3.1 upwards)
//
class DynamicR : public JetDefinition::Plugin {  // FastJet plugin implementing the Dynamic-R jet algorithm
public:
  /// Main constructor for the DynamicR class.  
  ///
  DynamicR (double radius, DynamicRJetAlgorithm jet_algorithm_in=DRAK_algorithm,
		  RecombinationScheme jet_recombination_scheme=E_scheme){  // Constructor with radius, algorithm choice, and recombination scheme
    _radius  = radius;                                             // Store jet radius parameter
    _radius2 = radius*radius;                                      // Store squared radius for repeated use
    _jet_algorithm = jet_algorithm_in;                             // Store selected DynamicR algorithm type
    _jet_recombiner = JetDefinition::DefaultRecombiner(jet_recombination_scheme);  // Initialize default recombiner
  }

  /// copy constructor
  DynamicR (const DynamicR & plugin) { 
    *this = plugin; 
  }

  // the things that are required by base class
  virtual std::string description () const;              // Return a textual description of the plugin/algorithm
  virtual void run_clustering(ClusterSequence &) const;  // Main clustering routine called by FastJet

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _radius;}             // Return the jet radius parameter
  virtual DynamicRJetAlgorithm  algorithm() const {return _jet_algorithm;}  // Return the selected algorithm type

private:
  double _radius, _radius2;      // Jet radius and its square

  DynamicRJetAlgorithm _jet_algorithm;  // Chosen DynamicR algorithm

  JetDefinition::DefaultRecombiner _jet_recombiner;  // Recombiner used for merging pseudojets

  static bool _first_time;       // Static flag, will be used for one-time banner printing or setup

  template<class BJ> void actual_run(ClusterSequence &) const;  // Internal actual implementation of the clustering run

  virtual double calculate_Rd(double , double , double ) const;
  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;  // Print a one-time informational banner about the plugin
};

// Convenience wrapper class for Dynamic-R anti-kt
class DRAK : public DynamicR {
  public:
    DRAK(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicR(radius, DRAK_algorithm, recombination_scheme) {}  // Construct anti-kt-like DynamicR plugin
};

// Convenience wrapper class for Dynamic-R Cambridge/Aachen
class DRCA : public DynamicR {
  public:
    DRCA(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicR(radius, DRCA_algorithm, recombination_scheme) {}  // Construct CA-like DynamicR plugin
};

// Convenience wrapper class for Dynamic-R kt
class DRKT : public DynamicR {
  public:
    DRKT(double radius, RecombinationScheme recombination_scheme=E_scheme) : DynamicR(radius, DRKT_algorithm, recombination_scheme) {}  // Construct kt-like DynamicR plugin
};
} // namespace contrib 

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_DYNAMICR_HH__
