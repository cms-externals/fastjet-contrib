// $Id: GHSAlgo.hh 1507 2025-01-30 10:33:36Z gstagn $
//
// Copyright (c) 2025, Rhorry Gauld, Alexander Huss, Giovanni Stagnitto
//
// based on initial version by Fabrizio Caola, Radoslaw Grabarczyk,
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
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

#ifndef __FASTJET_CONTRIB_GHSPLUGIN_HH__
#define __FASTJET_CONTRIB_GHSPLUGIN_HH__

#include <fastjet/internal/base.hh>

#include <iostream>

#ifndef __FJC_FLAVINFO_USEFJCORE__
#include "fastjet/NNH.hh"
#include "fastjet/PseudoJet.hh"
#endif

#include "fastjet/contrib/FlavInfo.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  /// given a list of base-jets (before applying a hardness cut) in
  /// the event (jets_base), return the jets with "dressed" flavour information
  ///
  /// @param jets_base: the input base-jets; the full list of event
  /// particles and associated flavours will be deduced from
  /// ClusterSequence associated with these jets.
  ///
  /// @param ptcut: this parameter applies a hardness cut on the base-jets
  /// and should be chosen to match the fiducial jet definition
  ///
  /// @param alpha: power of (ktmax/ktmin) used in the flavour-kt distance
  ///
  /// @param omega: relative weighting of rapidity separation
  ///
  /// @returns the list of hard jets dressed with their flavour
  std::vector<PseudoJet> run_GHS(
      const std::vector<PseudoJet> &jets_base, double ptcut, double alpha = 1.,
      double omega = 2.,
      const FlavRecombiner &flav_recombiner = FlavRecombiner());

  /// dress a given set of jets with the flavour of a given set of flavoured
  /// inputs (e.g. partons or heavy-flavour hadrons that were not themselves
  /// clustered into the jets), with the jet association chosen by the caller.
  /// run_GHS is this function applied to the hard jets and to the particles
  /// of their ClusterSequence, associated through PseudoJet::is_inside.
  ///
  /// @param jets: the jets to be dressed (e.g. the hard jets of the event);
  /// their user_index is preserved in the output
  ///
  /// @param flavoured_inputs: the flavour carriers, each with a FlavHistory
  /// user info and a user_index of 0 (not associated with any jet) or
  /// -(i+1) (associated with jets[i])
  ///
  /// @param alpha: power of (ktmax/ktmin) used in the flavour-kt distance
  ///
  /// @param omega: relative weighting of rapidity separation
  ///
  /// @returns the jets, in the input order, dressed with their flavour
  std::vector<PseudoJet> dress_jets_GHS(
      const std::vector<PseudoJet> &jets,
      const std::vector<PseudoJet> &flavoured_inputs, double alpha = 1.,
      double omega = 2.,
      const FlavRecombiner &flav_recombiner = FlavRecombiner());

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_GHSPLUGIN_HH__
