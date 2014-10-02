//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2013
//  David Krohn, Jesse Thaler, and Lian-Tao Wang
//
//  $Id$
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

#ifndef __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__
#define __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__

#include <fastjet/internal/base.hh>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/LimitedWarning.hh>

#include <map>
#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
   
   ////////
   //
   //  Core VR Code
   //
   ////////
   
   // In version 1.1, this replaces CoreJetAlgorithm.
   // This acts like any fastjet plugin since it implements run_clustering
   class VariableRPlugin : public JetDefinition::Plugin {

   public:
      // Type of clustering
      enum ClusterType {
         CALIKE,
         KTLIKE,
         AKTLIKE
      };
      
      // Constructor that sets VR algorithm parameters
      // rho = mass scale for effective radius (i.e. R ~ rho/pT)
      // min_r = minimum jet radius
      // max_r = maximum jet radius
      // clust_type = whether to use CA-like, kT-like, or anti-kT-like distance measure
      // precluster = whether to use optional kT subjets (of size min_r) for preclustering
      // (precluster = true is much faster)
      VariableRPlugin(double rho, double min_r, double max_r, ClusterType clust_type, bool precluster = false);
      
      // virtual function from JetDefinition::Plugin that implements the actual VR algorithm
      void run_clustering(fastjet::ClusterSequence & cs) const;
      
      // information string
      virtual string description() const;
      
      // TODO:  have this return a non-trivial answer.
      virtual double R() const;
      
   private:
      
      // parameters to define VR algorithm
      double _rho2, _min_r2, _max_r2;
      ClusterType _clust_type;
      
      // For preclustering, can use kT algorithm to make subclusters of size min_r
      bool _precluster;
      JetDefinition _pre_jet_def;
      
      // helper function to apply kT preclustering if desired
      void Precluster(ClusterSequence & cs, set<int>& unmerged_jets) const;
      
      // Helper struct to store two jets and a distance measure
      struct JetDistancePair{
         int j1,j2;
         double distance;
      };
      
      // Helper comparitor class for comparing JetDistancePairs
      class CompareJetDistancePair {
      public:
         CompareJetDistancePair(){};
         bool operator() (const JetDistancePair & lhs, const JetDistancePair &rhs) const {
            return (lhs.distance > rhs.distance);
         }
      };
      
      // helper function to merge jet with beam
      void MergeJetWithBeam(ClusterSequence & clust_seq, JetDistancePair & jdp, set<int>& unmerged_jets) const;
      
      // helper function to merge two jets into a new pseudojet
      void MergeJets(ClusterSequence & clust_seq,
                     JetDistancePair & jdp,
                     priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue,
                     set<int>& unmerged_jets) const;
      
      // use ClusterMode to determine jet-jet and jet-beam distance
      inline double GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2) const;
      inline double GetJBDistanceMeasure(const PseudoJet& jet) const;
      
      // helper function to establish measures
      void SetupDistanceMeasures(ClusterSequence & clust_seq,
                                 vector<JetDistancePair> &jet_vec,
                                 set<int>& unmerged_jets) const;

   };
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__
