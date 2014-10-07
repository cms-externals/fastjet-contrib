//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2013
//  David Krohn, Jesse Thaler, and Lian-Tao Wang
//
//  $Id: VariableR.cc 596 2014-04-16 22:15:15Z jthaler $
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

#include "VariableRPlugin.hh"

#include <cstdio>
#include "math.h"
#include <iomanip>
#include <cmath>
#include <map>
#include <sstream>
#include <queue>


#define PI 3.14159265
#define RPARAM 1.0
#define LARGE_NUMBER 10e20

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
   
   // Default constructor.  Just sets parameters
   VariableRPlugin::VariableRPlugin(double rho, double min_r, double max_r, ClusterType clust_type, bool precluster)
   :_rho2(rho*rho), _min_r2(min_r*min_r), _max_r2(max_r*max_r), _clust_type(clust_type), _precluster(precluster),
    _pre_jet_def(kt_algorithm, min_r) // at the moment, only option for preclustering is kT
   {
      // Added errors for user input.
      if (min_r < 0.0) throw Error("VariableRPlugin: Minimum radius must be positive.");
      if (precluster && min_r == 0.0) throw Error("VariableRPlugin: To apply preclustering, minimum radius must be non-zero.");
      if (max_r < 0.0) throw Error("VariableRPlugin: Maximum radius must be positive.");
      if (min_r > max_r) throw Error("VariableRPlugin: Minimum radius must be bigger than or equal to maximum radius.");
   }
   
   // precluster into kT subjets if desired.
   void VariableRPlugin::Precluster(ClusterSequence & cs, set<int>& unmerged_jets) const {
      int cntr(0);
      for(vector<PseudoJet>::const_iterator it = cs.jets().begin(); it != cs.jets().end(); it++)
         unmerged_jets.insert(unmerged_jets.end(), cntr++);
      
      // Make preclusters
      ClusterSequence pre_cs(cs.jets(), _pre_jet_def);
      vector<PseudoJet> preclustered_jets = pre_cs.inclusive_jets();
      vector<int> particle_jet_indices = pre_cs.particle_jet_indices(preclustered_jets);
      
      // This code take preclustered objects and puts them into the ClusterSequence tree of VR
      // TODO:  Figure out if there is a better way to do this step.
      for(int i = 0 ; i < (int)preclustered_jets.size(); i++){
         queue<int> constit_indices;
         for(int j = 0 ; j < (int)particle_jet_indices.size(); j++)
            if(particle_jet_indices[j] == i)
               constit_indices.push(j);
         
         int final_jet;
         while(constit_indices.size() > 1){
            int indx1 = constit_indices.front();
            unmerged_jets.erase(indx1);
            constit_indices.pop();
            int indx2 = constit_indices.front();
            unmerged_jets.erase(indx2);
            constit_indices.pop();
            cs.plugin_record_ij_recombination(indx1, indx2, 0., final_jet);
            constit_indices.push(final_jet);
            unmerged_jets.insert(unmerged_jets.end(), final_jet);
         }
      }
   }
   
   // Implements VR alorithm (virtual function from JetDefinition::Plugin)
   void VariableRPlugin::run_clustering(ClusterSequence & cs) const {
      set<int> unmerged_jets;
      
      if(_precluster){ // do kT preclustering
         assert(_min_r2 > 0.);  // since this is (min_r)^2, this should never happen
         Precluster(cs, unmerged_jets);
      } else // make a list of the unmerged jets
         for(int i = 0 ; i < (int)cs.jets().size(); i++)
            unmerged_jets.insert(unmerged_jets.end(), i);
      
      // priority_queue is part of the C++ standard library.
      // The objects being sorted are JetDistancePairs
      // The comparison function is just asking who has the smallest distance
      // The distances are set initially by SetupDistanceMeasures and then updated by the while loop below.
      vector<JetDistancePair> jet_vec;
      SetupDistanceMeasures(cs, jet_vec, unmerged_jets);
      priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > jet_queue(jet_vec.begin(),jet_vec.end());
      
      // go through the jet_queue until empty
      while(!jet_queue.empty()){
         
         // find the closest pair
         JetDistancePair jdpair = jet_queue.top();
         jet_queue.pop();
         
         // Rebuild the jet_queue
         // DK - the 1.5 below was just found empirically
         // JDT - this code is safe, but I'm not 100% sure why it is needed.
         // It rebuilds the jet_queue instead of letting the below functions do the hard work.
         if(jet_queue.size() > 50 && jet_queue.size() > 1.5*unmerged_jets.size()*unmerged_jets.size()){
            jet_vec.clear();
            SetupDistanceMeasures(cs, jet_vec, unmerged_jets);
            priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > tmp_jet_queue(jet_vec.begin(),jet_vec.end());
            swap(jet_queue,tmp_jet_queue);
         }
         
         // make sure not merged
         if((unmerged_jets.find(jdpair.j1) == unmerged_jets.end()) || (jdpair.j2 != -1 && unmerged_jets.find(jdpair.j2) == unmerged_jets.end()))
            continue;
         
         
         if(jdpair.j2 == -1) // If closest distance is to beam, then merge with beam
            MergeJetWithBeam(cs, jdpair, unmerged_jets);
         else // Otherwise, merge jets back together
            MergeJets(cs, jdpair, jet_queue, unmerged_jets);
      }
   }
   
   // Description of algorithm, including parameters
   string VariableRPlugin::description() const{
      stringstream myStream("");
      
      myStream << "Variable R (0903.0392), ";
      
      switch (_clust_type) {
         case AKTLIKE:
            myStream << "AKT";
            break;
         case CALIKE:
            myStream << "CA";
            break;
         case KTLIKE:
            myStream << "KT";
            break;
      }
      
      myStream << fixed << setprecision(1) << ", rho=" << sqrt(_rho2);
      myStream << ", min_r=" << sqrt(_min_r2);
      myStream << ", max_r=" << sqrt(_max_r2);
      myStream << (_precluster ? ", with precluster" : "");
      
      return myStream.str();
   }
   
   // TODO:  have this return something sensible (not sure if this is possible)
   double VariableRPlugin::R() const{
      return -1.;
   }
   
   // Add final jet to clust_seq.  No need to update jet_queue, since this jet has already been deleted.
   void VariableRPlugin::MergeJetWithBeam(ClusterSequence & clust_seq, JetDistancePair & jdp, set<int>& unmerged_jets) const{
      clust_seq.plugin_record_iB_recombination(jdp.j1, jdp.distance);
      unmerged_jets.erase(jdp.j1);
   }
   
   
   // Add jet merging to clust_seq.  Here, we need to update the priority_queue since some of the measures have changes.
   void VariableRPlugin::MergeJets(ClusterSequence & clust_seq,
                      JetDistancePair & jdp,
                      priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue,
                      set<int>& unmerged_jets) const{
      
      
      int new_jet_num;
      clust_seq.plugin_record_ij_recombination(jdp.j1,jdp.j2,jdp.distance,new_jet_num);

      unmerged_jets.erase(jdp.j1);
      unmerged_jets.erase(jdp.j2);
      
      // take the resulting jet and recompute all distances with it
      for(set<int>::iterator it = unmerged_jets.begin(); it != unmerged_jets.end(); it++){
         JetDistancePair jpair;
         jpair.j1 = new_jet_num;
         jpair.j2 = (*it);
         jpair.distance = GetJJDistanceMeasure(clust_seq.jets()[*it],clust_seq.jets()[new_jet_num]);
         jet_queue.push(jpair);
      }
      unmerged_jets.insert(unmerged_jets.end(), new_jet_num);
      
      // also add the new distance to beam
      JetDistancePair jpair;
      jpair.j1  = new_jet_num;
      jpair.j2 = -1; // -1 is for the beam
      jpair.distance = GetJBDistanceMeasure(clust_seq.jets()[new_jet_num]);
      jet_queue.push(jpair);
   }
   
   // Initial distance setup
   void VariableRPlugin::SetupDistanceMeasures(ClusterSequence & clust_seq,
                                  vector<JetDistancePair> &jet_vec,
                                  set<int> & unmerged_jets
                                  ) const{
      JetDistancePair jpair;
      
      // Add the jet-jet distances
      for(set<int>::iterator it1 = unmerged_jets.begin(); it1 != unmerged_jets.end(); it1++){
         for(set<int>::iterator it2 = it1; it2 != unmerged_jets.end(); it2++){
            if((*it1) != (*it2)){
               jpair.j1 = (*it1);
               jpair.j2 = (*it2);
               jpair.distance = GetJJDistanceMeasure(clust_seq.jets()[*it1],clust_seq.jets()[*it2]);
               jet_vec.push_back(jpair);
            }
         }
         // Add the jet-beam distances, and set initial merge info
         jpair.j1  = (*it1);
         jpair.j2 = -1; // -1 is for the beam
         jpair.distance = GetJBDistanceMeasure(clust_seq.jets()[*it1]);
         jet_vec.push_back(jpair);
      }
   }
   
   // get the dij between two jets
   // Different measures for AKTLIKE, CALIKE, and KTLIKE
   // TODO:  add arbitrary p weighting
   double VariableRPlugin::GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2) const {
      double ret;
      switch(_clust_type){
         case AKTLIKE:
            ret = min(1./j1.perp2(), 1./j2.perp2());
            break;
         case CALIKE :
            ret = 1.;
            break;
         case KTLIKE:
            ret = min(j1.perp2(), j2.perp2());
            break;
         default:
            assert(false);
      }
      
      ret *= j1.squared_distance(j2);
      return ret;
   }
   
   // jet diB between jet and beam
   // Different measures for AKTLIKE, CALIKE, and KTLIKE
   // TODO:  add arbitrary p weighting
   double VariableRPlugin::GetJBDistanceMeasure(const PseudoJet& jet) const{
      switch(_clust_type){
         case AKTLIKE:
         {
            if((_rho2 / jet.perp2()) < _min_r2)
               return _min_r2/jet.perp2();
            if((_rho2 / jet.perp2()) > _max_r2)
               return _max_r2/jet.perp2();
            return _rho2/jet.perp2()/jet.perp2();
         }
         case CALIKE :
         {
            if((_rho2 / jet.perp2()) < _min_r2)
               return _min_r2;
            if((_rho2 / jet.perp2()) > _max_r2)
               return _max_r2;
            return _rho2 /jet.perp2();
         }
            break;
         case KTLIKE:
         {
            if((_rho2 / jet.perp2()) < _min_r2)
               return _min_r2*jet.perp2();
            if((_rho2 / jet.perp2()) > _max_r2)
               return _max_r2*jet.perp2();
            return _rho2;
         }
         default:
            assert(false);
      }
   }
   
   
   
} // namespace contrib

FASTJET_END_NAMESPACE
