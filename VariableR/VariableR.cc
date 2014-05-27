//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2013
//  David Krohn, Jesse Thaler, and Lian-Tao Wang
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

#include "VariableR.hh"

#include <cstdio>
#include "math.h"
#include <cmath>
#include <map>
#include <sstream>
#include <queue>

#define PI 3.14159265
#define RPARAM 1.0
#define LARGE_NUMBER 10e20

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {


////////
//
//  Helper Functions and Types
//
////////

//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const fastjet::ClusterSequence & clust_seq, 
		 const vector<fastjet::PseudoJet> & jets) {

  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  

  // label the columns
  printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", 
	 "phi", "pt","m","e", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = clust_seq.constituents(sorted_jets[i]).size();
    printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].perp(),sorted_jets[i].m(),sorted_jets[i].e(), n_constituents);
  }
}


// greater than comparator
bool geq (int i,int j) { return (i>j); }


////////
//
//  CoreJetAlgorithm
//
////////



string CoreJetAlgorithm::description () const {
  ostringstream desc;
  
  desc << "CoreJetAlgorithm base class";

  return desc.str();
}
	
void CoreJetAlgorithm::SetupDistanceMeasures(ClusterSequence & clust_seq, priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue, map<int,MergeType> & jet_merged) const{
  int num_init_jets = clust_seq.jets().size();
  
  PTVectorAndMap ptvm = GetPTVectorAndMap(jet_merged,clust_seq);

  // Add the jet-jet distances
  for(int i = 0 ; i < num_init_jets-1 ; i++){
    PseudoJet j1 = clust_seq.jets()[i];
    for(int j = i+1 ; j < num_init_jets ; j++){
      JetDistancePair jpair;
      jpair.j1 = i;
      jpair.j2 = j;
      jpair.distance = GetJJDistanceMeasureFull(j1,clust_seq.jets()[j],RPARAM,ptvm);
      jet_queue.push(jpair);
    }
  }

  // Add the jet-beam distances, and set initial merge info
  for(int i = 0 ; i < num_init_jets ; i++){
    JetDistancePair jpair;
    jpair.j1  = i;
    jpair.j2 = -1; // -1 is for the beam
    jpair.distance = GetJBDistanceMeasureFull(clust_seq.jets()[i],ptvm);
    jet_queue.push(jpair);
  }
}

void CoreJetAlgorithm::run_clustering(ClusterSequence & clust_seq) const {
  run_smart_clustering(clust_seq);
}

void CoreJetAlgorithm::run_smart_clustering(ClusterSequence & clust_seq) const {
  map<int,MergeType> jet_merged;
  for(unsigned int i = 0 ; i < clust_seq.jets().size() ; i++)
    jet_merged[i] = Unmerged;

  priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > jet_queue;
  
  SetupDistanceMeasures(clust_seq,jet_queue, jet_merged);
   
  while(!jet_queue.empty()){
    JetDistancePair jdpair = jet_queue.top();
    jet_queue.pop();

    // make sure not merged
    if(jet_merged[jdpair.j1]!=Unmerged || jet_merged[jdpair.j2] != Unmerged)
      continue;
    
    if(jdpair.j2 == -1)
      MergeJetWithBeam(clust_seq,jdpair,jet_merged);
    else 
      MergeJets(clust_seq,jdpair,jet_queue,jet_merged);
  }
}

// get the dij between two jets
double CoreJetAlgorithm::GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const {
  return 0.0;
}

double CoreJetAlgorithm::GetAngleBetweenJets(const PseudoJet& j1, const PseudoJet& j2) const {
  ThreeVector tv1 = GetThreeVector(j1);
  ThreeVector tv2 = GetThreeVector(j2);

  NormalizeThreeVector(tv1);
  NormalizeThreeVector(tv2);

  double dot = Dot3D(tv1,tv2);

  if(dot>1)
    dot = 1.0;
  if(dot<-1)
    dot = -1.0;

  double angle = acos(dot);

  if(std::isnan(angle))
    angle = 2*PI;

  return angle;
}

double CoreJetAlgorithm::GetDeltaRBetweenJets(const PseudoJet& j1, const PseudoJet& j2) const {
  return std::pow(j1.plain_distance(j2),0.5);
}

ThreeVector CoreJetAlgorithm::GetThreeVector(const PseudoJet& jet) const {
  ThreeVector tv;
  tv.x = jet.px();
  tv.y = jet.py();
  tv.z = jet.pz();

  return tv;
}

void CoreJetAlgorithm::PrintJet(const PseudoJet& jet) const {
  cout << "(E,pt,eta,phi,m)=(" << jet.e()<<", "<<jet.perp()<<", "<<jet.eta() <<", " << jet.m() <<", " << jet.phi()<<  ")" << endl;
}

void CoreJetAlgorithm::NormalizeThreeVector(ThreeVector & tv) const {
  double mag = std::pow(std::pow(tv.x,2.0)+std::pow(tv.y,2.0)+std::pow(tv.z,2.0),0.5);
  tv.x *= 1/mag;
  tv.y *= 1/mag;
  tv.z *= 1/mag;
}
 
double CoreJetAlgorithm::Dot3D(ThreeVector &tv1, ThreeVector& tv2) const {
  return tv1.x * tv2.x +  tv1.y * tv2.y +  tv1.z * tv2.z;
}

double CoreJetAlgorithm::GetJBDistanceMeasure(const PseudoJet& jet) const{
  return 0.0;
}

double CoreJetAlgorithm::GetJJDistanceMeasureFull(const PseudoJet& j1, const PseudoJet& j2, double R, PTVectorAndMap & ptvm) const{
  return GetJJDistanceMeasure(j1,j2,R);
}

double CoreJetAlgorithm::GetJBDistanceMeasureFull(const PseudoJet& j1,PTVectorAndMap & ptvm) const{
  return GetJBDistanceMeasure(j1);
}

void CoreJetAlgorithm::MergeJetWithBeam(ClusterSequence & clust_seq, JetDistancePair & jdp, map<int,MergeType>& jet_merged) const{
  clust_seq.plugin_record_iB_recombination(jdp.j1,jdp.distance);
  jet_merged[jdp.j1] = MergedWithBeam;
}  

void CoreJetAlgorithm::MergeJets(ClusterSequence & clust_seq,JetDistancePair & jdp,priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue, map<int,MergeType>& jet_merged) const{
int new_jet_num;

  PTVectorAndMap ptvm = GetPTVectorAndMap(jet_merged,clust_seq);
  clust_seq.plugin_record_ij_recombination(jdp.j1,jdp.j2,jdp.distance,new_jet_num);

  jet_merged[jdp.j1] = MergedWithJet;
  jet_merged[jdp.j2] = MergedWithJet;
  jet_merged[new_jet_num] = Unmerged;
  
  // take the resulting jet and recompute all distances with it
  PseudoJet new_jet = clust_seq.jets()[new_jet_num];
  for(unsigned int j = 0 ; j < clust_seq.jets().size() ; j++){
    
    if((int) j != new_jet_num && jet_merged[j]==Unmerged){
      JetDistancePair jpair;
      jpair.j1  = new_jet_num;
      jpair.j2 = j;
      jpair.distance = GetJJDistanceMeasureFull(clust_seq.jets()[j],new_jet,RPARAM,ptvm);
      jet_queue.push(jpair);
    }
  }
  // also add the new distance to beam
  JetDistancePair jpair;
  jpair.j1  = new_jet_num;
  jpair.j2 = -1; // -1 is for the beam
  jpair.distance = GetJBDistanceMeasureFull(clust_seq.jets()[new_jet_num], ptvm);
  jet_queue.push(jpair);
}

/*
  This information is needed for some algorithms
*/
PTVectorAndMap CoreJetAlgorithm::GetPTVectorAndMap(map<int,MergeType>& jet_merged,ClusterSequence & clust_seq) const{
  PTVectorAndMap ptvm;
  ptvm.pts_sorted = false;
  ptvm.boost_vector_set = false;
  ptvm.clust_seq = &clust_seq;
  
  if(ShouldCompilePTInfo())
    for(unsigned int i = 0 ; i < clust_seq.jets().size(); i++)
      if(jet_merged[i]!=MergedWithJet){
	double pt = clust_seq.jets()[i].perp();
	ptvm.pts.push_back(pt);
	ptvm.map_pt_to_indx[pt] = i;
      }
  return ptvm;
}

bool CoreJetAlgorithm::ShouldCompilePTInfo() const{
  return false;
}

void CoreJetAlgorithm::SortPTVectorAndMap(PTVectorAndMap &ptvm, int num_vects)const{
  if(!ptvm.pts_sorted){
    partial_sort (ptvm.pts.begin(), ptvm.pts.begin()+num_vects, ptvm.pts.end(),geq);
    ptvm.pts_sorted = true;
  }
}


////////
//
//  AKTVR
//
////////

AKTVR::AKTVR(double rho, double max_r)
  :my_rho(rho), my_max_r(max_r)
{
}

string AKTVR::description () const {
  ostringstream desc;
  
  desc << "AKTVR";

  return desc.str();
}

// get the dij between two jets
double AKTVR::GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const{
  double ptm2_j1 = std::pow(j1.perp(),-2.0);
  double ptm2_j2 = std::pow(j2.perp(),-2.0);
  double min_ptm2 = std::min(ptm2_j1,ptm2_j2);
  
  double dr = GetDeltaRBetweenJets(j1,j2);

  // want a maximum effective delta-R, so don't cluster if 
  // delta-r beyond the specified value
  
  if(dr > my_max_r)
    return LARGE_NUMBER;

  return min_ptm2 * std::pow(dr,2.0);
}

// get the jet-beam distance
double AKTVR::GetJBDistanceMeasure(const PseudoJet& jet) const{
  double ptm4 = std::pow(jet.perp(),-4.0);

  return std::pow(my_rho,2.0) * ptm4;
}

////////
//
//  CAVR
//
////////


CAVR::CAVR(double rho, double max_r)
  :my_rho(rho), my_max_r(max_r)
{
}

string CAVR::description () const {
  ostringstream desc;
  
  desc << "CAVR";

  return desc.str();
}

// get the dij between two jets
double CAVR::GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const{
  double dr = GetDeltaRBetweenJets(j1,j2);

  // want a maximum effective delta-R, so don't cluster if 
  // delta-r beyond the specified value
  
  if(dr > my_max_r)
    return LARGE_NUMBER;

  return std::pow(dr,2.0);
}

// get the jet-beam distance
double CAVR::GetJBDistanceMeasure(const PseudoJet& jet) const{
  double ptm2 = std::pow(jet.perp(),-2.0);
  return std::pow(my_rho,2.0) * ptm2;
}


} // namespace contrib

FASTJET_END_NAMESPACE
