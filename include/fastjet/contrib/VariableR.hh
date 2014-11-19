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

#ifndef __FASTJET_CONTRIB_VARIABLER_HH__
#define __FASTJET_CONTRIB_VARIABLER_HH__

#include <fastjet/internal/base.hh>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <map>
#include <queue>

using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

////////
//
//  Helper Functions and Types
//
////////

void print_jets (const fastjet::ClusterSequence &, 
                 const vector<fastjet::PseudoJet> &);

// Store two jets and a distance measure
struct JetDistancePair{
  int j1,j2;
  double distance;
};

struct PTVectorAndMap{
  bool pts_sorted, boost_vector_set;
  PseudoJet boost_vector;
  vector<double> pts;
  map<double,int> map_pt_to_indx;
  ClusterSequence * clust_seq;
};

struct ThreeVector{
  double x,y,z;
};

enum MergeType {Unmerged, MergedWithBeam, MergedWithJet};

// Class for comparing JetDistancePairs
class CompareJetDistancePair {
public:
  CompareJetDistancePair(){};
  bool operator() (const JetDistancePair & lhs, const JetDistancePair &rhs) const {
    return (lhs.distance > rhs.distance);
  }
};

////////
//
//  CoreJetAlgorithm
//
////////

class CoreJetAlgorithm : public JetDefinition::Plugin {
public:
  CoreJetAlgorithm (){};

  virtual bool supports_ghosted_passive_areas() const {return true;}
  virtual void set_ghost_separation_scale(double scale) const {};

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;                      
  virtual double R() const {return 0.0;}

private:

  PTVectorAndMap GetPTVectorAndMap(map<int,MergeType>& jet_merged,ClusterSequence & clust_seq) const;
  void run_smart_clustering(ClusterSequence &) const;                      
  virtual double GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const;
  virtual double GetJBDistanceMeasure(const PseudoJet& j1) const;                   
  virtual double GetJJDistanceMeasureFull(const PseudoJet& j1, const PseudoJet& j2, double R, PTVectorAndMap & ptvm) const;
  virtual double GetJBDistanceMeasureFull(const PseudoJet& j1, PTVectorAndMap & ptvm) const;
  virtual bool ShouldCompilePTInfo() const;
protected:
  void SortPTVectorAndMap(PTVectorAndMap &ptvm, int num_vects)const;
  void SetupDistanceMeasures(ClusterSequence & clust_seq, priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair >& jet_queue,   map<int,MergeType> & jet_merged) const;
  void MergeJetWithBeam(ClusterSequence & clust_seq, JetDistancePair & jdp, map<int,MergeType>& jet_merged) const;
  void MergeJets(ClusterSequence & clust_seq,JetDistancePair & jdp,priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue, map<int,MergeType>& jet_merged) const;
  void NormalizeThreeVector(ThreeVector & tv) const;
  ThreeVector GetThreeVector(const PseudoJet& jet) const;
  double Dot3D(ThreeVector &tv1, ThreeVector& tv2) const;
  double GetAngleBetweenJets(const PseudoJet& j1, const PseudoJet& j2) const;
  double GetDeltaRBetweenJets(const PseudoJet& j1, const PseudoJet& j2) const;
  void PrintJet(const PseudoJet& jet) const;
};

////////
//
//  AKTVR
//
////////

class AKTVR : public CoreJetAlgorithm{
public:
  AKTVR (double rho, double max_r);
  
  virtual std::string description () const;
  virtual double R() const {return 0.0;}
private:
  double my_rho, my_max_r;

  double GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const;
  double GetJBDistanceMeasure(const PseudoJet& j1) const;
};

////////
//
//  CAVR
//
////////

class CAVR : public CoreJetAlgorithm{
public:
  CAVR (double rho, double max_r);
  
  virtual std::string description () const;
  virtual double R() const {return 0.0;}
private:
  double my_rho, my_max_r;

  double GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2, double R) const;
  double GetJBDistanceMeasure(const PseudoJet& j1) const;
};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VARIABLER_HH__
