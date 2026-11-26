//--------------------------------------------------------------------
//
// file: DynamicR.cc
//
//----------------------------------------------------------------------
// $Id: DynamicR.cc 1567 2026-06-24 12:19:44Z tousik $
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

// fastjet functionalities
#include "fastjet/ClusterSequence.hh"   // FastJet cluster sequence
#include "fastjet/NNH.hh"               // Nearest-neighbour helper used during clustering
#include "fastjet/PseudoJet.hh"         // FastJet PseudoJet
#include "fastjet/contrib/DynamicR.hh"  // Declaration of the DynamicR class

// other stuff
#include <list>     // std::list
#include <memory>   // std::unique_ptr
#include <cmath>    // sqrt and other math utilities
#include <vector>   // std::vector
#include <sstream>  // std::ostringstream

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

//----------------------------------------------------------------------
// Lightweight jet wrapper for the Dynamic-R anti-kt-like metric
//----------------------------------------------------------------------
// This class provides the minimal interface required by NNH:
//   - init(...)
//   - distance(...)
//   - beam_distance(...)
// It stores a copy of the jet and the extra information needed to compute
// the clustering and beam distances for the DRAK metric.
class DRAKBriefJet {
public:
  // Initialise from a PseudoJet and plugin-wide global information.
  // Also extract stored DynamicRJetInfo, if available, in order to use the
  // running mean/mean-squared radius modifier information of previously
  // clustered evolving proto-jets.
  void init(const PseudoJet & jet, const DynamicRGlobalInfo * globalInfo){
    _jeti = jet;               // store jet object
    _pti2 = jet.pt2();         // store transverse momentum squared
    if (jet.has_user_info<DynamicRJetInfo>() ) {
      // retrieve accumulated mean radius information from previous clustering
      _avg_R = jet.user_info<DynamicRJetInfo>().mean_R();
      // variance canculated from mean-squared and mean
      _sd_R  = jet.user_info<DynamicRJetInfo>().ms_R() - _avg_R * _avg_R;
      // standard deviation (preserve sign in case of numerical negative values)
      _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    }
    else {
      // for original input particles there is no stored radius information
      _avg_R = 0.0;
      _sd_R  = 0.0;
    }
    _r = globalInfo->R0();      // initial reference radius parameter
  }

  // Pairwise distance for the anti-kt-like Dynamic-R metric:
  // dij = DeltaR^2 / max(pt_i^2, pt_j^2)
  double distance(const DRAKBriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti); // geometric distance in rapidity-phi plane
    double pt2max = jet->_pti2;                          // start with other jet pt^2
    if (pt2max < _pti2) pt2max = _pti2;                 // take the maximum pt^2
    return distance2/pt2max;
  }

  // Beam distance for this metric.
  // The effective radius is modified by (_r + _sd_R).
  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num/_pti2; 
    return diB; 
  }

private:
  PseudoJet _jeti;                  // stored jet
  double _pti2, _avg_R, _sd_R, _r;  // pt^2, mean R, s.d. in R, reference radius

};

//----------------------------------------------------------------------
// Lightweight jet wrapper for the Dynamic-R Cambridge/Aachen-like metric
//----------------------------------------------------------------------
class DRCABriefJet {
public:
  // Initialise the wrapper from the PseudoJet and plugin settings.
  void init(const PseudoJet & jet, const DynamicRGlobalInfo * globalInfo){
    _jeti = jet;               // store jet object
    if (jet.has_user_info<DynamicRJetInfo>()) {
      // retrieve accumulated mean/mean-squared radius information if available
      _avg_R = jet.user_info<DynamicRJetInfo>().mean_R();
      _sd_R  = jet.user_info<DynamicRJetInfo>().ms_R() - _avg_R * _avg_R;
      _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    }
    else {
      // no previous user info for input particles
      _avg_R = 0.0;
      _sd_R  = 0.0;
    }
    _r = globalInfo->R0();      // base radius parameter
  }

  // Pairwise distance for the CA-like Dynamic-R metric:
  // dij = DeltaR^2
  double distance(const DRCABriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    return distance2;
  }

  // Beam distance for the CA-like metric.
  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num; 
    return diB; 
  }

private:
  PseudoJet _jeti;               // stored jet
  double _avg_R, _sd_R, _r;      // mean R, s.d. in R, reference radius

};

//----------------------------------------------------------------------
// Lightweight jet wrapper for the Dynamic-R kt-like metric
//----------------------------------------------------------------------
class DRKTBriefJet {
public:
  // Initialise from the jet and plugin-wide global information.
  void init(const PseudoJet & jet, const DynamicRGlobalInfo * globalInfo){
    _jeti = jet;               // store jet object
    _pti2 = jet.pt2();         // store transverse momentum squared
    if (jet.has_user_info<DynamicRJetInfo>()) {
      // retrieve accumulated mean/mean-squared radius information if available
      _avg_R = jet.user_info<DynamicRJetInfo>().mean_R();
      _sd_R  = jet.user_info<DynamicRJetInfo>().ms_R() - _avg_R * _avg_R;
      _sd_R  = (_sd_R > 0.0) ? sqrt(_sd_R) : -sqrt(-_sd_R);
    }
    else {
      // original input particles start with zero stored radius information
      _avg_R = 0.0;
      _sd_R  = 0.0;
    }
    _r = globalInfo->R0();      // base radius parameter
  }

  // Pairwise distance for the kt-like Dynamic-R metric:
  // dij = DeltaR^2 * min(pt_i^2, pt_j^2)
  double distance(const DRKTBriefJet * jet) const {
    double distance2 = jet->_jeti.plain_distance(_jeti);
    double pt2min = jet->_pti2;                          // start with other jet pt^2
    if (pt2min > _pti2) pt2min = _pti2;                 // take the minimum pt^2
    return distance2*pt2min;
  }

  // Beam distance for the kt-like metric.
  double beam_distance() const {
    double num = (_r + _sd_R);
    double diB = num * num * _pti2; 
    return diB; 
  }

private:
  PseudoJet _jeti;                  // stored jet
  double _pti2, _avg_R, _sd_R, _r;  // pt^2, mean R, s.d. in R, reference radius

};


//------------------------------------------------------------------
// implementation of the DynamicR plugin
//------------------------------------------------------------------

// Static flag to ensure the banner is printed only once.
bool DynamicR::_first_time = true;

// Return a human-readable description of the plugin configuration.
string DynamicR::description () const {
  ostringstream desc;
  switch(_jet_algorithm) {
    case DRKT_algorithm:
      // kt-like dynamic radius algorithm
      desc << "Dynamic Radius (KT) Jet algorithm with R0 = " << R();
      break;
    case DRCA_algorithm:
      // Cambridge/Aachen-like dynamic radius algorithm
      desc << "Dynamic Radius (CA) Jet algorithm with R0 = " << R();
      break;
    case DRAK_algorithm:
      // anti-kt-like dynamic radius algorithm
      desc << "Dynamic Radius (AK) Jet algorithm with R0 = " << R();
      break;
    default:
      // safety check for unsupported algorithm choice
      throw Error("unrecognized jet_algorithm");
  }

  return desc.str();
}



// Dispatch the clustering to the appropriate brief-jet metric class
// depending on the selected dynamic-radius algorithm.
void DynamicR::run_clustering(ClusterSequence & clust_seq) const {

  switch(_jet_algorithm) {
    case DRKT_algorithm:
      actual_run<DRKTBriefJet>(clust_seq); // kt-like metric
      break;
    case DRCA_algorithm:
      actual_run<DRCABriefJet>(clust_seq); // CA-like metric
      break;
    case DRAK_algorithm:
      actual_run<DRAKBriefJet>(clust_seq); // anti-kt-like metric
      break;
    default:
      throw Error("unrecognized jet_algorithm");
  }
}

// Main clustering routine, templated on the brief-jet type BJ.
// BJ defines how dij and diB are computed.
template<class BJ> void DynamicR::actual_run(ClusterSequence & clust_seq) const {

  int njets = clust_seq.jets().size(); // current number of active jets
  int insize = njets;                  // number of original input particles/jets
  PseudoJet newjet, jet_i, jet_j;      // temporary jet holders
  double mean_DR, ms_DR, wt;           // accumulated radius statistics for merged jet
  double mean_DR1, ms_DR1;             // temporary accumulators for nested loops
  double jet_i_pt, jet_j_pt;           // transverse momentum of the two jets
  double DR2ij;                        // pairwise distance quantities
  double pti, ptj;                     // constituent transverse momenta

  double radius = R();                 // base radius parameter
  DynamicRGlobalInfo globalInfo(radius);         // helper object passed to NNH

  // Build the nearest-neighbour helper on the current jets using the
  // metric encoded by BJ.
  NNH<BJ, DynamicRGlobalInfo> nn(clust_seq.jets(), &globalInfo);

  // Standard FastJet plugin clustering loop:
  // repeatedly merge the pair with smallest dij or remove a jet to beam.
  while (njets > 0) {
    int i, j, k;                       // i,j = pair to merge; k = new jet index
    double dij = nn.dij_min(i, j);     // smallest pair/beam distance
    if (j >= 0) {
      // A pairwise recombination is to be performed.
      jet_i  = clust_seq.jets()[i];
      jet_j  = clust_seq.jets()[j];
      
      if (i < insize && j < insize) {
        // Case 1:
        // both jets are original input particles.
        // The merged-jet mean and mean-squared radius modifiers are computed directly from their
        // separation and pt weights.
        jet_i_pt = jet_i.pt();
        jet_j_pt = jet_j.pt();
        ms_DR = jet_i.plain_distance(jet_j);
        mean_DR = sqrt(ms_DR);
        wt = jet_i_pt*jet_j_pt;
      }
      else if ( i < insize && j >= insize) {
        // Case 2:
        // i is an original particle, j is a previously merged jet.
        // Recompute cross terms between particle i and all constituents of j,
        // then combine with the stored DynamicRJetInfo of jet j.
        jet_i_pt = jet_i.pt();
        mean_DR = 0;
        ms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_j_constituents = jet_j.constituents();
        for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++){
          DR2ij = jet_i.plain_distance(jet_j_constituents[jt]);
          ptj  = jet_j_constituents[jt].pt();
          ms_DR  += DR2ij*ptj;
          mean_DR += sqrt(DR2ij)*ptj;
          wt      += ptj;
        }
        mean_DR = mean_DR*jet_i_pt;
        ms_DR  = ms_DR *jet_i_pt;
        wt      = wt     *jet_i_pt;

        // Add previously stored intra-jet information from jet_j and normalize.
        wt      = wt + jet_j.user_info<DynamicRJetInfo>().wt();
        mean_DR = (mean_DR + jet_j.user_info<DynamicRJetInfo>().mean_R()*jet_j.user_info<DynamicRJetInfo>().wt())/wt;
        ms_DR  = (ms_DR + jet_j.user_info<DynamicRJetInfo>().ms_R()*jet_j.user_info<DynamicRJetInfo>().wt())/wt;
      }
      else if ( i >= insize && j < insize) {
        // Case 3:
        // i is a previously merged jet, j is an original particle.
        // Same logic as Case 2, with i and j interchanged.
        jet_j_pt = jet_j.pt();
        mean_DR = 0;
        ms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_i_constituents = jet_i.constituents();
        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          DR2ij = jet_j.plain_distance(jet_i_constituents[it]);
          pti  = jet_i_constituents[it].pt();
          ms_DR  += DR2ij*pti;
          mean_DR += sqrt(DR2ij)*pti;
          wt      += pti;
        }
        mean_DR = mean_DR*jet_j_pt;
        ms_DR  = ms_DR *jet_j_pt;
        wt      = wt     *jet_j_pt;

        // Add previously stored intra-jet information from jet_i and normalize.
        wt      = wt + jet_i.user_info<DynamicRJetInfo>().wt();
        mean_DR = (mean_DR + jet_i.user_info<DynamicRJetInfo>().mean_R()*jet_i.user_info<DynamicRJetInfo>().wt())/wt;
        ms_DR  = (ms_DR + jet_i.user_info<DynamicRJetInfo>().ms_R()*jet_i.user_info<DynamicRJetInfo>().wt())/wt;
      }
      else  {
        // Case 4:
        // both i and j are previously merged jets.
        // Compute all cross terms between their constituents, then add the
        // stored internal radius information of each jet.
        mean_DR = 0;
        ms_DR  = 0;
        wt      = 0;
        vector<PseudoJet> jet_i_constituents = jet_i.constituents();
        vector<PseudoJet> jet_j_constituents = jet_j.constituents();
        vector<double> ptjs, ptis;
        ptjs.clear();
        ptis.clear();

        double wti = 0;
        double wtj = 0;

        // Cache constituent pt values of jet_i and total pt sum.
        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          pti  = jet_i_constituents[it].pt();
          wti = wti + pti;
          ptis.push_back(pti);
        }

        // Cache constituent pt values of jet_j and total pt sum.
        for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++){
          ptj  = jet_j_constituents[jt].pt();
          wtj = wtj + ptj;
          ptjs.push_back(ptj);
        }

        // Compute weighted cross terms between constituents of jet_i and jet_j.
        for (unsigned int it=0; it < jet_i_constituents.size(); it++){
          mean_DR1 = 0;
          ms_DR1  = 0;
          for (unsigned int jt=0; jt < jet_j_constituents.size(); jt++) {
            DR2ij = jet_i_constituents[it].plain_distance(jet_j_constituents[jt]);
            ms_DR1  += DR2ij*ptjs[jt];
            mean_DR1 += sqrt(DR2ij)*ptjs[jt];
          }
          mean_DR = mean_DR + mean_DR1*ptis[it];
          ms_DR  = ms_DR  + ms_DR1 *ptis[it];
        }

        // Add internal information already stored in the two merged jets.
        wt = wti*wtj + jet_i.user_info<DynamicRJetInfo>().wt() + jet_j.user_info<DynamicRJetInfo>().wt();
        mean_DR = (mean_DR + jet_i.user_info<DynamicRJetInfo>().mean_R()*jet_i.user_info<DynamicRJetInfo>().wt() + jet_j.user_info<DynamicRJetInfo>().mean_R()*jet_j.user_info<DynamicRJetInfo>().wt())/wt;
        ms_DR  = (ms_DR + jet_i.user_info<DynamicRJetInfo>().ms_R()*jet_i.user_info<DynamicRJetInfo>().wt() + jet_j.user_info<DynamicRJetInfo>().ms_R()*jet_j.user_info<DynamicRJetInfo>().wt())/wt;
      }

      // Recombine the two jets using FastJet's default PseudoJet addition.
      newjet = jet_i + jet_j;

      // Calculate the dynamic radius at this stage
      double Rd = calculate_Rd(mean_DR, ms_DR, wt);

      // Attach updated radius statistics to the merged jet.
      std::unique_ptr<DynamicRJetInfo> userInfo(new DynamicRJetInfo(mean_DR, ms_DR, wt));
      // Also attach updated dynamic radius at this state
      userInfo->set_Rd(Rd);
      newjet.set_user_info(userInfo.release());

      // Record the recombination in the cluster sequence and update the
      // nearest-neighbour structure with the new jet.
      clust_seq.plugin_record_ij_recombination(i, j, dij, newjet, k);
      nn.merge_jets(i, j, clust_seq.jets()[k], k);
    }
    else {
      // No partner jet: jet i is recombined with the beam.
      clust_seq.plugin_record_iB_recombination(i, dij);
      nn.remove_jet(i);
    }
    njets--;
  }
    
}

// Calculate the dynamic radius given mean and mean-square radius modifier
double DynamicR::calculate_Rd(double mean_DR, double ms_DR, double wt) const {
  double R0 = R();
  double Rd;
  Rd = ms_DR - mean_DR*mean_DR; // variance
  Rd = (Rd > 0.0) ? R0 + sqrt(Rd) : R0 - sqrt(-Rd); // final radius Rd = R0 + sd_i
  return Rd;
}

// print a banner for reference to the 3rd-party code
void DynamicR::_print_banner(ostream *ostr) const{
  // Print the banner only once.
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  // Banner with citation information.
  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the Dynamic Radius Jet Clustering Algorithm.            " << endl;
  (*ostr) << "# This implementation is built within the FastJet3 framework.             " << endl;
  (*ostr) << "# Please cite the FastJet3 references when using this software.           " << endl;
  (*ostr) << "#                                                                         " << endl;
  (*ostr) << "# If you use this algorithm, please cite:                                 " << endl;
  (*ostr) << "# B. Mukhopadhyaya, T. Samui, R. K. Singh                                 " << endl;
  (*ostr) << "# Dynamic Radius Jet Clustering Algorithm                                 " << endl;
  (*ostr) << "# JHEP 04 (2023) 019                                                      " << endl;
  (*ostr) << "# DOI: 10.1007/JHEP04(2023)019                                            " << endl;
  (*ostr) << "# arXiv:2301.13074                                                        " << endl;
  (*ostr) << "# Repository: https://github.com/tousiksamui/DynamicRJetAlgorithm /       " << endl;
  (*ostr) << "#                                                                         " << endl;
  (*ostr) << "# Additional studies using this algorithm are listed in the README file.  " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

} // namespace contrib
FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
