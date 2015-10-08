//  EnergyCorrelator Package
//  Questions/Comments?  larkoski@mit.edu gavin.salam@cern.ch jthaler@jthaler.net
//
//  Copyright (c) 2013
//  Andrew Larkoski, Gavin Salam, and Jesse Thaler
//
//  $Id: EnergyCorrelator.cc 759 2014-11-20 23:00:34Z larkoski $
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

#include "EnergyCorrelator.hh"
#include <sstream>
#include <limits>
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

double EnergyCorrelator::result(const PseudoJet& jet) const {

   // if jet does not have constituents, throw error
   if (!jet.has_constituents()) throw Error("EnergyCorrelator called on jet with no constituents.");

   // get N = 0 case out of the way
   if (_N == 0) return 1.0;

   // find constituents 
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   double answer = 0.0;

   // take care of N = 1 case.
   if (_N == 1) {
      for (unsigned int i = 0; i < particles.size(); i++) {
         answer += energy(particles[i]);
      }
      return answer;
   }

   double half_beta = _beta/2.0;

   // take care of N = 2 case.
   if (_N == 2) {
      for (unsigned int i = 0; i < particles.size(); i++) {
         for (unsigned int j = i + 1; j < particles.size(); j++) { //note offset by one so that angle is never called on identical pairs
            answer += energy(particles[i])
                      * energy(particles[j])
                      * pow(angleSquared(particles[i],particles[j]), half_beta);
         }
      }   
      return answer;
   }
   

   // if N > 4, then throw error
   if (_N > 5) {
      throw Error("EnergyCorrelator is only hard coded for N = 0,1,2,3,4,5");
   }

   
   // Now deal with N = 3,4,5.  Different options if storage array is used or not.  
   if (_strategy == storage_array) {
   
      // For N > 2, fill static storage array to save computation time.

      // Make energy storage
      std::vector<double> energyStore;
      energyStore.resize(particles.size());

      // Make angular storage
      std::vector< std::vector<double> > angleStore;
      angleStore.resize(particles.size());
      for (unsigned int i = 0; i < angleStore.size(); i++) {
         angleStore[i].resize(i);
      }
      
      // Fill storage with energy/angle information
      for (unsigned int i = 0; i < particles.size(); i++) {
         energyStore[i] = energy(particles[i]);
         for (unsigned int j = 0; j < i; j++) {
           angleStore[i][j] = pow(angleSquared(particles[i],particles[j]), half_beta);
         }
      }

      // now do recursion
      if (_N == 3) {
         for (unsigned int i = 2; i < particles.size(); i++) {
            for (unsigned int j = 1; j < i; j++) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (unsigned int k = 0; k < j; k++) {
                  answer += ans_ij
                            * energyStore[k]
                            * angleStore[i][k]
                            * angleStore[j][k];
               }
            }
         } 
      } else if (_N == 4) {
         for (unsigned int i = 3; i < particles.size(); i++) {
            for (unsigned int j = 2; j < i; j++) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (unsigned int k = 1; k < j; k++) {
                  double ans_ijk = ans_ij
                                 * energyStore[k]
                                 * angleStore[i][k]
                                 * angleStore[j][k];
                  for (unsigned int l = 0; l < k; l++) {
                     answer += ans_ijk
                                       * energyStore[l]
                                       * angleStore[i][l]
                                       * angleStore[j][l]
                                       * angleStore[k][l];
                  }
               }
            }
         } 
      } else if (_N == 5) {
         for (unsigned int i = 4; i < particles.size(); i++) {
            for (unsigned int j = 3; j < i; j++) {
               double ans_ij = energyStore[i]
                               * energyStore[j]
                               * angleStore[i][j];
               for (unsigned int k = 2; k < j; k++) {
                  double ans_ijk = ans_ij
                                 * energyStore[k]
                                 * angleStore[i][k]
                                 * angleStore[j][k];
                  for (unsigned int l = 1; l < k; l++) {
                     double ans_ijkl = ans_ijk
                                       * energyStore[l]
                                       * angleStore[i][l]
                                       * angleStore[j][l]
                                       * angleStore[k][l];
                     for (unsigned int m = 0; m < l; m++) {
                        answer += ans_ijkl
                                  * energyStore[m]
                                  * angleStore[i][m]
                                  * angleStore[j][m]
                                  * angleStore[k][m]
                                  * angleStore[l][m];
                     }                  
                  }
               }
            }
         } 
      
      } else {
         assert(_N <= 5);
      }

   } else if (_strategy == slow) {
      if (_N == 3) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energy(particles[i])
                               * energy(particles[j])
                               * pow(angleSquared(particles[i],particles[j]), half_beta);
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  answer += ans_ij
                            * energy(particles[k])
                            * pow(angleSquared(particles[i],particles[k]), half_beta)
                            * pow(angleSquared(particles[j],particles[k]), half_beta);
               }
            }
         }
      } else if (_N == 4) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energy(particles[i])
                               * energy(particles[j])
                               * pow(angleSquared(particles[i],particles[j]), half_beta);
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  double ans_ijk = ans_ij
                                   * energy(particles[k])
                                   * pow(angleSquared(particles[i],particles[k]), half_beta)
                                   * pow(angleSquared(particles[j],particles[k]), half_beta);
                  for (unsigned int l = k + 1; l < particles.size(); l++) {
                     answer += ans_ijk
                               * energy(particles[l])
                               * pow(angleSquared(particles[i],particles[l]), half_beta)
                               * pow(angleSquared(particles[j],particles[l]), half_beta)
                               * pow(angleSquared(particles[k],particles[l]), half_beta);
                  }
               }
            }
         }
      } else if (_N == 5) {
         for (unsigned int i = 0; i < particles.size(); i++) {
            for (unsigned int j = i + 1; j < particles.size(); j++) {
               double ans_ij = energy(particles[i])
                               * energy(particles[j])
                               * pow(angleSquared(particles[i],particles[j]), half_beta);
               for (unsigned int k = j + 1; k < particles.size(); k++) {
                  double ans_ijk = ans_ij
                                   * energy(particles[k])
                                   * pow(angleSquared(particles[i],particles[k]), half_beta)
                                   * pow(angleSquared(particles[j],particles[k]), half_beta);
                  for (unsigned int l = k + 1; l < particles.size(); l++) {
                     double ans_ijkl = ans_ijk
                               * energy(particles[l])
                               * pow(angleSquared(particles[i],particles[l]), half_beta)
                               * pow(angleSquared(particles[j],particles[l]), half_beta)
                               * pow(angleSquared(particles[k],particles[l]), half_beta);
                     for (unsigned int m = l + 1; m < particles.size(); m++) {
                        answer += ans_ijkl
                                  * energy(particles[m])
                                  * pow(angleSquared(particles[i],particles[m]), half_beta)
                                  * pow(angleSquared(particles[j],particles[m]), half_beta)
                                  * pow(angleSquared(particles[k],particles[m]), half_beta)
                                  * pow(angleSquared(particles[l],particles[m]), half_beta);
                     }
                  }
               }
            }
         }
      } else {
         assert(_N <= 5);
      } 
   } else {
      assert(_strategy == slow || _strategy == storage_array);
   }
   
   return answer;
}

double EnergyCorrelator::energy(const PseudoJet& jet) const {
   if (_measure == pt_R) {
      return jet.perp();
   }  else if (_measure == E_theta) {
      return jet.e();
   } else {
      assert(_measure==pt_R || _measure==E_theta);
      return std::numeric_limits<double>::quiet_NaN();
   }
}

double EnergyCorrelator::angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const {
   if (_measure == pt_R) {
      return jet1.squared_distance(jet2);
   } else if (_measure == E_theta) {
      // doesn't seem to be a fastjet built in for this
      double dot = jet1.px()*jet2.px() + jet1.py()*jet2.py() + jet1.pz()*jet2.pz();
      double norm1 = sqrt(jet1.px()*jet1.px() + jet1.py()*jet1.py() + jet1.pz()*jet1.pz());
      double norm2 = sqrt(jet2.px()*jet2.px() + jet2.py()*jet2.py() + jet2.pz()*jet2.pz());
      
      double costheta = dot/(norm1 * norm2);
      if (costheta > 1.0) costheta = 1.0; // Need to handle case of numerical overflow
      double theta = acos(costheta);
      return theta*theta;    
        
   } else {
      assert(_measure==pt_R || _measure==E_theta);
      return std::numeric_limits<double>::quiet_NaN();
   }
}

string EnergyCorrelator::description_parameters() const {
  ostringstream oss;
  oss << "N=" << _N << ", beta=" << _beta;

  if      (_measure == pt_R)    oss << ", pt_R measure";
  else if (_measure == E_theta) oss << ", E_theta measure";
  else throw Error("unrecognized measure");

  if      (_strategy == slow)   oss << " and 'slow' strategy";
  else if (_strategy == storage_array)   oss << " and 'storage_array' strategy";
  else throw Error("unrecognized strategy");

  return oss.str();
}

string EnergyCorrelator::description_no_N() const {
  ostringstream oss;
  oss << "beta=" << _beta;

  if      (_measure == pt_R)    oss << ", pt_R measure";
  else if (_measure == E_theta) oss << ", E_theta measure";
  else throw Error("unrecognized measure");

  if      (_strategy == slow)   oss << " and 'slow' strategy";
  else if (_strategy == storage_array)   oss << " and 'storage_array' strategy";
  else throw Error("unrecognized strategy");

  return oss.str();
}

string EnergyCorrelator::description() const {
  ostringstream oss;
  oss << "Energy Correlator ECF(N,beta) for ";
  oss << description_parameters();
  return oss.str();
}

string EnergyCorrelatorRatio::description() const {
  ostringstream oss;
  oss << "Energy Correlator ratio ECF(N+1,beta)/ECF(N,beta) for ";
  oss << EnergyCorrelator(_N,_beta,_measure,_strategy).description_parameters();
  return oss.str();
}

string EnergyCorrelatorDoubleRatio::description() const {
  ostringstream oss;
  oss << "Energy Correlator double ratio ECF(N-1,beta)ECF(N+1,beta)/ECF(N,beta)^2 for ";
  oss << EnergyCorrelator(_N,_beta,_measure,_strategy).description_parameters();
  return oss.str();
}

string EnergyCorrelatorD2::description() const {
  ostringstream oss;
  oss << "Energy Correlator observable D2 ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3 for ";
  oss << EnergyCorrelator(3,_beta,_measure,_strategy).description_no_N();
  return oss.str();
}

string EnergyCorrelatorC1::description() const {
  ostringstream oss;
  oss << "Energy Correlator observable C1 ECF(2,beta)/ECF(1,beta)^2 for ";
  oss << EnergyCorrelator(2,_beta,_measure,_strategy).description_no_N();
  return oss.str();
}

string EnergyCorrelatorC2::description() const {
  ostringstream oss;
  oss << "Energy Correlator observable C2 ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2 for ";
  oss << EnergyCorrelator(3,_beta,_measure,_strategy).description_no_N();
  return oss.str();
}





} // namespace contrib

FASTJET_END_NAMESPACE
