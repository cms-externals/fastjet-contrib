// ConstituentSubtractor package                                                                                                                       
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz                   
                                                                                                                                                      
// Copyright (c) 2014-, Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner                                                                   
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
//                                                                                                                                                     // You should have received a copy of the GNU General Public License                                                                                   
// along with this code. If not, see <http://www.gnu.org/licenses/>.                                                                                   
//----------------------------------------------------------------------                                                                               


#include "RescalingClasses.hh"
 
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh                                                                                    


namespace contrib{




  // BackgroundRescalingYPhi
  BackgroundRescalingYPhi::BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2){
    _v2=v2;
    _v3=v3;
    _v4=v4;
    _psi=psi;
    _a1=a1;
    _sigma1=sigma1;
    _a2=a2;
    _sigma2=sigma2;
    _use_rap=true;
    _use_phi=true;
  }

  void BackgroundRescalingYPhi::use_rap_term(bool use_rap){
    _use_rap=use_rap;
  }

  void BackgroundRescalingYPhi::use_phi_term(bool use_phi){
    _use_phi=use_phi;
  }

  double BackgroundRescalingYPhi::result(const PseudoJet & particle) const {
    double phi_term=1;
    if (_use_phi){
      double phi=particle.phi();
      phi_term=1 + 2*_v2*_v2*cos(2*(phi-_psi)) + 2*_v3*_v3*cos(3*(phi-_psi)) +  2*_v4*_v4*cos(4*(phi-_psi));
    }
    double rap_term=1;
    if (_use_rap){
      double y=particle.rap();
      rap_term=_a1*exp(-y*y/(2*_sigma1*_sigma1)) + _a2*exp(-y*y/(2*_sigma2*_sigma2));
    }

    return phi_term*rap_term;
  }

} // namespace contrib                                                                                                                                  


FASTJET_END_NAMESPACE



