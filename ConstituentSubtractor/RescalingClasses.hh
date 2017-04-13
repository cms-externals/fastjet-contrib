// ConstituentSubtractor package                                                                                                                        
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz                   
//                                                                                                                                                     
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
//                                                                                                                                                     
// You should have received a copy of the GNU General Public License                                                                                   
// along with this code. If not, see <http://www.gnu.org/licenses/>.                                                                                   
//----------------------------------------------------------------------   

#ifndef __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__
#define __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__



#include <fastjet/FunctionOfPseudoJet.hh>
#include <iostream>
//#include "TH1.h" 

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh                                                                                     

namespace contrib{


  template<class T>
  class BackgroundRescalingYFromRoot : public FunctionOfPseudoJet<double> {
  public:
    /// construct a background rescaling function using ROOT TH1 histogram bin contents
    BackgroundRescalingYFromRoot(): _hist(0) {}
    BackgroundRescalingYFromRoot(T* hist=0) {_hist = hist;}

    // return the rescaling factor associated with this jet  
    virtual double result(const PseudoJet & particle) const {
      if (!_hist){
	std::cout << "histogram for rescaling not defined!!!" << std::endl;
	throw;
      }
      double y = particle.rap();
      int bin=_hist->FindBin(y);
      return _hist->GetBinContent(bin);
    }

  private:
    T* _hist;
  };



class BackgroundRescalingYPhi : public FunctionOfPseudoJet<double> {
public:
  ///  Construct background rescaling function in rapidity and azimuth using this parameterization:

  ///  f(y,phi) = phi_term(phi) * rap_term(y)
  ///  where

  ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))
  ///  with four parameters v2, v3, v4, and psi.

  ///  rap_term(y) = a1*exp(-pow(y,2)/(2*sigma1^2)) + a2*exp(-pow(y,2)/(2*sigma2^2))
  ///  with four parameters sigma1, sigma2, a1, and a2. 

  ///  This function is used to rescale the background which is subtracted such that one can correctly account
  ///  for the modulation of the UE due to rapidity dependence of the particle production
  ///  and/or due to the modulation in the azimuthal angle which is characteristic for heavy ion collisions.
  ///  The overall normalization of function f is arbitrary since it divides out in the calculation of position dependent rho (background is first demodulated to obtain unbiased position independent rho, and then it is modulated to obtain position dependent rho, see fastjet classes GridMedianBackgroundEstimator and JetMedianBackgroundEstimator for detailed calculation).

  BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
  BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2);

  void use_rap_term(bool use_rap);
  void use_phi_term(bool use_phi);

  /// return the rescaling factor associated with this jet  
  virtual double result(const PseudoJet & jet) const;
private:
  double _v2, _v3, _v4, _psi, _a1, _sigma1, _a2, _sigma2;
  bool _use_rap, _use_phi;
};









} // namespace contrib                                                                                                                                  

FASTJET_END_NAMESPACE


#endif   //__FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__
