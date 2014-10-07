// $Id: example.hh 587 2014-04-06 13:44:24Z berta $
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

#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "ConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/ConstituentSubtractor.hh


// \class JetWidth
// computes the jet width used as one shape in the example
class JetWidth : public fastjet::FunctionOfPseudoJet<double>{
public:

  // action of the function
  double result(const fastjet::PseudoJet &jet) const{
    if (!jet.has_constituents()){
      return -0.1;
    }    
    double width = 1e-6;
    double ptSum = 0;
    
    if (jet.constituents().size() < 2) return width;

    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    for (unsigned int iconst=0; iconst<constituents.size(); iconst++){
      fastjet::PseudoJet cons=constituents.at(iconst);
      double dR = sqrt(cons.squared_distance(jet));
      double pt = cons.pt();
      width += dR * pt;
      ptSum += pt;
    }

    return width/ptSum;
  }

};
