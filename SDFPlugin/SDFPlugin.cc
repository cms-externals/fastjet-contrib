// -*- C++ -*-
//
// Copyright (c) -,
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

#include "fastjet/contrib/SDFPlugin.hh"
#include "fastjet/contrib/FlavInfo.hh"
#include <istream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
namespace contrib{

//using fastjet::contrib::FlavInfo;
//using fastjet::contrib::FlavHistory;

void SDFlavourCalc::operator()(fastjet::PseudoJet& jet) {
  fastjet::ClusterSequence cs(jet.constituents(), p_plugin.get());
  cs.exclusive_jets_up_to(1);
  std::vector<fastjet::PseudoJet> sdjets = p_sd->operator()(cs.exclusive_jets_up_to(1));
  if(sdjets.size() != 1) jet.set_user_info(new fastjet::contrib::FlavHistory(fastjet::contrib::FlavInfo()));
  std::vector<fastjet::PseudoJet> sdj = sdjets[0].constituents();
  fastjet::contrib::FlavInfo flavj = fastjet::contrib::FlavHistory::current_flavour_of(sdj[0]);
  for(size_t i=1; i < sdj.size(); i++){
    flavj = flavj+fastjet::contrib::FlavHistory::current_flavour_of(sdj[i]);
  }
  // if(jet.constituents().size() != sdj.size()) {
  //   std::cout<<jet.constituents().size()<<" "<<sdj.size()<<"\n";
  //   for(auto p: jet.constituents()) std::cout<<FlavHistory::current_flavour_of(p).description()<<"\n";
  //   std::cout<<"---------\n";
  //   for(auto p: sdj) std::cout<<FlavHistory::current_flavour_of(p).description()<<"\n";
  //   std::cout<<"**********\n";
  // }
  jet.set_user_info(new fastjet::contrib::FlavHistory(flavj));
}

void SDFlavourCalc::operator()(std::vector<fastjet::PseudoJet>& jets) {
  for(fastjet::PseudoJet& jet: jets) {
    operator()(jet);
  }
}

} // namespace contrib

FASTJET_END_NAMESPACE
