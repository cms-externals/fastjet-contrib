//  QCDAwarePlugin Package
//  Questions/Comments?  abuckley@cern.ch, cpollard@cern.ch
//
//  Copyright (c) 2014-2025
//  Andy Buckley, Chris Pollard, Donatas Zaripovas, Xinyuan Tan
//
// $Id: QCDAwarePlugin.hh 1533 2026-03-03 22:30:36Z buckley $
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

#ifndef __FASTJET_CONTRIB_QCDAWAREPLUGIN_HH__
#define __FASTJET_CONTRIB_QCDAWAREPLUGIN_HH__

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <queue>
#include <string>
#include <vector>
#include "DistanceMeasure.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {
    namespace QCDAwarePlugin {

        /**
         * \brief the PJDist struct: simple storage of distance
         *        information between two pseudojets
         */

        struct PJDist {
            double dist;
            int pj1;
            int pj2;
        };


        /**
         * \brief the QCDAwarePlugin class performs QCD-aware jet
         *        clustering given a particular DistanceMeasure
         *
         *  the QCDAwarePlugin class performs the QCD-aware jet
         *  clustering given a particular DistanceMeasure.
         */
        //------------------------------------------------------------------------
        // @todo
        // this class isn't really QCDAware; the distance measure is...
        class QCDAwarePlugin : public JetDefinition::Plugin {

            public:
                // User still owns the pointer to dm after using this
                // constructor.
                QCDAwarePlugin(const DistanceMeasure *dm, 
                        bool use_couplings=false, 
                        double coupling_power=-2.0, 
                        int running_coupling_order_alpha_s=0, 
                        double alpha_s = 0.1181, 
                        int running_coupling_order_alpha_em=0, 
                        double alpha_em = 1.0/128.92)
                    : _dm(dm),
                      _use_couplings(use_couplings),
                      _coupling_power(coupling_power),
                      _running_coupling_order_alpha_s{running_coupling_order_alpha_s},
                      _alpha_s{alpha_s},
                      _running_coupling_order_alpha_em{running_coupling_order_alpha_em},
                      _alpha_em(alpha_em) {}

                /// default destructor
                // we don't delete _dm here because it is owned by the
                // user.
                virtual ~QCDAwarePlugin() {}

                void run_clustering(fastjet::ClusterSequence& cs) const;

                std::string description() const;

                double R() const;

                // setters for run-time configuration
                QCDAwarePlugin& set_use_couplings(bool u) { _use_couplings = u; return *this; }
                QCDAwarePlugin& set_coupling_power(double p) { _coupling_power = p; return *this; }
                QCDAwarePlugin& set_alpha_s(double as) { _alpha_s = as; return *this; }
                QCDAwarePlugin& set_alpha_em(double a) { _alpha_em = a; return *this; }
                QCDAwarePlugin& set_running_coupling_order_alpha_s(int o) { _running_coupling_order_alpha_s = o; return *this; }
                QCDAwarePlugin& set_running_coupling_order_alpha_em(int o) { _running_coupling_order_alpha_em = o; return *this; }
                QCDAwarePlugin& set_enable_qcd(bool b) { _enable_qcd = b; return *this; }
                QCDAwarePlugin& set_enable_qed(bool b) { _enable_qed = b; return *this; }


            private:
                const DistanceMeasure *_dm;
                bool _use_couplings;
                double _coupling_power;
                int _running_coupling_order_alpha_s;
                double _alpha_s;
                int _running_coupling_order_alpha_em;
                double _alpha_em;

                bool _enable_qcd{true};
                bool _enable_qed{true};

                void insert_pj(ClusterSequence &cs,
                        std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                        unsigned int iJet,
                        std::vector<bool>& ismerged) const;

                void merge_iB(ClusterSequence &cs,
                        const PJDist& dist,
                        std::vector<bool>& ismerged) const;

                void merge_ij(ClusterSequence &cs,
                        std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                        const PJDist& dist,
                        std::vector<bool>& ismerged) const;

                // returns zero if p and q aren't allowed to combine.
                // returns the combined pid otherwise.
                std::pair<int, double> flavor_sum(const fastjet::PseudoJet& p, const fastjet::PseudoJet& q) const;
        };


    } // QCDAware
} // contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_QCDAWARE_HH__
