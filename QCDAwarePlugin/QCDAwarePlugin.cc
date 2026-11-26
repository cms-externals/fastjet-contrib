//  QCDAwarePlugin Package
//  Questions/Comments?  abuckley@cern.ch, cpollard@cern.ch
//
//  Copyright (c) 2014-2025
//  Andy Buckley, Chris Pollard, Donatas Zaripovas, Xinyuan Tan
//
// $Id: QCDAwarePlugin.cc 1533 2026-03-03 22:30:36Z buckley $
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

#include "QCDAwarePlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace fastjet;


namespace contrib {
    namespace QCDAwarePlugin {

        /// @name QCD and EM couplings and colour factors
        /// @{
        static constexpr double C_F       = 4.0/3.0;
        static constexpr double C_A       = 3.0;
        static constexpr double T_R       = 1.0/2.0;
        static constexpr double m_Z       = 91.1876;
        /// @}

        /// Define an unbound pseudojet-distance comparator
        bool operator > (const PJDist& pj1, const PJDist& pj2) {
            return pj1.dist > pj2.dist;
        }

        /// Consts for the running couplings
        ///
        /// @todo Make overrideable? Surely too fine a detail
        static const map<int,double> quark_masses = {
            {1, 0.00216},  // up
            {2, 0.00470},  // down
            {3, 0.0935},   // strange
            {4, 1.273},    // charm
            {5, 4.183},    // bottom
            {6, 172.5}     // top
        };

        // helper functions in anon namespace for inaccessibility
        namespace {

            /// @brief Helper function to extract the PDG ID from the user_index of a PseudoJet.
            ///
            /// Assumes the user_index was encoded as: user_index = 1000000 * object_id + pdg_id,
            /// where:
            ///   - object_id is an integer identifying the original source object (e.g. Rivet particle index, etc.)
            ///   - pdg_id is the particle’s PDG Monte Carlo ID code.
            ///
            /// The object_id can come from any external framework or internal labeling strategy.
            /// This approach allows multiplexing of all conceivable SM clustering flavours with
            /// an index to re-identify the constituents with their originals in the calling code.
            inline int getPID(const fastjet::PseudoJet& p) {
                return (p.user_index() > 0 ? p.user_index() % 1000000 : p.user_index() % (-1000000));
            }

            /// Check if the particle is a quark or antiquark
            inline bool isQuark(const fastjet::PseudoJet& p) {
                return p != 0 && abs(getPID(p)) <= 6;
            }

            /// Check if the particle is a gluon
            inline bool isGluon(const fastjet::PseudoJet& p) {
                return getPID(p) == 21;
            }

            /// Check if the particle is a photon
            inline bool isPhoton(const fastjet::PseudoJet& p) {
                return getPID(p) == 22;
            }

            /// Check if the particle is a charged lepton
            inline bool isLepton(const fastjet::PseudoJet& p) {
                int abspid = abs(getPID(p));
                return abspid == 11 || abspid == 13 || abspid == 15;
            }

            /// Number of active flavours at scale Q^2, for alpha_s running
            inline int numFlavors(double Q2) {
                int nf = 0;
                for (const auto& kv : quark_masses) {
                    double m = kv.second;
                    if (Q2 > m*m) ++nf;
                }
                return nf;
            }
            

            /// 5-flavour Lambda calculation from given alpha_s(MZ)
            inline double lambda5(double alpha_s) {
                double beta0 = (33.0 - 2.0*5) / (12.0 * M_PI);
                return m_Z * exp(-1.0 / (2.0 * beta0 * alpha_s));
            }

            /// 4-flavour Lambda calculation from given alpha_s(MZ)
            ///
            /// @brief Take alpha_s(MZ) as an argument, and chain 5->4 internally
            inline double lambda4(double lambda5) {
                double mb = quark_masses.at(5);
                double beta0_5 = (33.0 - 2.0*5) / (12.0 * M_PI);
                double beta0_4 = (33.0 - 2.0*4) / (12.0 * M_PI);
                double L5 = log(mb*mb / (lambda5*lambda5));
                double L4 = (beta0_5 / beta0_4) * L5;
                return mb * exp(-L4/2.0);
            }

            /// Generic Lambda calculation
            inline double lambda(int nf, double alpha_s) {
                double val_lambda5 = lambda5(alpha_s);
                double val_lambda4 = lambda4(val_lambda5);
                return (nf > 4) ? val_lambda5 : val_lambda4;
            }

            /// Running coupling
            inline double alphaS(double Q2, int order, double alpha_s) {
                if (order == 0) return alpha_s;
                if (order == 1) {
                    int nf = numFlavors(Q2);
                    double Lambda = lambda(nf, alpha_s);
                    double beta0 = (33.0 - 2.0*nf) / (12.0 * M_PI);
                    double L = log(Q2 / (Lambda*Lambda));
                    return 1.0 / (beta0 * L);
                }
                else {
                    throw std::invalid_argument("alphaS: only orders 0 (fixed) and 1 (one-loop) are supported");
                }
            }

            inline double alphaEM(double Q2, int order, double alpha_em) {
                if (order == 0) return alpha_em;
                else if (order == 1) {
                    double z_f = 60/9;
                    double m_Z2 = m_Z * m_Z;
                    return alpha_em / (1 - (alpha_em * z_f / (3 * M_PI)) * log(Q2/m_Z2));
                }
                else {
                    throw std::invalid_argument("alphaEM: only orders 0 (fixed) and 1 (one-loop) are supported");
                }
            }

        }


        void QCDAwarePlugin::insert_pj(ClusterSequence &cs,
                priority_queue<PJDist, vector<PJDist>, greater<PJDist> >& pjds,
                unsigned int iJet,
                vector<bool>& ismerged) const {

            const PseudoJet& ijet = cs.jets()[iJet];

            for (unsigned int jJet = 0; jJet < iJet; jJet++) {
                // don't calculate distances for already-merged pjs
                if (ismerged[jJet])
                    continue;

                const PseudoJet& jjet = cs.jets()[jJet];

                PJDist pjd;
                pjd.pj1 = iJet;
                pjd.pj2 = jJet;

                pair<int,double> res = flavor_sum(ijet, jjet);
                int c = res.first;
                double factor = res.second;

                if (c == 0) {
                    pjd.dist = DBL_MAX;
                }
                   
                else {
                    double deltaR = ijet.delta_R(jjet);
                    if (deltaR > _dm->R()) {
                        pjd.dist = DBL_MAX;
                    }
                    else{
                        double couplings = _use_couplings ? (pow(factor, _coupling_power)) : 1.0;
                        pjd.dist = couplings * _dm->dij(ijet, jjet);
                    }
                }

                pjds.push(pjd);

            }

            // calculate the beam distance
            PJDist pjd;
            pjd.pj1 = iJet;
            pjd.pj2 = -1;
            double diB0 = _dm->diB(ijet);
            double bf = 1.0;
            if (_use_couplings) {
                double alphas = alphaS(ijet.perp2(), _running_coupling_order_alpha_s, _alpha_s);
                double alpha = alphaEM(ijet.perp2(), _running_coupling_order_alpha_em, _alpha_em);
                if (isQuark(ijet) || isGluon(ijet))        bf = pow(alphas, _coupling_power);
                else if (isLepton(ijet) || isPhoton(ijet)) bf = pow(alpha, _coupling_power);
            }
            pjd.dist = bf * diB0;
            pjds.push(pjd);

            ismerged.push_back(false);

            return;
        }


        void QCDAwarePlugin::merge_iB(ClusterSequence &cs,
                const PJDist& pjd,
                std::vector<bool>& ismerged) const {

            cs.plugin_record_iB_recombination(pjd.pj1, pjd.dist);

            ismerged[pjd.pj1] = true;

            return;
        }

        void QCDAwarePlugin::merge_ij(ClusterSequence &cs,
                std::priority_queue<PJDist, std::vector<PJDist>, std::greater<PJDist> >& pjds,
                const PJDist& pjd,
                std::vector<bool>& ismerged) const {

            // mark both old pjs as merged
            ismerged[pjd.pj1] = true;
            ismerged[pjd.pj2] = true;

            const PseudoJet& pj1 = cs.jets()[pjd.pj1];
            const PseudoJet& pj2 = cs.jets()[pjd.pj2];
            PseudoJet pj3 = pj1 + pj2;

            int c = flavor_sum(pj1, pj2).first;
            if (c == 0) {
                cout << "ERROR: attempting to merge pseudojets with pdgids "
                    << getPID(pj1) << " and " << getPID(pj2)
                    << ", which is not allowed. This will probably break." << endl;
                pj3.set_user_index(-999);
            } else
                pj3.set_user_index(c);

            int newidx;
            cs.plugin_record_ij_recombination(pjd.pj1, pjd.pj2, pjd.dist, pj3, newidx);

            insert_pj(cs, pjds, newidx, ismerged);

            return;
        }

        /// Compute the resulting flavour of a clustering
        std::pair<int,double> QCDAwarePlugin::flavor_sum(const fastjet::PseudoJet& p, const fastjet::PseudoJet& q) const {

            // Set colour factor and couplings
            double f;

            // use the geometric-mean for Q^2 = p.perp() * q.perp()
            double Q2 = p.perp() * q.perp();
            double alphas = alphaS(Q2, _running_coupling_order_alpha_s, _alpha_s);
            double alpha = alphaEM(Q2, _running_coupling_order_alpha_em, _alpha_em);

            if (_enable_qcd) {
                // a quark can cluster with a gluon.
                if ( isQuark(p) && isGluon(q) ) {
                    f =  C_F * alphas;
                    return std::make_pair(p.user_index(), f);
                }
                else if ( isGluon(p) && isQuark(q)) {
                    f = C_F * alphas;
                    return std::make_pair(q.user_index(), f);
                }

                // gluons can cluster.
                else if (isGluon(p) && isGluon(q)) {
                    f = C_A * alphas;
                    return std::make_pair(21, f);
                }

                // same-flavor quarks and anti-quarks can cluster.
                else if (isQuark(p) && isQuark(q) && (getPID(p) + getPID(q) == 0)) {
                    f = T_R * alphas;
                    return std::make_pair(21, f);
                }
            }
            
            if (_enable_qed) {
                // a quark can cluster with a photon.
                if ( isQuark(p) && isPhoton(q) ) {
                    f = alpha;
                    return std::make_pair(p.user_index(), f);
                }
                else if ( isPhoton(p) && isQuark(q) ) {
                    f = alpha;
                    return std::make_pair(q.user_index(), f);
                }

                // leptons and photons can cluster.
                else if (isLepton(p) && isPhoton(q)) {
                    f = alpha;
                    return std::make_pair(p.user_index(), f);
                }
                else if (isPhoton(p) && isLepton(q)) {
                    f = alpha;
                    return std::make_pair(q.user_index(), f);
                }

                // lepton anti-lepton pairs can cluster.
                else if (isLepton(p) && isLepton(q) && (getPID(p) + getPID(q) == 0)) {
                    f = alpha;
                    return std::make_pair(22, f);
                }
            }

            return std::make_pair(0, 0.0);
        }



        void QCDAwarePlugin::run_clustering(ClusterSequence& cs) const {

            vector<bool> ismerged;

            priority_queue<PJDist, vector<PJDist>, greater<PJDist> > pjds;
            for (unsigned int iJet = 0; iJet < cs.jets().size(); iJet++)
                insert_pj(cs, pjds, iJet, ismerged);

            while (!pjds.empty()) {
                PJDist pjd = pjds.top();
                pjds.pop();

                // check for already merged pj1
                if (ismerged[pjd.pj1])
                    continue;

                // check for the beam
                if (pjd.pj2 < 0) {
                    merge_iB(cs, pjd, ismerged);
                    continue;
                }

                // check for already merged pj2
                if (ismerged[pjd.pj2])
                    continue;

                merge_ij(cs, pjds, pjd, ismerged);
            }

            return;
        }

        string QCDAwarePlugin::description() const {
            stringstream ss;
            ss << "QCDAwarePlugin jet algorithm with R = " << R() <<
                " and " << _dm->algname() << " distance measure";
            return ss.str();
        }

        double QCDAwarePlugin::R() const {
            return _dm->R();
        }



    } // QCDAware
} // contrib


FASTJET_END_NAMESPACE
