#ifndef __DISTANCEMEASURE_HH__
#define __DISTANCEMEASURE_HH__

#include <float.h>
#include <algorithm>
#include "fastjet/PseudoJet.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
    namespace QCDAwarePlugin {

        /**
         * the DistanceMeasure abstract base class calculates the
         * distance between two PseudoJets. the default implemented
         * DistanceMeasures include kt, C/A, and anti-kt
         */

        class DistanceMeasure {
            public:
                virtual double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const = 0;
                virtual double diB(const fastjet::PseudoJet& pji) const = 0;
                virtual double R() const = 0;
                virtual std::string algname() const = 0;

                virtual ~DistanceMeasure() {};
        };


        class KtMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                KtMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    double drbyR2 = pji.squared_distance(pjj) / (_R * _R);
                    return std::min(pji.perp2(), pjj.perp2()) * drbyR2;
                }

                inline double diB(const fastjet::PseudoJet& pji) const {
                    return pji.perp2();
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "kt";
                }
        };


        class AntiKtMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                AntiKtMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    double drbyR2 = pji.squared_distance(pjj) / (_R * _R);
                    return 1.0 / std::max(pji.perp2(), pjj.perp2()) * drbyR2;
                }

                inline double diB(const fastjet::PseudoJet& pji) const {
                    return 1.0 / pji.perp2();
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "anti-kt";
                }
        };


        class CAMeasure : public DistanceMeasure {
            private:
                double _R;

            public:
                CAMeasure(double R) : _R(R) {}

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    return pji.squared_distance(pjj) / (_R * _R);
                }


                inline double diB(const fastjet::PseudoJet& pji) const {
                    return 1.0;
                }

                virtual double R() const {
                    return _R;
                }

                std::string algname() const {
                    return "Cambridge-Aachen";
                }
        };

        class FlavourKtMeasure : public DistanceMeasure {
            private:
                double _R;
                double _alpha;

                /// Helper function to determine if the PseudoJet corresponds to a quark.
                /// This assumes the user_index was encoded as: user_index = 1000000 * object_id + pdg_id,
                /// where:
                ///   - object_id is an integer identifying the original source object (e.g. Rivet particle index, etc.)
                ///   - pdg_id is the PDG Monte Carlo particle ID code.
                inline bool isFlav(const fastjet::PseudoJet& p) const {
                    int pdgid = (p.user_index() > 0 ? p.user_index() % 1000000 : p.user_index() % (-1000000));
                    pdgid = std::abs(pdgid);
                    return (pdgid >= 1 && pdgid <= 6);
                }
                
            public:
                FlavourKtMeasure(double R, double alpha=2.0) : _R(R), _alpha(alpha) {
                    if (_alpha < 0.0 || _alpha > 2.0) {
                        throw std::invalid_argument(
                            "FlavourKtMeasure: alpha exponent must be in [0,2], got " 
                            + std::to_string(_alpha));
                    }
                }

                inline double dij(const fastjet::PseudoJet& pji, const fastjet::PseudoJet& pjj) const {
                    double drbyR2 = pji.squared_distance(pjj) / (_R * _R);
                    // see arXiv:hep-ph/0601139 for definitions
                    bool flav = isFlav(pji.perp2() < pjj.perp2() ? pji : pjj);  // get softer particles
                    double kt2 = flav ? std::pow(std::max(pji.perp(), pjj.perp()), _alpha) * std::pow(std::min(pji.perp(), pjj.perp()), 2.0 - _alpha)
                                      : std::min(pji.perp2(), pjj.perp2());
                    return kt2 * drbyR2;
                }

                inline double diB(const fastjet::PseudoJet& pji) const {
                    return pji.perp2();
                }

                virtual double R() const {
                    return _R; 
                }

                std::string algname() const {
                    return "Flavour-kt"; 
                }
        };

    } // QCDAware
} // contrib

FASTJET_END_NAMESPACE

#endif
