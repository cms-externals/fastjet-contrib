#ifndef SDFLAV_DEF
#define SDFLAV_DEF
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/Recluster.hh"
#include <memory>

FASTJET_BEGIN_NAMESPACE // defined in fastjet/internal/base.hh
namespace contrib{

/**
 * @class SDFlavourCalc
 * @brief A utility class applying SoftDrop algorithm with JADE reclustering
 * as required by the SDFPlugin algorithm.
 */


class SDFlavourCalc {
public:

  /** Constructor for SDFlavourCalc
   *
   * Initializes the SoftDrop algorithm with specified parameters for beta, zcut and R.
   * The constructor also sets up a reclusterign scheme using the JADE plugin.
   *
   * @param beta SoftDrop beta parameter, controlling the angular exponent. Must be strictly positive. Default is 2.
   * @param zcut SoftDrop zcut parameter, defining the energy fraction threshold. Default is 0.1.
   * @param R The jet rdius parameter. Default is 0.4.
   */
  SDFlavourCalc(const double beta = 2,
                const double zcut = 0.1,
                const double R = 0.4) : p_sd(new fastjet::contrib::SoftDrop(beta, zcut,
                                                                            fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure::scalar_z,
                                                                            R,
                                                                            std::numeric_limits<double>::infinity(),
                                                                            fastjet::contrib::RecursiveSymmetryCutBase::RecursionChoice::larger_pt,
                                                                            0)),
                                        p_plugin(new fastjet::JadePlugin()),
                                        p_recluster(new fastjet::Recluster(fastjet::JetDefinition(p_plugin.get()))){
    p_sd->set_reclustering(true,p_recluster.get());
  }

  /// Unique pointer to the SoftDrop instance
  std::unique_ptr<fastjet::contrib::SoftDrop> p_sd;
  /// Unique pointer to the JADE plugin for reclustering
  std::unique_ptr<fastjet::JetDefinition::Plugin> p_plugin;
  /// Unique pointer to the Recluster instance
  std::unique_ptr<fastjet::Recluster> p_recluster;


  /**
   * @brief Applies the SoftDrop procedure to a single jet.
   *
   * @param jets A reference to a PseudoJet to process.
   */
  void operator()(fastjet::PseudoJet& jet);

  /**
   * @brief Applies the SoftDrop procedure to a vector of jets.
   *
   * @param jets A vector of PseudoJets to process.
   */
  void operator()(std::vector<fastjet::PseudoJet>& jets);

};

} // namespace contrib 

FASTJET_END_NAMESPACE

#endif
