#include <iostream>
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/SDFPlugin.hh" 
#include "fastjet/contrib/FlavInfo.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(int iargc, char **argv){

  // give user control over printout (mainly relevant for make check)
  // usage: "./example-SDF [nevmax [njetmax]] < ../data/pythia8_Zq_vshort.dat"
  unsigned int nevmax = 2;
  unsigned int njetmax = 1;
  if (iargc > 1) nevmax  = stoi(argv[1]);
  if (iargc > 2) njetmax = stoi(argv[2]);

  // print banner for FastJet at the start, so it doesn't mix
  // into the other output
  ClusterSequence::print_banner(); 

  // we start with a base jet definition (should be either
  // antikt_algorithm or cambridge_algorithm, or their e+e- variants)
  JetDefinition base_jet_def(antikt_algorithm, 0.4);
  // enable it to track flavours (default is net flavour)
  FlavRecombiner flav_recombiner;
  base_jet_def.set_recombiner(&flav_recombiner);

  SDFlavourCalc sdFlavCalc;

  // loop over some number of events
  int n_events = 10;
  for (int iev = 0; iev < n_events && iev < nevmax; iev++) {
    cout << "\n#---------------------------------------------------------------\n";
    cout << "# read event " << iev;

    // read in input particles: see that routine for info 
    // on how to set up the PseudoJets with flavour information
    vector<PseudoJet> event;
    read_event(event);
    cout<< " with " << event.size() << " particles" << endl;
    // run the jet clustering with the base jet definition and the
    // IFNPlugin-based jet definition
    vector<PseudoJet> base_jets = base_jet_def(event);
    cout<<"now do the sd..."<<endl;
    vector<PseudoJet> sdflav_jets = base_jet_def(event);
    sdFlavCalc(sdflav_jets);

    // // ----------------------------------------------------
    // // loop over the two leading jets and print out their properties
    for (unsigned int ijet = 0; ijet < base_jets.size() && ijet < njetmax; ijet++) {
      // first print out the original anti-kt jets and the IFN jets
      const auto & base_jet = base_jets[ijet];
      const auto & sdflav_jet  = sdflav_jets [ijet];
      cout << endl;
      cout << "base jet " << ijet << ": ";
      cout << "pt=" << base_jet.pt() << " rap=" << base_jet.rap() << " phi=" << base_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(base_jet).description() << endl;
      cout << "SD flav jet  " << ijet << ": ";
      cout << "pt=" << sdflav_jet.pt() << " rap=" << sdflav_jet.rap() << " phi=" << sdflav_jet.phi();
      cout << ", flav = " << FlavHistory::current_flavour_of(sdflav_jet).description() << endl;

    }
  }

  return 0;
}

// read in input particles and set up PseudoJets with flavour information
void read_event(vector<PseudoJet> &event){  
    // read in the input particles and their PDG IDs
    string line;
    double px, py, pz, E;
    int    pdg_id;
    event.resize(0);
    while(getline(cin,line)) {
      if(line[0] == '#') continue;

      istringstream iss(line);
      iss >> px >> py >> pz >> E >> pdg_id;
      // create a fastjet::PseudoJet with these components and put it onto
      // back of the input_particles vector
      PseudoJet p(px,py,pz,E);

      // assign information about flavour (will be deleted automatically)
      p.set_user_info(new FlavHistory(pdg_id));
      event.push_back(p);

      if (cin.peek() == '\n' || cin.peek() == EOF) {
        getline(cin,line);
        break;
      }
    }
}
