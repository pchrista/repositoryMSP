// Authors: Ioannis Sidiras, Paul Veen, Rik Spijkers

// C++ libraries
#include <iostream>
#include <vector>
//ROOT libraries
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TAttMarker.h"
#include "TLegend.h"

#define PI 3.14159265
using namespace std;
using namespace TMath;

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	// Returns delta phi in range -pi/2 3pi/2
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
}
	
bool IsStrange(Int_t particlepdg) {
	// Checks if given pdg code is that of a (anti)strange hadron, returns True/False
	// This also find ssbar mesons such as the phi
	Int_t pdg = std::abs(particlepdg);
	pdg /= 10; // get rid of the last digit, is not important for quark content
	if (pdg % 10 == 3) return true; // 3rd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 2nd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 1st quark
	return false;
}
	
void ssbar_SpectraAnalysis(TString inputFiles, TString outfilename) {
  // Define the TChain
  TChain *ch1 = new TChain("tree");
  ch1->Add(inputFiles);
  
  TFile *output = new TFile(outfilename,"RECREATE");
  
  // Now we define vectors that carry the information at event level.
  vector<Int_t>* vID = 0;
  vector<Double_t>* vPt = 0;
  vector<Double_t>* vPhi = 0;
  vector<Double_t>* vStatus = 0;
  vector<Double_t>* vEta = 0;
  vector<Double_t>* vMother1 = 0;
  vector<Double_t>* vMotherID = 0;

  // Setting up chain branch addresses to the vectors defined above
  ch1->SetBranchAddress("ID",&vID);
  ch1->SetBranchAddress("PT",&vPt);
  ch1->SetBranchAddress("PHI",&vPhi);
  ch1->SetBranchAddress("ETA",&vEta);
  ch1->SetBranchAddress("STATUS",&vStatus);
  ch1->SetBranchAddress("MOTHER",&vMother1);
  ch1->SetBranchAddress("MOTHERID",&vMotherID);
  
  // Definition of variables
  Int_t aID, pID;
  Double_t pPt, pPhi, pStatus, pEta, pMotherID; //For trigger
  Double_t aPt, aPhi, aStatus, aEta, aMotherID; //For associate
  int sign = 1;
  int nTrigger = 0;
  
  // Each vector is an event number of events analyzed is the total number of vectors
  int nEvents = ch1->GetEntries();
  cout<<"The number of events for this analysis is: "<<nEvents<<endl;
	
  //===========================================================//
  // Definition of produced histograms
  // 1D histograms
  TH1D *fHistPDG = new TH1D("fHistPt\DG", "PDG code; pdg code; Counts", 11999, -5999.5, 5999.5);
  TH1D *fHistPt = new TH1D("fHistPt", "Transverse Momentum; p_{T} GeV/c; Counts", 100, 0, 50);
  //===========================================================//
  
  //===========================================================//
  // Event Loop
  for(int iEvent = 0; iEvent < nEvents; iEvent++) {
    ch1->GetEntry(iEvent);
    int nparticles = vID->size();
    for(int ipart = 0; ipart < nparticles; ipart++){
      pID = (*vID)[ipart];
      pPhi = (*vPhi)[ipart];
      pPt = (*vPt)[ipart];
      pStatus = (*vStatus)[ipart];
      pEta =(*vEta)[ipart];
      pMotherID = (*vMotherID)[ipart];
      
      fHistPDG->Fill(pID);
      fHistPt->Fill(pPt);
    } // Particle loop
  } // Event loop

  output->Write();
  output->Close();
  
  cout << "File: " << outfilename << " has been created!" << endl;
}

