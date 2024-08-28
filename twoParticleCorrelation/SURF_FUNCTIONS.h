#ifndef SURF_FUNCTIONS_H
#define SURF_FUNCTIONS_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THelix.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <string>
//#include <Pythia8/Pythia.h> 

// Daughter Particle Calculations
Float_t calculatePt(Float_t x, Float_t y);
Float_t calculateTheta(Float_t x, Float_t y, Float_t z);
Float_t calculateEta(Float_t theta);

// Calculate components of the four-vector 
Float_t calculatePx(Float_t pT, Float_t phi);
Float_t calculatePy(Float_t pT, Float_t phi);
Float_t calculatePz(Float_t pT, Float_t eta);
Float_t calculateEnergy(Float_t px, Float_t py, Float_t pz, Float_t m0);

//  Jet Calculations
Float_t calculateJetPt(std::vector<Float_t> pT, std::vector<Float_t> phi);
Float_t calculateJetEta(std::vector<Float_t> pT, std::vector<Float_t> eta, std::vector<Float_t> phi);
Float_t calculateJetPhi(std::vector<Float_t> pT, std::vector<Float_t> phi);

// Transformation Functions
std::vector<TLorentzVector> boostFunction(TLorentzVector jet, std::vector<TLorentzVector> particles);

// Event display functions
float dR(float eta1, float phi1, float eta2, float phi2);
void formatHelix(THelix *h, float pz, float pt, int color = 0, bool doWTA = false);

#endif