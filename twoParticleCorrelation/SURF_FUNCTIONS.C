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
//#include "TPythia8.h"

#include <iostream>
#include <vector>
#include <string>
//#include "Pythia8/Pythia.h"
#include "SURF_FUNCTIONS.h"

// Daughter Particle Calculations
Float_t calculatePt(Float_t px, Float_t py) {return sqrt(px*px + py*py);}

Float_t calculateTheta(Float_t px, Float_t py, Float_t pz) {
    Float_t r = sqrt(px*px + py*py + pz*pz);
    return acos(pz/r);
}

Float_t calculateEta(Float_t theta) {return -log(tan(theta/2));}

Float_t calculatePhi(Float_t px, Float_t py) {
    if (px > 0 && py > 0) {return atan(py/px);} // Quadrant 1
    else if (px < 0 && py > 0) {return atan(py/px) + TMath::Pi();} // Quadrant 2
    else if (px < 0 && py < 0) {return atan(py/px) - TMath::Pi();} // Quadrant 3
    else if (px > 0 && py < 0) {return atan(py/px);} // Quadrant 4
    else {return TMath::ATan2(py, px);} // Deals with the case if y or x is 0
}

// Calculate components of the four-momentum vector 
Float_t calculatePx(Float_t pT, Float_t phi) {return pT * cos(phi);}
Float_t calculatePy(Float_t pT, Float_t phi) {return pT * sin(phi);}
Float_t calculatePz(Float_t pT, Float_t eta) {return pT * sinh(eta);}
Float_t calculateEnergy(Float_t px, Float_t py, Float_t pz, Float_t m0) {return sqrt(m0*m0 + px*px + py*py + pz*pz);}

// Jet Calculations
Float_t calculateJetPt(std::vector<Float_t> pT, std::vector<Float_t> phi) {
    Float_t pxSum, pySum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pxSum += calculatePx(pT[i], phi[i]);
        pySum += calculatePy(pT[i], phi[i]);
    }
    return sqrt((pxSum*pxSum) + (pySum+pySum));
} 

Float_t calculateJetEta(std::vector<Float_t> pT, std::vector<Float_t> eta, std::vector<Float_t> phi) {
    Float_t pxSum, pySum, pzSum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pxSum += calculatePx(pT[i], phi[i]);
        pySum += calculatePy(pT[i], phi[i]);
        pzSum += calculatePz(pT[i], eta[i]);
    }
    return calculateEta(calculateTheta(pxSum, pySum, pzSum));
}

Float_t calculateJetPhi(std::vector<Float_t> pT, std::vector<Float_t> phi) {
    Float_t pxSum, pySum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pxSum += calculatePx(pT[i], phi[i]);
        pySum += calculatePy(pT[i], phi[i]);
    }
    if (pxSum > 0 && pySum > 0) {return atan(pySum/pxSum);} // Quadrant 1
    else if (pxSum < 0 && pySum > 0) {return atan(pySum/pxSum) + TMath::Pi();} // Quadrant 2
    else if (pxSum < 0 && pySum < 0) {return atan(pySum/pxSum) - TMath::Pi();} // Quadrant 3
    else if (pxSum > 0 && pySum < 0) {return atan(pySum/pxSum);} // Quadrant 4
    else {return TMath::ATan2(pySum, pxSum);} // Deals with that case if y or x is 0
    /* Will return a value from 0 to 360 degrees
    if (pxSum < 0 && pySum < 0) {return (atan(pySum/pxSum) + TMath::Pi());} // Quadrant 3 
    else if (pxSum < 0 && pySum > 0) {return (atan(pySum/pxSum) + TMath::Pi());} // Quadrant 2
    else if (pxSum > 0 && pySum < 0) {return (atan(pySum/pxSum) + (2*TMath::Pi()));} // Quadrant 4
    else {return atan(pySum/pxSum);}; // Quadrant 1 */
}

Float_t calculateJetPx(std::vector<Float_t> pT, std::vector<Float_t> phi) {
    Float_t pxSum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pxSum += calculatePx(pT[i], phi[i]);
    }
    return pxSum;
}

Float_t calculateJetPy(std::vector<Float_t> pT, std::vector<Float_t> phi) {
    Float_t pySum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pySum += calculatePy(pT[i], phi[i]);
    }
    return pySum;
}

Float_t calculateJetPz(std::vector<Float_t> pT, std::vector<Float_t> eta) {
    Float_t pzSum;
    for (Int_t i = 0; i < pT.size(); i++) {
        pzSum += calculatePz(pT[i], eta[i]);
    }
    return pzSum;
}

// Transformation Functions
// Returns a vector (or TLorentzVectors) with all the boosted particles
std::vector<TLorentzVector> boostFunction(TLorentzVector jet, std::vector<TLorentzVector> particles) {
    TVector3 jetBoostVector = jet.BoostVector(); // BoostVector() returns a TVector3 object
    std::vector<TLorentzVector> boostedParticles;
    for (Int_t i = 0; i < particles.size(); i++) {
        particles[i].Boost(-jetBoostVector);
        boostedParticles.push_back(particles[i]);
    }
    return boostedParticles;
}

// Returns an individual boosted particle as a TLorentzVector object
TLorentzVector boostParticle(TLorentzVector jet, TLorentzVector particle) {
    TVector3 jetBoostVector = jet.BoostVector(); // BoostVector() returns a TVector3 object
    particle.Boost(-jetBoostVector);
    return particle;

}

// Event display functions
// Calculates the distance between two particles in the eta-phi plane (AKA angular distance, tells how much obj moving in same dir.)
float dR(float eta1, float phi1, float eta2, float phi2) { 
  return TMath::Power(TMath::Power(TMath::ACos(TMath::Cos(phi1-phi2)), 2) + TMath::Power(eta1-eta2, 2), 0.5);
}

// Format the helix appearance
void formatHelix(THelix *h, float pz, float pt, int color, bool doWTA) {
  
  // Setting line attributes
  if (color==0) h->SetLineColor(kGreen);
  if (color==1) h->SetLineColor(kWhite);
  if (color==2) h->SetLineColor(kMagenta);
  if (color==3) h->SetLineColor(kRed);
  if (color==4) h->SetLineColor(kYellow);
  if (color==5) h->SetLineColor(kCyan-4);
  if (color==-1) h->SetLineColor(kBlack);
  h->SetLineWidth(1);
  
  // Adding bounds to the helix depending on pt and pz
  float rangeBound = 1.0;

  if (!doWTA && pt<1.0 && TMath::Abs(pz)<0.5) {rangeBound = 2 * TMath::Abs(pz);}
  if (!doWTA && pt<1.5 && TMath::Abs(pz)<0.5) {rangeBound = TMath::Abs(pz);}
  if (!doWTA && pt>1.5) {rangeBound = rangeBound*pz/pt*10;}
  if (!doWTA && pt>5) {rangeBound = 0.2;}

  h->SetRange(0, rangeBound);

  if (pz < 0) {h->SetRange(-rangeBound, 0);} 
  if (pz > 0) {h->SetRange(0, 1, kHelixZ);}
  else {h->SetRange(-1, 0, kHelixZ);}

  //h->SetRange(-2,2,kLabX);
  //h->SetRange(-2,2,kLabY);
}


