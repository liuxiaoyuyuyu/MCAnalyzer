#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TAttLine.h"
#include "THelix.h"
#include "TView.h"
#include "TRandom.h"
#include "TAttPad.h"
#include "TMath.h"
#include "TVector3.h"
#include "TView3D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include <iostream>
#include <string>
#include <vector>

#include "SURF_FUNCTIONS.h"
#include "coordinateTools.h"
R__LOAD_LIBRARY(SURF_FUNCTIONS_C);

// Global variables
std::string title = "#mu^{+}#mu^{#minus} (10,000 GeV)";
std::string fullPathDir = "/Users/xl155/Dropbox/FutureColliders/mupmum_10000GeV.root";
//std::string fullPathDir = "/Users/katherinexie/SURF2024/futureCollidersData/mupmum_10000GeV.root";
Float_t etaCutOff = 2.5;
Int_t multCutOff = 194;
TFile *f = new TFile(fullPathDir.c_str(), "read"); // Opening file

// Creating TTreeReader object and linking branches
TTreeReader* reader = new TTreeReader("trackTree");

// Setup branches for other particles
TTreeReaderValue<std::vector<Float_t>> pxBranch(*reader, "px");
TTreeReaderValue<std::vector<Float_t>> pyBranch(*reader, "py");
TTreeReaderValue<std::vector<Float_t>> pzBranch(*reader, "pz");
TTreeReaderValue<std::vector<Int_t>> chgBranch(*reader, "chg");
TTreeReaderValue<std::vector<Int_t>> pidBranch(*reader, "pid");

// Function that returns the distribution of the simple signal function 
TH2F createSimpleSignalDist(std::vector<Int_t> multiplicityVector, std::vector<Int_t> weightVector) {
    
    std::string signalTitle = "Normalized Signal Distribution for " + title;

    // Histogram for the signal distribution
    //TH2F hSignal("hSignal", signalTitle.c_str(), 100, -5, 5, 100, -TMath::Pi(), TMath::Pi());
    const float         EtaBW = 0.3; 
    const float         PhiBW = TMath::Pi()/16;
    TH2F hSignal("hSignal", signalTitle.c_str(),41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);

    // ***** EVENT LOOP *****
    reader->Restart(); // Restarting event loop
    Int_t numSelectedEvents = 0;
    Int_t jetCounter = 0; 
    Int_t numTrigg = 0;

    while (reader->Next()) {

        // Check to see if the event is in the multiplicity bin
        Int_t eventIndex = reader->GetCurrentEntry();
        if (multiplicityVector[eventIndex] < multCutOff) {continue;}
        numSelectedEvents++;
        //std::cout << entryIndex << std::endl;

        // ***** PARTICLE LOOP ******
        // Particle 1 Loop
        for (Int_t i = 0; i < pxBranch->size()-1; i++) {

            // Checking if particle 1 is charged
            if ((*chgBranch)[i] == 0) {continue;}

            // Caluclating eta and phi for particle 1
            Float_t eta1 = calculateEta(calculateTheta((*pxBranch)[i], (*pyBranch)[i], (*pzBranch)[i]));
            if (fabs(eta1) > etaCutOff) {continue;} // Checking eta range
            Float_t phi1 = calculatePhi((*pxBranch)[i], (*pyBranch)[i]);
            numTrigg++;
            //std::cout << "Trigger Particle " << i << ", eta: " << eta1 << ", phi: " << phi1 << std::endl;

            // Particle 2 Loop
            for (Int_t j = i + 1; j < pxBranch->size(); j++) {

                // Checking if particle 2 is charged
                if ((*chgBranch)[j] == 0) {continue;}

                // Calculating eta and phi for particle 2
                Float_t eta2 = calculateEta(calculateTheta((*pxBranch)[j], (*pyBranch)[j], (*pzBranch)[j]));
                if (fabs(eta2) > etaCutOff) {continue;} // Checking eta range
                Float_t phi2 = calculatePhi((*pxBranch)[j], (*pyBranch)[j]);

                // Calculating delta eta and delta phi
                Float_t deltaEta = eta2 - eta1;
                Float_t deltaPhi = TMath::ACos(TMath::Cos(phi2-phi1));

                // Filling the histograms multiple times due to symmetries
                hSignal.Fill(deltaEta, deltaPhi);
                hSignal.Fill(deltaEta, -deltaPhi);
                hSignal.Fill(-deltaEta, -deltaPhi);
                hSignal.Fill(-deltaEta, deltaPhi);
                hSignal.Fill(deltaEta, 2*TMath::Pi()-deltaPhi);
                hSignal.Fill(-deltaEta, 2*TMath::Pi()-deltaPhi);

                //hSignal.Fill(deltaEta, deltaPhi, 1.0/(weightVector[eventIndex]));
                //hSignal.Fill(deltaEta, -deltaPhi, 1.0/(weightVector[eventIndex]));
                   
                //hSignal.Fill(deltaEta, deltaPhi, weightVector[eventIndex]);
            }
        } 
    } 
    std::cout << "Number of selected events for the signal distribution: " << numSelectedEvents << std::endl;

    // ***** NORMALIZATION ******
    //hSignal.Scale(1.0/(reader->GetEntries()));
    hSignal.Scale(1.0/numTrigg); 

    // ***** HISTOGRAM CUSTOMIZATION ******
    hSignal.SetXTitle("#Delta#eta");
    hSignal.SetYTitle("#Delta#phi");
    hSignal.SetZTitle("S(#Delta#eta, #Delta#phi)");
    
    hSignal.SetTitleOffset(1.5, "X");
    hSignal.SetTitleOffset(1.5, "Y");
    hSignal.SetTitleOffset(1.3, "Z");

    hSignal.SetTitleFont(132, "T");
    hSignal.SetTitleFont(132, "XYZ");

    hSignal.SetLabelFont(132, "T");
    hSignal.SetLabelFont(132, "XYZ");

    hSignal.SetStats(0);

    return hSignal;
}

// Function that creates the eta-phi distribution for the pseudoparticle mixing technique
TH2F createEtaPhiDist(std::vector<Int_t> multiplicityVector) {

    // First intialize the eta-phi distribution for all particles
    TH2F etaPhiDist("etaPhiDist", "(#eta, #phi) Distribution for all Particles", 100, -etaCutOff, etaCutOff, 100, -TMath::Pi(), TMath::Pi());
     
    reader->Restart(); // Restarting event loop

    while (reader->Next()) {

        // Check to see if the event is in the multiplicity bin
        Int_t entryIndex = reader->GetCurrentEntry();
        if (multiplicityVector[entryIndex] < multCutOff) {continue;}

        for (Int_t i = 0; i < pxBranch->size(); i++) {

            // Checking if particle is charged
            if ((*chgBranch)[i] == 0) {continue;}

            // Caluclating eta and phi for particle 
            Float_t eta = calculateEta(calculateTheta((*pxBranch)[i], (*pyBranch)[i], (*pzBranch)[i]));
            if (fabs(eta) > etaCutOff) {continue;} // Checking eta range
            Float_t phi = calculatePhi((*pxBranch)[i], (*pyBranch)[i]);

            etaPhiDist.Fill(eta, phi);
        }
    }

    // ***** CUSTOMIZATION *****
    etaPhiDist.SetXTitle("#eta");
    etaPhiDist.SetYTitle("#phi");
    etaPhiDist.SetZTitle("Entries");
    etaPhiDist.SetTitleOffset(1.5, "X");
    etaPhiDist.SetTitleOffset(1.5, "Y");
    etaPhiDist.SetTitleOffset(1.2, "Z");

    etaPhiDist.SetTitleFont(132, "T");
    etaPhiDist.SetTitleFont(132, "X");
    etaPhiDist.SetTitleFont(132, "Y");
    etaPhiDist.SetTitleFont(132, "Z");
    
    return etaPhiDist;
}


// Function that returns the distribution of the simple mixed-event background function
TH2F createBackgroundDist(std::vector<Int_t> multiplicityVector, Int_t numMixFactor, TH2F hSignal, TH2F etaPhiDist) {

    std::string backgroundTitle = "Background Distribution for " + title;

    // Histogram for the background distribution
    //TH2F hBackground("hBackground", backgroundTitle.c_str(), 100, -5, 5, 100, -TMath::Pi(), TMath::Pi());
    
    const float         EtaBW = 0.3; 
    const float         PhiBW = TMath::Pi()/16;
    TH2F hBackground("hBackground", backgroundTitle.c_str(),41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);

    // ***** PSEUDOPARTICLE MIXING *****
    // Calculating the number of pseudoparticles samples:
    double numSigEntries = hSignal.GetEntries();
    //Int_t numPseudo = numSigEntries * 8 * numMixFactor;
    //numPseudo *= 4;
    //numPseudo *= 2;
    double numPseudo = (1 + floor(sqrt(1+(4*2*(double)numMixFactor*numSigEntries)))) / 2; // Not sure why floor is added after the sqrt (used in Austin's code)
    //cout<<numSigEntries<<endl;
    //cout<<numMixFactor*numSigEntries<<endl;
    //cout<<2*numMixFactor*numSigEntries<<endl;
    //cout<<(double)(8*numMixFactor*numSigEntries)<<endl;
    // (Manually calculated)
    //Int_t numPseudo = 31190;
    std::cout << "Number of entries in signal distribution: " << numSigEntries << std::endl;
    std::cout << "Number of pseudoparticles to select: " << numPseudo << std::endl;

    for (Int_t i = 0; i < numPseudo-1; i++) {
        Double_t selectedEta1, selectedPhi1;
        etaPhiDist.GetRandom2(selectedEta1, selectedPhi1);
        
        for (Int_t j = i+1; j < numPseudo; j++) {
            Double_t selectedEta2, selectedPhi2;
            etaPhiDist.GetRandom2(selectedEta2, selectedPhi2);
            
            // Calculating delta eta and delta phi
            Float_t deltaEta = selectedEta2 - selectedEta1;
            Float_t deltaPhi = TMath::ACos(TMath::Cos(selectedPhi2-selectedPhi1));

            //std::cout << i << ": (" << deltaEta << ", " << deltaPhi << ")" << std::endl;

            // Filling the histograms twice times due to symmetries
            hBackground.Fill(deltaEta, deltaPhi);
            hBackground.Fill(deltaEta, -deltaPhi);
            hBackground.Fill(-deltaEta, -deltaPhi);
            hBackground.Fill(-deltaEta, deltaPhi);
            hBackground.Fill(deltaEta, 2*TMath::Pi()-deltaPhi);
            hBackground.Fill(-deltaEta, 2*TMath::Pi()-deltaPhi);

            ///hBackground.Fill(deltaEta, deltaPhi);
        }
    } 

    // ***** NORMALIZATION *****
    //hBackground.Scale(1.0/numPseudo);

    // ***** HISTOGRAM CUSTOMIZATION *****
    hBackground.SetXTitle("#Delta#eta");
    hBackground.SetYTitle("#Delta#phi");
    hBackground.SetZTitle("B(#Delta#eta, #Delta#phi)");

    hBackground.SetTitleOffset(1.5, "X");
    hBackground.SetTitleOffset(1.5, "Y");
    hBackground.SetTitleOffset(1.3, "Z");

    hBackground.SetTitleFont(132, "T");
    hBackground.SetTitleFont(132, "XYZ");

    hBackground.SetLabelFont(132, "T");
    hBackground.SetLabelFont(132, "XYZ");

    hBackground.SetStats(0);

    return hBackground; 
} 

void simpleTwoParticleCorr() {

    // Finding the multiplicities of each event and storing it in a vector
    std::vector<Int_t> multVec;
    reader->Restart(); // Ensuring event loop starts from beginning
    while (reader->Next()) {
        Int_t multiplicity = 0; // Counter for N_ch
        for (Int_t i = 0; i < pxBranch->size(); i++) {
            if ((*chgBranch)[i] != 0) {multiplicity++;}
        }
        multVec.push_back(multiplicity);
        //std::cout << "Index " << reader->GetCurrentEntry() << ": " << multiplicity << std::endl;
    }

    // Finding the weights (num. trigger particles) for the signal distribution
    std::vector<Int_t> weightVec; // Stores the weights (number of trigger particles)
    reader->Restart(); // Ensuring event loop starts from beginning
    while (reader->Next()) {
        Int_t numTrigg = 0;
        // Don't apply the N_ch cut because we want to keep the size of weightVec the same as the total number of events
        if (multVec[reader->GetCurrentEntry()] < multCutOff) {
            weightVec.push_back(numTrigg);
            continue;
        }
        for (Int_t i = 0; i < pxBranch->size(); i++) {
            if ((*chgBranch)[i] == 0) {continue;}
            Float_t eta1 = calculateEta(calculateTheta((*pxBranch)[i], (*pyBranch)[i], (*pzBranch)[i]));
            if (fabs(eta1) > etaCutOff) {continue;} // Checking eta range
            numTrigg++;
        }
        weightVec.push_back(numTrigg);
    }
    
    //std::cout << "Weight vector size: " << weightVec.size() << std::endl;
    //std::cout << "Multiplicity vector size: " << multVec.size() << std::endl;

    TFile *fout = new TFile("simpleTwoParticleCorr.root", "recreate"); // Creating output file

    // Creating canvas for the signal histogram
    TCanvas *cSignal = new TCanvas("cSignal", "Canvas for the Signal Distribution", 800, 600);
    TH2F simpleSignalHist = createSimpleSignalDist(multVec, weightVec);

    cSignal->cd();
    simpleSignalHist.Draw("SURF1");
    cSignal->Write();

    // Testing eta/phi distribution
    TCanvas *cEtaPhi = new TCanvas("cEtaPhi", "Canvas for the Eta-Phi Distribution", 800, 600);
    TH2F etaPhiHist = createEtaPhiDist(multVec);
    
    cEtaPhi->cd();
    etaPhiHist.Draw("SURF1");
    cEtaPhi->Write();

    // Creating canvas for the background histogram
    TCanvas *cBackground = new TCanvas("cBackground", "Canvas for the Background Distribution", 800, 600);
    TH2F simpleBackgroundHist = createBackgroundDist(multVec, 10, simpleSignalHist, etaPhiHist);

    cBackground->cd();
    simpleBackgroundHist.Draw("SURF1");
    cBackground->Write();

    // Creating canvas for the corrected signal distribution
    TCanvas *cCorrected = new TCanvas("cCorrected", "Canvas for the Corrected Signal Distribution", 800, 600);
    TH2F correctedHist = simpleSignalHist;

    correctedHist.Divide(&simpleBackgroundHist);
    std::cout << "B(0,0): " << simpleBackgroundHist.GetBinContent(simpleBackgroundHist.FindBin(0,0)) << std::endl;
    correctedHist.Scale(simpleBackgroundHist.GetBinContent(simpleBackgroundHist.FindBin(0,0)));

    cCorrected->cd();
    correctedHist.Draw("SURF1");

    // ***** HISTOGRAM CUSTOMIZATION ***** //
    std::string correctedTitle = "Corrected Signal Distribution for " + title;
    correctedHist.SetTitle(correctedTitle.c_str());
    correctedHist.GetZaxis()->SetTitle("C(#Delta#eta, #Delta#phi)");

    correctedHist.SetTitleOffset(1.1, "Z");
    correctedHist.SetTitleFont(132, "T");
    correctedHist.SetTitleFont(132, "XYZ");
    correctedHist.SetLabelFont(132, "T");
    correctedHist.SetLabelFont(132, "XYZ");
    //correctedHist.SetMaximum(2);

    cCorrected->Write(); 

    // Creating canvas for the projected delta phi histgram
    TCanvas *cProjection = new TCanvas("cProjection", "Canvas for the Projected Distributions", 800, 800);
    cProjection->Divide(1, 2);

    cProjection->cd(1);
    gPad->SetGrid();
    TH1D *projectedSignalHist = simpleSignalHist.ProjectionY("projectedSignalHist", 1, -1);
    projectedSignalHist->SetLineWidth(2);
    projectedSignalHist->SetLineColor(kBlue);
    projectedSignalHist->Draw("HIST L");

    cProjection->cd(2);
    gPad->SetGrid();
    TH1D *projectedBackgroundHist = simpleBackgroundHist.ProjectionY("projectedBackgroundHist", 1, -1);
    projectedBackgroundHist->SetMinimum(0);
    projectedBackgroundHist->SetMaximum(6000000);
    projectedBackgroundHist->SetLineWidth(2);
    projectedBackgroundHist->SetLineColor(kBlue);
    projectedBackgroundHist->Draw("HIST L SAME");

    cProjection->Write(); 

    delete cSignal;
    delete cEtaPhi;
    delete cBackground;
    delete cCorrected;
    delete cProjection; 
    f->Close();
    fout->Close();
}