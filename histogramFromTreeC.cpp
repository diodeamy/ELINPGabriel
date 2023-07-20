#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include <cmath>

void histogramFromTreeC() {
    
    // 
    TTree *tree;
    UShort_t detectorId;
    UShort_t energy;
    TBranch *b_detector, *b_energy; 
    
    //
    TFile *digitizer_1 = TFile::Open("data/deuterated-digitizer_1.root");
    digitizer_1->GetObject("events", tree);

    //
    tree->SetBranchAddress("detectorId", &detectorId, &b_detector);
    tree->SetBranchAddress("energy", &energy, &b_energy);

    //
    Long64_t nentries = tree->GetEntries();
    TH1D *hist1 = new TH1D("hist1", "energy", pow(2,14), 0, pow(2,14));

    TCanvas *c1 = new TCanvas("c1", "Histograms", 200, 10, 700, 900);
    c1-> Draw();    
    hist1->Draw();

    for (int jentry=0; jentry<nentries; jentry++) {
        tree->GetEntry(jentry);
        if (detectorId == 1) {
            hist1->Fill(energy);
        }
    }
}