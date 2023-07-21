// ## 3meV?
// ## extragem arii din zona 
// ## aliniem un anumit peak which gets shifted left and right because of shifts in temperatue
// ## n fisiere trebuie alineate cateva peak-uri
// ## faci un grafic energii prim fisier si un al doilea (E1 vs E2) care os a iti dea exact o linie/polinom cu care sa faci un fit;
// ## scriptul o sa trb ssa faca asta pt fiecare alt grafic 
// ## TSpectrum (poate te ajuta sa gasesti peak-uri) "ROOT::TSpectrum in Python" 
// ## poate merge cu niste hard-coded ranges intre care nu se shiftuiesc peakurile
// ## poti sa faci alta chestie pt ultimul: mergi din dreapta pana gasesti primul peak
// ## MakeFile files for easy compilation (" ex g++ *.cpp -std=c++17")
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <cmath>
#include "TFile.h"
#include "TTree.h"

void energyCalibration() {
    TFile *f1 = TFile::Open("data/deuterated-digitizer_1.root");
    TFile *f2 = TFile::Open("data/deuterated-digitizer_2.root");
    TFile *f3 = TFile::Open("data/deuterated-digitizer_3.root");

    TTree *tree1, *tree2, *tree3;
    UShort_t detectorId;
    UShort_t energy;
    Int_t npeaks = 10;
    const Int_t binnum = pow(2, 14);
    auto binmin = 0;
    auto binmax = binnum;

    f1->GetObject("events", tree1);
    f2->GetObject("events", tree2);
    f3->GetObject("events", tree3);

    tree1->SetBranchAddress("detectorId", &detectorId);
    tree1->SetBranchAddress("energy", &energy);
    
    tree2->SetBranchAddress("detectorId", &detectorId);
    tree2->SetBranchAddress("energy", &energy);

    tree3->SetBranchAddress("detectorId", &detectorId);
    tree3->SetBranchAddress("energy", &energy);
    
    auto nentries1 = tree1->GetEntries();

    TH1D *h1 = new TH1D("h1", "energy", binnum, binmin, binmax);
    
    for (int jentry=0; jentry<nentries1; jentry++) {
        tree1->GetEntry(jentry);
        if (detectorId == 0) {
            h1->Fill(energy);
        }
    }
    h1->Draw();

    TSpectrum *s1 = new TSpectrum();
    Int_t pfound1 = s1->Search(h1, 2.7,"", 0.05);     
    printf("Found %d candidate peaks to fit \n", pfound1);

    //attempts to find peaks in certain areas of the graph using TSpectrum::Search()
}