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

void energyCalibration() {
    TFile *f1 = TFile::Open("data/deuterated-digitizer_1.root");
    TFile *f2 = TFile::Open("data/deuterated-digitizer_2.root");
    TFile *f3 = TFile::Open("data/deuterated-digitizer_3.root");

    TTree *tree1, *tree2, *tree3;
    UShort_t *detectorId;
    UShort_t *energy;


    f1->GetObject("events", tree1);
    f2->GetObject("events", tree2);
    f3->GetObject("events", tree3);

    tree1->SetBranchAddress("detectorId", &detectorId);
    tree1->SetBranchAddress("energy", &energy);
    
    tree2->SetBranchAddress("detectorId", &detectorId);
    tree2->SetBranchAddress("energy", &energy);

}