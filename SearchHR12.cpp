#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include <iostream>
using namespace std;

TH1F *getHistogram(string fileName);
vector<Double_t> getPeaks(TH1F *hist, Int_t peaknr, Int_t xmin, Int_t xmax);
vector<Double_t> getFitParams(vector<Double_t> v1, vector<Double_t> v2);
vector<Double_t> getHistPeaks(string name);
TTree makefriendtr(string fileName);

void calibrationFinal()
{
   auto peaks1 = getHistPeaks("data/deuterated-digitizer_1.root");
   auto peaks2 = getHistPeaks("data/deuterated-digitizer_2.root");
}

vector<Double_t> getHistPeaks(std::string name)
{
   auto hist = getHistogram(name);

   auto peaks1 = getPeaks(hist, 3, 800, 2500);
   auto peaks2 = getPeaks(hist, 6, 10000, 14000);
   peaks1.insert(peaks1.end(), peaks2.begin(), peaks2.end());

   return peaks1;
}

void SearchHR12()
{
}

TH1F *getHistogram(string fileName)
{
   UShort_t energy, detectorId;

   const char *filenamechar = fileName.c_str();
   TFile *f1 = TFile::Open(filenamechar);

   auto tree = (TTree *)f1->Get("events");

   tree->SetBranchAddress("energy", &energy);
   tree->SetBranchAddress("detectorId", &detectorId);

   Long64_t nentries = tree->GetEntries();

   TH1F *h_1 = new TH1F("hist1", "energy", pow(2, 14), 0, pow(2, 14));

   for (int jentry = 0; jentry < nentries; jentry++)
   {
      tree->GetEntry(jentry);
      if (detectorId == 0)
      {
         h_1->Fill(energy);
      }
   }
   return h_1;
}

vector<Double_t> getPeaks(TH1F *hist, Int_t peaknr, Int_t xmin, Int_t xmax)
{
   TSpectrum *s_1 = new TSpectrum(peaknr);
   hist->GetXaxis()->SetRangeUser(xmin, xmax);
   s_1->Search(hist, 50, "", 0.0005);

   Double_t *xpeaks = s_1->GetPositionX();

   auto found = s_1->GetNPeaks();

   std::vector<Double_t> v(xpeaks, xpeaks + found);
   return v;
}

vector<Double_t> getFitParams(vector<Double_t> v1, vector<Double_t> v2)
{
   std::sort(v1.data(), v1.data() + v1.size());
   std::sort(v2.data(), v2.data() + v2.size());

   auto xPosPlot = new TGraph(v1.size(), v1.data(), v2.data());
   xPosPlot->Draw("AC*");

   auto linearFit = new TF1("linearFit", "[0]*x+[1]", v1[0], v1[v1.size()]);
   xPosPlot->Fit(linearFit, "R");
   auto slope = linearFit->GetParameter(0);
   auto intercept = linearFit->GetParameter(1);

   vector<Double_t> v(slope, intercept);
   return v;
}
TTree makefriendtr(string fileName, string file2Name)
// should make a friend tree of the file to be calibrated (so not digitizer_1 but digitizer_3)
// and add a new branch to hold the calibrated energy
{
   const char *file2namechar = file2Name.c_str();
   const char *filenamechar = fileName.c_str();
   TFile *f1 = TFile::Open(filenamechar);
   auto T = (TTree *)f1->Get("events");

   TFile *ff = new TFile("treefriend.root", "recreate");
   TTree *TF = T->CloneTree(0);
   TF->SetName("TF");

   UShort_t energy, detectorId;
   TF->SetBranchAddress("energy", &energy);
   TF->SetBranchAddress("detectorId", &detectorId);

   if (detectorId == 0)
   {
      TF->GetBranch("events")->SetFile("eventsd0c.root");
      TF->CopyEntries(T);
   }

   TF->Print();
   TF->Write();
   delete ff;

   auto parameters = getFitParams(getHistPeaks(fileName), getHistPeaks(file2Name));

   Long64_t nentries = TF->GetEntries();
   for (int jentry = 0; jentry < nentries; jentry++)
   {
      TF->GetEntry(jentry);
      jentry *parameters[0] + parameters[1];
   }

   for (int jentry = 0; jentry < 100; jentry++)
   {
      cout << TF->GetEntry(jentry) << endl;
      cout << T->GetEntry(jentry) << endl;
   }
   return {};
}
