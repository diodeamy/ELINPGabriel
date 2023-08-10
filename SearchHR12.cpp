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
vector<Double_t> calibratePeaks(vector<Double_t> v1, vector<Double_t> v2);
vector<Double_t> getHistPeaks(string name);

void calibrationFinal()
{
   auto peaks1 = getHistPeaks("data/deuterated-digitizer_1.root");
   auto peaks2 = getHistPeaks("data/deuterated-digitizer_2.root");

   calibratePeaks(peaks1, peaks2);
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
   TTree *tree, *tree_n;
   UShort_t detectorId, detectorId_n;
   UShort_t energy, energy_n;
   Double_t fPositionX_1[100], fPositionX_2[100], fPositionX_1_n[100], fPositionX_2_n[100];
   Double_t fPositionY_1[100], fPositionY_2[100], fPositionY_1_n[100], fPositionY_2_n[100];
   Double_t refPositionX[18], refPositionX_n[18];
   Double_t refPositionY[18], refPositionY_n[18];
   Int_t fNPeaks = 0;
   Int_t i, nfound_1, nfound_2, nfound_1_n, nfound_2_n, bin, bin_n;
   const Int_t nbins = 4000;
   Double_t xmin_1 = 0, xmin_2 = 10500; // xmin_1 is for the first 4 peaks and xmin_2 for the last 5
   Double_t xmax_1 = nbins, xmax_2 = nbins + xmin_2;
   Double_t a, a_n;
   Double_t source[nbins], dest_1[nbins], dest_2[nbins], dest_1_n[nbins], dest_2_n[nbins];
   Double_t sigma_1 = 12, sigma_2 = 19, threshold_1 = 0.8, threshold_2 = 0.9;
   Int_t deconiterations_1 = 4, deconiterations_2 = 2, averWindow_1 = 2, averWindow_2 = 7;

   // Dataset 1 - the standard

   TFile *f1 = TFile::Open("data/deuterated-digitizer_1.root");
   f1->GetObject("events", tree);

   tree->SetBranchAddress("energy", &energy);
   tree->SetBranchAddress("detectorId", &detectorId);

   Long64_t nentries = tree->GetEntries();

   // first part stuff
   TH1F *h_1 = new TH1F("hist1", "energy", nbins, xmin_1, xmax_1);

   for (int jentry = 0; jentry < nentries; jentry++)
   {
      tree->GetEntry(jentry);
      if (detectorId == 0)
      {
         h_1->Fill(energy);
      }
   }

   h_1->SetTitle("High resolution peak searching, number of iterations = 4");
   h_1->GetXaxis()->SetRange(1, nbins);
   TH1F *d_1 = new TH1F("d_1", "", nbins, xmin_1, xmax_1);
   h_1->Draw("L");
   for (i = 0; i < nbins; i++)
      source[i] = h_1->GetBinContent(i + 1);

   TSpectrum *s_1 = new TSpectrum();

   nfound_1 = s_1->SearchHighRes(source, dest_1, nbins, sigma_1, threshold_1, kTRUE, deconiterations_1, kTRUE, averWindow_1);
   Double_t *xpeaks_1 = s_1->GetPositionX();
   for (i = 0; i < nfound_1; i++)
   {
      a = xpeaks_1[i];
      bin = 1 + Int_t(a + 0.5);
      fPositionX_1[i] = h_1->GetBinCenter(bin);
      fPositionY_1[i] = h_1->GetBinContent(bin);
   }

   TPolyMarker *pm_1 = (TPolyMarker *)h_1->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1)
   {
      h_1->GetListOfFunctions()->Remove(pm_1);
      delete pm_1;
   }
   pm_1 = new TPolyMarker(nfound_1, fPositionX_1, fPositionY_1);
   h_1->GetListOfFunctions()->Add(pm_1);
   pm_1->SetMarkerStyle(23);
   pm_1->SetMarkerColor(kRed);
   pm_1->SetMarkerSize(1.3);

   for (i = 0; i < nbins; i++)
      d_1->SetBinContent(i + 1, dest_1[i]);
   d_1->SetLineColor(kRed);
   d_1->Draw("SAME");

   // second part stuff//
   TH1F *h_2 = new TH1F("hist2", "energy", nbins, xmin_2, xmax_2);

   for (int jentry = 0; jentry < nentries; jentry++)
   {
      tree->GetEntry(jentry);
      if (detectorId == 0)
      {
         h_2->Fill(energy);
      }
   }

   h_2->SetTitle("High resolution peak searching, number of iterations = 2");
   h_2->GetXaxis()->SetRange(1, nbins);
   printf("hello");
   TH1F *d_2 = new TH1F("d_2", "", nbins, xmin_2, xmax_2);
   h_2->Draw("L");
   for (i = 0; i < nbins; i++)
      source[i] = h_2->GetBinContent(i + 1);

   TSpectrum *s_2 = new TSpectrum();

   nfound_2 = s_2->SearchHighRes(source, dest_2, nbins, sigma_2, threshold_2, kTRUE, deconiterations_2, kTRUE, averWindow_2);
   Double_t *xpeaks_2 = s_2->GetPositionX();
   for (i = 0; i < nfound_2; i++)
   {
      a = xpeaks_2[i];
      bin = 1 + Int_t(a + 0.5);
      fPositionX_2[i] = h_2->GetBinCenter(bin);
      fPositionY_2[i] = h_2->GetBinContent(bin);
   }

   TPolyMarker *pm_2 = (TPolyMarker *)h_2->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2)
   {
      h_2->GetListOfFunctions()->Remove(pm_2);
      delete pm_2;
   }
   pm_2 = new TPolyMarker(nfound_2, fPositionX_2, fPositionY_2);
   h_2->GetListOfFunctions()->Add(pm_2);
   pm_2->SetMarkerStyle(23);
   pm_2->SetMarkerColor(kRed);
   pm_2->SetMarkerSize(1.3);

   for (i = 0; i < nbins; i++)
      d_2->SetBinContent(i + 1, dest_2[i]);
   d_2->SetLineColor(kRed);
   d_2->Draw("SAME");

   // printf("Reference dataset : Found %d candidate peaks in first part\n and %d in second part \n",nfound_1, nfound_2);
   // for (i=0;i<nfound_1;i++) printf("posx= %f, posy= %f\n",fPositionX_1[i], fPositionY_1[i]);
   // for (i=0;i<nfound_2;i++) printf("posx = %f, posy = %f\n", fPositionX_2[i], fPositionY_2[i]);

   // make single array of the found possitions
   for (i = 0; i < nfound_1; i++)
   {
      refPositionX[i] = fPositionX_1[i];
      refPositionY[i] = fPositionY_1[i];
   }
   for (i = 0; i < nfound_2; i++)
   {
      refPositionX[i + 4] = fPositionX_2[i];
      refPositionY[i + 4] = fPositionY_2[i];
   }
   // for(i=0; i<9;i++) printf("The positions of all peaks is: \n posx = %f    posy = %f \n", refPositionX[i], refPositionY[i]);
   printf("The positions for the 9 peaks of the reference dataset have been found. \n Calibrating other datasetss...\n");

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   // now that we have the "standard" (i.e. the reference peakpoints for calibration) we can import any other file
   // from the folder and proceed to try and calibrate them

   // Some dataset n

   TFile *fn = TFile::Open("data/deuterated-digitizer_3.root");
   fn->GetObject("events", tree_n);

   tree_n->SetBranchAddress("energy", &energy_n);
   tree_n->SetBranchAddress("detectorId", &detectorId_n);

   Long64_t nentries_n = tree_n->GetEntries();

   // first part stuff
   TH1F *h_1_n = new TH1F("hist1_n", "energy", nbins, xmin_1, xmax_1);

   for (int jentry = 0; jentry < nentries_n; jentry++)
   {
      tree_n->GetEntry(jentry);
      if (detectorId_n == 0)
      {
         h_1_n->Fill(energy_n);
      }
   }

   h_1_n->SetTitle("High resolution peak searching, number of iterations = 4");
   h_1_n->GetXaxis()->SetRange(1, nbins);
   TH1F *d_1_n = new TH1F("d_1_n", "", nbins, xmin_1, xmax_1);
   h_1_n->Draw("L");
   for (i = 0; i < nbins; i++)
      source[i] = h_1_n->GetBinContent(i + 1);

   TSpectrum *s_1_n = new TSpectrum();

   nfound_1_n = s_1_n->SearchHighRes(source, dest_1_n, nbins, sigma_1, threshold_1, kTRUE, deconiterations_1, kTRUE, averWindow_1);
   Double_t *xpeaks_1_n = s_1_n->GetPositionX();
   for (i = 0; i < nfound_1_n; i++)
   {
      a_n = xpeaks_1_n[i];
      bin_n = 1 + Int_t(a_n + 0.5);
      fPositionX_1_n[i] = h_1_n->GetBinCenter(bin_n);
      fPositionY_1_n[i] = h_1_n->GetBinContent(bin_n);
   }

   TPolyMarker *pm_1_n = (TPolyMarker *)h_1_n->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1_n)
   {
      h_1_n->GetListOfFunctions()->Remove(pm_1_n);
      delete pm_1_n;
   }
   pm_1_n = new TPolyMarker(nfound_1_n, fPositionX_1_n, fPositionY_1_n);
   h_1_n->GetListOfFunctions()->Add(pm_1_n);
   pm_1_n->SetMarkerStyle(23);
   pm_1_n->SetMarkerColor(kRed);
   pm_1_n->SetMarkerSize(1.3);

   for (i = 0; i < nbins; i++)
      d_1_n->SetBinContent(i + 1, dest_1_n[i]);
   d_1_n->SetLineColor(kRed);
   d_1_n->Draw("SAME");

   // second part stuff//
   TH1F *h_2_n = new TH1F("hist2_n", "energy", nbins, xmin_2, xmax_2);

   for (int jentry = 0; jentry < nentries_n; jentry++)
   {
      tree_n->GetEntry(jentry);
      if (detectorId_n == 0)
      {
         h_2_n->Fill(energy_n);
      }
   }

   h_2_n->SetTitle("High resolution peak searching, number of iterations = 2");
   h_2_n->GetXaxis()->SetRange(1, nbins);
   printf("hello");
   TH1F *d_2_n = new TH1F("d_2_n", "", nbins, xmin_2, xmax_2);
   h_2_n->Draw("L");
   for (i = 0; i < nbins; i++)
      source[i] = h_2_n->GetBinContent(i + 1);

   TSpectrum *s_2_n = new TSpectrum();

   nfound_2_n = s_2_n->SearchHighRes(source, dest_2_n, nbins, sigma_2, threshold_2, kTRUE, deconiterations_2, kTRUE, averWindow_2);
   Double_t *xpeaks_2_n = s_2_n->GetPositionX();
   for (i = 0; i < nfound_2_n; i++)
   {
      a_n = xpeaks_2_n[i];
      bin_n = 1 + Int_t(a_n + 0.5);
      fPositionX_2_n[i] = h_2_n->GetBinCenter(bin_n);
      fPositionY_2_n[i] = h_2_n->GetBinContent(bin_n);
   }

   TPolyMarker *pm_2_n = (TPolyMarker *)h_2_n->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2_n)
   {
      h_2_n->GetListOfFunctions()->Remove(pm_2_n);
      delete pm_2_n;
   }
   pm_2_n = new TPolyMarker(nfound_2_n, fPositionX_2_n, fPositionY_2_n);
   h_2_n->GetListOfFunctions()->Add(pm_2_n);
   pm_2_n->SetMarkerStyle(23);
   pm_2_n->SetMarkerColor(kRed);
   pm_2_n->SetMarkerSize(1.3);

   for (i = 0; i < nbins; i++)
      d_2_n->SetBinContent(i + 1, dest_2_n[i]);
   d_2_n->SetLineColor(kRed);
   d_2_n->Draw("SAME");

   // printf("Dataset n: Found %d candidate peaks in first part\n and %d in second part \n",nfound_1_n, nfound_2_n);
   // for( i=0;i<nfound_1_n;i++) printf("posx= %f, posy= %f\n",fPositionX_1_n[i], fPositionY_1_n[i]);
   // for (i=0;i<nfound_2_n;i++) printf("posx = %f, posy = %f\n", fPositionX_2_n[i], fPositionY_2_n[i]);

   for (i = 0; i < nfound_1_n; i++)
   {
      refPositionX_n[i] = fPositionX_1_n[i];
      refPositionY_n[i] = fPositionY_1_n[i];
   }
   for (i = 0; i < nfound_2_n; i++)
   {
      refPositionX_n[i + 4] = fPositionX_2_n[i];
      refPositionY_n[i + 4] = fPositionY_2_n[i];
   }
   // for(i=0; i<9;i++) printf("The positions of all peaks is: \n posx = %f    posy = %f \n", refPositionX_n[i], refPositionY_n[i]);

   // LESSGOOOO, now we can finally plot peak vs peak and try to fit a function through the graph for a conversion

   double x[100], y[100];
   for (int i = 0; i < 9; i++)
   {
      x[i] = refPositionX[i];
      y[i] = refPositionX_n[i];
   }
   std::sort(x, x + 10);
   std::sort(y, y + 10);

   auto g = new TGraph(10, x, y);
   g->SetTitle("Peak calibration ;X: reference peaks;Y: dataset n peaks");
   g->Draw("AC*");

   // TF1* linearFit = new TF1("linearFit", "[0]*x + [1]", x[0], x[9]);
   // g->Fit(linearFit, "R");

   // double a_linear = linearFit->GetParameter(0);
   // double b_linear = linearFit->GetParameter(1);
   // double a_error_linear = linearFit->GetParError(0);
   // double b_error_linear = linearFit->GetParError(1);

   // std::cout << "Linear Fit Parameters:" << std::endl;
   // std::cout << "a = " << a_linear << " +/- " << a_error_linear << std::endl;
   // std::cout << "b = " << b_linear << " +/- " << b_error_linear << std::endl;

   // TF1* parabolicFit = new TF1("parabolicFit", "[0]*x*x + [1]*x + [2]", x[0], x[9]);
   // g->Fit(parabolicFit, "R");

   // double a_parabolic = parabolicFit->GetParameter(0);
   // double b_parabolic = parabolicFit->GetParameter(1);
   // double c_parabolic = parabolicFit->GetParameter(2);
   // double a_error_parabolic = parabolicFit->GetParError(0);
   // double b_error_parabolic = parabolicFit->GetParError(1);
   // double c_error_parabolic = parabolicFit->GetParError(2);

   // std::cout << "Parabolic Fit Parameters:" << std::endl;
   // std::cout << "a = " << a_parabolic << " +/- " << a_error_parabolic << std::endl;
   // std::cout << "b = " << b_parabolic << " +/- " << b_error_parabolic << std::endl;
   // std::cout << "c = " << c_parabolic << " +/- " << c_error_parabolic << std::endl;

   // i need to aply the linear/parabolic function to the datapoints in dataset n in order to bring each peak more or less
   // where the peak of the reference dataset is

   // for detectorId == 0 we need to apply the parameters to dataset_n at energy[i] for all i in dataset

   //    for (i=0; i<tree_n->GetEntries(); i++);
   //    tree_n->GetEntry(i);

   //    if (detectorId ==0) {
   //       energy_n[i] = appliedFit(energy_n[i]);
   //    }

   // }

   // void appliedFit(UShort energy_n, Int_t fitType) {
   //    if (fitType ==0) {
   //       energy_n = `
   //    }
   //    else if (fitType ==1) {

   //    }
}

TH1F *getHistogram(string fileName)
{
   UShort_t energy, detectorId;

   const char *filenamechar = fileName.c_str();
   TFile *f1 = TFile::Open(filenamechar); // fix this! TFIle does not take anything but const char and we need to pass each dataset

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

vector<Double_t> calibratePeaks(vector<Double_t> v1, vector<Double_t> v2)
{
   /*auto xPosPlot = TGraph(); // 10, 0, v1.size()
   for (Int_t i = 0; i < v1.size(); i++)
   {
      xPosPlot.AddPoint(v1[i], v2[i]);
   }
   xPosPlot.Draw("AC"); // AC*
   return {};
   */

   for (Int_t i = 0; i < v1.size(); i++)
      cout << v1[i] << endl;
   for (Int_t i = 0; i < v2.size(); i++)
      cout << v2[i] << endl;
   return {};
}