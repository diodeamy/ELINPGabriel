#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TTree.h"
#include <iostream>


void SearchHR12() {
   TTree *tree, *tree_n;
   UShort_t detectorId, detectorId_n;
   UShort_t energy, energy_n;
   Double_t fPositionX_1[100], fPositionX_2[100], fPositionX_1_n[100], fPositionX_2_n[100];
   Double_t fPositionY_1[100], fPositionY_2[100], fPositionY_1_n[100], fPositionY_2_n[100];
   Double_t refPositionX[18];
   Double_t refPositionY[18];
   Int_t fNPeaks = 0;
   Int_t i, nfound_1, nfound_2, nfound_1_n, nfound_2_n, bin, bin_n;
   const Int_t nbins = 4000;
   Double_t xmin_1 = 0, xmin_2 = 10500; //xmin_1 is for the first 4 peaks and xmin_2 for the last 5
   Double_t xmax_1= nbins, xmax_2= nbins+xmin_2; 
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

   //first part stuff
   TH1F *h_1 = new TH1F("hist1", "energy", nbins, xmin_1, xmax_1);
    
   for (int jentry=0; jentry<nentries; jentry++) {
      tree->GetEntry(jentry);
      if (detectorId == 0) {
         h_1->Fill(energy);
        }
    }

   h_1->SetTitle("High resolution peak searching, number of iterations = 4");
   h_1->GetXaxis()->SetRange(1,nbins);
   TH1F *d_1 = new TH1F("d_1","",nbins,xmin_1,xmax_1);
   h_1->Draw("L");
   for (i = 0; i < nbins; i++) source[i]=h_1->GetBinContent(i + 1);
 
   TSpectrum *s_1 = new TSpectrum();
 
   nfound_1 = s_1->SearchHighRes(source, dest_1, nbins, sigma_1, threshold_1, kTRUE, deconiterations_1, kTRUE, averWindow_1);
   Double_t *xpeaks_1 = s_1->GetPositionX();
   for (i = 0; i < nfound_1; i++) {
      a=xpeaks_1[i];
      bin = 1 + Int_t(a + 0.5);
      fPositionX_1[i] = h_1->GetBinCenter(bin);
      fPositionY_1[i] = h_1->GetBinContent(bin);
   }
 
   TPolyMarker * pm_1 = (TPolyMarker*)h_1->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1) {
      h_1->GetListOfFunctions()->Remove(pm_1);
      delete pm_1;
   }
   pm_1 = new TPolyMarker(nfound_1, fPositionX_1, fPositionY_1);
   h_1->GetListOfFunctions()->Add(pm_1);
   pm_1->SetMarkerStyle(23);
   pm_1->SetMarkerColor(kRed);
   pm_1->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_1->SetBinContent(i + 1,dest_1[i]);
   d_1->SetLineColor(kRed);
   d_1->Draw("SAME");
   
   //second part stuff//
   TH1F *h_2 = new TH1F("hist2", "energy", nbins, xmin_2, xmax_2);

   for (int jentry=0; jentry<nentries; jentry++) {
      tree->GetEntry(jentry);
      if (detectorId == 0) {
         h_2->Fill(energy);
        }
    }

   h_2->SetTitle("High resolution peak searching, number of iterations = 2");
   h_2->GetXaxis()->SetRange(1,nbins);
   printf("hello");
   TH1F *d_2 = new TH1F("d_2","",nbins,xmin_2,xmax_2);
   h_2->Draw("L");
   for (i = 0; i < nbins; i++) source[i]=h_2->GetBinContent(i + 1);


   TSpectrum *s_2 = new TSpectrum();
 
   nfound_2 = s_2->SearchHighRes(source, dest_2, nbins, sigma_2, threshold_2, kTRUE, deconiterations_2, kTRUE, averWindow_2);
   Double_t *xpeaks_2 = s_2->GetPositionX();
   for (i = 0; i < nfound_2; i++) {
      a=xpeaks_2[i];
      bin = 1 + Int_t(a + 0.5);
      fPositionX_2[i] = h_2->GetBinCenter(bin);
      fPositionY_2[i] = h_2->GetBinContent(bin);
   }
 
   TPolyMarker * pm_2 = (TPolyMarker*)h_2->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2) {
      h_2->GetListOfFunctions()->Remove(pm_2);
      delete pm_2;
   }
   pm_2 = new TPolyMarker(nfound_2, fPositionX_2, fPositionY_2);
   h_2->GetListOfFunctions()->Add(pm_2);
   pm_2->SetMarkerStyle(23);
   pm_2->SetMarkerColor(kRed);
   pm_2->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_2->SetBinContent(i + 1,dest_2[i]);
   d_2->SetLineColor(kRed);
   d_2->Draw("SAME");
 
   printf("Reference dataset : Found %d candidate peaks in first part\n and %d in second part \n",nfound_1, nfound_2);
   for (i=0;i<nfound_1;i++) printf("posx= %f, posy= %f\n",fPositionX_1[i], fPositionY_1[i]);
   for (i=0;i<nfound_2;i++) printf("posx = %f, posy = %f\n", fPositionX_2[i], fPositionY_2[i]);

   
   // make single array of the found possitions
   for (i=0;i<nfound_1;i++) {
      refPositionX[i] = fPositionX_1[i];
      refPositionY[i] = fPositionY_1[i];
   }
   for (i=0;i<nfound_2;i++) {
      refPositionX[i] = fPositionX_2[i];
      refPositionY[i] = fPositionY_2[i];
   }
   printf("The positions for the 9 peaks of the reference dataset have been found. \n Calibrating other datasetss...\n");
   for(i=0; i<17;i++) printf("The positions of all peaks is: \n posx = %f    posy = %f \n", refPositionX[i], refPositionY[i]);
   
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

   //first part stuff
   TH1F *h_1_n = new TH1F("hist1_n", "energy", nbins, xmin_1, xmax_1);
    
   for (int jentry=0; jentry<nentries_n; jentry++) {
      tree_n->GetEntry(jentry);
      if (detectorId_n == 0) {
         h_1_n->Fill(energy_n);
        }
    }

   h_1_n->SetTitle("High resolution peak searching, number of iterations = 4");
   h_1_n->GetXaxis()->SetRange(1,nbins);
   TH1F *d_1_n = new TH1F("d_1_n","",nbins,xmin_1,xmax_1);
   h_1_n->Draw("L");
   for (i = 0; i < nbins; i++) source[i]=h_1_n->GetBinContent(i + 1);
 
   TSpectrum *s_1_n = new TSpectrum();
 
   nfound_1_n = s_1_n->SearchHighRes(source, dest_1_n, nbins, sigma_1, threshold_1, kTRUE, deconiterations_1, kTRUE, averWindow_1);
   Double_t *xpeaks_1_n = s_1_n->GetPositionX();
   for (i = 0; i < nfound_1_n; i++) {
      a_n=xpeaks_1_n[i];
      bin_n = 1 + Int_t(a_n + 0.5);
      fPositionX_1_n[i] = h_1_n->GetBinCenter(bin_n);
      fPositionY_1_n[i] = h_1_n->GetBinContent(bin_n);
   }
 
   TPolyMarker * pm_1_n = (TPolyMarker*)h_1_n->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1_n) {
      h_1_n->GetListOfFunctions()->Remove(pm_1_n);
      delete pm_1_n;
   }
   pm_1_n = new TPolyMarker(nfound_1_n, fPositionX_1_n, fPositionY_1_n);
   h_1_n->GetListOfFunctions()->Add(pm_1_n);
   pm_1_n->SetMarkerStyle(23);
   pm_1_n->SetMarkerColor(kRed);
   pm_1_n->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_1_n->SetBinContent(i + 1,dest_1_n[i]);
   d_1_n->SetLineColor(kRed);
   d_1_n->Draw("SAME");
   
   //second part stuff//
   TH1F *h_2_n = new TH1F("hist2_n", "energy", nbins, xmin_2, xmax_2);

   for (int jentry=0; jentry<nentries_n; jentry++) {
      tree_n->GetEntry(jentry);
      if (detectorId_n == 0) {
         h_2_n->Fill(energy_n);
        }
    }

   h_2_n->SetTitle("High resolution peak searching, number of iterations = 2");
   h_2_n->GetXaxis()->SetRange(1,nbins);
   printf("hello");
   TH1F *d_2_n = new TH1F("d_2_n","",nbins,xmin_2,xmax_2);
   h_2_n->Draw("L");
   for (i = 0; i < nbins; i++) source[i]=h_2_n->GetBinContent(i + 1);


   TSpectrum *s_2_n = new TSpectrum();
 
   nfound_2_n = s_2_n->SearchHighRes(source, dest_2_n, nbins, sigma_2, threshold_2, kTRUE, deconiterations_2, kTRUE, averWindow_2);
   Double_t *xpeaks_2_n = s_2_n->GetPositionX();
   for (i = 0; i < nfound_2_n; i++) {
      a_n=xpeaks_2_n[i];
      bin_n = 1 + Int_t(a_n + 0.5);
      fPositionX_2_n[i] = h_2_n->GetBinCenter(bin_n);
      fPositionY_2_n[i] = h_2_n->GetBinContent(bin_n);
   }
 
   TPolyMarker * pm_2_n = (TPolyMarker*)h_2_n->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2_n) {
      h_2_n->GetListOfFunctions()->Remove(pm_2_n);
      delete pm_2_n;
   }
   pm_2_n = new TPolyMarker(nfound_2_n, fPositionX_2_n, fPositionY_2_n);
   h_2_n->GetListOfFunctions()->Add(pm_2_n);
   pm_2_n->SetMarkerStyle(23);
   pm_2_n->SetMarkerColor(kRed);
   pm_2_n->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_2_n->SetBinContent(i + 1,dest_2_n[i]);
   d_2_n->SetLineColor(kRed);
   d_2_n->Draw("SAME");
 
   printf("Dataset n: Found %d candidate peaks in first part\n and %d in second part \n",nfound_1_n, nfound_2_n);
   for( i=0;i<nfound_1_n;i++) printf("posx= %f, posy= %f\n",fPositionX_1_n[i], fPositionY_1_n[i]);
   for (i=0;i<nfound_2_n;i++) printf("posx = %f, posy = %f\n", fPositionX_2_n[i], fPositionY_2_n[i]);
   

   // LESSGOOOO, now we put all these peaks

}