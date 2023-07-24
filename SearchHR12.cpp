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
   Double_t fPositionX[100];
   Double_t fPositionY[100];
   Double_t refPositionX[9];
   Double_t refPositionY[9];
   Int_t fNPeaks = 0;
   Int_t i,nfound_1, nfound_2,bin;
   const Int_t nbins = 4000;
   Double_t xmin_1 = 0, xmin_2 = 10500; //xmin_1 is for the first 4 peaks and xmin_2 for the last 5
   Double_t xmax_1= nbins, xmax_2= nbins+xmin_2; 
   Double_t a;
   Double_t source[nbins], dest_1[nbins], dest_2[nbins];
   Double_t sigma_1 = 12, sigma_2 = 19, threshold_1 = 0.8, threshold_2 = 0.9;
   Int_t deconiterations_1 = 4, deconiterations_2 = 2, averWindow_1 = 2, averWindow_2 = 7;

   // Dataset 1 - the standard
   {
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
      fPositionX[i] = h_1->GetBinCenter(bin);
      fPositionY[i] = h_1->GetBinContent(bin);
   }
 
   TPolyMarker * pm_1 = (TPolyMarker*)h_1->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1) {
      h_1->GetListOfFunctions()->Remove(pm_1);
      delete pm_1;
   }
   pm_1 = new TPolyMarker(nfound_1, fPositionX, fPositionY);
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
      fPositionX[i] = h_2->GetBinCenter(bin);
      fPositionY[i] = h_2->GetBinContent(bin);
   }
 
   TPolyMarker * pm_2 = (TPolyMarker*)h_2->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2) {
      h_2->GetListOfFunctions()->Remove(pm_2);
      delete pm_2;
   }
   pm_2 = new TPolyMarker(nfound_2, fPositionX, fPositionY);
   h_2->GetListOfFunctions()->Add(pm_2);
   pm_2->SetMarkerStyle(23);
   pm_2->SetMarkerColor(kRed);
   pm_2->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_2->SetBinContent(i + 1,dest_2[i]);
   d_2->SetLineColor(kRed);
   d_2->Draw("SAME");
 
   printf("Reference dataset : Found %d candidate peaks in first part\n and %d in second part \n",nfound_1, nfound_2);
   for( i=0;i<nfound_1;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
   for (i=0;i<nfound_2;i++) printf("posx = %f, posy = %f\n", fPositionX[i], fPositionY[i]);

   
   // make single array of the found possitions
   for (i=0;i<nfound_1;i++) {
      refPositionX[i] = fPositionX[i];
      refPositionY[i] = fPositionY[i];
   }
   for (i=0;i<nfound_2;i++) {
      refPositionX[i] = fPositionX[i];
      refPositionY[i] = fPositionY[i];
   }
   printf("The positions for the 9 peaks of the reference dataset have been found. \n Calibrating other datasetss...");
   
   delete h_1;
   delete h_2;
   delete s_1;
   delete s_2;
   delete pm_1;
   delete pm_2;
   delete d_1;
   delete d_2;
   }
   
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   
   // now that we have the "standard" (i.e. the reference peakpoints for calibration) we can import any other file 
   // from the folder and proceed to try and calibrate them


   // Some dataset n
   {

   TFile *fn = TFile::Open("data/deuterated-digitizer_3.root");
   fn->GetObject("events", tree_n);

   tree_n->SetBranchAddress("energy", &energy_n);
   tree_n->SetBranchAddress("detectorId", &detectorId_n);

   Long64_t nentries = tree_n->GetEntries();

   //first part stuff
   TH1F *h_1 = new TH1F("hist1", "energy", nbins, xmin_1, xmax_1);
    
   for (int jentry=0; jentry<nentries; jentry++) {
      tree_n->GetEntry(jentry);
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
      fPositionX[i] = h_1->GetBinCenter(bin);
      fPositionY[i] = h_1->GetBinContent(bin);
   }
 
   TPolyMarker * pm_1 = (TPolyMarker*)h_1->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_1) {
      h_1->GetListOfFunctions()->Remove(pm_1);
      delete pm_1;
   }
   pm_1 = new TPolyMarker(nfound_1, fPositionX, fPositionY);
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
      tree_n->GetEntry(jentry);
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
      fPositionX[i] = h_2->GetBinCenter(bin);
      fPositionY[i] = h_2->GetBinContent(bin);
   }
 
   TPolyMarker * pm_2 = (TPolyMarker*)h_2->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm_2) {
      h_2->GetListOfFunctions()->Remove(pm_2);
      delete pm_2;
   }
   pm_2 = new TPolyMarker(nfound_2, fPositionX, fPositionY);
   h_2->GetListOfFunctions()->Add(pm_2);
   pm_2->SetMarkerStyle(23);
   pm_2->SetMarkerColor(kRed);
   pm_2->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d_2->SetBinContent(i + 1,dest_2[i]);
   d_2->SetLineColor(kRed);
   d_2->Draw("SAME");
 
   printf("Dataset n: Found %d candidate peaks in first part\n and %d in second part \n",nfound_1, nfound_2);
   for( i=0;i<nfound_1;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
   for (i=0;i<nfound_2;i++) printf("posx = %f, posy = %f\n", fPositionX[i], fPositionY[i]);
   }

}