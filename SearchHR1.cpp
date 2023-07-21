#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TTree.h"


void SearchHR1() {
   TTree *tree;
   UShort_t detectorId;
   UShort_t energy;
   Double_t fPositionX[100];
   Double_t fPositionY[100];
   Int_t fNPeaks = 0;
   Int_t i,nfound,bin;
   const Int_t nbins = 16384;
   Double_t xmin     = 0;
   Double_t xmax     = nbins;
   Double_t a;
   Double_t source[nbins], dest[nbins];
   Double_t sigma = 10, threshold = 0.78;
   Int_t deconiterations = 4, averWindow = 2;
   // gROOT->ForceStyle();
 
   // import and pointing
   // TString dir  = gROOT->GetTutorialDir();
   // TString file = dir+"/spectrum/TSpectrum.root";
   // TFile *f     = new TFile(file.Data());
   // TH1F *h = (TH1F*) f->Get("back2");

   //
   // TString file1 = ;
   TFile *f1 = TFile::Open("data/deuterated-digitizer_1.root");
   f1->GetObject("events", tree);

   tree->SetBranchAddress("energy", &energy);
   tree->SetBranchAddress("detectorId", &detectorId); 

   Long64_t nentries = tree->GetEntries();
   TH1F *h = new TH1F("hist1", "energy", nbins, xmin, xmax);
   
   // TH1F *h = (TH1F*) f1->GetObject("events", tree);  //this imports a histogtam from tree; we don't have one
   for (int jentry=0; jentry<nentries; jentry++) {
      tree->GetEntry(jentry);
      if (detectorId == 0) {
         h->Fill(energy);
        }
    }

   //

   h->SetTitle("High resolution peak searching, number of iterations = 3");
   h->GetXaxis()->SetRange(1,nbins);
   printf("hello");
   TH1F *d = new TH1F("d","",nbins,xmin,xmax);
   h->Draw("L");
 
   for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
 
   h->Draw("L");
 
   TSpectrum *s = new TSpectrum();
 
   nfound = s->SearchHighRes(source, dest, nbins, sigma, threshold, kTRUE, deconiterations, kTRUE, averWindow);
   Double_t *xpeaks = s->GetPositionX();
   for (i = 0; i < nfound; i++) {
      a=xpeaks[i];
      bin = 1 + Int_t(a + 0.5);
      fPositionX[i] = h->GetBinCenter(bin);
      fPositionY[i] = h->GetBinContent(bin);
   }
 
   TPolyMarker * pm = (TPolyMarker*)h->GetListOfFunctions()->FindObject("TPolyMarker");
   if (pm) {
      h->GetListOfFunctions()->Remove(pm);
      delete pm;
   }
   pm = new TPolyMarker(nfound, fPositionX, fPositionY);
   h->GetListOfFunctions()->Add(pm);
   pm->SetMarkerStyle(23);
   pm->SetMarkerColor(kRed);
   pm->SetMarkerSize(1.3);
 
   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,dest[i]);
   d->SetLineColor(kRed);
   d->Draw("SAME");
 
   printf("Found %d candidate peaks\n",nfound);
   for( i=0;i<nfound;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
}