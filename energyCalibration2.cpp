    
    
    Int_t i, bin;
    Double_t binmax = binnum; 
    Double_t source[binnum], dest[binnum];
    Double_t fPositionX[100];
    Double_t fPositionY[100];
    Double_t a;

    Int_t pfound1 = s1->SearchHighRes(source, dest, binnum, 8, 2, kTRUE, 3, kTRUE, 3);
    Double_t *xpeaks = s1->GetPositionX();
    for (i = 0; i < pfound1; i++) {
        a=xpeaks[i];
        bin = 1 + Int_t(a + 0.5);
        fPositionX[i] = h1->GetBinCenter(bin);
        fPositionY[i] = h1->GetBinContent(bin);
    }
    TPolyMarker * pm = (TPolyMarker*)h1->GetListOfFunctions()->FindObject("TPolyMarker");
    if (pm) {
        h1->GetListOfFunctions()->Remove(pm);
        delete pm;
    }

    pm = new TPolyMarker(pfound1, fPositionX, fPositionY);
    h1->GetListOfFunctions()->Add(pm);
    pm->SetMarkerStyle(23);
    pm->SetMarkerColor(kRed);
    pm->SetMarkerSize(1.3);
 

    TH1F *d = new TH1F("d","",binnum,binmin,binmax);
    for (i = 0; i < binnum; i++) d->SetBinContent(i + 1,dest[i]);
    d->SetLineColor(kRed);
    d->Draw("SAME");
 
    printf("Found %d candidate peaks\n",pfound1);
    for( i=0;i<pfound1;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
}