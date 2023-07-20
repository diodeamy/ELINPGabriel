from ROOT import TFile, TH1D, gPad, TCanvas

file = TFile("/data/deuterated-digitizer_1.root")

tree = file.Get("events")

hist1 = TH1D("hist1", "hist1", 2**14, 0, 2**14)
hist2 = TH1D("hist2", "hist2", 2**14, 0, 2**14)

cnt = 0

for ev in tree:
    if ev.detectorId == 0:
        
        hist1.Fill(ev.energy)

        cnt += 1

        if cnt > 100000:
            break

c1 = TCanvas("c1", "c1", 800, 600)
hist1.Draw("hist")

c2 = TCanvas("c2", "c2", 800, 600)
hist2.Draw("hist")

gPad.WaitPrimitive("ggg")

