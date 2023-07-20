import ROOT
import uproot
import matplotlib.pyplot as plt
import numpy as np
from ROOT import gROOT

gROOT.Reset()

myFile = ROOT.TFile.Open("data/deuterated-digitizer_1.root")
myTree = myFile.Get('events')

energy = np.zeros(1, dtype=np.float64)
detector = np.zeros(1, dtype=np.ulonglong)

myTree.SetBranchAddress('energy', energy)
myTree.SetBranchAddress('detectorId', detector)



histograms = {}

binnum = 2**14
binmin = 0
binmax = 2**14

ROOT.TH1.AddDirectory(False)
for entry, i in range(myTree.GetEntries()//1000, detector:
    myTree.GetEntry(entry)
    
    if detector[0] not in histograms:
        histograms[detector[0]] = ROOT.TH1F(f"hist_{detector[0]}", f"Energy Distribution - Detector {detector[0]}", binnum, binmin, binmax)
        
    histograms[detector[0]].Fill(energy[0])


if "canvas" in ROOT.gROOT.GetListOfCanvases():
    ROOT.gROOT.GetListOfCanvases().FindObject("canvas").Close()

canvas = ROOT.TCanvas("canvas", "Histograms", 800, 600)
canvas.Divide(5, 1)

pad_index = 1
for detector_id, histogram in histograms.items():
    canvas.cd(pad_index)  # Activate the pad at the current index
    histogram.Draw()
    histogram.GetXaxis().SetTitle("Energy")
    histogram.GetYaxis().SetTitle("Frequency")
    pad_index += 1

canvas.Update()
canvas.Draw()

myFile.Close()