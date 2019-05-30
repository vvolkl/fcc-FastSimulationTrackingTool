#include <TString.h>
#include <TROOT.h>

void LoadAll(TString dname)
{
gROOT->Reset();
TString Action = ".L SolGeom" + dname + ".cxx+";
gROOT->ProcessLine(Action);
gROOT->ProcessLine(".L SolTrack.cxx+");
gROOT->ProcessLine(".L SolGridCov.cxx+");
gROOT->ProcessLine(".L ObsTrk.cxx+");
gROOT->ProcessLine(".L CompRes.c+");
}
