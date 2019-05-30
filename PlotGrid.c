#include <TMath.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include "SolGeom.h"
#include "SolTrack.h"

//
void PlotGrid()
{
	//
	// File title encoding/decoding
	//
	const Int_t nPar = 5;
	const Int_t nAng = 8;
	TString sPar[nPar] = { "d0", "Phi", "Pt", "Theta", "Z" };
	TString sAng[nAng] = { "10deg", "20deg", "30deg", "40deg",
						   "45deg", "60deg", "75deg", "90deg" };
	//
	// Output files
	TFile *OutFile = new TFile("./dump/graphs.root", "RECREATE");
	//
	// graphs
	TGraphErrors *gr[nPar][nAng];
	//
	// Loop on all files
	for (Int_t ip = 0; ip < nPar; ip++)
	{
		for (Int_t ja = 0; ja < nAng; ja++)
		{
			TString fname = sPar[ip] + "-"+ sAng[ja];
			TString gname = "g" + fname;
//	Read file
			float pt ; float spt; 
			float res; float sres;
			char strng[100]; int nbytes = 100;
			Double_t  x[200]; Double_t  y[200];
			Double_t sx[200]; Double_t sy[200];
			TString name = "./dump/" + fname;
			FILE *fdata = fopen(name, "r");
			if (!fdata)
			{
				cout << "SolGeom::GeoRead - can't open input file" << endl;
				return;
			}
			Int_t ndat = 0;
			while (fgets(strng, nbytes, fdata) != NULL)
			{
				//cout << strng;
				char hash[2] = "#";
				if (strncmp(strng,hash,1) == 0)cout << "Processing file: " << fname << endl;
				else
				{
					int status = sscanf(strng, "%g %g %g %g",
						&pt, &spt, &res, &sres);
					if (pt >= 1.0)
					{ 
						 x[ndat]  = pt;  y[ndat] = res;
						sx[ndat] = spt; sy[ndat] = sres;
						ndat++; 
					}
				}
			}
			gr[ip][ja] = new TGraphErrors(ndat, x, y, sx, sy);
			gr[ip][ja]->SetName(gname);
			if(OutFile->IsOpen())gr[ip][ja]->Write();
		}
	}
	OutFile->Close();
}