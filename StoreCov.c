#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <TTree.h>
#include <TString.h>
#include "SolGeom.h"
#include "SolTrack.h"
#include "SolGridCov.h"
//
//
//
void StoreCov(TString fname)
{
	//
	// Initialize geometry
	// 
	Bool_t Res = kTRUE;	// Resolution effects enabled
	Bool_t MS  = kTRUE;	// Multiple scattering effects enabled
	SolGeom *G;				// Initialize geometry
	const Int_t nDet = 9;
	Bool_t OK[nDet] = {		// Enable selected parts of the detector 
		1,					// Beam pipe
		1,					// Inner VTX pixel layers
		1,					// Outer VTX layers
		1,					// Drift chamber
		1,					// Barrel Si wrapper
		1,					// Barrel pre-shower
		1,					// Forw. VTX pixel layers
		1,					// Forw. Si wrapper
		1 };				// Forw. pre-shower
	G = new SolGeom(OK);	// Geometry with selected detectors
	// Write covariance grid to root file
	SolGridCov *GC = new SolGridCov();
	GC->Write(fname, G);
}
//
void test(Double_t ptr, Double_t angr)
{
	//StoreCov();
	//
	// Test interpolation over the grid
	//
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA-BASE.root");
	//
	// Compare to direct calculation
	//
	SolGeom *G = new SolGeom();
	Bool_t Res = kTRUE; Bool_t MS = kTRUE;
	Double_t th = TMath::Pi()*angr / 180.; 
	Double_t x[3]; Double_t p[3];
	x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
	p[0] = ptr; p[1] = 0; p[2] = ptr / TMath::Tan(th);
	SolTrack *tr = new 	SolTrack(x, p, G);	// Initialize track
	tr->CovCalc(Res, MS);				// Calculate covariance
	TMatrixDSym Cdirect = tr->Cov();
	cout << "Direct calculation: pt = "<< ptr
		<< ", ang = "<<angr<<endl; Cdirect.Print();
	TMatrixDSym Cv = GC->GetCovL(ptr, angr);	// Get covariance interpolation
	cout << "Interpolated matrix" << endl; Cv.Print();
	TMatrixDSym Diff = Cv- Cdirect;		// Get covariance difference
	cout << "Difference between interpolation and direct calculation:" << endl;
	//Diff.Print();
	TMatrixDSym DiffN(5);
	for (Int_t i = 0; i < 5; i++)
	{
		for (Int_t j = 0; j < 5; j++)
		{
			DiffN(i, j) = Diff(i, j) / (0.5*(Cv(i, j) + Cdirect(i, j)));
		}
	}
	DiffN.Print();
	//
}
