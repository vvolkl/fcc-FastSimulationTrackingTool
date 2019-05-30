
///
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// \macro_code
///
/// \author Andreas Morsch, modifications to write root file by F. Bedeschi

#include "TSystem.h"
#include "TH1F.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

void WW(Int_t nev = 100, TString fname = "PyOut.root", Int_t ndeb = 1)
{
	// Load libraries
	gSystem->Load("libEG");
	gSystem->Load("libEGPythia8");
	// Histograms
	TH1F* etaMu = new TH1F("etaMu", "Pseudorapidity", 120, -3., 3.);
	TH1F* ptMu  = new TH1F("ptMu", "pt", 100, 0., 100.);
	TH1F* mInv  = new TH1F("mInv", "Invariant mass", 100, 0., 200.);
	// Setup File and Tree for output
	TFile *fout = new TFile(fname, "RECREATE");
	TTree *T = new TTree("tID", "Pythia particles");
	// Array of particles
	TClonesArray* particles = new TClonesArray("TParticle", 1000);
	// Link to TTree
	T->Branch("Particles", "TClonesArray", &particles, 64000, 0);
	// Create pythia8 object
	TPythia8* pythia8 = new TPythia8();
	// Configure
	pythia8->ReadConfigFile("config_ee_WW.cmd");
	//   
	// Initialize
	pythia8->Initialize(11 /* e- */, -11 /* e+ */, 240 /* GeV */);
	//
	//	Relevant particle codes
	//
	Int_t Muon = 13;		// Muon PDG code
	Int_t Z0 = 23;			// Z0   PDG code
	Double_t Mz = TDatabasePDG::Instance()->GetParticle(Z0)->Mass();
	Int_t Nsel = 0;			// Number of selected events
	// Event loop
	for (Int_t iev = 0; iev < nev; iev++) {
		pythia8->GenerateEvent();
		if (iev < ndeb) pythia8->EventListing();
		pythia8->ImportParticles(particles, "All");
		Int_t np = particles->GetEntriesFast();
		// Particle loop
		Int_t Nmu = 0;
		Int_t NmuZ = 0;
		Int_t Zmth = 0;
		const Int_t MaxMu = 20;
		TParticle* MuPart[MaxMu]; 
		for (Int_t ip = 0; ip < np; ip++) {
			TParticle* part = (TParticle*)particles->At(ip);
			Int_t ist = part->GetStatusCode();
			// Positive codes are final particles.
			if (ist <= 0) continue;
			Int_t pdg = part->GetPdgCode();
			Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
			if (charge == 0.) continue;
			Float_t eta = part->Eta();
			Float_t pt = part->Pt();

			if (TMath::Abs(pdg) == Muon) etaMu->Fill(eta);
			if (TMath::Abs(pdg) == Muon && pt > 0.) ptMu->Fill(pt);

			if (TMath::Abs(pdg) == Muon && Nmu < MaxMu)
			{
				MuPart[Nmu] = part;
				Nmu++;
			}
		}		// end particle loop
		//
		// Check if there are at least 2 muon pairs consistent with Z
		// Write out event
		Double_t RangeCut = 60.0;
		if (Nmu > 1)
		{
			Int_t Nz = 0;
			for (Int_t j1 = 0; j1 < Nmu - 1; j1++)		// Loop over all muon pair combinations
			{
				TParticle *part1 = MuPart[j1];
				Int_t pdg1 = part1->GetPdgCode();
				Float_t Q1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Charge() / 3.0;
				for (Int_t j2 = j1 + 1; j2 < Nmu; j2++)
				{
					TParticle *part2 = MuPart[j2];
					Int_t pdg2 = part2->GetPdgCode();
					Float_t Q2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Charge() / 3.0;
					if (Q1*Q2 < 0)						// Require opposite charge
					{
						TLorentzVector p1(0., 0., 0., 0.);
						p1.SetPxPyPzE(part1->Px(), part1->Py(), part1->Pz(), part1->Energy());
						TLorentzVector p2(0., 0., 0., 0.);
						p2.SetPxPyPzE(part2->Px(), part2->Py(), part2->Pz(), part2->Energy());
						TLorentzVector ptot = p1 + p2;	// 4-momentum sum of 2 muons
						Double_t Minv = ptot.M();
						mInv->Fill(Minv);
						if (TMath::Abs(Minv - Mz) < RangeCut) Nz++;
					}
				}
			}
			if (Nz > 0)	// Write out event if at least one candidate withing range
			{
				Nsel++;
				T->Fill();
				fout->Write();
			}
		}

	}	// End event loop
	//
	pythia8->PrintStatistics();
	//
	fout->Close();
	delete fout;
	//
	TCanvas* c1 = new TCanvas("c1", "Pythia8 Z->mu mu", 800, 800);
	c1->Divide(2, 2);
	c1->cd(1);
	//etaMu->Scale(5./Float_t(nev));
	etaMu->SetXTitle("#eta");
	etaMu->Draw();
	etaMu->SetYTitle("dN/d#eta");

	c1->cd(2);
	//gPad->SetLogy();
	//ptMu->Scale(5./Float_t(nev));
	ptMu->SetXTitle("p_{t} [GeV/c]");
	ptMu->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");
	ptMu->Draw();

	c1->cd(3);
	mInv->SetXTitle("M_{inv} [GeV/c^{2}]");
	mInv->SetYTitle("Di-muon invariant mass [GeV/c^{2}");
	mInv->Draw();

	cout << "Total number of generated events:" << nev << ", Selected: " << Nsel << endl;
}
