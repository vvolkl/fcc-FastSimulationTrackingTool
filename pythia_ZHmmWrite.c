/// \file
/// \ingroup tutorial_pythia
/// pythia8 basic example
///
/// to run, do:
///
/// ~~~{.cpp}
///  root > .x pythia8.C
/// ~~~
///
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// \macro_code
///
/// \author Andreas Morsch, modifications to write root file be F. Bedeschi

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
#include <iostream>

void ZHmm(Int_t nev  = 100, TString fname="PyOut.root",Int_t ndeb = 1)
{
// Load libraries
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");
// Histograms
   TH1F* etaMu = new TH1F("etaMu", "Pseudorapidity", 120, -3., 3.);
   TH1F* ptMu  = new TH1F("ptMu",  "pt",             100,   0., 100.);
// Setup File and Tree for output
   TFile *fout = new TFile(fname,"RECREATE");
   TTree *T    = new TTree("tID","Pythia particles");
// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 1000);
// Link to TTree
   T->Branch("Particles","TClonesArray",&particles,64000,0);
// Create pythia8 object
   TPythia8* pythia8 = new TPythia8();

// Configure
   pythia8->ReadConfigFile("config_ee_zh_zmumu.cmd");
//   
// Initialize

   pythia8->Initialize(11 /* e- */, -11 /* e+ */, 240 /* GeV */);

// Event loop
   for (Int_t iev = 0; iev < nev; iev++) {
      pythia8->GenerateEvent();
      if (iev < ndeb) pythia8->EventListing();
      pythia8->ImportParticles(particles,"All");
      T->Fill();
      fout->Write();
      Int_t np = particles->GetEntriesFast();
// Particle loop
      Int_t Nmu = 0;
      Int_t NmuZ = 0;
      Int_t Zmth = 0;
      for (Int_t ip = 0; ip < np; ip++) {
         TParticle* part = (TParticle*) particles->At(ip);
         Int_t ist = part->GetStatusCode();
         // Positive codes are final particles.
         if (ist <= 0) continue;
         Int_t pdg = part->GetPdgCode();
         Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
         if (charge == 0.) continue;
         Float_t eta = part->Eta();
         Float_t pt  = part->Pt();

         //if (TMath::Abs(pdg) == 13) etaMu->Fill(eta);
         //if (TMath::Abs(pdg) == 13 && pt > 0.) ptMu->Fill(pt);

	 if(TMath::Abs(pdg) == 13)
	 {
	   Nmu++;
	   Int_t Mother =  part->GetFirstMother();
	   while(Mother > 0)
	   {
	     TParticle* mpart = (TParticle*)particles->At(Mother);
	     Int_t mpdg = mpart->GetPdgCode();
	     if(mpdg == 23)
	     {
	       Zmth++;
	       etaMu->Fill(eta);
	       ptMu->Fill(pt);
	       break;
	     }
	     Mother =  mpart->GetFirstMother();
	   }
	 }

      }
      if(Nmu == 0)
      {
        cout<<"No muons found in event!!!"<<endl;
        pythia8->EventListing();
      }
      else if(Nmu != Zmth)
      {
        //cout<<"At least one muon not from Z"<<endl;
	//pythia8->EventListing();
      }
   }

   pythia8->PrintStatistics();

   fout->Close();
   delete fout;

   TCanvas* c1 = new TCanvas("c1","Pythia8 Z->mu mu",800,800);
   c1->Divide(1, 2);
   c1->cd(1);
   //etaMu->Scale(5./Float_t(nev));
   etaMu->Draw();
   etaMu->SetXTitle("#eta");
   etaMu->SetYTitle("dN/d#eta");

   c1->cd(2);
   //gPad->SetLogy();
   //ptMu->Scale(5./Float_t(nev));
   ptMu->Draw();
   ptMu->SetXTitle("p_{t} [GeV/c]");
   ptMu->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");
 }
