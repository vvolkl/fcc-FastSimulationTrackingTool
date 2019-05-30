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
/// \author Franco Bedeschi

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"

void ZmmWrite(Int_t nev  = 100, TString fname = "PyOut.root", Int_t ndeb = 1)
{
// Load libraries
   gSystem->Load("libEG");
   gSystem->Load("libEGPythia8");
// Histograms
   TH1F* etaMu = new TH1F("etaMu", "Pseudorapidity", 120, -3., 3.);
   TH1F* ptMu  = new TH1F("ptMu",  "pt",             100,   0., 50.);
   TH1F* Minv  = new TH1F("Minv",  "di-muon invariant mass",100,80.,100.);
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
   pythia8->ReadConfigFile("config_ee_z_mumu.cmd");
//
// Initialize

   pythia8->Initialize(11 /* e- */, -11 /* e+ */, 91.187 /* GeV */);

// Event loop
   for (Int_t iev = 0; iev < nev; iev++) {
      pythia8->GenerateEvent();
      if (iev < ndeb) pythia8->EventListing();
      pythia8->ImportParticles(particles,"All");
      T->Fill();
      fout->Write();
      Int_t np = particles->GetEntriesFast();
// Particle loop
      TVectorD p1(4);
      TVectorD p2(4);
      Int_t Nmu = 0;
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

         etaMu->Fill(eta);
         if (pt > 0.) ptMu->Fill(pt);
         if(TMath::Abs(pdg) == 13)
         {
                Nmu++;
                if(Nmu == 1)
                {
                        p1(0) = part->Px();
                        p1(1) = part->Py();
                        p1(2) = part->Pz();
                        p1(3) = part->Energy();
                }
                if(Nmu == 2)
                {
                        p2(0) = part->Px();
                        p2(1) = part->Py();
                        p2(2) = part->Pz();
                        p2(3) = part->Energy();
                }

         }
      }
      cout<<"P1 = "<<p1(0)<<", "<<p1(1)<<", "<<p1(2)<<", "<<p1(3)<<endl;
      cout<<"P1 = "<<p2(0)<<", "<<p2(1)<<", "<<p2(2)<<", "<<p2(3)<<endl;

      Float_t Mi = TMath::Sqrt(pow(p1(3)+p2(3),2)-pow(p1(0)+p2(0),2)
                                                 -pow(p1(1)+p2(1),2)
                                                 -pow(p1(2)+p2(2),2));
      cout<<"Invariant Mass is: "<<Mi<<" GeV"<<endl;
      Minv->Fill(Mi);
   }

   pythia8->PrintStatistics();

   fout->Close();
   delete fout;

   TCanvas* c1 = new TCanvas("c1","Pythia8 Z->mu mu",800,800);
   c1->Divide(2, 2);
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

   c1->cd(3);
   // Invariant mass plot
   Minv->SetXTitle("#m^{+}#m^{-} invariant mass (GeV)");
   Minv->Draw();
 }
