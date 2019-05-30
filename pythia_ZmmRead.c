
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// Description:
/// Read in Z-> mu mu events and calculate invariant mass
///
/// \author: Franco Bedeschi

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <iostream>
#include "SolGeom.h"
#include "SolGridCov.h"
#include "ObsTrk.h"

void ZmmRead(Int_t nev  = 100, TString fname = "PyOut.root", Int_t ndeb = 1)
{
// Initialize Geometry
	SolGeom *G = new SolGeom();		// Initialize geometry
	Double_t Bfield = G->B();		// Get B field in Tesla
// Initialize tracking resolution
	SolGridCov *GC = new SolGridCov();
	GC->Read("Cov.root");			// Read in covariance array
// Histograms
   TH1F* etaMu = new TH1F("etaMu", "Pseudorapidity", 120, -3., 3.);
   TH1F* ptMu  = new TH1F("ptMu",  "pt",             100,   0., 50.);
   TH1F* gMinv = new TH1F("Minv", "Generated di-muon invariant mass", 100, 85., 95.);
   TH1F* oMinv = new TH1F("Minv", "Observed  di-muon invariant mass", 100, 85., 95.);
// Setup File and Tree for input
   TFile *fin = new TFile(fname,"READ");
   TTree *T = (TTree*)fin->Get("tID");
// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 1000);
// Link to TTree
   T->SetBranchAddress("Particles", &particles);
   Int_t MaxEv = T->GetEntries();
   Int_t Nevents = TMath::Min(nev, MaxEv);
// Event loop
   for (Int_t iev = 0; iev < Nevents; iev++) {
	  Int_t nRead = T->GetEntry(iev);
	  if (nRead <= 0) continue;
      Int_t np = particles->GetEntriesFast();
// Particle loop
      TLorentzVector p1(0.,0.,0.,0.);
	  TLorentzVector p2(0., 0., 0., 0.);
	  Double_t Q1 = 0; Double_t M1 = 0;
	  Double_t Q2 = 0; Double_t M2 = 0;
      Int_t Nmu = 0;
      for (Int_t ip = 0; ip < np; ip++) 
	  {
         TParticle* part = (TParticle*) particles->At(ip);
         Int_t ist = part->GetStatusCode();
         // Positive codes are final particles.
         if (ist <= 0) continue;
         Int_t pdg = part->GetPdgCode();
		 Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3.0;
		 Float_t mass   = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
         if (charge == 0.) continue;
		 Int_t Mother1 = part->GetFirstMother();
		 TParticle* mpart1 = (TParticle*)particles->At(Mother1);
		 Int_t mpdg1 = mpart1->GetPdgCode();
		 Int_t Mother2 = part->GetSecondMother();
		 Int_t mpdg2 = 0;
		 if (Mother2 >= 0)
		 {
			 TParticle* mpart2 = (TParticle*)particles->At(Mother2);
			 Int_t mpdg2 = mpart2->GetPdgCode();
		 }
         Float_t eta = part->Eta();
         Float_t pt  = part->Pt();
		 if (mpdg1 == 23){
			 etaMu->Fill(eta);
			 if (pt > 0.) ptMu->Fill(pt);
		 }
         if(TMath::Abs(pdg) == 13 && mpdg1 == 23)		// Only muons from Z0
         {
                Nmu++;
				if (Nmu == 1)
				{
					Q1 = charge; M1 = mass;
					p1.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
				}
				if (Nmu == 2)
				{
					Q2 = charge; M2 = mass;
					p2.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
				}
         }
      }
	  TLorentzVector ptot = p1 + p2;
      Float_t Mi = ptot.M();			// Generated mu mu invariant mass
      gMinv->Fill(Mi);
	  //
	  // Now apply track resolution effects
	  //
	  TVector3 tX(0.0, 0.0, 0.0);	// Set origin to (0,0,0)
	  TVector3 tP1 = p1.Vect();		// Get generated momenta
	  TVector3 tP2 = p2.Vect();	
	  // Acceptance above 10 degrees and pt>1 GeV
	  Double_t AngMin = 10.0;
	  Double_t th1 = TMath::ACos(TMath::Abs(tP1.CosTheta()))*180. / TMath::Pi();
	  Double_t th2 = TMath::ACos(TMath::Abs(tP2.CosTheta()))*180. / TMath::Pi();
	  if (th1 < AngMin || th2 < AngMin)continue;
	  if (tP1.Pt() < 1. || tP2.Pt() < 1.)continue;	// Pt cuts
	  //
	  ObsTrk *Tr1 = new ObsTrk(tX, tP1, Q1, Bfield, GC);	// Apply track resolution
	  ObsTrk *Tr2 = new ObsTrk(tX, tP2, Q2, Bfield, GC);
	  TVector3 obsP1 = Tr1->GetObsP();				// Get smeared momenta
	  TVector3 obsP2 = Tr2->GetObsP();
	  Double_t E1 = TMath::Sqrt(M1*M1 + obsP1.Mag2());	// Smeared energies
	  Double_t E2 = TMath::Sqrt(M2*M2 + obsP2.Mag2());
	  TLorentzVector oP1(obsP1, E1);	// Fill smeared Lorentz vectors
	  TLorentzVector oP2(obsP2, E2);
	  TLorentzVector oPtot = oP1 + oP2;	// Total momentum 4-vector
	  Float_t MiObs = oPtot.M();		// Invariant mass
	  oMinv->Fill(MiObs);				// Fill histogram
   }
   //
   fin->Close();
   delete fin;
   //
   // Plot histograms
   TCanvas* c1 = new TCanvas("c1","Pythia8 Z->mu mu",10,10,800,800);
   c1->Divide(2, 2);
   c1->cd(1);
   //etaMu->Scale(5./Float_t(nev));
   etaMu->Draw();
   etaMu->SetXTitle("#eta");
   etaMu->SetYTitle("dN/d#eta");
   //
   c1->cd(2);
   //gPad->SetLogy();
   //ptMu->Scale(5./Float_t(nev));
   ptMu->Draw();
   ptMu->SetXTitle("p_{t} [GeV/c]");
   ptMu->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");
   //
   c1->cd(3);
   // Invariant mass plot
   gMinv->SetXTitle("#m^{+}#m^{-} invariant mass (GeV)");
   gMinv->Draw();
   //oMinv->SetXTitle("#m^{+}#m^{-} invariant mass (GeV)");
   //oMinv->SetLineColor(kRed);
   //oMinv->Draw("SAME");
   c1->cd(4);
   // Invariant mass plot
   oMinv->SetXTitle("#m^{+}#m^{-} invariant mass (GeV)");
   oMinv->Draw();
 }
