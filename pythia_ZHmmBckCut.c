
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// Description:
/// Read in HZ-> mu mu events and calculate recoil mass
///
/// \author: Franco Bedeschi

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
//#include "SolGeom.h"
#include "SolGridCov.h"
#include "ObsTrk.h"

void ZHmmBckCut(Int_t nev  = 100, TString fname = "ZHmm.root", Int_t ndeb = 1)
{
//
//	Relevant particle codes
//
	Int_t Muon = 13;		// Muon PDG code
	Int_t Z0   = 23;		// Z0   PDG code
	Double_t Mz = TDatabasePDG::Instance()->GetParticle(Z0)->Mass();
// Initialize B field
	//SolGeom *G = new SolGeom();	// Initialize geometry
	//Double_t Bfield = G->B();		// Get B field in Tesla
	Double_t Bfield = 2.0;			// Set B field in Tesla
// Initialize tracking resolution
	SolGridCov *GCid = new SolGridCov();
	GCid->Read("CovIDEA-BASE.root");			// Read in covariance array
	SolGridCov *GCcl = new SolGridCov();
	GCcl->Read("CovCLD.root");			// Read in covariance array
// Histograms
   TH1D* etaMu  = new TH1D("etaMu",  "Pseudorapidity",                             120,  -3.,   3.);
   TH1D* ptMu   = new TH1D("ptMu",   "pt",                                         100,   0., 100.);
   TH1D* gMinv  = new TH1D("Minv",   "Generated di-muon invariant mass",           100,  70., 110.);
   TH1D* oMinv_id  = new TH1D("Minv_id", "IDEA: Observed  di-muon invariant mass", 100, 70., 110.);
   TH1D* oMinv_cl  = new TH1D("Minv_cl", "CLD:  Observed  di-muon invariant mass", 100, 70., 110.);
   TH1D *gHrec  = new TH1D("gHrec",  "Generated Higgs recoil mass",                200, 50., 150.);
   TH1D *eHrec  = new TH1D("eHrec",  "Generated Higgs recoil mass (E smear)",      200, 50., 150.);
   TH1D *oHrec_id  = new TH1D("oHrec_id", "IDEA: Observed  Higgs recoil mass (No E smear)",   200, 50., 150.);
   TH1D *oHrec_cl  = new TH1D("oHrec_cl", "CLD:  Observed  Higgs recoil mass (No E smear)", 200, 50., 150.);
   TH1D *eoHrec_id = new TH1D("eoHrec_id", "IDEA: Observed  Higgs recoil mass (With E smear)", 200, 50., 150.);
   TH1D *eoHrec_cl = new TH1D("eoHrec_cl", "CLD: Observed  Higgs recoil mass (With E smear)", 200, 50., 150.);
   TH1D *gCos = new TH1D("gCos", "Generated cos of muon opening angle", 100, -1.0, 1.0);
   TH1D *oCos_id = new TH1D("oCos_id", "IDEA: cos of muon opening angle", 100, -1.0, 1.0);
   TH1D *oCos_cl = new TH1D("oCos_cl", "CLD:  cos of muon opening angle", 100, -1.0, 1.0);
   TH2D *gpCos = new TH2D("gpCos", "Generated Cos vs pair p", 100, 0., 100., 50, -1.0, 1.0);

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
   for (Int_t iev = 0; iev < Nevents; iev++)	// Main event loop
   {
	  Int_t nRead = T->GetEntry(iev);
	  if (nRead <= 0) continue;
      Int_t np = particles->GetEntriesFast();
// Particle loop
	  TLorentzVector p1(0., 0., 0., 0.);
	  TLorentzVector p2(0., 0., 0., 0.);
	  Double_t Q1 = 0; Double_t M1 = 0;
	  Double_t Q2 = 0; Double_t M2 = 0;
      Int_t Nmu = 0;
	  Int_t NmuZ = 0;
	  Int_t Zmth = 0;
	  const Int_t MaxMu = 20;
	  TParticle* MuPart[MaxMu];
      for (Int_t ip = 0; ip < np; ip++)			// Main loop on event particles
	  {
         TParticle* part = (TParticle*) particles->At(ip);
         Int_t ist = part->GetStatusCode();
         // Positive codes are final particles.
         if (ist <= 0) continue;
         Int_t pdg = part->GetPdgCode();
		 Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3.0;
		 Float_t mass   = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
         if (charge == 0.) continue;
		 // Select muons pairs closest to Z0
		 Float_t eta = part->Eta();
		 Float_t pt = part->Pt();

		 if (TMath::Abs(pdg) == Muon) etaMu->Fill(eta);
		 if (TMath::Abs(pdg) == Muon && pt > 0.) ptMu->Fill(pt);

		 if (TMath::Abs(pdg) == Muon && Nmu < MaxMu)
		 {
			 MuPart[Nmu] = part;
			 Nmu++;
		 }
      }			// End main loop on particles in the event
	  //
	  Double_t RangeCut = 20.0;
	  Int_t Nz = 0;
	  if (Nmu >1)	// Found at least two muons. See if close enough to Z and select best pair
	  {
		  Double_t diffM = 1000;
		  for (Int_t j1 = 0; j1 < Nmu - 1; j1++)		// Loop over all muon pair combinations
		  {
			  TParticle *part1 = MuPart[j1];
			  Int_t pdg1 = part1->GetPdgCode();
			  Q1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Charge() / 3.0;
			  for (Int_t j2 = j1 + 1; j2 < Nmu; j2++)
			  {
				  TParticle *part2 = MuPart[j2];
				  Int_t pdg2 = part2->GetPdgCode();
				  Q2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Charge() / 3.0;
				  if (Q1*Q2 < 0)						// Require opposite charge
				  {
					  TLorentzVector p_1(0., 0., 0., 0.);
					  p_1.SetPxPyPzE(part1->Px(), part1->Py(), part1->Pz(), part1->Energy());
					  TLorentzVector p_2(0., 0., 0., 0.);
					  p_2.SetPxPyPzE(part2->Px(), part2->Py(), part2->Pz(), part2->Energy());
					  TLorentzVector p_tot = p_1 + p_2;	// 4-momentum sum of 2 muons
					  Double_t Minv = p_tot.M();
					  Double_t diff = TMath::Abs(Minv - Mz);
					  if (diff < RangeCut)
					  {
						  Nz++;
						  if (diff < diffM)
						  {
							  diffM = diff;
							  p1 = p_1; p2 = p_2;
						  }
					  }
				  }
			  }
		  }
		  if (Nz > 0)
		  {
			  TLorentzVector ptot = p1 + p2;
			  Double_t Mi = ptot.M();						// Generated mu mu invariant mass
			  gMinv->Fill(Mi);								// Fill generated invariant mass histogram
			  Double_t rs = 240;	Double_t sEoE = 1.92e-3 / TMath::Sqrt(2);		// CoM energy and relative energy spread
			  Double_t s = rs * rs;
			  Double_t rss = gRandom->Gaus(rs, rs*sEoE);		// Smear CoM energy for energy spread
			  Double_t ss = rss * rss;
			  Double_t gHrc = TMath::Sqrt(s + Mi * Mi - 2 * rs *ptot.E());	// Generated recoil mass
			  Double_t eHrc = TMath::Sqrt(ss + Mi * Mi - 2 * rss*ptot.E());	// Generated recoil mass w/ E spread
			  gHrec->Fill(gHrc);		// Fill generated recoil mass histogram
			  //
			  // Now apply track resolution effects
			  //
			  TVector3 tP1 = p1.Vect();		// Get generated momentum first  muon
			  TVector3 tP2 = p2.Vect();		// Get generated momentum second muon
			  Double_t gCs = tP1.Dot(tP2) / (tP1.Mag()*tP2.Mag()); // Cos of muon opening angle
			  gCos->Fill(gCs);
			  TVector3 tptot = ptot.Vect();
			  gpCos->Fill(tptot.Mag(), gCs);
			  // Acceptance above 10 degrees and pt>1 GeV
			  Double_t AngMin = 10.0;
			  Double_t th1 = TMath::ACos(TMath::Abs(tP1.CosTheta()))*180. / TMath::Pi();
			  Double_t th2 = TMath::ACos(TMath::Abs(tP2.CosTheta()))*180. / TMath::Pi();
			  //
			  if (th1 < AngMin || th2 < AngMin)continue;		// Angular acceptance cuts
			  if (tP1.Pt() < 1. || tP2.Pt() < 1.)continue;		// Pt cuts
			  //
			  // Fill IDEA Histograms
			  Double_t CsMin = -1.0;
			  Double_t CsMax = -0.4;
			  Double_t pMin = 35.0;
			  Double_t pMax = 60.0;
			  TVector3 tX(0.0, 0.0, 0.0);						// Set origin to (0,0,0)
			  ObsTrk *Tr1id = new ObsTrk(tX, tP1, Q1, Bfield, GCid);	// Apply track resolution
			  ObsTrk *Tr2id = new ObsTrk(tX, tP2, Q2, Bfield, GCid);
			  TVector3 obsP1id = Tr1id->GetObsP();					// Get observed momenta
			  TVector3 obsP2id = Tr2id->GetObsP();           
			  Double_t oCsid = obsP1id.Dot(obsP2id) / 
				                (obsP1id.Mag()*obsP2id.Mag()); // Cos of muon opening angle
			  oCos_id->Fill(oCsid);
			  TVector3 pPairId = obsP1id + obsP2id;
			  if (oCsid > CsMin && oCsid < CsMax 
				  && pPairId.Mag()>pMin && pPairId.Mag()<pMax)
			  {
				  Double_t E1id = TMath::Sqrt(M1*M1 + obsP1id.Mag2());	// Observed energies
				  Double_t E2id = TMath::Sqrt(M2*M2 + obsP2id.Mag2());
				  TLorentzVector oP1id(obsP1id, E1id);						// Fill observed Lorentz vectors
				  TLorentzVector oP2id(obsP2id, E2id);
				  TLorentzVector oPtotid = oP1id + oP2id;					// Total momentum 4-vector
				  Float_t MiObsid = oPtotid.M();							// Invariant mass
				  oMinv_id->Fill(MiObsid);								// Fill invariant mass histogram
				  Double_t oHrcid = TMath::Sqrt(s + MiObsid * MiObsid - 2 * rs *oPtotid.E());	// Observed recoil
				  Double_t eoHrcid = TMath::Sqrt(ss + MiObsid * MiObsid - 2 * rss*oPtotid.E());	// Observed recoil w/ E spread 
				  eHrec->Fill(eHrc);									// Fill recoil histograms
				  oHrec_id->Fill(oHrcid);
				  eoHrec_id->Fill(eoHrcid);
			  }
			  //
			  // Fill CLD Histograms
			  ObsTrk *Tr1cl = new ObsTrk(tX, tP1, Q1, Bfield, GCcl);	// Apply track resolution
			  ObsTrk *Tr2cl = new ObsTrk(tX, tP2, Q2, Bfield, GCcl);
			  TVector3 obsP1cl = Tr1cl->GetObsP();					// Get observed momenta
			  TVector3 obsP2cl = Tr2cl->GetObsP();
			  Double_t oCscl = obsP1cl.Dot(obsP2cl) /
				  (obsP1cl.Mag()*obsP2cl.Mag()); // Cos of muon opening angle
			  oCos_cl->Fill(oCscl);
			  TVector3 pPairCl = obsP1cl + obsP2cl;
			  if (oCsid > CsMin && oCsid < CsMax
				  && pPairCl.Mag()>pMin && pPairCl.Mag() < pMax)
			  {
				  Double_t E1cl = TMath::Sqrt(M1*M1 + obsP1cl.Mag2());	// Observed energies
				  Double_t E2cl = TMath::Sqrt(M2*M2 + obsP2cl.Mag2());
				  TLorentzVector oP1cl(obsP1cl, E1cl);						// Fill observed Lorentz vectors
				  TLorentzVector oP2cl(obsP2cl, E2cl);
				  TLorentzVector oPtotcl = oP1cl + oP2cl;					// Total momentum 4-vector
				  Float_t MiObscl = oPtotcl.M();							// Invariant mass
				  oMinv_cl->Fill(MiObscl);								// Fill invariant mass histogram
				  Double_t oHrccl = TMath::Sqrt(s + MiObscl * MiObscl - 2 * rs *oPtotcl.E());	// Observed recoil
				  Double_t eoHrccl = TMath::Sqrt(ss + MiObscl * MiObscl - 2 * rss*oPtotcl.E());	// Observed recoil w/ E spread 
				  // Fill recoil histograms
				  oHrec_cl->Fill(oHrccl);
				  eoHrec_cl->Fill(eoHrccl);
			  }
		  }
	  }
   }				// End event loop
   //
   fin->Close();
   delete fin;
   //
   // Plot histograms
   TCanvas* c1 = new TCanvas("c1","Pythia8 HZ (Z{#rightarrow}#mu #mu",10,10,800,800);
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
   gMinv->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV)");
   gMinv->Draw();
   c1->cd(4);
   // Invariant mass plot
   oMinv_id->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV)");
   oMinv_id->SetLineColor(kBlue);
   oMinv_id->SetTitle("Observed di-muon invariant mass");
   oMinv_id->Draw();
   oMinv_cl->SetLineColor(kRed);
   oMinv_cl->Draw("SAME");
   TLegend *lgZ = new TLegend(0.1, 0.9, 0.3, 0.70);
   lgZ->AddEntry(oMinv_id, "IDEA");
   lgZ->AddEntry(oMinv_cl, "CLD");
   lgZ->Draw();
   //
   TCanvas* c2 = new TCanvas("c2", "Pythia8 HZ (Z{#rightarrow}#mu #mu)", 50, 50, 800, 500);
   c2->Divide(2, 2);
   c2->cd(1);
   gHrec->SetXTitle("Recoil mass (GeV)");
   gHrec->Draw();
   c2->cd(2);
   eHrec->SetTitle("Higgs recoil mass with 0.136% beam spread");
   eHrec->SetXTitle("Recoil mass (GeV)");
   eHrec->Draw();
   c2->cd(3);
   oHrec_id->SetTitle("Higgs recoil mass - No Energy spread");
   oHrec_id->SetXTitle("Recoil mass (GeV)");
   oHrec_id->SetLineColor(kBlue);
   oHrec_id->Draw();
   oHrec_cl ->SetTitle("Higgs recoil mass - No Energy spread");
   oHrec_cl	->SetXTitle("Recoil mass (GeV)");
   oHrec_cl->SetLineColor(kRed);
   oHrec_cl->Draw("SAME");
   TLegend *lg = new TLegend(0.5, 0.9, 0.7, 0.70);
   TString LgTitle = "Detectors:";
   lg->SetHeader(LgTitle);
   lg->AddEntry(oHrec_id, "IDEA", "L");
   lg->AddEntry(oHrec_cl, "CLD", "L");
   lg->Draw();
   c2->cd(4);
   //gHrec->Draw();
   eHrec->SetLineColor(kBlack);
   eHrec->SetTitle("Higgs recoil mass with 0.136% beam spread");
   eHrec->SetXTitle("Recoil mass (GeV)");
   eHrec->Draw();
   eoHrec_id->SetTitle("Higgs recoil mass with 0.136% beam spread");
   eoHrec_id->SetXTitle("Recoil mass (GeV)");
   eoHrec_id->SetLineColor(kBlue);
   //eoHrec_id->SetMaximum(900);
   eoHrec_id->Draw("SAME");
   eoHrec_cl->SetTitle("Higgs recoil mass with 0.1% beam spread");
   eoHrec_cl->SetXTitle("Recoil mass (GeV)");
   eoHrec_cl->SetLineColor(kRed);
   eoHrec_cl->Draw("SAME");
   TLegend *lg1 = new TLegend(0.5, 0.9, 0.7, 0.70);
   lg1->SetHeader(LgTitle);
   lg1->AddEntry(eHrec, "Beam only", "L");
   lg1->AddEntry(eoHrec_id, "IDEA", "L");
   lg1->AddEntry(eoHrec_cl, "CLD", "L");
   lg1->Draw();
   TCanvas *c3 = new TCanvas("c3", "Opening angle", 70, 70, 500, 500);
   c3->Divide(2, 2);
   c3->cd(1);
   gCos->Draw();
   c3->cd(2);
   oCos_id->Draw();
   c3->cd(3);
   oCos_cl->Draw();
   c3->cd(4);
   gpCos->Draw();
   //
   // Write out histograms
   TString name = fname;
   Int_t idot = name.Index("."); 
   Int_t iend = name.Length() - 1;
   name.Remove(idot, iend);
   TString hfname = "HfCut_" + name + ".root";
   TFile *hfil = new TFile(hfname, "RECREATE");
   oHrec_id->Write();
   oHrec_cl->Write();
   eHrec->Write();
   eoHrec_id->Write();
   eoHrec_cl->Write();
   gCos->Write();
   oCos_id->Write();
   oCos_cl->Write();
   gpCos->Write();
   hfil->Close();

 }
