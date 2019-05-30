#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <iostream>

//
void PlotRecoil()
{
	TFile *HstF_HZ = new TFile("Hf_ZHmm.root", "READ");
	Double_t scale_HZ = 3.28;
	TFile *HstF_ZZ = new TFile("Hf_ZZ.root", "READ");
	Double_t scale_ZZ = 68.1;
	TFile *HstF_WW = new TFile("Hf_WW.root", "READ");
	Double_t scale_WW = 81.9;
	//Double_t scale_WW = 8.19;
	//
	// Plot recoil with no energy spread
	// A. IDEA
	Int_t rebin = 2;
	TH1D *Wrec_id = (TH1D *)HstF_WW->Get("oHrec_id");
	Wrec_id->SetFillColor(kGreen);
	Wrec_id->Scale(scale_WW);
	Wrec_id->Rebin(rebin);
	TH1D *Zrec_id = (TH1D *)HstF_ZZ->Get("oHrec_id");
	Zrec_id->SetFillColor(kBlue);
	Zrec_id->Scale(scale_ZZ);
	Zrec_id->Rebin(rebin);
	TH1D *Hrec_id = (TH1D *)HstF_HZ->Get("oHrec_id");
	Hrec_id->SetFillColor(kRed);
	Hrec_id->Scale(scale_HZ);
	Hrec_id->Rebin(rebin);
	THStack *hRec_id = new THStack("hRec_id", "IDEA: Higgs recoil (No energy spread)");
	hRec_id->Add(Wrec_id);
	hRec_id->Add(Zrec_id);
	hRec_id->Add(Hrec_id);
	THStack *hRec_id1 = new THStack("hRec_id1", "IDEA: Higgs recoil (No energy spread)");
	hRec_id1->Add(Zrec_id);
	hRec_id1->Add(Hrec_id);
	TCanvas *CNBSA = new TCanvas("CNBSA", "IDEA: Higgs recoil (no energy spread)", 10, 10, 800, 500);
	CNBSA->Divide(1,2);
	Double_t bin_size = rebin * 0.5;
	char Ylabel[20];
	sprintf(Ylabel, "Events/%1.1f GeV", bin_size);
	//cout << "Ylabel= " << Ylabel << endl;
	CNBSA->cd(1);
	TLegend *lga1 = new TLegend(0.1, 0.9, 0.3, 0.7);
	lga1->AddEntry(Hrec_id, "HZ", "f");
	lga1->AddEntry(Zrec_id, "ZZ", "f");
	lga1->AddEntry(Wrec_id, "WW", "f");
	hRec_id->Draw();
	hRec_id->GetYaxis()->SetTitle(Ylabel);
	hRec_id->GetXaxis()->SetTitle("Recoil mass (GeV)");
	lga1->Draw();
	CNBSA->Modified();
	CNBSA->cd(2);
	TLegend *lga2 = new TLegend(0.1, 0.9, 0.3, 0.7);
	lga2->AddEntry(Hrec_id, "HZ", "f");
	lga2->AddEntry(Zrec_id, "ZZ", "f");
	hRec_id1->Draw();
	hRec_id1->GetYaxis()->SetTitle(Ylabel);
	hRec_id1->GetXaxis()->SetTitle("Recoil mass (GeV)");
	lga2->Draw();
	CNBSA->Modified();
	// B. CLD
	TH1D *Wrec_cl = (TH1D *)HstF_WW->Get("oHrec_cl");
	Wrec_cl->SetFillColor(kGreen);
	Wrec_cl->Scale(scale_WW);
	Wrec_cl->Rebin(rebin);
	TH1D *Zrec_cl = (TH1D *)HstF_ZZ->Get("oHrec_cl");
	Zrec_cl->SetFillColor(kBlue);
	Zrec_cl->Scale(scale_ZZ);
	Zrec_cl->Rebin(rebin);
	TH1D *Hrec_cl = (TH1D *)HstF_HZ->Get("oHrec_cl");
	Hrec_cl->SetFillColor(kRed);
	Hrec_cl->Scale(scale_HZ);
	Hrec_cl->Rebin(rebin);
	THStack *hRec_cl = new THStack("hRec_cl", "CLD: Higgs recoil (No energy spread)");
	hRec_cl->Add(Wrec_cl);
	hRec_cl->Add(Zrec_cl);
	hRec_cl->Add(Hrec_cl);
	THStack *hRec_cl1 = new THStack("hRec_cl1", "CLD: Higgs recoil (No energy spread)");
	hRec_cl1->Add(Zrec_cl);
	hRec_cl1->Add(Hrec_cl);
	TCanvas *CNBSB = new TCanvas("CNBSB", "CLD: Higgs recoil (no energy spread)", 30, 30, 800, 500);
	CNBSB->Divide(1, 2);
	CNBSB->cd(1);
	TLegend *lgb1 = new TLegend(0.1, 0.9, 0.3, 0.7);
	lgb1->AddEntry(Hrec_cl, "HZ", "f");
	lgb1->AddEntry(Zrec_cl, "ZZ", "f");
	lgb1->AddEntry(Wrec_cl, "WW", "f");
	hRec_cl->Draw();
	hRec_cl->GetYaxis()->SetTitle(Ylabel);
	hRec_cl->GetXaxis()->SetTitle("Recoil mass (GeV)");
	lgb1->Draw();
	CNBSA->Modified();
	CNBSB->cd(2);
	TLegend *lgb2 = new TLegend(0.1, 0.9, 0.3, 0.7);
	lgb2->AddEntry(Hrec_cl, "HZ", "f");
	lgb2->AddEntry(Zrec_cl, "ZZ", "f");
	hRec_cl1->Draw();
	hRec_cl1->GetYaxis()->SetTitle(Ylabel);
	hRec_cl1->GetXaxis()->SetTitle("Recoil mass (GeV)");
	lgb2->Draw();
	CNBSA->Modified();
	// Plot recoil with energy spread
	// A. IDEA

	TH1D *eWrec_id = (TH1D *)HstF_WW->Get("eoHrec_id");
	eWrec_id->SetFillColor(kGreen);
	eWrec_id->Scale(scale_WW);
	eWrec_id->Rebin(rebin);
	TH1D *eZrec_id = (TH1D *)HstF_ZZ->Get("eoHrec_id");
	eZrec_id->SetFillColor(kBlue);
	eZrec_id->Scale(scale_ZZ);
	eZrec_id->Rebin(rebin);
	TH1D *eHrec_id = (TH1D *)HstF_HZ->Get("eoHrec_id");
	eHrec_id->SetFillColor(kRed);
	eHrec_id->Scale(scale_HZ);
	eHrec_id->Rebin(rebin);
	THStack *ehRec_id = new THStack("ehRec_id", "IDEA: Higgs recoil #Delta E/E = .136%");
	ehRec_id->Add(eWrec_id);
	ehRec_id->Add(eZrec_id);
	ehRec_id->Add(eHrec_id);
	THStack *ehRec_id1 = new THStack("ehRec_id1", "IDEA: Higgs recoil #Delta E/E = .136% ");
	ehRec_id1->Add(eZrec_id);
	ehRec_id1->Add(eHrec_id);
	TCanvas *eCNBSA = new TCanvas("eCNBSA", "IDEA: Higgs recoil #Delta E/E = .136%", 10, 10, 800, 500);
	eCNBSA->Divide(1, 2);
	eCNBSA->cd(1);
	TLegend *elga1 = new TLegend(0.1, 0.9, 0.3, 0.7);
	elga1->AddEntry(eHrec_id, "HZ", "f");
	elga1->AddEntry(eZrec_id, "ZZ", "f");
	elga1->AddEntry(eWrec_id, "WW", "f");
	ehRec_id->Draw();
	ehRec_id->GetYaxis()->SetTitle(Ylabel);
	ehRec_id->GetXaxis()->SetTitle("Recoil mass (GeV)");
	elga1->Draw();
	eCNBSA->Modified();
	eCNBSA->cd(2);
	TLegend *elga2 = new TLegend(0.1, 0.9, 0.3, 0.7);
	elga2->AddEntry(Hrec_id, "HZ", "f");
	elga2->AddEntry(Zrec_id, "ZZ", "f");
	ehRec_id1->Draw();
	ehRec_id1->GetYaxis()->SetTitle(Ylabel);
	ehRec_id1->GetXaxis()->SetTitle("Recoil mass (GeV)");
	elga2->Draw();
	eCNBSA->Modified();
	// B. CLD

	TH1D *eWrec_cl = (TH1D *)HstF_WW->Get("eoHrec_cl");
	eWrec_cl->SetFillColor(kGreen);
	eWrec_cl->Scale(scale_WW);
	eWrec_cl->Rebin(rebin);
	TH1D *eZrec_cl = (TH1D *)HstF_ZZ->Get("eoHrec_cl");
	eZrec_cl->SetFillColor(kBlue);
	eZrec_cl->Scale(scale_ZZ);
	eZrec_cl->Rebin(rebin);
	TH1D *eHrec_cl = (TH1D *)HstF_HZ->Get("eoHrec_cl");
	eHrec_cl->SetFillColor(kRed);
	eHrec_cl->Scale(scale_HZ);
	eHrec_cl->Rebin(rebin);
	THStack *ehRec_cl = new THStack("ehRec_cl", "CLD: Higgs recoil #Delta E/E = .136%");
	ehRec_cl->Add(eWrec_cl);
	ehRec_cl->Add(eZrec_cl);
	ehRec_cl->Add(eHrec_cl);
	THStack *ehRec_cl1 = new THStack("ehRec_cl1", "CLD: Higgs recoil #Delta E/E = .136%");
	ehRec_cl1->Add(eZrec_cl);
	ehRec_cl1->Add(eHrec_cl);
	TCanvas *eCNBSB = new TCanvas("eCNBSB", "CLD: Higgs recoil #Delta E/E = .136%", 30, 30, 800, 500);
	eCNBSB->Divide(1, 2);
	eCNBSB->cd(1);
	TLegend *elgb1 = new TLegend(0.1, 0.9, 0.3, 0.7);
	elgb1->AddEntry(eHrec_cl, "HZ", "f");
	elgb1->AddEntry(eZrec_cl, "ZZ", "f");
	elgb1->AddEntry(eWrec_cl, "WW", "f");
	ehRec_cl->Draw();
	ehRec_cl->GetYaxis()->SetTitle(Ylabel);
	ehRec_cl->GetXaxis()->SetTitle("Recoil mass (GeV)");
	elgb1->Draw();
	eCNBSA->Modified();
	eCNBSB->cd(2);
	TLegend *elgb2 = new TLegend(0.1, 0.9, 0.3, 0.7);
	elgb2->AddEntry(eHrec_cl, "HZ", "f");
	elgb2->AddEntry(eZrec_cl, "ZZ", "f");
	ehRec_cl1->Draw();
	ehRec_cl1->GetYaxis()->SetTitle(Ylabel);
	ehRec_cl1->GetXaxis()->SetTitle("Recoil mass (GeV)");
	elgb2->Draw();
	eCNBSA->Modified();
	//
	//HstF_HZ->Close();
	//HstF_ZZ->Close();
	//HstF_WW->Close();
}