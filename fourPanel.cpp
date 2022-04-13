{
TCanvas*c1=new TCanvas("c1","c1",400,800);
c1->Divide(1,4);

gStyle->SetTitleX(.1);
gStyle->SetTitleW(.4);
gStyle->SetTitleY(.97);
gStyle->SetTitleH(.08);
 gStyle->SetTitleW(.3);
 gStyle->SetOptStat(0);
 gStyle->SetPalette(kRainBow);

  heff->Draw("LEGO2");
  gPad->SetPhi(120);
  heff->GetXaxis()->SetRangeUser(-35,35);
  heff->GetXaxis()->CenterTitle(1);
  heff->GetYaxis()->CenterTitle(1);
  heff->GetZaxis()->CenterTitle(1);
  heff->GetZaxis()->SetTitle("Effective Solid Angle (msr)");
  heff->GetXaxis()->SetTitleSize(.048);
  heff->GetYaxis()->SetTitleSize(.048);
  heff->GetZaxis()->SetTitleSize(.048);
  heff->GetXaxis()->SetTitleOffset(1.7);
  heff->GetYaxis()->SetTitleOffset(1.5);
  heff->GetZaxis()->SetTitleOffset(1.04);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.1);
  heff->Scale(1000.);
  heff->Draw("LEGO2");

c1->cd(1);
TFile *_file0 = TFile::Open("acc30bins_shms_21deg3p3.root");
heff->SetTitle("10 cm Cryo Target");
//.x drawAcc.cpp;
  heff->Draw("LEGO2");
  gPad->SetPhi(120);
  heff->GetXaxis()->SetRangeUser(-35,35);
  heff->GetXaxis()->CenterTitle(1);
  heff->GetYaxis()->CenterTitle(1);
  heff->GetZaxis()->CenterTitle(1);
  heff->GetZaxis()->SetTitle("Effective Solid Angle (msr)");
  heff->GetXaxis()->SetTitleSize(.048);
  heff->GetYaxis()->SetTitleSize(.048);
  heff->GetZaxis()->SetTitleSize(.048);
  heff->GetXaxis()->SetTitleOffset(1.7);
  heff->GetYaxis()->SetTitleOffset(1.5);
  heff->GetZaxis()->SetTitleOffset(1.04);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.1);
  heff->Scale(1000.);
  heff->Draw("LEGO2");
c1->cd(2);
TFile *_file1 = TFile::Open("acc30bins_shms_c21deg3p3.root");
heff->SetTitle("Carbon foil z=0cm");
//.x drawAcc.cpp;
  heff->Draw("LEGO2");
  gPad->SetPhi(120);
  heff->GetXaxis()->SetRangeUser(-35,35);
  heff->GetXaxis()->CenterTitle(1);
  heff->GetYaxis()->CenterTitle(1);
  heff->GetZaxis()->CenterTitle(1);
  heff->GetZaxis()->SetTitle("Effective Solid Angle (msr)");
  heff->GetXaxis()->SetTitleSize(.048);
  heff->GetYaxis()->SetTitleSize(.048);
  heff->GetZaxis()->SetTitleSize(.048);
  heff->GetXaxis()->SetTitleOffset(1.7);
  heff->GetYaxis()->SetTitleOffset(1.5);
  heff->GetZaxis()->SetTitleOffset(1.04);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.1);
  heff->Scale(1000.);
  heff->Draw("LEGO2");
c1->cd(3);
TFile *_file2 = TFile::Open("acc30bins_shms_ald21deg3p3.root");
heff->SetTitle("Al Dummy z=5cm");
//.x drawAcc.cpp;
  heff->Draw("LEGO2");
  gPad->SetPhi(120);
  heff->GetXaxis()->SetRangeUser(-35,35);
  heff->GetXaxis()->CenterTitle(1);
  heff->GetYaxis()->CenterTitle(1);
  heff->GetZaxis()->CenterTitle(1);
  heff->GetZaxis()->SetTitle("Effective Solid Angle (msr)");
  heff->GetXaxis()->SetTitleSize(.048);
  heff->GetYaxis()->SetTitleSize(.048);
  heff->GetZaxis()->SetTitleSize(.048);
  heff->GetXaxis()->SetTitleOffset(1.7);
  heff->GetYaxis()->SetTitleOffset(1.5);
  heff->GetZaxis()->SetTitleOffset(1.04);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.1);
  heff->Scale(1000.);
  heff->Draw("LEGO2");
c1->cd(4);
TFile *_file3 = TFile::Open("acc30bins_shms_alu21deg3p3.root");
//.x drawAcc.cpp;
  heff->Draw("LEGO2");
heff->SetTitle("Al Dummy z=-5cm");

  gPad->SetPhi(120);
  heff->GetXaxis()->SetRangeUser(-35,35);
  heff->GetXaxis()->CenterTitle(1);
  heff->GetYaxis()->CenterTitle(1);
  heff->GetZaxis()->CenterTitle(1);
  heff->GetZaxis()->SetTitle("Effective Solid Angle (msr)");
  heff->GetXaxis()->SetTitleSize(.048);
  heff->GetYaxis()->SetTitleSize(.048);
  heff->GetZaxis()->SetTitleSize(.048);
  heff->GetXaxis()->SetTitleOffset(1.7);
  heff->GetYaxis()->SetTitleOffset(1.5);
  heff->GetZaxis()->SetTitleOffset(1.04);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.1);
  heff->Scale(1000.);
  heff->Draw("LEGO2");
c1->SaveAs("accDiffTargets.pdf")
}
