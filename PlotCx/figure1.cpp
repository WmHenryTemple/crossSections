#include "plot_cx_f1.cpp"
//#include "plot_cx_raffaella.cpp"

void figure1(double angle=21){
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  //  gStyle->SetGridStyle(3);
  //  gStyle->SetTitleX(.5);
  //  gStyle->SetTitleW(.7);
  gStyle->SetOptTitle(0);
  //  gStyle->SetBorderSize(1);
  TCanvas *c1=new TCanvas("c1","Figure 1",700,700);
  gPad->SetTopMargin(.09);
  gPad->SetLeftMargin(.12);
  gPad->SetBottomMargin(.09);
  gPad->SetRightMargin(.12);

  c1->SetGrayscale(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  gROOT->ProcessLine("plot_cx(\"r\",\"21\",\"hms\",\"pass400\",.2,.95,0.54,1.13,.85,1.15,\"xb\",1)");

  //  l->DeleteEntry();

  TLegend* l= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_hms");
  //  l1->Delete();
  l->SetNColumns(2);
  l->SetX1(.15);
  l->SetX2(.81);
  l->SetY1(.67);
  l->SetY2(.88);
  l->Draw();

  TPaveText *t1=new TPaveText(.2,.823,.8,.905,"NDC");

  //  t1->AddText("#theta_{c} = 25.0#circ");
  //  t1->AddText("4.41 < Q^{2} (GeV^{2}) < 9.91");
  t1->AddText("3.39 < Q^{2} (GeV^{2}) < 8.25"); //#theta_{c} = 21.0#circ");

  t1->SetBorderSize(0);
  t1->SetFillColor(kWhite);
  //  t1->Draw();

  TGraphErrors *gs=(TGraphErrors*)gPad->GetListOfPrimitives()->FindObject("grcx1");
  gROOT->ProcessLine("plot_cx(\"r\",\"21\",\"shms\",\"pass400\",.2,.95,0.48,1.1,.85,1.15,\"xb\",1)");
  //  TLegend* l0= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_shms");
  //  l0->Delete();


  c1->SaveAs("figure1.pdf");
  //  return;
}
