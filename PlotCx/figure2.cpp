
#include "plot_cx.cpp"
//#include "plot_cx_raffaella.cpp"

void figure2(double angle=21){
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
  TCanvas *c1=new TCanvas("c1","Figure 2",800,500);
  c1->SetGrayscale(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();

  Float_t small = 1e-5;
  Float_t small2 = .05;
  c1->Divide(2,2,small,small);

  //void plot_cx(string target="d", string angle="39", string spec="shms",float "xb"min=0.0, float "xb"max=1.0, float cxmin=0.0, float cxmax = 50.0, float ratiomin = 0.7, float ratiomax =1.3 ){


  c1->cd(1);
  gPad->SetBottomMargin(small);
  gPad->SetRightMargin(small);
  gPad->SetTopMargin(.1);
  gPad->SetLeftMargin(.15);
  gROOT->ProcessLine("plot_cx(\"r\",\"25\",\"shms\",\"pass400\",.25,.95,0.57,1.14,.85,1.15,\"xb\",1)");
  TLegend* l= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_shms");
  l->Delete();
  TPaveText *t1=new TPaveText(.22,.8,.93,.89,"NDC");
  //  t1->AddText("#theta_{c} = 25.0#circ");
  //  t1->AddText("4.41 < Q^{2} (GeV^{2}) < 9.91");
  t1->AddText("4.41 < Q^{2} (GeV^{2}) < 9.91; #theta_{c} = 25.0#circ");

  t1->SetBorderSize(0);
  t1->SetFillColor(kWhite);
  t1->Draw();

  c1->cd(2);
  gPad->SetTopMargin(.1);
  gPad->SetBottomMargin(small);
  gPad->SetRightMargin(.15);
  gPad->SetLeftMargin(small);
  gROOT->ProcessLine("plot_cx(\"r\",\"29\",\"shms\",\"pass400\",.25,.95,0.57,1.14,.85,1.15,\"xb\",1)");
  TLegend* l2= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_shms");
  l2->Delete();
  TPaveText *t2=new TPaveText(.05,.8,.8,.89,"NDC");
  //  t2->AddText("#theta_{c} = 29.0#circ");
  t2->AddText("4.73 < Q^{2} (GeV^{2}) < 11.06; #theta_{c} = 29.0#circ");
  //  t2->AddText("4.73 < Q^{2} (GeV^{2}) < 11.06");
  t2->SetBorderSize(0);
  t2->SetFillColor(kWhite);
  t2->Draw();

  c1->cd(3);
  gPad->SetTopMargin(small);
  gPad->SetRightMargin(small);
  gPad->SetLeftMargin(0.15);
  gROOT->ProcessLine("plot_cx(\"r\",\"33\",\"shms\",\"pass400\",.25,.95,0.57,1.14,.85,1.15,\"xb\",1)");
  TLegend* l3= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_shms");
  l3->Delete();
  TPaveText *t3=new TPaveText(.2,.9,.95 ,0.99,"NDC");
  //  t3->AddText("#theta_{c} = 33.0#circ");
  t3->AddText("5.13 < Q^{2} (GeV^{2}) < 12.09; #theta_{c} = 33.0#circ");
  //  t3->AddText("5.13 < Q^{2} (GeV^{2}) < 12.09");
  t3->SetBorderSize(0);
  t3->SetFillColor(kWhite);
  t3->Draw();
  c1->cd(4);
  gPad->SetTopMargin(small);
  gPad->SetRightMargin(.15);
  gPad->SetLeftMargin(small);
  gROOT->ProcessLine("plot_cx(\"r\",\"39\",\"shms\",\"pass400\",.25,.95,0.57,1.14,.85,1.15,\"xb\",1)");
  TLegend* l4= (TLegend*)gPad->GetListOfPrimitives()->FindObject("leg_shms");
  //  l4->SetNColumns(2);
  l4->SetX1(.15);
  l4->SetX2(.7);
  l4->SetY1(.59);
  l4->SetY2(.89);
  /* 
     l4->SetX1(.25);
  l4->SetX2(.6);
  l4->SetY1(.55);
  l4->SetY2(.89);
  */
  l4->Draw();
  TPaveText *t4=new TPaveText(.05,.9,.8,.99,"NDC");
  //  t4->AddText("#theta_{c} = 39.0#circ");
  t4->AddText("5.46 < Q^{2} (GeV^{2}) < 13.40; #theta_{c} = 39.0#circ");
  //  t4->AddText("5.46 < Q^{2} (GeV^{2}) < 13.40");
  t4->SetBorderSize(0);
t4->SetFillColor(kWhite);
  t4->Draw();

  c1->SaveAs("figure2.pdf");
  //  return;
}
