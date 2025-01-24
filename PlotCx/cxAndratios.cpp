#include "plot_cxAr.cpp"
//#include "plot_cx_raffaella.cpp"

void cxAndratios(double angle=59){
  //  gROOT->ForceStyle();
  //  gStyle->SetGridStyle(3);
  //  gStyle->SetTitleX(.5);
  //  gStyle->SetTitleW(.7);

  TCanvas *c1=new TCanvas();
  c1->Divide(3,1);
  c1->SetGrid();
  //void plot_cx(string target="d", string angle="39", string spec="shms",float "xb"min=0.0, float "xb"max=1.0, float cxmin=0.0, float cxmax = 50.0, float ratiomin = 0.7, float ratiomax =1.3 ){
  if(angle==21)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"21\",\"shms\",\"pass400\",2,6,-.1,4.,.9,1.1,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"21\",\"shms\",\"pass400\",2,6,-.1,8.,.9,1.1,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"21\",\"shms\",\"pass400\",2,6,0.6,1.1,.9,1.1,\"xb\")");
    }
  if(angle==20)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"21\",\"hms\",\"pass400\",2.,6.,-.1,3.,.9,1.1,\"xb\")");
      c1->cd(2);
	gROOT->ProcessLine("plot_cx(\"d\",\"21\",\"hms\",\"pass400\",2.,6.,-.1,6.,.9,1.1,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"21\",\"hms\",\"pass400\",2.,6.,0.5,1.1,.9,1.1,\"xb\")");
    }

  if(angle==25)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"25\",\"shms\",\"pass400\",2,5.5,-.1,2.,.85,1.15,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"25\",\"shms\",\"pass400\",2,5.5,-.1,4.,.85,1.15,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"25\",\"shms\",\"pass400\",2,5.5,0.6,1.1,.85,1.15,\"xb\")");
    }

  if(angle==29)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"29\",\"shms\",\"pass400\",1.5,4.5,-.1,1.2,.85,1.15,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"29\",\"shms\",\"pass400\",1.5,4.5,-.1,2.4,.85,1.15,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"29\",\"shms\",\"pass400\",1.5,4.5,0.6,1.1,.85,1.15,\"xb\")");
    }

  if(angle==33)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"33\",\"shms\",\"pass400\",1,4.,-.1,.8,.85,1.15,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"33\",\"shms\",\"pass400\",1,4.,-.1,1.6,.85,1.15,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"33\",\"shms\",\"pass400\",1,4.,0.6,1.1,.85,1.15,\"xb\")");
    }

  if(angle==39)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"39\",\"shms\",\"pass400\",1.,3.2,-.1,.5,.85,1.15,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"39\",\"shms\",\"pass400\",1.,3.2,-.1,1.,.85,1.15,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"39\",\"shms\",\"pass400\",1.,3.2,0.6,1.1,.85,1.15,\"xb\")");
    }

  /*

  if(angle==59)
    {
      c1->cd(1);
      gROOT->ProcessLine("plot_cx(\"h\",\"59\",\"hms\",\"pass327\",.2,1.,-.1,.5,.85,1.15,\"xb\")");
      c1->cd(2);
      gROOT->ProcessLine("plot_cx(\"d\",\"59\",\"hms\",\"pass327\",.2,1.,-.1,1.,.85,1.15,\"xb\")");
      c1->cd(3);
      gROOT->ProcessLine("plot_cx(\"r\",\"59\",\"hms\",\"pass327\",.2,1.,0.6,1.1,.85,1.15,\"xb\")");
    }
  */
  c1->SaveAs(Form("shms%2.0fdeg_400.ep.pdf",angle));
  return;
}
