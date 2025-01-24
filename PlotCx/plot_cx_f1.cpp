#include "extractCS.cpp"
#include "src/whitlow/readIt.cpp"
#include "src/formatHist.cpp"
#include "src/xsecTable.cpp"
//#include "src/getRadCorr.cpp"
//#include "src/readReport.cpp"
#include "src/getMomList.cpp"
#include "src/getCSB.cpp"
#include "src/getCJ.cpp"
#include "src/getJAM.cpp"
#include "src/getKull.cpp"
#include "src/getKP.cpp"
#include "jmu/getJmu.cpp"

void plot_cx(string target="d", string angle="39", string spec="shms", string pass="pass105", float xbmin=0.0, float xbmax=1.0, float cxmin=0.0, float cxmax = 50.0, float ratiomin = 0.7, float ratiomax =1.3, string xaxis = "xb", int numPad=2 ){

 Float_t delta, ratio, err, ep, model_k, modeld, modelh, xmin, xmax;
 string pset[5];
 float psetf[5];
 getMomList(angle,pset,psetf,spec);
  xmin=psetf[0]*(.9);
  xmax=psetf[3]*(1.22);
    xmin=.85*2.7;
  //  xmin=1.22*5.1;
  // if(spec=="hms")xmax=psetf[4]*(1.1);

  ////////////////////////////////////////
  ////////   MODEL   /////////////////////
  ////////////////////////////////////////

 vector <float> mx;
 vector <float> my_kph;  //kp hybrid
 vector <float> my_kpd;  //kp dis 
 vector <float> my;    //CJ15
 double wmin,wmax;

 Double_t thetac=getAngle(angle,spec);
 TGraph2D *grd=getRadCorrW2("d",1,"shms");  
 TGraph2D *grh=getRadCorrW2("h",1,"shms");  
  for (Int_t i=0;i<1000;i++)
   {
     ep=xmin+i/1000.*(xmax-xmin);
     Double_t ebeam=10.602;
     Double_t sin2 = sin(thetac/2/180*TMath::Pi())*sin(thetac/2/180*TMath::Pi());
     Double_t nu= ebeam-ep;
     Double_t q2 = 4.*ep*ebeam*sin2;
     Float_t mp = .9382723;
     Double_t xb=q2/2./mp/nu;
     Double_t w2= mp*mp + 2.*mp*nu-q2;
     modeld=grd->Interpolate(w2,thetac);  //<<"\t"<<
     modelh=grh->Interpolate(w2,thetac);  //<<"\t"<<
     if(i==0)wmax=w2+.5;
     if(i==999)wmin=w2-.1;
     if(modelh!=0 && modeld!=0 && xb < 0.94){ 

       if(xaxis=="xb")mx.push_back(xb);
       if(xaxis=="w2")mx.push_back(w2);
       if(xaxis=="ep")mx.push_back(ep);

       if(target=="r")my.push_back(modeld/modelh/2.);
       if(target=="h")my.push_back(modelh);
       if(target=="d")my.push_back(modeld);
     }
   }
 Int_t pts=mx.size();
 TGraph *gm=new TGraph(pts,&mx[0],&my[0]);
 gm->SetName("gm");
 TGraph *gm_cj=getCJ(angle);  //at SHMS angles
 TGraph *gm_k=getKP("dis",angle);  
 TGraph *gm_k1=getKP("hybrid",angle);  
 TGraph *gm_k2=getKP("akp17",angle);  
 TGraphErrors *gm_jam=getJAM(thetac);  

 if(xaxis=="xb") gm->GetXaxis()->SetTitle("x");
 if(xaxis=="ep") gm->GetXaxis()->SetTitle("E' (GeV)");
 if(xaxis=="w2") gm->GetXaxis()->SetTitle("W2 (GeV^2)");
 gm->GetXaxis()->SetLimits(xbmin,xbmax);
 if(xaxis=="ep")gm->GetXaxis()->SetLimits(xmin-.1,xmax+.1);
 if(xaxis=="w2")gm->GetXaxis()->SetLimits(wmin,wmax);

 gm->SetLineColor(kBlue);
 gm_cj->SetLineColor(kRed);
 gm_k1->SetLineColor(kGreen+3);
 gm_jam->SetLineColor(kOrange-6);
 gm_k2->SetLineColor(9);
 
 gm->SetLineStyle(2);
 // gm_cj->SetLineStyle(2);
 gm_k->SetLineStyle(3);
 gm_k1->SetLineStyle(4);
 gm_k2->SetLineStyle(6);
 gm_jam->SetLineStyle(3);

 int lw=2;
 if (numPad==2)lw=1;

 gm->SetLineWidth(lw);
 gm_cj->SetLineWidth(lw);
 gm_k->SetLineWidth(lw);
 gm_k1->SetLineWidth(lw);
 gm_k2->SetLineWidth(lw);
 gm_jam->SetLineWidth(lw);


 if(target=="h"){
   gm->GetYaxis()->SetTitle("d#sigma/d#OmegadE' (nb/GeV/sr)");
   gm->SetTitle(Form("%s Hydrogen %sdeg",spec.c_str(),angle.c_str()));
   gm->GetYaxis()->SetRangeUser(cxmin,cxmax);
 }
 if(target=="d"){
   gm->GetYaxis()->SetTitle("d#sigma/d#OmegadE' (nb/GeV/sr)");
   gm->SetTitle(Form("%s Deuterium %sdeg",spec.c_str(),angle.c_str()));
   gm->GetYaxis()->SetRangeUser(cxmin,cxmax);
 }
if(target=="r"){
  gm->GetYaxis()->SetTitle("#sigma_{D}/#sigma_{H}");
  gm->SetTitle(Form("%s D/H/2 %sdeg",spec.c_str(),angle.c_str()));
  gm->GetYaxis()->SetRangeUser(cxmin,cxmax);
 }
// gm->GetYaxis()->SetTitle("");
// gm->GetXaxis()->SetTitle("");


 gm->GetYaxis()->SetTitleSize(.055);
  gm->GetXaxis()->SetTitleSize(.065);
  gm->GetXaxis()->SetTitleOffset(.7);
 gm->GetXaxis()->CenterTitle();
 gm->GetYaxis()->CenterTitle();
 // gm->GetXaxis()->SetNdivisions(-507 );
 // gm->GetYaxis()->SetNdivisions(-507 );
 // gm->GetXaxis()->SetLabelSize(0.06);// size of axis value font
 // TCanvas *c1=new TCanvas();
 TPad *pad1, *pad2;
 if(numPad==2){
   pad1 = new TPad("pad1","pad1",0,0.4,1,1);
   pad2 = new TPad("pad2","pad2",0,0.02,1,0.38);
   //  gStyle->SetGridStyle(1);
   // pad1->SetGrid();
   // pad2->SetGrid();
   pad1->SetBottomMargin(0.1);
   pad1->SetLeftMargin(0.13);
   pad2->SetLeftMargin(0.13);
   pad1->SetBorderMode(1);
   pad2->SetTopMargin(0.1);
   pad2->SetBottomMargin(0.17);
   pad2->SetBorderMode(0);
   pad1->Draw();
   pad2->Draw();
   pad1->SetGrid();
   pad2->SetGrid();
   pad1->cd();
 }
 if(spec=="hms")gm->Draw("la");
 // if(target=="r" && spec=="hms" && xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}
 if(target=="r" && spec=="hms" && xaxis=="xb"){gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");gm_jam->Draw("same");}


  ////////////////////////////////////////
  ////////   Data   /////////////////////
  ////////////////////////////////////////
 TGraphErrors *gwhit=readIt();

  double mrthn_x[7]={.195, .225, .255, .285, .315, .345, .375};
  double mrthn_w[7]={12.2,11.7,11.3,10.9,10.5,10.1,9.63};
  if(xaxis=="w2"){
  for(int i=0;i<7;i++)mrthn_x[i]=mrthn_w[i];
  }
  double mrthn_y[7]={1.725, 1.697, 1.674, 1.656, 1.629, 1.588, 1.544};
  double mrthn_ye[7]={.015, .014, .014, .016, .016, .016, .023};
  double mrthn_corr[7]={1.00, .998, .998, 1.00, 1.003, 1.006};
  for(int i=0;i<7;i++)mrthn_y[i]=mrthn_y[i]/2*mrthn_corr[i];
  for(int i=0;i<7;i++)mrthn_ye[i]=mrthn_ye[i]/2*mrthn_corr[i];
  TGraphErrors *mrthn=new TGraphErrors(7,mrthn_x,mrthn_y,0,mrthn_ye);
  mrthn->SetMarkerStyle(21);
  mrthn->SetMarkerSize(0.5);
  mrthn->SetMarkerColor(kGreen+2);
  //  if(target=="r"&& angle=="21")mrthn->Draw("sp");
  //  if(spec=="hms")gwhit->Draw("sp");
  /////////////////////////////////////////////////////////////////////

  TGraphErrors *grcx1=extractCS(spec,target,angle,pset[0],1,pass,xaxis);
  grcx1->SetName("grcx1");
  TGraphErrors *grcx2=extractCS(spec,target,angle,pset[1],1,pass,xaxis);
  TGraphErrors *grcx3=extractCS(spec,target,angle,pset[2],1,pass,xaxis);
  TGraphErrors *grcx4=extractCS(spec,target,angle,pset[3],1,pass,xaxis);
  TGraphErrors *grcx5;
  if(spec=="hms")grcx5=extractCS(spec,target,angle,pset[4],1,pass,xaxis);

  TGraphErrors *grcx1a=extractCS(spec,target,angle,pset[0],3,pass,xaxis);
  grcx1->SetName("grcx1");
  TGraphErrors *grcx2a=extractCS(spec,target,angle,pset[1],3,pass,xaxis);
  TGraphErrors *grcx3a=extractCS(spec,target,angle,pset[2],3,pass,xaxis);
  TGraphErrors *grcx4a=extractCS(spec,target,angle,pset[3],3,pass,xaxis);
  TGraphErrors *grcx5a;
  if(spec=="hms")grcx5a=extractCS(spec,target,angle,pset[4],3,pass,xaxis);


  TGraphErrors *grsys1=extractCS(spec,target,angle,pset[0],2,pass,xaxis);
  grsys1->SetName("grsys1");
  TGraphErrors *grsys2=extractCS(spec,target,angle,pset[1],2,pass,xaxis);
  grsys2->SetName("grsys2"); 
  TGraphErrors *grsys3=extractCS(spec,target,angle,pset[2],2,pass,xaxis);
  grsys3->SetName("grsys3");
  TGraphErrors *grsys4=extractCS(spec,target,angle,pset[3],2,pass,xaxis);
  grsys4->SetName("grsys4");
  TGraphErrors *grsys5;
  if(spec=="hms")grsys5=extractCS(spec,target,angle,pset[4],2,pass,xaxis); 

  TGraphErrors *gr1=extractCS(spec,target,angle,pset[0],0,pass,xaxis);
  TGraphErrors *gr2=extractCS(spec,target,angle,pset[1],0,pass,xaxis);
  TGraphErrors *gr3=extractCS(spec,target,angle,pset[2],0,pass,xaxis);
  TGraphErrors *gr4=extractCS(spec,target,angle,pset[3],0,pass,xaxis);
  TGraphErrors *gr5;
  if(spec=="hms")gr5=extractCS(spec,target,angle,pset[4],0,pass,xaxis);
  /*
  ///////////////////////////////////////////////////////
  string kin=spec+"_"+target+angle+"deg";
  xsecTable(grcx1a, grcx1, grsys1, gr1, thetac, kin+pset[0], target);
  xsecTable(grcx2a, grcx2, grsys2, gr2, thetac,  kin+pset[1], target);
  xsecTable(grcx3a, grcx3, grsys3, gr3, thetac,  kin+pset[2], target);
  xsecTable(grcx4a, grcx4, grsys4, gr4, thetac,  kin+pset[3], target);
  if(spec=="hms")xsecTable(grcx5a, grcx5, grsys5, gr5, thetac, kin+pset[4], target);
  */ 
 ////////////////////////////////////////////////////////
  cout << "I got out of the table routine"<<endl;
  TGraphErrors *jmu=getJmu(target,angle,grd,grh,1);

 jmu->SetMarkerStyle(33);

 grcx1->SetMarkerStyle(33);
 grcx2->SetMarkerStyle(33);
 grcx3->SetMarkerStyle(33);
 grcx4->SetMarkerStyle(33);
 if(spec=="hms")grcx5->SetMarkerStyle(33);

 jmu->SetMarkerSize(1.);
 grcx1->SetMarkerSize(1.);
 grcx2->SetMarkerSize(1.);
 grcx3->SetMarkerSize(1.);
 grcx4->SetMarkerSize(1.);
 if(spec=="hms")grcx5->SetMarkerSize(1.);

 grcx1a->SetMarkerStyle(33);
 grcx2a->SetMarkerStyle(33);
 grcx3a->SetMarkerStyle(33);
 grcx4a->SetMarkerStyle(33);
 if(spec=="hms")grcx5a->SetMarkerStyle(33);

 grcx1a->SetMarkerSize(1.);
 grcx2a->SetMarkerSize(1.);
 grcx3a->SetMarkerSize(1.);
 grcx4a->SetMarkerSize(1.);
 if(spec=="hms")grcx5a->SetMarkerSize(1.);

 jmu->SetMarkerColor(kBlack);
 // grcx1->SetMarkerColor(kRed);
 // grcx2->SetMarkerColor(kBlue);
 // grcx3->SetMarkerColor(kGreen+3);
 // grcx4->SetMarkerColor(kViolet);

 grsys1->SetFillColorAlpha(kGray+2,.5);
 grsys2->SetFillColorAlpha(kGray+2,.5);
 grsys3->SetFillColorAlpha(kGray+2,.5);
 grsys4->SetFillColorAlpha(kGray+2,.5);
 
 grsys1->SetFillStyle(3254);
 grsys2->SetFillStyle(3254);
 grsys3->SetFillStyle(3254);
 grsys4->SetFillStyle(3254);


 if(spec=="hms"){
   grcx1->SetMarkerColor(kPink+7);
   grcx2->SetMarkerColor(kPink+7);
   grcx3->SetMarkerColor(kPink+7);
   grcx4->SetMarkerColor(kPink+7);
   grcx5->SetMarkerColor(kPink+7);
   grcx1->SetLineColor(kPink+7);
   grcx2->SetLineColor(kPink+7);
   grcx3->SetLineColor(kPink+7);
   grcx4->SetLineColor(kPink+7);
   grcx5->SetLineColor(kPink+7);
   grcx1->SetMarkerStyle(20);
   grcx2->SetMarkerStyle(20);
   grcx3->SetMarkerStyle(20);
   grcx4->SetMarkerStyle(20);
   grcx5->SetMarkerStyle(20);

   grcx1a->SetMarkerColor(kPink+7);
   grcx2a->SetMarkerColor(kPink+7);
   grcx3a->SetMarkerColor(kPink+7);
   grcx4a->SetMarkerColor(kPink+7);
   grcx5a->SetMarkerColor(kPink+7);
   grcx1a->SetLineColor(kPink+7);
   grcx2a->SetLineColor(kPink+7);
   grcx3a->SetLineColor(kPink+7);
   grcx4a->SetLineColor(kPink+7);
   grcx5a->SetLineColor(kPink+7);
   grcx1a->SetMarkerStyle(20);
   grcx2a->SetMarkerStyle(20);
   grcx3a->SetMarkerStyle(20);
   grcx4a->SetMarkerStyle(20);
   grcx5a->SetMarkerStyle(20);

   grsys1->SetFillColorAlpha(kPink+7,.5);
   grsys2->SetFillColorAlpha(kPink+7,.5);
   grsys3->SetFillColorAlpha(kPink+7,.5);
   grsys4->SetFillColorAlpha(kPink+7,.5);
   grsys5->SetFillColorAlpha(kPink+7,.5);
   grsys1->SetFillStyle(3245);
   grsys2->SetFillStyle(3245);
   grsys3->SetFillStyle(3245);
   grsys4->SetFillStyle(3245);
   grsys5->SetFillStyle(3245);
 }



 // jmu->Draw("sp");
 // if(spec=="shms"){
   grsys1->Draw("E3 same");
   grsys2->Draw("E3 same");
   grsys3->Draw("E3 same");
   grsys4->Draw("E3 same");
   // }
   if(spec=="hms") grsys5->Draw("E3 same");
 grcx1->Draw("sp");
 grcx2->Draw("sp");
 grcx3->Draw("sp");
 grcx4->Draw("sp");
 if(spec=="hms") grcx5->Draw("sp");

 grcx1a->Draw("sp");
 grcx2a->Draw("sp");
 grcx3a->Draw("sp");
 grcx4a->Draw("sp");
 if(spec=="hms") grcx5a->Draw("sp");
 

  
 //  gm->Draw("same");

 // if(target=="r"&&spec=="hms"&&xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}
 // if(target=="r"&&spec=="hms"&&xaxis=="xb"){gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}
 if(target=="r" && spec=="hms" && xaxis=="xb"){gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");gm_jam->Draw("same");}
 if(target=="r" && spec=="shms" && xaxis=="xb"){gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");gm_jam->Draw("same");gm->Draw("same");}  //just added
 if(numPad==2)
   {
     pad2->cd();
     TGraphErrors *gr1=extractCS(spec,target,angle,pset[0],0,pass,xaxis);
     TGraphErrors *gr2=extractCS(spec,target,angle,pset[1],0,pass,xaxis);
     TGraphErrors *gr3=extractCS(spec,target,angle,pset[2],0,pass,xaxis);
     TGraphErrors *gr4=extractCS(spec,target,angle,pset[3],0,pass,xaxis);
     TGraphErrors *gr5;
     if(spec=="hms")gr5=extractCS(spec,target,angle,pset[4],0,pass,xaxis);
     TGraphErrors *jmur=getJmu(target,angle,grd,grh,0);
     
     //gStyle->SetTitleFontSize(.1);
     if(xaxis=="xb") gr1->GetXaxis()->SetTitle("Bjorken X");
     if(xaxis=="ep") gr1->GetXaxis()->SetTitle("E' (GeV)");
     if(xaxis=="w2") gr1->GetXaxis()->SetTitle("W2 (GeV^2)");
     gr1->GetXaxis()->SetTitleSize(.06);
     gr1->GetXaxis()->CenterTitle();
     gr1->SetTitle("Data/F1F221 v0.995");
     gr1->GetYaxis()->SetRangeUser(ratiomin,ratiomax);
     //gr1->GetXaxis()->SetRangeUser(xbmin,xbmax);
     if(xaxis=="xb") gr1->GetXaxis()->SetLimits(xbmin,xbmax);
     if(xaxis=="ep") gr1->GetXaxis()->SetLimits(xmin-.1,xmax+.1);
     if(xaxis=="w2") gr1->GetXaxis()->SetLimits(wmin,wmax);
     gr1->GetXaxis()->SetLabelSize(.06);
     gr1->GetYaxis()->SetLabelSize(.06);
     gr1->GetXaxis()->SetNdivisions(-5);
 // gr1->GetYaxis()->SetNdivisions(5);
     jmur->SetMarkerStyle(33);
     gr1->SetMarkerStyle(33);
     gr2->SetMarkerStyle(33);
     gr3->SetMarkerStyle(33);
     gr4->SetMarkerStyle(33);
     if(spec=="hms")gr5->SetMarkerStyle(33);

     jmur->SetMarkerSize(1.);
     gr1->SetMarkerSize(1.);
     gr2->SetMarkerSize(1.);
     gr3->SetMarkerSize(1.);
     gr4->SetMarkerSize(1.);
     if(spec=="hms")gr5->SetMarkerSize(1.);

     jmur->SetMarkerColor(kBlack);
     gr1->SetMarkerColor(kRed);
     gr2->SetMarkerColor(kOrange);
     gr3->SetMarkerColor(kGreen+3);
     gr4->SetMarkerStyle(kViolet);
     if(spec=="hms")gr5->SetMarkerColor(kBlack);

     gr1->Draw("ap");
     // jmur->Draw("sp");
     gr1->Draw("sp");
     gr2->Draw("sp");
     gr3->Draw("sp");
     gr4->Draw("sp");
     if(spec=="hms")gr5->Draw("sp");
     pad1->cd();
   }
 TPaveText *t1=new TPaveText(.15,.18,.25,.2,"NDC");



 t1->AddText(pass.c_str());


 // t1->AddText("Casey's pass4");
 
 // t1->Draw();

 TLegend *leg;
 if(xaxis!="w2")
   {
     leg=new TLegend(.35,.6,.65,.9);
   }
 else leg=new TLegend(.17,.55,.58,.92);
 leg->SetBorderSize(1);
 // leg->AddEntry(gm,"F1F221","l");
 if(xaxis=="xb")
   {
     //     leg->AddEntry(gm_k,"KP DIS","l");
   }
 if(spec=="shms")leg->AddEntry(grcx1,"Data (E12-10-002)","p");
 if(spec=="hms"){
   double dum[3]={1,2,3};
   TGraphErrors *g=new TGraphErrors(3,dum,dum,dum,dum);
   g->SetMarkerStyle(33);
   g->SetMarkerSize(1.0);
   g->SetMarkerColor(kBlack);
   leg->AddEntry(gm_k1,"KP HYBRID","l");
   leg->AddEntry(g,"SHMS Hall C, E12-10-002","p");
   leg->AddEntry(gm_k2,"AKP17","l");
   leg->AddEntry(grcx1,"HMS Hall C, E12-10-002","p");
   //   leg->AddEntry(mrthn,"Hall A, MARATHON","p");
   leg->AddEntry(gm_cj,"CJ15","l");
   //   leg->AddEntry(gwhit,"SLAC, (Whitlow)","p");
   leg->AddEntry(gm,"F1F221","l");
   leg->AddEntry(gm_jam,"JAM","l");




 }
 // leg->AddEntry(grcx1,Form("E'=%2.2f GeV/c",psetf[0]),"p");
 // leg->AddEntry(grcx2,Form("E'=%2.2f GeV/c",psetf[1]),"p");
 // leg->AddEntry(grcx3,Form("E'=%2.2f GeV/c",psetf[2]),"p");
 // leg->AddEntry(grcx4,Form("E'=%2.2f GeV/c",psetf[3]),"p");
 // if(spec=="hms") leg->AddEntry(grcx5,Form("E'=%2.2f GeV/c",psetf[4]),"p");

 // leg->AddEntry(jmu,"JMU","p");
 leg->SetName(Form("leg_%s",spec.c_str()));
 if(target=="r"&&spec=="hms") leg->Draw();
 /// c1->SaveAs(Form("cs_%s%sdeg.pdf",target.c_str(),angle.c_str()));
  return;
}
