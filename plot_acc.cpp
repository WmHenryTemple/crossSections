#include "ExtractCSfromAcc.cpp"
//#include "extractCS.cpp"
#include "src/formatHist.cpp"
#include "src/xsecTable.cpp"
//#include "src/getRadCorr.cpp"
//#include "src/readReport.cpp"
#include "src/getMomList.cpp"
#include "src/getCSB.cpp"
#include "src/getCJ.cpp"
#include "src/getKull.cpp"
#include "src/getKP.cpp"
#include "jmu/getJmu.cpp"

void plot_cx(string target="d", string angle="39", string spec="shms", string pass="pass105", float xbmin=0.0, float xbmax=1.0, float cxmin=0.0, float cxmax = 50.0, float ratiomin = 0.7, float ratiomax =1.3, string xaxis = "xb", int numPad=2 ){
  gStyle->SetTitleY(1.);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleX(.2);
  gStyle->SetTitleW(.6);
 Float_t delta, ratio, err, ep, model_k, modeld, modelh, xmin, xmax;
 string pset[5];
 float psetf[5];
 getMomList(angle,pset,psetf,spec);
 xmin=psetf[0]*(.9);
 xmax=psetf[3]*(1.22);
 if(spec=="hms")xmax=psetf[4]*(1.1);

  ////////////////////////////////////////
  ////////   MODEL   /////////////////////
  ////////////////////////////////////////

 vector <float> mx;
 vector <float> my_kph;  //kp hybrid
 vector <float> my_kpd;  //kp dis 
 vector <float> my;    //CJ15
 double wmin,wmax;

 Double_t thetac=getAngle(angle,spec);
  string version="v996t2";
  cout << "Going to use version "<<version<<" for "<<pass<<endl;
  TGraph2D *grh=getRadCorrW2("h",1,spec,version);  
  grh->SetName("grh");
  TGraph2D *grd=getRadCorrW2("d",1,spec,version);  
  grd->SetName("grd");

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
 if(xaxis=="xb") gm->GetXaxis()->SetTitle("Bjorken x");
 if(xaxis=="ep") gm->GetXaxis()->SetTitle("E' (GeV)");
 if(xaxis=="w2") gm->GetXaxis()->SetTitle("W2 (GeV^2)");
 gm->GetXaxis()->SetLimits(xbmin,xbmax);
 if(xaxis=="ep")gm->GetXaxis()->SetLimits(xmin-.1,xmax+.1);
 if(xaxis=="w2")gm->GetXaxis()->SetLimits(wmin,wmax);

 gm->SetLineColor(kBlue);
 gm_cj->SetLineColor(kRed);
 gm_k->SetLineColor(kGreen+3);
 gm_k1->SetLineColor(kOrange-6);
 gm_k2->SetLineColor(kMagenta);
 
 // gm->SetLineStyle(1);
 gm_cj->SetLineStyle(2);
 gm_k->SetLineStyle(3);
 gm_k1->SetLineStyle(4);
 gm_k2->SetLineStyle(6);
 int lw=2;
 if (numPad==2)lw=1;

 gm->SetLineWidth(lw);
 gm_cj->SetLineWidth(lw);
 gm_k->SetLineWidth(lw);
 gm_k1->SetLineWidth(lw);
 gm_k2->SetLineWidth(lw);


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
 gm->GetXaxis()->CenterTitle();
 gm->GetYaxis()->CenterTitle();
 // gm->GetXaxis()->SetNdivisions(-507 );
 // gm->GetYaxis()->SetNdivisions(-507 );
 // gm->GetXaxis()->SetLabelSize(0.06);// size of axis value font
 // TCanvas *c1=new TCanvas();
 TPad *pad1, *pad2;
 if(numPad==2){
   pad1 = new TPad("pad1","pad1",0,0.52,1,1);
   pad2 = new TPad("pad2","pad2",0,0.02,1,0.5);
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
 // if(spec=="shms")
gm->Draw("la");
// if(target=="r" && spec=="shms" && xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}


  ////////////////////////////////////////
  ////////   Data   /////////////////////
  ////////////////////////////////////////
  double mrthn_x[7]={.195, .225, .255, .285, .315, .345, .375};
  double mrthn_w[7]={12.2,11.7,11.3,10.9,10.5,10.1,9.63};
  if(xaxis=="w2"){
  for(int i=0;i<7;i++)mrthn_x[i]=mrthn_w[i];
  }
  double mrthn_y[7]={1.725, 1.697, 1.674, 1.656, 1.629, 1.588, 1.544};
  double mrthn_ye[7]={.015, .014, .014, .016, .016, .016, .023};
  for(int i=0;i<7;i++)mrthn_y[i]=mrthn_y[i]/2;
  for(int i=0;i<7;i++)mrthn_ye[i]=mrthn_ye[i]/2;
  TGraphErrors *mrthn=new TGraphErrors(7,mrthn_x,mrthn_y,0,mrthn_ye);
  mrthn->SetMarkerStyle(21);
  mrthn->SetMarkerSize(1);
  mrthn->SetMarkerColor(kMagenta+2);
  //  if(target=="r"&& angle=="21")mrthn->Draw("sp");
  /////////////////////////////////////////////////////////////////////
  cout << spec <<"\t"<< target <<"\t"<< angle <<"\t"<< pset[0]<<"\t"<<1<<"\t"<<pass<<"\t"<<xaxis<<endl;

  TGraphErrors *grcx1a=ExtractCSfromAcc(spec,target,angle,pset[0],3,pass,xaxis);
  grcx1a->SetName("grcx1");
  TGraphErrors *grcx2a=ExtractCSfromAcc(spec,target,angle,pset[1],3,pass,xaxis);
  TGraphErrors *grcx3a=ExtractCSfromAcc(spec,target,angle,pset[2],3,pass,xaxis);
  TGraphErrors *grcx4a=ExtractCSfromAcc(spec,target,angle,pset[3],3,pass,xaxis);
  TGraphErrors *grcx5a;
  if(spec=="hms")grcx5a=ExtractCSfromAcc(spec,target,angle,pset[4],3,pass,xaxis);
  
  TGraphErrors *gr1=ExtractCSfromAcc(spec,target,angle,pset[0],0,pass,xaxis);
  TGraphErrors *gr2=ExtractCSfromAcc(spec,target,angle,pset[1],0,pass,xaxis);
  TGraphErrors *gr3=ExtractCSfromAcc(spec,target,angle,pset[2],0,pass,xaxis);
  TGraphErrors *gr4=ExtractCSfromAcc(spec,target,angle,pset[3],0,pass,xaxis);
  TGraphErrors *gr5;
  if(spec=="hms")gr5=ExtractCSfromAcc(spec,target,angle,pset[4],0,pass,xaxis);

  TGraphErrors *grMC1=extractCS(spec,target,angle,pset[0],0,"pass308",xaxis);
  TGraphErrors *grMC2=extractCS(spec,target,angle,pset[1],0,"pass308",xaxis);
  TGraphErrors *grMC3=extractCS(spec,target,angle,pset[2],0,"pass308",xaxis);
  TGraphErrors *grMC4=extractCS(spec,target,angle,pset[3],0,"pass308",xaxis);
  TGraphErrors *grMC5;
  if(spec=="hms")grMC5=extractCS(spec,target,angle,pset[4],0,"pass308",xaxis);


  ///////////////////////////////////////////////////////
  string kin=spec+"_"+target+angle+"deg";
  TGraphErrors *jmu=getJmu(target,angle,grd,grh,1);

 jmu->SetMarkerStyle(33);

 grcx1a->SetMarkerStyle(27);
 grcx2a->SetMarkerStyle(27);
 grcx3a->SetMarkerStyle(27);
 grcx4a->SetMarkerStyle(27);
 if(spec=="hms")grcx5a->SetMarkerStyle(27);

 grcx1a->SetMarkerSize(.6);
 grcx2a->SetMarkerSize(.6);
 grcx3a->SetMarkerSize(.6);
 grcx4a->SetMarkerSize(.6);
 if(spec=="hms")grcx5a->SetMarkerSize(.6);

 jmu->SetMarkerColor(kBlack);

 grcx1a->SetMarkerColor(kRed);
 grcx2a->SetMarkerColor(kBlue);
 grcx3a->SetMarkerColor(kGreen+3);
 grcx4a->SetMarkerColor(kViolet);

 if(spec=="hms")grcx5a->SetMarkerColor(kRed);

 grcx1a->Draw("sp");
 grcx2a->Draw("sp");
 grcx3a->Draw("sp");
 grcx4a->Draw("sp");
 if(spec=="hms") grcx5a->Draw("sp");

 grcx1a->Draw("sp");
 grcx2a->Draw("sp");
 grcx3a->Draw("sp");
 grcx4a->Draw("sp");
 if(spec=="hms") grcx5a->Draw("sp");
 

  
 //  gm->Draw("same");

 // if(target=="r"&&spec=="hms"&&xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}
 // if(target=="r"&&spec=="hms"&&xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}
 if(numPad==2)
   {
     pad2->cd();

     
     //gStyle->SetTitleFontSize(.1);
     if(xaxis=="xb") gr1->GetXaxis()->SetTitle("Bjorken X");
     if(xaxis=="ep") gr1->GetXaxis()->SetTitle("E' (GeV)");
     if(xaxis=="w2") gr1->GetXaxis()->SetTitle("W2 (GeV^2)");
     gr1->GetXaxis()->SetTitleSize(.06);
     gr1->GetXaxis()->CenterTitle();
     gr1->SetTitle("Data/F1F221 v0.996t2");
     gr1->GetYaxis()->SetRangeUser(ratiomin,ratiomax);
     //gr1->GetXaxis()->SetRangeUser(xbmin,xbmax);
     if(xaxis=="xb") gr1->GetXaxis()->SetLimits(xbmin,xbmax);
     if(xaxis=="ep") gr1->GetXaxis()->SetLimits(xmin-.1,xmax+.1);
     if(xaxis=="w2") gr1->GetXaxis()->SetLimits(wmin,wmax);
     gr1->GetXaxis()->SetLabelSize(.06);
     gr1->GetYaxis()->SetLabelSize(.06);
     gr1->GetXaxis()->SetNdivisions(-5);
 // gr1->GetYaxis()->SetNdivisions(5);
     gr1->SetMarkerStyle(27);
     gr2->SetMarkerStyle(27);
     gr3->SetMarkerStyle(27);
     gr4->SetMarkerStyle(27);
     if(spec=="hms")gr5->SetMarkerStyle(27);

     gr1->SetMarkerSize(.6);
     gr2->SetMarkerSize(.6);
     gr3->SetMarkerSize(.6);
     gr4->SetMarkerSize(.6);
     if(spec=="hms")gr5->SetMarkerSize(.6);

     gr1->SetMarkerColor(kRed);
     gr2->SetMarkerColor(kBlue);
     gr3->SetMarkerColor(kGreen+3);
     gr4->SetMarkerColor(kViolet);
     if(spec=="hms")gr5->SetMarkerColor(kBlack);

     grMC1->SetMarkerStyle(4);
     grMC2->SetMarkerStyle(4);
     grMC3->SetMarkerStyle(4);
     grMC4->SetMarkerStyle(4);
     if(spec=="hms")grMC5->SetMarkerStyle(4);

     grMC1->SetMarkerSize(.6);
     grMC2->SetMarkerSize(.6);
     grMC3->SetMarkerSize(.6);
     grMC4->SetMarkerSize(.6);
     if(spec=="hms")grMC5->SetMarkerSize(.6);

     gr1->Draw("ap");
     gr1->Draw("sp");
     gr2->Draw("sp");
     gr3->Draw("sp");
     gr4->Draw("sp");
     if(spec=="hms")gr5->Draw("sp");

     grMC1->Draw("sp");
     grMC2->Draw("sp");
     grMC3->Draw("sp");
     grMC4->Draw("sp");
     if(spec=="hms")grMC5->Draw("sp");

     gr1->Draw("sp");
     gr2->Draw("sp");
     gr3->Draw("sp");
     gr4->Draw("sp");
     if(spec=="hms")gr5->Draw("sp");
     pad1->cd();
   }
 // if(target=="r"&&spec=="shms"&&xaxis=="xb"){gm_k->Draw("same");gm_cj->Draw("same");gm_k1->Draw("same");gm_k2->Draw("same");}

 TPaveText *t1=new TPaveText(.3,.8,.7,.9,"NDC");



 t1->AddText(pass.c_str());
 // if(target=="r")t1->AddText("D/H ratio");
 // if(target=="d")t1->AddText("Duet. xsec");
 // t1->AddText(Form("%sdeg",angle.c_str()));
 // t1->AddText("Casey's pass4");
 
 t1->Draw();

 TLegend *leg;
 if(xaxis!="w2")
   {
     leg=new TLegend(.35,.6,.65,.9);
   }
 else leg=new TLegend(.15,.5,.55,.9);
 leg->SetBorderSize(1);
 leg->AddEntry(gm,"F1F221","l");
 if(xaxis=="xb")
   {
   }
 if(spec=="shms")leg->AddEntry(grcx1a,"SHMS Data (E12-10-002)","p");
 if(spec=="hms"){
   double dum[3]={1,2,3};
   TGraphErrors *g=new TGraphErrors(3,dum,dum,dum,dum);
   g->SetMarkerStyle(27);
   g->SetMarkerSize(.6);
   g->SetMarkerColor(kBlack);
   leg->AddEntry(grcx1a,"HMS Hall C, E12-10-002","p");
   //   leg->AddEntry(g,"SHMS Hall C, E12-10-002","p");
   //   leg->AddEntry(mrthn,"Hall A, MARATHON","p");
 }
 // leg->AddEntry(grcx1,Form("E'=%2.2f GeV/c",psetf[0]),"p");
 // leg->AddEntry(grcx2,Form("E'=%2.2f GeV/c",psetf[1]),"p");
 // leg->AddEntry(grcx3,Form("E'=%2.2f GeV/c",psetf[2]),"p");
 // leg->AddEntry(grcx4,Form("E'=%2.2f GeV/c",psetf[3]),"p");
 // if(spec=="hms") leg->AddEntry(grcx5,Form("E'=%2.2f GeV/c",psetf[4]),"p");

 // leg->AddEntry(jmu,"JMU","p");
 leg->SetName(Form("leg_%s",spec.c_str()));
 // if(target=="r") leg->Draw();
 /// c1->SaveAs(Form("cs_%s%sdeg.pdf",target.c_str(),angle.c_str()));
 t1->Draw();
  return;
}
