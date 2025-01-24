#include "src/readReport.cpp"
#include "src/getCharge.cpp"
#include "src/getHMSCharge.cpp"

void fixRange(TH1F *h){
  double min=0;
  double max=0;
  bool foundmin=0;
  bool foundmax=0;
  bool foundpeak=0;
  for(int i=0; i < h->GetNbinsX(); i++) {
    double val = h ->GetBinContent(i);
    if(val==h->GetMaximum())foundpeak=1;
    if(!foundmin && val>0){min=h->GetBinCenter(i);foundmin=1;}
    if(foundpeak==1 && foundmin && !foundmax && val==0){max=h->GetBinCenter(i);foundmax=1;}
  }
  cout<<"Found Min and max to be \t"<<min<<"  "<<max<<endl;
  Double_t fact=h->GetBinWidth(0);
  h->GetXaxis()->SetRangeUser(min-3*fact,max+3*fact);
  return;
}

void ratios(string tgt="h",string angle="21", string mom="3p3",string spec="shms",   bool rebin=true){
  //  gStyle->SetOptStat(0);
  string kin=tgt+angle+"deg"+mom; 
  cout << "Kinematic is : "<<kin<<endl;
  Int_t DRAW=1;
  bool scale_em=false;
  bool simpleScale=false;
  rebin=false;

  // histos
  // hm* MC
  // hd* Data
  // he* Dummy (empty)
  // hs* Data - dummy

  //*****MC Histograms*****
  TFile *fm=new TFile(Form("mcWtOut/pass100/%s_mcWt%s.root",spec.c_str(),kin.c_str()));
  if(!fm->IsOpen())return;
  TH1F *hmd=(TH1F*)fm->Get("delWt");
  TH1F *hmy=(TH1F*)fm->Get("yWt");
  TH1F *hmxp=(TH1F*)fm->Get("xpWt");
  TH1F *hmyp=(TH1F*)fm->Get("ypWt");
  TH1F *hmw2=(TH1F*)fm->Get("w2Wt2");
  //****Data Histograms***** 
  TFile *fd=new TFile(Form("dataYieldOut/pass100/%s_dataYield_%s.root",spec.c_str(),kin.c_str()));
  if(!fd->IsOpen())return;
  TH1F *hdd=(TH1F*)fd->Get("hdd");
  TH1F *hdd_stat=(TH1F*)hdd->Clone();
  hdd_stat->SetName("hdd_stat");
  TH1F *herr=(TH1F*)fd->Get("herr_pt2pt");
  TH1F *hdy=(TH1F*)fd->Get("hyd");
  TH1F *hdxp=(TH1F*)fd->Get("hxpd");
  TH1F *hdyp=(TH1F*)fd->Get("hypd");
  TH1F *hdw2=(TH1F*)fd->Get("hw2d_calc2");
  TH1F *hdw2_stat=(TH1F*)hdw2->Clone();
  hdw2_stat->SetName("hdw2_stat");
  TH1F *hsys=(TH1F*)fd->Get("herr_global");
  TH1F *hsysR=(TH1F*)fd->Get("herr_globalR");
  TH1F *herr_live=(TH1F*)fd->Get("herr_live");
  TH1F *herr_boil=(TH1F*)fd->Get("herr_boil");
  TH1F *hBoilCorr=(TH1F*)fd->Get("hBoilCorr");
  herr->Divide(hdd);
  hsys->Divide(hdd);
  hsysR->Divide(hdd);
  herr_live->Divide(hdd);
  herr_boil->Divide(hdd);
  hBoilCorr->Divide(hdd);
  double densityCorr=hBoilCorr->GetBinContent(30);
  cout << " The boiling factor is "<<densityCorr<<endl;
  if(rebin)
    {
      hdd->Rebin(3);
      hdd_stat->Rebin(3);
      hmd->Rebin(3);
      hdw2->Rebin(3);
      hdw2_stat->Rebin(3);
      hmw2->Rebin(3);
      herr->Rebin(3);
      hsys->Rebin(3);
      hsysR->Rebin(3);
      herr_live->Rebin(3);
      herr_boil->Rebin(3);
      hBoilCorr->Rebin(3);
      herr->Scale(1/3.);
      hsys->Scale(1/3.);
      hsysR->Scale(1/3.);
      herr_live->Scale(1/3.);
      herr_boil->Scale(1/3.);
      hBoilCorr->Scale(1/3.);
      densityCorr=hBoilCorr->GetBinContent(10);
    }
  //=============================
  // add point to point errors to statisical
  // (boiling, tracking, trigger error)
  //=============================
  int nbins=hdd->GetNbinsX();
  for(int i=0;i<nbins;i++){
    double content=hdd->GetBinContent(i);
    double statErr=hdd->GetBinError(i)/content;
    double ptErr=herr->GetBinContent(i);
    double total=sqrt(statErr*statErr+ptErr*ptErr)*content;
    if(content!=0){
      hdd->SetBinError(i,total);
      //     cout<< "Stat: "<<statErr<<"Pt2pt: "<<ptErr<<"Total: "<<total<<endl;
    } 
  }
  nbins=hdw2->GetNbinsX();
  for(int i=0;i<nbins;i++){
    double content=hdw2->GetBinContent(i);
    double ptErr=herr->GetBinContent(30);
    if(rebin)ptErr=herr->GetBinContent(10);
    double statErr=hdw2->GetBinError(i)/content;
    double total=sqrt(statErr*statErr+ptErr*ptErr)*content;
    if(content!=0){
      hdw2->SetBinError(i,total);
      //     cout<< "Stat: "<<statErr<<"Pt2pt: "<<ptErr<<"Total: "<<total<<endl;
    } 
  }

  //=============================
  // Charge Normalize data yields
  //=============================
  Float_t charge=0;
  if(spec=="shms")charge=getCharge(tgt,angle,mom);
  if(spec=="hms")charge=getHMSCharge(tgt,angle,mom);
  cout << "The Charge is: "<<charge<<endl;
  //charge=charge*densityCorr;   //this was to remove the density correction
  hdd->Scale(1./charge);
  hdd_stat->Scale(1./charge);
  hdy->Scale(1./charge);
  hdxp->Scale(1./charge);
  hdyp->Scale(1./charge);
  hdw2->Scale(1./charge);
  hdw2_stat->Scale(1./charge);

  // Clone data histos to subtract dummy
  TH1F *hsd=(TH1F*)hdd->Clone();
  TH1F *hsd_stat=(TH1F*)hdd_stat->Clone();
  TH1F *hsy=(TH1F*)hdy->Clone();
  TH1F *hsxp=(TH1F*)hdxp->Clone();
  TH1F *hsyp=(TH1F*)hdyp->Clone();
  TH1F *hsw2=(TH1F*)hdw2->Clone();
  TH1F *hsw2_stat=(TH1F*)hdw2->Clone();
  hsd->SetName("hsd");
  hsy->SetName("hsy");
  hsxp->SetName("hsxp");
  hsyp->SetName("hsyp");
  hsw2->SetName("hsw2");
  hsw2_stat->SetName("hsw2_stat");
  hsd_stat->SetName("hsd_stat");
  cout << "hsd->GetBinContent(25)"<<hsd->GetBinContent(25)<<endl;
  cout << "Charge normalized data"<<endl;
  cout <<"d data="<<hsd->Integral()<<endl;
  cout <<"y data="<<hsy->Integral()<<endl;
  cout <<"yp data="<<hsyp->Integral()<<endl;
  cout <<"xp data="<<hsxp->Integral()<<endl;
  cout <<"w2 data="<<hsw2->Integral()<<endl;

  //=============================
  // Dummy normalized data yields
  //=============================
  TH1F *hed, *hey, *hexp, *heyp, *hew2, *heerr, *hed_stat, *hew2_stat; 
  if(tgt!="c"){
    // Dummy charge
    Float_t charged=0;
    if(spec=="hms")
      {
	if(kin=="h21deg3p3"||kin=="d21deg3p3"){
	  charged+=readReport(1641,"BCM4C charge");
	}
	if(kin=="h21deg4p0"||kin=="d21deg4p0"){
	  charged+=readReport(1628,"BCM4C charge");
	  charged+=readReport(1629,"BCM4C charge");
	}
	if(kin=="h21deg4p5"||kin=="d21deg4p5"){
	  charged+=readReport(1609,"BCM4C charge");
	  charged+=readReport(1610,"BCM4C charge");
	}
	if(kin=="h21deg5p1"||kin=="d21deg5p1"){
	  charged+=readReport(1596,"BCM4C charge");
	}
	if(kin=="h21deg5p7"||kin=="d21deg5p7"){
	  charged+=readReport(1573,"BCM4C charge");
	}
      }
    if(spec=="shms")charged=getCharge("al",angle,mom);
    if(spec=="hms")charged=getHMSCharge("al",angle,mom);
    // Get dummy histos
    string dummyFile="al"+angle+"deg"+mom+"_"+tgt; 
    if(simpleScale)dummyFile="al"+angle+"deg"+mom+"_no"; 
    TFile *fdum=new TFile(Form("dataYieldOut/pass100/%s_dataYield_%s.root",spec.c_str(),dummyFile.c_str()));
    if(!fdum->IsOpen())return;
    //*****Dummy Histos *****
    hed=(TH1F*)fdum->Get("hdd");
    hed_stat=(TH1F*)hed->Clone();
    //     heerr=(TH1F*)fdum->Get("herr");
    heerr=(TH1F*)fdum->Get("herr_pt2pt");
    hey=(TH1F*)fdum->Get("hyd");
    hexp=(TH1F*)fdum->Get("hxpd");
    heyp=(TH1F*)fdum->Get("hypd");
    hew2=(TH1F*)fdum->Get("hw2d_calc2");
    hew2_stat=(TH1F*)hew2->Clone();
    hed->SetName("hed");
    hed_stat->SetName("hed_stat");
    heerr->SetName("heerr");
    hey->SetName("hey");
    hexp->SetName("hexp");
    heyp->SetName("heyp");
    hew2->SetName("hew2");
    hew2_stat->SetName("hew2_stat");
    if(rebin){
      hed->Rebin(3);
      hed_stat->Rebin(3);
      heerr->Rebin(3);
      heerr->Scale(1/3.);
      hew2->Rebin(3);
      hew2_stat->Rebin(3);
    }
    //=============================
    // add point to point errors to statisical
    //=============================
    
    heerr->Divide(hed);
    for(int i=0;i<nbins;i++){
      double content=hed->GetBinContent(i);
      double statErr=hed->GetBinError(i)/content;
      double ptErr=heerr->GetBinContent(i);
      double total=sqrt(statErr*statErr+ptErr*ptErr)*content;
      if(content!=0){
	//	 cout<< "Stat: "<<statErr<<"Pt2pt: "<<ptErr<<"Total: "<<total<<endl;
	hed->SetBinError(i,total);
      }
    }
    cout<<"Did I get here 2"<<endl;
    nbins=hew2->GetNbinsX();
    for(int i=0;i<nbins;i++){
      double content=hew2->GetBinContent(i);
      double ptErr=heerr->GetBinContent(30);
      if(rebin)ptErr=heerr->GetBinContent(10);
      double statErr=hew2->GetBinError(i)/content;
      double total=sqrt(statErr*statErr+ptErr*ptErr)*content;
      if(content!=0){
	hew2->SetBinError(i,total);
	//     cout<< "Stat: "<<statErr<<"Pt2pt: "<<ptErr<<"Total: "<<total<<endl;
      } 
    }
    hed->Scale(1./charged);
    hed_stat->Scale(1./charged);
    hey->Scale(1./charged);
    hexp->Scale(1./charged);
    heyp->Scale(1./charged);
    hew2->Scale(1./charged);
    hew2_stat->Scale(1./charged);
    cout << "Before scaling dummy hdd="<<  hdd->GetBinContent(30) << "  and hed= ";
    cout << hed->GetBinContent(30) <<"  and hsd= " << hsd->GetBinContent(30) <<endl;
    if(simpleScale)
      {
	if(tgt=="h"){   
	  hed->Scale(1/3.789);
	  hey->Scale(1/3.789);
	  hexp->Scale(1/3.789);
	  heyp->Scale(1/3.789);
	  hew2->Scale(1/3.789);
	}
	if(tgt=="d"){
	  hed->Scale(1/4.063);
	  hey->Scale(1/4.063);
	  hexp->Scale(1/4.063);
	  heyp->Scale(1/4.063);
	  hew2->Scale(1/4.063);
	}
      }
    cout << "Before dummy subtraction hdd="<<  hdd->GetBinContent(30) << "  and hed= ";
    cout << hed->GetBinContent(30) <<"  and hsd= " << hsd->GetBinContent(30) <<endl;
    // Dummy subtraction
    /*
    hsd->Add(hed,-1);
    hsd_stat->Add(hed_stat,-1);
    hsy->Add(hey,-1);
    hsxp->Add(hexp,-1);
    hsyp->Add(heyp,-1);
    hsw2->Add(hew2,-1);
    hsw2_stat->Add(hew2,-1);
    cout << "After dummy subtraction hdd="<<  hdd->GetBinContent(30) << "  and hed= ";
    cout << hed->GetBinContent(30) <<"  and hsd= " << hsd->GetBinContent(30) <<endl;
    */
    cout<<"After subtracting dummy"<<endl;
    cout <<"d data="<<hsd->Integral()<<endl;
    cout <<"y data="<<hsy->Integral()<<endl;
    cout <<"yp data="<<hsyp->Integral()<<endl;
    cout <<"xp data="<<hsxp->Integral()<<endl;
    cout <<"w2 data="<<hsw2->Integral()<<endl;
    // Boiling Correction
    hsd->Scale(1/densityCorr);
    hsd_stat->Scale(1/densityCorr); 
    hsy->Scale(1/densityCorr);
    hsyp->Scale(1/densityCorr);
    hsxp->Scale(1/densityCorr); 
    hsw2->Scale(1/densityCorr); 
    hsw2_stat->Scale(1/densityCorr); 
    cout<<"After applying boiling correction"<<endl;
    cout <<"d data="<<hsd->Integral()<<endl;
    cout <<"y data="<<hsy->Integral()<<endl;
    cout <<"yp data="<<hsyp->Integral()<<endl;
    cout <<"xp data="<<hsxp->Integral()<<endl; 
    cout <<"w2 data="<<hsw2->Integral()<<endl;
  }
  cout << "hsd->GetBinContent(25)"<<hsd->GetBinContent(25)<<endl;
  TH1F *hrd=(TH1F*)hsd->Clone();
  TH1F *hrd_stat=(TH1F*)hsd_stat->Clone();
  hrd->SetName("hrd");
  hrd_stat->SetName("hrd_stat");
  //  RATIO DATA/MC
  hrd->Divide(hmd);
  hrd_stat->Divide(hmd);
  cout << "hrd->GetBinContent(25)"<<hrd->GetBinContent(25)<<endl;  
  //  Foramatting
  hrd->SetTitle("Ratio Data/MC");
  hmd->SetTitle("dp/p (%)");
  hmy->SetTitle("Y_{tar} (cm)");
  hmxp->SetTitle("X'_{tar} (mr)");
  hmyp->SetTitle("Y'_{tar} (mr)");
  hmw2->SetTitle("W2 (GeV^{2})");

  hdd->SetLineColor(kBlack);
  hdy->SetLineColor(kBlack);
  hdxp->SetLineColor(kBlack);
  hdyp->SetLineColor(kBlack);
  hdw2->SetLineColor(kBlack);

  hmd->SetLineColor(kRed);
  hmy->SetLineColor(kRed);
  hmxp->SetLineColor(kRed);
  hmyp->SetLineColor(kRed);
  hmw2->SetLineColor(kRed);

  hmd->SetFillStyle(3002);
  hmd->SetFillColor(kRed);
  hmy->SetFillStyle(3002);
  hmy->SetFillColor(kRed);
  hmxp->SetFillStyle(3002);
  hmxp->SetFillColor(kRed);
  hmyp->SetFillStyle(3002);
  hmyp->SetFillColor(kRed);
  hmw2->SetFillStyle(3002);
  hmw2->SetFillColor(kRed);

  hsd->SetMarkerStyle(21);
  hsy->SetMarkerStyle(21);
  hsxp->SetMarkerStyle(21);
  hsyp->SetMarkerStyle(21);
  hsw2->SetMarkerStyle(21);

  hsd->SetMarkerSize(.5);
  hsy->SetMarkerSize(.5);
  hsxp->SetMarkerSize(.5);
  hsyp->SetMarkerSize(.5);
  hsw2->SetMarkerSize(.5);

  hsd->SetLineColor(kBlue);
  hsy->SetLineColor(kBlue);
  hsxp->SetLineColor(kBlue);
  hsyp->SetLineColor(kBlue);
  hsw2->SetLineColor(kBlue);

  hsd->SetMarkerColor(kBlue);
  hsy->SetMarkerColor(kBlue);
  hsxp->SetMarkerColor(kBlue);
  hsyp->SetMarkerColor(kBlue);
  hsw2->SetMarkerColor(kBlue);

  if(tgt!="c")
    {
      hed->SetLineColor(kMagenta);
      hey->SetLineColor(kMagenta);
      hexp->SetLineColor(kMagenta);
      heyp->SetLineColor(kMagenta);
      hew2->SetLineColor(kMagenta);
     
      hed->SetFillStyle(3003);
      hed->SetFillColor(kMagenta);
      hey->SetFillStyle(3003);
      hey->SetFillColor(kMagenta);
      hexp->SetFillStyle(3003);
      hexp->SetFillColor(kMagenta);
      heyp->SetFillStyle(3003);
      heyp->SetFillColor(kMagenta);
      hew2->SetFillStyle(3003);
      hew2->SetFillColor(kMagenta);
    }
  hrd->GetYaxis()->SetRangeUser(0.7,1.3); 
  // hrd->GetYaxis()->SetRangeUser(0.8,1.20); 
  // hsw2->GetXaxis()->SetRangeUser(-5,15); 
  // hsw2->GetXaxis()->SetRangeUser(0,15); 
  // hsw2->GetXaxis()->SetRangeUser(8,15); 
  hsd->GetXaxis()->SetRangeUser(-11,23); 
  hrd->GetXaxis()->SetRangeUser(-11,23); 
  if(spec=="hms")
    {
      hsd->GetXaxis()->SetRangeUser(-6.5,9.5); 
      hrd->GetXaxis()->SetRangeUser(-6.5,9.5); 
    }
  hsxp->GetXaxis()->SetRangeUser(-100,100); 
  hsyp->GetXaxis()->SetRangeUser(-50,50); 
  hmy->GetXaxis()->SetRangeUser(-6,6); 
  if(tgt=="c") hmy->GetXaxis()->SetRangeUser(-1.5,1.5); 
  // hsy->GetXaxis()->SetRangeUser(-5,5); 

  cout << "After formatting"<<endl;
  cout <<"d data="<<hsd->Integral()<<endl;
  cout <<"y data="<<hsy->Integral()<<endl;
  cout <<"yp data="<<hsyp->Integral()<<endl;
  cout <<"xp data="<<hsxp->Integral()<<endl;
  cout <<"w2 data="<<hsw2->Integral()<<endl;
  cout <<"d mc="<<hmd->Integral()<<endl;
  cout <<"y mc="<<hmy->Integral()<<endl;
  cout <<"yp mc="<<hmyp->Integral()<<endl;
  cout <<"xp mc="<<hmxp->Integral()<<endl;
  cout <<"w2 data="<<hmw2->Integral()<<endl;

  cout <<"d data/mc="<<hsd->Integral()/hmd->Integral()<<endl;
  cout <<"y data/mc="<<hsy->Integral()/hmy->Integral()<<endl;
  cout <<"yp data/mc="<<hsyp->Integral()/hmyp->Integral()<<endl;
  cout <<"xp data/mc="<<hsxp->Integral()/hmxp->Integral()<<endl;
  cout <<"w2 data/mc="<<hsw2->Integral()/hmw2->Integral()<<endl;

  Double_t fact=hsd->Integral()/hmd->Integral();
  // fact=1;
 
  if (DRAW==1)
    {
      if(scale_em)
	{
	  hsd->Scale(1/fact);
	  hsy->Scale(1/fact);
	  hsxp->Scale(1/fact);
	  hsyp->Scale(1/fact);
	  //	 hsw2->Scale(1/fact);
	}
      // Drawing
      TCanvas *c1=new TCanvas("c1","c1",600,900);
      c1->Divide(2,3);
      c1->cd(4);
      hsd->Draw("E");
      hmd->Draw("HIST E same");
      if(tgt!="c")hed->Draw("HIST E same");
     

      TLegend *leg=new TLegend(.3,.1,.7,.3);

      //     TLegend *leg=new TLegend(.3,.6,.7,.9);  standard

      //     TLegend *leg=new TLegend(.12,.62,.5,.85);

      //     TLegend *leg=new TLegend(.5,.6,.9,.9);
      leg->AddEntry(hmd,"MC","lfp");
      //     leg->AddEntry(hdd,"Data","lfp");
      leg->AddEntry(hsd,"Data - dummy","lfp");
      if(tgt!="c")leg->AddEntry(hed,"Dummy","lfp");
      //     leg->SetHeader(Form("Scale data/mc=%2.4f",fact));
      c1->cd(6);
      hrd->Draw();
      gStyle->SetOptFit(1);
      hrd->Fit("pol1");
      //     fixRange(hrd);
      //     hrd->Draw();
      leg->Draw("same");
      c1->cd(1);
      hmy->Draw("HIST E");
      //     fixRange(hsy);
      if(tgt!="c")hey->Draw("same HIST E");

      hsy->Draw("E same");
      c1->cd(3);
      hsxp->Draw("E");
      //          fixRange(hsxp);
      if(tgt!="c")hexp->Draw("same HIST E");
      hmxp->Draw("same HIST E");

      c1->cd(2);
      hsyp->Draw("E");
      //     fixRange(hsyp);
      if(tgt!="c")heyp->Draw("same HIST E");
      hmyp->Draw("same HIST E");

      c1->cd(5);
      hsw2->Draw("E");
      fixRange(hsw2);
      if(tgt!="c")hew2->Draw("same HIST E");
      hmw2->Draw("same HIST E");

      //     TPaveText *pt=new TPaveText(0.38,.3,.58,.55,"NDC");
      TPaveText *pt=new TPaveText(0.25,.45,.45,.55,"NDC");
      pt->AddText(kin.c_str());
      //     pt->AddText("Std; z target 3mm upstream");
      //     pt->AddText("Std; Spec. X' -3mr");
      //     pt->AddText("no_offset ROOTfiles");
      //     pt->SetFillColor(20);
      c1->cd(4);pt->Draw("BR");
      c1->SaveAs(Form("ratiosOut/pass401/%s_ratios%s.pdf",spec.c_str(),kin.c_str()));
      /*
      //    Figure for write up
      TCanvas *c2=new TCanvas("c2","c2",1200,600);
      c2->Divide(2,1);
      c2->cd(1);
      hsd->GetXaxis()->SetTitle("Delta");
      hsd->SetTitle("Charge Normalized Yield");
      hrd->GetXaxis()->SetTitle("Delta");
      hsd->Draw("E");
      hsd->GetYaxis()->SetRangeUser(0,1.6);
      hdd->Draw("E same");
      hmd->Draw("HIST E same");
      leg->Draw("same");
      if(tgt!="c")hed->Draw("HIST E same");
      pt->Draw("BR");
      c2->cd(2);
      hrd->Draw();
      c2->SaveAs("deltaHistos.pdf");
      */
    }


  TFile *oFile=new TFile(Form("ratiosOut/pass401/%s_ratios%s.root",spec.c_str(),kin.c_str()),"RECREATE");

  hdd->Write("hdd");
  hdd_stat->Write("hdd_stat");
  hdy->Write("hdy");
  hdxp->Write("hdxp");
  hdyp->Write("hdyp");
  hdw2->Write("hdw2");
  hdw2_stat->Write("hdw2_stat");
  hsys->Write("hsys");
  hsysR->Write("hsysR");
  herr_live->Write("herr_live");
  herr_boil->Write("herr_boil");
  if(tgt!="c")
    {
      hed->Write("hed");
      hed_stat->Write("hed_stat");
      hey->Write("hey");
      hexp->Write("hexp");
      heyp->Write("heyp");
      hew2->Write("hew2");
      hew2_stat->Write("hew2_stat");
    }

  hsd->Write("hsd");
  hsd_stat->Write("hsd_stat");
  hsy->Write("hsy");
  hsxp->Write("hsxp");
  hsyp->Write("hsyp");
  hsw2->Write("hsw2");
  hsw2->Write("hsw2_stat");

  hmd->Write("hmd");
  hmy->Write("hmy");
  hmxp->Write("hmxp");
  hmyp->Write("hmyp");
  hmw2->Write("hmw2");

  hrd->Write("hrd");
  hrd_stat->Write("hrd_stat");
  oFile->Close();
  if(tgt!="c")cout << "The very end     hdd="<<  hdd->GetBinContent(30) << "  and hed= " << hed->GetBinContent(30) <<"  and hsd= " << hsd->GetBinContent(30) <<endl;
  /* 
     cout<<fixed<<setprecision(4);
     cout <<endl<<"Delta Histigrams Integrals"<<endl;
     cout <<"Only trust this if DRAw=0"<<endl;
     cout << "RawData\t\tDummyt\t\tDummySub\t\tMC"<<endl;
     cout << hdd->Integral()<<"\t";
     cout << hed->Integral()<<"\t";
     cout << hsd->Integral()<<"\t";
     cout << hmd->Integral()<<endl;
  */
  return;
}
