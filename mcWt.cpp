#include "TH1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include "src/readReport.cpp"
//#include "src/getRadCorr.cpp"
#include "src/getRadCorrW2.cpp"
#include "src/getCerEff.cpp"
#include "src/getMom.cpp"
#include "src/getCharge.cpp"
#include "src/collCut.cpp"
#include "src/fidCut.cpp"

using namespace std;

//  Weights mc-single-arm output with a given model
//  Calculates the MC scale factor 
//  Outputs histograms of ytar, yp, xp, w2, and delta 
//  Output txt file with luminocities, charge, factors etc

//  Things to check
// * mc generation range (dxp, dyp, delup, down)
// * 
void mcWt(string tgt="h",string angle="21", string mom="3p3", string spec="shms"){

  string kin=tgt+angle+"deg"+mom;
  string kin_al="h"+angle+"deg"+mom;
  cout <<" The Kinematic is " << kin<<endl;

  Double_t pOffset=0.;
  Double_t mp = .9382723;
  Double_t mp2 = mp*mp;
  //  Grids with model
  TGraph2D *gr;   //born
  TGraph2D *gr2;  //born/rad
  TGraph2D *gr3;  //rad 
  Double_t dxp, dyp, delup, deldown, delCutHi, delCutLo, xpCut, ypCut, yCut;
  Double_t target, ebeam, hsec, thetac, tgtDenLen, tgtMass, offset,beamTheta;
  Double_t charge=0;
  ofstream oFile;
  ofstream outFile;
  outFile.open(Form("pass401_mcWt_%s.txt",spec.c_str()),ios::app | ios::out );


  if(spec=="shms")
    {
      dxp=60;//phi in MC xp_tar
      dyp=65;
      delup=25;
      deldown=-15;
      //match data cuts
      delCutLo = -10.;
      delCutHi = +22.;
      xpCut= 0.10;
      ypCut= 0.10;
      yCut= 10.;
    }

  if(spec=="hms")
    {
      //      dxp=110;//phi
      //            dyp=55;
      //      delup=15;
      //      deldown=-15;
      //deb
      //      dxp=100;
      //      dyp=50;
      //      delup=12;
      //      deldown=-12;
      //casey
      dxp=110;
      dyp=65;
      delup=15;
      deldown=-15;      

//match data cuts
      delCutLo = -6.;
      delCutHi = +9.;
      xpCut= 0.1;
      ypCut= 0.05;
      yCut= 10.;
    }
  //kin err
  ebeam=10.602;//*(1.0-0.001);
  //  charge=getCharge(tgt,angle,mom);
  //  Want counts/uA
  charge=1.;

  hsec=getMom(kin,spec);
  if(tgt=="alu"||tgt=="ald"){hsec=getMom(kin_al,spec);cout << "Alum "<< kin_al << endl;}
  cout << "The central momentum is "<<hsec<<endl;
  if(spec=="hms"&&hsec<5.5){  
    offset = -0.000276*pow(hsec,3) + 0.002585*pow(hsec,2) - 0.008697*hsec+1.0064;
    hsec=hsec*offset;
  }
  //kin err
  //  hsec=hsec*(1.-0.001);
  cout << "The central momentum is "<<hsec<<endl;
  if(angle=="39")thetac=38.975;
  if(angle=="33")thetac=32.975;
  if(angle=="29")thetac=28.99;
  if(angle=="25")thetac=24.98;
  if(angle=="21"&&spec=="shms")thetac=21.035;
  if(angle=="21"&&spec=="hms")thetac=20.995;

  cout << "The central angle is "<<thetac<<endl;
  if(angle=="59")thetac=58.98;

  if(tgt=="alu"||tgt=="ald")
    {
      gr=getRadCorrW2(tgt.c_str(),1,spec.c_str(),"v996t2");  //born
      gr->SetName("gr");
      gr2=getRadCorrW2(tgt.c_str(),2,spec.c_str(),"v996t2");  //born/rad
      gr2->SetName("gr2");
      gr3=getRadCorrW2(tgt.c_str(),3,spec.c_str(),"v996t2");  //rad
      gr3->SetName("gr3");
      tgtDenLen=.18155;
      tgtMass=26.98;
    }

  if(tgt=="c")
    {
      gr=getRadCorrW2("c",1,spec.c_str(),"v996t2");  //born
      gr->SetName("gr");
      gr2=getRadCorrW2("c",2,spec.c_str(),"v996t2");  //born/rad
      gr2->SetName("gr2");
      gr3=getRadCorrW2("c",3,spec.c_str(),"v996t2");  //rad
      gr3->SetName("gr3");
      tgtDenLen=2.19*(.5244/2.19);
      tgtMass=12.011;
    }

  if(tgt=="d")
    {
      gr=getRadCorrW2("d",1,spec.c_str(),"v996t2");  
      gr->SetName("gr");
      gr2=getRadCorrW2("d",2,spec.c_str(),"v996t2");  
      gr2->SetName("gr2");
      gr3=getRadCorrW2("d",3,spec.c_str(),"v996t2");  
      gr3->SetName("gr3");
      tgtDenLen=.167*9.9682*.996;
      tgtMass=2.014;
    }

  if(tgt=="h")
    {
      gr=getRadCorrW2("h",1,spec.c_str(),"v996t2");  
      gr->SetName("gr");
      gr2=getRadCorrW2("h",2,spec.c_str(),"v996t2");  
      gr2->SetName("gr2");
      gr3=getRadCorrW2("h",3,spec.c_str(),"v996t2");  
      gr3->SetName("gr3");
      tgtDenLen=0.0723*9.9659*.996;
      tgtMass=1.0079;
    }
  ////****  Beam is not along z ******
  cout << "Correcting central angle for beam angle"<<endl;
  cout << "Before: "<<thetac<<endl;
  //kin err
  //  beamTheta=0.00045+0.00025; //shooting beam right .45mr
  beamTheta=0.00045; //shooting beam right .45mr
  if(spec=="hms")beamTheta*=-1;
  thetac+=beamTheta*180./TMath::Pi();
  cout << "After: "<<thetac<<endl;
  //////////////////////////////////////
  Double_t thetacrad=thetac*TMath::Pi()/180;
  Double_t lumdata=tgtDenLen*6.022137e-10/tgtMass*charge/1.602177e-13; // targets/nb
  outFile <<kin<<"\t"<< ebeam <<"\t"<< thetac<<"\t"<<hsec <<"\t" << lumdata <<"\t";
  Float_t xfoc, yfoc, dxdz, dydz, ztarini, ytarini, delini, xptarini, yptarini;
  Float_t zrec, ytarrec, delrec, yptarrec, xptarrec, xtarini, xstop, ystop, fail_id;
  TString fmc = "mc/casey/"+spec+"_mc_"+kin+".root";
  //  fmc="/u/group/shms/wmhenry/mc-single-arm/worksim/"+spec+"_mc_"+kin+".root";
  //  if(tgt=="alu"||tgt=="ald")fmc="mc/aruni/"+spec+"_mc_"+kin+".root";
  if(tgt=="alu"||tgt=="ald")fmc="mc/deb/"+spec+"_mc_"+kin+".root";

  if(tgt=="alu")fmc="mc/deb/"+spec+"_"+angle+"deg_m"+mom+"_al_upstream.root";
  if(tgt=="ald")fmc="mc/deb/"+spec+"_"+angle+"deg_m"+mom+"_al_downstream.root";

  //  if(tgt=="alu")fmc="/w/hallc-scifs17exp/xem2/aruni/mc-reweight/output/mc-ntuples/mc145.root";
  //  if(tgt=="ald")fmc="/w/hallc-scifs17exp/xem2/aruni/mc-reweight/output/mc-ntuples/mc146.root";

  //fmc = "mc/testing/ELOG457/"+spec+"_mc_"+kin+".root";
  //  fmc = "mc/march2021/"+spec+"_"+kin+".root";
  // fmc = "mc/"+spec+"_"+kin+"_deb.root";
  // fmc = "mc/testing/"+spec+"_"+kin+".root";
  //  if(spec=="hms") fmc = "mc/rewind/"+spec+"_"+kin+".root";
  //  if(spec=="hms") fmc = "mc/testing/"+spec+"_"+kin+".root";
  //  if(spec=="hms") fmc = "/work/hallc/shms/wmhenry/mc-single-arm/worksim/"+spec+"_"+kin+".root";
  //    TString fmc = "mc/deb/shms_39deg_m1p3_h.root";
  //  TString fmc = "mc/abel/h39_1.3.root";
  TFile *fm=new TFile(fmc);
  fm->Print();


  TTree *trm;
  string arm;
  if(spec=="shms")
    {
      trm=(TTree*)fm->Get("h1411");
      //      trm=(TTree*)fm->Get("h9040");
      arm="p";
    }
  if(spec=="hms")
    {
      trm=(TTree*)fm->Get("h1");
      arm="h";
    }

  cout << "The MC files has "<<trm->GetEntries()<< " entries"<< endl;

  trm->SetBranchAddress(Form("%ssxfp",arm.c_str()), &xfoc);
  trm->SetBranchAddress(Form("%ssyfp",arm.c_str()), &yfoc);
  trm->SetBranchAddress(Form("%ssxpfp",arm.c_str()), &dxdz);
  trm->SetBranchAddress(Form("%ssypfp",arm.c_str()), &dydz);
  trm->SetBranchAddress(Form("%ssztari",arm.c_str()), &ztarini);
  trm->SetBranchAddress(Form("%ssytari",arm.c_str()), &ytarini);
  trm->SetBranchAddress(Form("%ssdeltai",arm.c_str()), &delini);
  trm->SetBranchAddress(Form("%ssyptari",arm.c_str()), &yptarini);
  trm->SetBranchAddress(Form("%ssxptari",arm.c_str()), &xptarini);
  trm->SetBranchAddress(Form("%ssztar",arm.c_str()), &zrec);
  trm->SetBranchAddress(Form("%ssytar",arm.c_str()), &ytarrec);
  trm->SetBranchAddress(Form("%ssdelta",arm.c_str()), &delrec);
  trm->SetBranchAddress(Form("%ssyptar",arm.c_str()), &yptarrec);
  trm->SetBranchAddress(Form("%ssxptar",arm.c_str()), &xptarrec);
  trm->SetBranchAddress(Form("%ssxtari",arm.c_str()), &xtarini);
  trm->SetBranchAddress("xsieve", &xstop);
  trm->SetBranchAddress("ysieve", &ystop);
  trm->SetBranchAddress("stop_id", &fail_id);

  TString fOut=Form("mcWtOut/pass401/%s_mcWt%s.root",spec.c_str(),kin.c_str());
  TFile *out=new TFile(fOut,"RECREATE");
  TTree *tree=new TTree("tree","Monte Carlo Weighted");
  cout << "opened two more files"<<endl;
  //Copy the original mc-single-arm branchs
  tree->Branch("xfoc",&xfoc);
  tree->Branch("yfoc",&yfoc);
  tree->Branch("dxdz",&dxdz);
  tree->Branch("dydz",&dydz);
  tree->Branch("ztarini",&ztarini);
  tree->Branch("ytarini",&ytarini);
  tree->Branch("delini",&delini);
  tree->Branch("yptarini",&yptarini);
  tree->Branch("xptarini",&xptarini);
  tree->Branch("zrec",&zrec);
  tree->Branch("ytarrec",&ytarrec);
  tree->Branch("delrec",&delrec);
  tree->Branch("yptarrec",&yptarrec);
  tree->Branch("xptarrec",&xptarrec);
  tree->Branch("xtarini",&xtarini);
  tree->Branch("xstop",&xstop);
  tree->Branch("ystop",&ystop);
  tree->Branch("fail_id",&fail_id);

  Float_t born_corr, rad, rci, hstheta, sigmac, q2, w2, csb_cx;
  Float_t hse, hsev, thetaini, sin2, nu, wt, wt_noCsb, xb, dt, phasespcor, phasespcorCos; 
  Int_t ngen=0;
  Int_t nacc=0;
  Float_t born=0;
  Double_t totalWt=0;
  Double_t sigave=0;
  // These will be the new branches
  tree->Branch("born",&born);
  tree->Branch("rad",&rad);
  tree->Branch("wt",&wt);
  tree->Branch("born_corr",&born_corr);
  tree->Branch("rci",&rci);
  tree->Branch("hsev",&hsev);
  tree->Branch("hse",&hse);
  tree->Branch("thetaini",&thetaini);
  tree->Branch("hstheta",&hstheta);
  tree->Branch("sigmac",&sigmac);
  tree->Branch("q2",&q2);
  tree->Branch("w2",&w2);
  tree->Branch("xb",&xb);
  tree->Branch("csb_cx",&csb_cx);
  tree->Branch("phaseSpaceCorr",&phasespcor);
  //  tree->Branch("phaseSpaceCorrCos",&phasespcorCos);

  //histos for comparisons
  TH1F *delWt=new TH1F("delWt","Monte Carlo Weighted delta",60,-30.,30.);
  TH1F *delWt_noCsb=new TH1F("delWt_noCsb","Monte Carlo Weighted delta no CSB",60,-30.,30.);  
  TH1F *hAvgTheta=new TH1F("hAvgTheta","Average Theta",60,-30.,30.);
  //  TH1F *delWt=new TH1F("delWt","Monte Carlo Weighted delta",15,-6.,9.);
  TH1F *xpWt=new TH1F("xpWt","Monte Carlo Weighted xp_tar",100,-100.,100.);
  TH1F *ypWt=new TH1F("ypWt","Monte Carlo Weighted yp_tar",100,-100.,100.);
  TH1F *yWt=new TH1F("yWt","Monte Carlo Weighted y_tar",334,-10,10);
  TH1F *w2Wt=new TH1F("w2Wt","Monte Carlo Weighted W2",375,-10,20);
  TH1F *w2Wt2=new TH1F("w2Wt2","Monte Carlo Weighted W2",720,-10,26);
  TH1F *w2Wt2_noCsb=new TH1F("w2Wt2_noCSB","Monte Carlo Weighted W2 no CSB",720,-10,26);  
  TH1F *xbWt=new TH1F("xbWt","Monte Carlo Weighted X_{B}",120,0,3.);
  //
  TH1F *centralBorn=new TH1F("avgBorn","#sigma_{born}(#delta,#theta_{central) ",32,-10.,22.);
  TH1F *htemp=new TH1F("htemp","Born for event with 9<#delta<10",100,0,10);
  TH2F *thetaBorn=new TH2F("thetaBorn","Born vs thetaini;#theta (rad);born (nb/GeV/sr)",50,.25,.50,250,0,100);
  TH2F *epBorn=new TH2F("epBorn","Born vs Eprime;E' (GeV);born (nb/GeV/sr)",50,0,8,250,0,100);
  TH2F *thetaRad=new TH2F("thetaRad","Rad vs thetaini;#theta (rad);born (nb/GeV/sr)",50,.25,.50,250,0,100);
  TH2F *epRad=new TH2F("epRad","Rad vs Eprime;E' (GeV);born (nb/GeV/sr)",50,0.,8.0,250,0,100);

  // Focal Plane Plots
  TH2F *mc_xVy=new TH2F("mc_xVy","x_fp vs y_fp; y_fp (cm); x_fp (cm)",100,-40.,40.0,100,-40.,40.);
  TH2F *mc_xpVyp=new TH2F("mc_xpVyp","xp_fp vs yp_fp; yp_fp (rad); xp_fp (rad)",100,-0.06,0.06,100,-0.1,0.1);
  TH2F *mc_xVxp=new TH2F("mc_xVxp","x_fp vs x_fp; xp_fp (rad); x_fp (cm)",100,-0.1,0.1,100,-40.,40.);
  TH2F *mc_ypVy=new TH2F("mc_ypVy","yp_fp vs y_fp; y_fp (cm); yp_fp (rad)",100,-40.,40.0,100,-0.06,0.06);

  TH2F *mc_yptarVytar=new TH2F("mc_yptarVytar","yp_tar vs y_tar; y_tar (cm); yp_tar (rad)",100,-6.,6.,100,-0.05,0.05);


  sin2 = sin(thetacrad/2.)*sin(thetacrad/2.);
  nu = ebeam - hsec;
  q2 = 4.*hsec*ebeam*sin2;
  w2= mp2 + 2.*mp*nu-q2;
  sigmac = gr->Interpolate(w2,thetacrad*180./TMath::Pi());// born
  TProfile* hp=new TProfile("hp","Profile of xB Born",100,0,1); 
  TProfile* deltaBornProf=new TProfile("deltaBornProf","Profile of delta Born",32,-10,22); 
  Int_t nEvents=trm->GetEntries();
  //  nEvents=5000;
  Int_t wtf=0;
  cout << "About to loop"<<endl;

  for (Int_t i=0; i<nEvents; i++)
    {

      wt=0;
      if(i%250000==0)cout<<i<<endl;
      trm->GetEntry(i);
      if(fail_id==0 ){//I forget why i need this condition
      //	cout << "Event:\t"<<i<<"\t";
      //	cout << "Before ytarrec:\t"<<ytarrec;
	//      ytarrec=ytarrec + gRandom->Gaus(0,.4);
	//	cout << "/tAfter ytarrec:\t"<<ytarrec<<endl;
      }

      //Calculate E' and apply offset
      hsev=hsec*(1. + delini/100.);
      hse=hsec*(1. + delrec/100.);
      //      hsev=hsev*(1+pOffset);
      //Calculate theta
      if(spec=="shms"){
	thetaini = acos(cos(thetacrad + yptarini)*cos(xptarini));
	hstheta = acos(cos(thetacrad + yptarrec)*cos(xptarrec));
      }
      if(spec=="hms"){
	thetaini = acos(cos(thetacrad - yptarini)*cos(xptarini));
	hstheta = acos(cos(thetacrad - yptarrec)*cos(xptarrec));
      }

      //call cross sections with initial quantities
      sin2 = sin(thetaini/2.)*sin(thetaini/2.);
      nu = ebeam - hsev;
      q2 = 4.*hsev*ebeam*sin2;
      w2= mp2 + 2.*mp*nu - q2;
      xb=q2/2./mp/nu;
      //            born = gr->Interpolate(hsev,thetaini*180./TMath::Pi());// born
      //    rci = gr2->Interpolate(hsev,thetaini*180./TMath::Pi());// born/rad
      //    rad = gr3->Interpolate(hsev,thetaini*180./TMath::Pi());// rad
      born = gr->Interpolate(w2,thetaini*180./TMath::Pi());// born
      rci = gr2->Interpolate(w2,thetaini*180./TMath::Pi());// born/rad
      rad = gr3->Interpolate(w2,thetaini*180./TMath::Pi());// rad
      //      cout << "born;rci;rad"<<"\t"<<born<<"\t"<<rci<<"\t"<<rad<<endl;
      //      cout << w2 <<"\t"<<thetaini*180./3.14<<"\t"<<born<<endl;
      thetaBorn->Fill(thetaini,born);
      thetaRad->Fill(thetaini,rad);
      epBorn->Fill(hsev,born);
      epRad->Fill(hsev,rad);
      hp->Fill(xb,born);
      phasespcor=pow(1+pow(xptarini,2)+pow(yptarini,2),-1.5);
      born_corr=born*phasespcor;

      //Add CSB
      csb_cx=0;
      if(spec=="shms"|| (spec=="hms"&&angle=="59") )
	{
	  Double_t p0=-2.09 * thetaini*180./TMath::Pi() +12.47;
	  Double_t p1=0.2 * thetaini*180./TMath::Pi() -0.6338;
	  csb_cx=exp(p0)*(exp(p1*(ebeam-hsev))-1.);

	  //***************************************************************
	  // Testing how much a difference using **************************
	  // separate csb_cx for lh2 and ld2 makes*************************
	  // From Gyang's studies *****************************************
	  //	  if(tgt=="d")csb_cx=exp(-41.69)*(exp(4.46*(10.602-hsev))-1.);//ld2
	  //	  if(tgt=="h")csb_cx=exp(-48.32)*(exp(5.18*(10.602-hsev))-1.);//ld2
	  //***************************************************************

	  if(tgt=="d")csb_cx=2*csb_cx;
	  if(tgt=="c")csb_cx=12*csb_cx*19.32/((890.4+769.1)/2); //need to add. (use rad length)
	  //	  if(tgt=="alu"||tgt=="ald")csb_cx=27*csb_cx*19.32/((890.4+769.1)/2); //need to add. (use rad length)
	  if(tgt=="alu"||tgt=="ald")csb_cx=0; //need to add. (use rad length)

	}

      //********************************
      //      csb_cx=0;
      //********************************

      //calculate W2 with reconstructed quantitites for ratios with MC
      sin2 = sin(hstheta/2.)*sin(hstheta/2.);
      nu = ebeam -hse;
      q2 = 4.*hse*ebeam*sin2;
      w2= mp2 + 2.*mp*nu - q2;

      //      if(i%25000==0)cout<<"hsev,thetaini,born,rci\t";
      //	if(born==0)
      //	  {
      //		wtf++;
      //		cout<<hsev<<"\t"<<thetaini*180./TMath::Pi()<<"\t"<<born<<"\t"<<rci<<"\t"<<rad<<"\t"<<wtf<<endl;
      // 	  }
      bool fid=fidCut(xfoc, yfoc, dxdz, dydz);
      bool coll=collCut(xptarrec, yptarrec, delrec, ytarrec);
      if(abs(xptarini)<dxp && abs(yptarini)<dyp && delini>deldown && delini<delup && born>0)
	{
	  //	  cout << "wt1" << wt <<endl;
	  wt = (born/rci+csb_cx)*phasespcor;
	  wt_noCsb=(born/rci)*phasespcor;
	  //	  cout << "wt  " << wt <<endl;
	  sigave+=born;
	  ngen++;
	  //	  wt=(rad)/phasespcor;
	}

      bool pi_thsh_cut=true;
      //      if(mom=="5p7" && spec=="hms" && tgt=="h") pi_thsh_cut = w2>1.2;
      //      if(mom=="5p1" && spec=="shms" && tgt=="h") pi_thsh_cut = w2>1.2;
      if(tgt=="h") pi_thsh_cut = w2>1.2;
      if(fail_id==0 && delrec<delCutHi && delrec >delCutLo && pi_thsh_cut)// && coll && fid)
	{
	  if(abs(xptarrec)<xpCut && abs(yptarrec)<ypCut && abs(ytarrec)<yCut)
	  //	  if(abs(xptarrec)<xpCut && abs(yptarrec)<ypCut && abs(ytarrec-1.5)<0.5)
	    {
	      if(wt==0)
	     {
	       wtf++;
	       //	       cout<<"WARNING    "<<w2<<"\t"<<thetaini*180./TMath::Pi()<<"\t"<<born<<"\t"<<rci<<"\t"<<rad<<"\t"<<wtf<<endl;
	       //	       wt=.1;
	       //		 wt=100;
	     }
	      //	      cout << "wt  " << wt <<"\t"<< born <<"\t"<< rci <<"\t"<< csb_cx <<"\t"<< phasespcor << endl;
	      
	   deltaBornProf->Fill(delrec,born);
	   delWt->Fill(delrec,wt);//delta weighted histogram
	   delWt_noCsb->Fill(delrec,wt_noCsb);//delta weighted histogram	   
	   hAvgTheta->Fill(delrec,hstheta*wt*180/TMath::Pi());
	   xpWt->Fill(xptarrec*1000.,wt);//mr
	   ypWt->Fill(yptarrec*1000.,wt);//mr
	   yWt->Fill(ytarrec,wt);
	   w2Wt->Fill(w2,wt);//W2 weighted histogram
	   w2Wt2->Fill(w2,wt);
	   w2Wt2->Fill(w2,wt_noCsb);	   
	   xbWt->Fill(xb,wt);
	   //focal plane variables
	   mc_xVy->Fill(yfoc,xfoc,wt);
	   mc_xpVyp->Fill(dydz,dxdz,wt);
	   mc_xVxp->Fill(dxdz,xfoc,wt);
	   mc_ypVy->Fill(yfoc,dydz,wt);
	   mc_yptarVytar->Fill(ytarrec,yptarrec,wt);

	   totalWt+=wt;
	   nacc++;
	   if(9.<delrec<10.)htemp->Fill(born);
	 }
       }
     tree->Fill();
    }
  
  for (Int_t i=0;i<32;i++){
      Double_t ene=hsec*(1. + (-9.5+i)/100.);
      sin2 = sin(thetacrad/2.)*sin(thetacrad/2.);
      nu = ebeam - ene;
      q2 = 4.*ene*ebeam*sin2;
      w2= mp2 + 2.*mp*nu-q2;
      born = gr->Interpolate(w2,thetacrad*180./TMath::Pi());// born
      centralBorn->Fill(-9.5+i,born);
      cout << "Delta: "<<-9.5+i<<"  born: "<<born<<endl;
  }
  TH1F *h100=(TH1F*)epBorn->ProfileX();
  TH1F *h101=(TH1F*)thetaBorn->ProfileX();
  h100->SetTitle("E' Born Profile");
  h101->SetTitle("Theta Born Profile");
  fm->Close();
  cout <<"from loop ngen= "<<ngen<<endl;
  // *************************************
  Double_t ngen2=ngen;
  ngen=nEvents;
  // *************************************

  cout <<"ngen= "<<ngen<<endl;
  cout << "Integrals after filling "<<nacc<<" events, del, yp, xp, y, w2"<<endl;
  cout << setprecision(8);
  cout << "Total Wt: "<<totalWt<< endl;
  cout << delWt->Integral() << endl;
  cout << ypWt->Integral() << endl;
  cout << xpWt->Integral() << endl;
  cout << yWt->Integral() << endl;
  cout << w2Wt->Integral() << endl;

  Double_t dep=(delup-deldown)/100.*hsec;
  Double_t phase_space=4.0*dxp*dyp*dep/1000.;
  Double_t lummc=1/phase_space*ngen*1000.;
  Double_t fract=lumdata*phase_space/ngen/1000.;
  cout << "rho*t: " << tgtDenLen <<endl;
  cout << "target mass: " << tgtMass <<endl;
  cout << "charge: " << charge <<endl;
  cout << "luminosity data: "<<lumdata<<"targets/cm^2"<<endl;

  cout << "Generation Limits, delta, xp, yp:" << dep << "\t" << dxp <<"\t" << dyp << endl;
  cout << "Generated Events" << ngen << endl;
  cout << "luminosity MC: "<<lummc<<"  targets/cm^2"<<endl;
  cout << "MC factor (lumdata/lumMC):  "<< fract <<endl;
  cout << "4*dxp*dyp: "<<4.0*dxp*dyp/1000./1000.<<endl;
  cout << "dep" << "\t"<< "phase_space" << "\t"<< "fract" << endl;
  cout << dep << "\t"<< phase_space << "\t"<< fract << endl;
  cout << "Norm: "<< fract << endl;
  cout << "Charge: "<< charge << endl;

  outFile << lummc << "\t"<<charge<<"\t"<<ngen2 <<"\t"<< nacc <<"\t"<< sigave/ngen <<"\t"<< fract/charge << "\t"<<totalWt/nacc<<"\t";

  delWt->Scale(fract/charge);
  hAvgTheta->Scale(fract/charge);
  ypWt->Scale(fract/charge);
  xpWt->Scale(fract/charge);
  yWt->Scale(fract/charge);
  w2Wt->Scale(fract/charge);
  w2Wt2->Scale(fract/charge);
  xbWt->Scale(fract/charge);
  mc_xVy->Scale(fract/charge);
  mc_xpVyp->Scale(fract/charge);
  mc_xVxp->Scale(fract/charge);
  mc_ypVy->Scale(fract/charge);

  cout << "Integrals after sacling, del, yp, xp, y, w2"<<endl;
  cout << delWt->Integral() << endl;
  cout << ypWt->Integral() << endl;
  cout << xpWt->Integral() << endl;
  cout << yWt->Integral() << endl;
  cout << w2Wt->Integral() << endl;
  //  cout <<"delWt->Integral: "<< delWt->Integral() <<endl;
  cout << "MC Yield: "<<totalWt*fract/charge<<endl; 
  outFile<<totalWt*fract/charge<<endl; 
  delWt->Sumw2();
  centralBorn->Write();
  htemp->Write();
  delWt->Write();
  hAvgTheta->Write();
  xpWt->Write();
  ypWt->Write();
  yWt->Write();
  w2Wt->Write();
  w2Wt2->Write();
  xbWt->Write();
  deltaBornProf->Write();
  thetaBorn->Write();
  epBorn->Write();
  thetaRad->Write();
  epRad->Write();
  mc_xVy->Write();
  mc_xpVyp->Write();
  mc_xVxp->Write();
  mc_ypVy->Write();
  mc_yptarVytar->Write();
  hp->Write();
  tree->Write();
  h100->Write();
  h101->Write();
  out->Close();
  return;
}

