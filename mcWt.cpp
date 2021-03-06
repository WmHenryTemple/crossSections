
#include <iostream>
#include <iomanip>
//#include "src/readReport.cpp"
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
void mcWt(string tgt="h",string angle="21", string mom="2p7", string spec="shms"){

  string kin=tgt+angle+"deg"+mom;
  cout <<" The Kinematic is " << kin<<endl;
  Double_t pOffset=0.;
  Double_t mp = .9382723;
  Double_t mp2 = mp*mp;
  //  Grids with model
  TGraph2D *gr;   //born
  TGraph2D *gr2;  //born/rad
  TGraph2D *gr3;  //rad 
  Double_t dxp, dyp, delup, deldown, delCutHi, delCutLo, xpCut, ypCut, yCut;
  Double_t target, ebeam, hsec, thetac, tgtDenLen, tgtMass, offset, offseti;
  Double_t charge=0;
  ofstream oFile;

  ofstream outFile;
  outFile.open(Form("test_mcWt_%s.txt",spec.c_str()),ios::app | ios::out );

  if(spec=="shms")
    {
      dxp=60;
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
      dxp=110;
      dyp=55;
      delup=15;
      deldown=-15;
      //deb
      dxp=100;
      dyp=50;
      delup=12;
      deldown=-12;
      //match data cuts
      delCutLo = -6.;
      delCutHi = +9.;
      xpCut= 0.1;
      ypCut= 0.05;
      yCut= 10.;
    }

  ebeam=10.602;
  //  charge=getCharge(tgt,angle,mom);
  //  Want counts/uA
  charge=1.;
  hsec=getMom(kin,spec);
  cout << "The central momentum is "<<hsec<<endl;
  if(spec=="hms"&&hsec<5.5){  
    offset = -0.000276*pow(hsec,3) + 0.002585*pow(hsec,2) - 0.008697*hsec+1.0064;
    hsec=hsec*offset;
  }
  cout << "The central momentum is "<<hsec<<endl;
  if(angle=="39")thetac=38.975;
  if(angle=="33")thetac=32.975;
  if(angle=="29")thetac=28.99;
  if(angle=="25")thetac=24.98;
  if(angle=="21"&&spec=="shms")thetac=21.035;
  if(angle=="21"&&spec=="hms")thetac=20.995;

  cout << "The central angle is "<<thetac<<endl;
  if(angle=="59")thetac=58.98;

  if(tgt=="c")
    {
      gr=getRadCorrW2("c",1,spec.c_str());  //born
      gr->SetName("gr");
      gr2=getRadCorrW2("c",2,spec.c_str());  //born/rad
      gr2->SetName("gr2");
      gr3=getRadCorrW2("c",3,spec.c_str());  //rad
      gr3->SetName("gr3");
      tgtDenLen=2.19*(.5244/2.19);
      tgtMass=12.011;
    }

  if(tgt=="d")
    {
      gr=getRadCorrW2("d",1,spec.c_str());  
      gr->SetName("gr");
      gr2=getRadCorrW2("d",2,spec.c_str());  
      gr2->SetName("gr2");
      gr3=getRadCorrW2("d",3,spec.c_str());  
      gr3->SetName("gr3");
      tgtDenLen=.167*9.9682*.996;
      tgtMass=2.014;
    }

  if(tgt=="h")
    {
      gr=getRadCorrW2("h",1,spec.c_str());  
      gr->SetName("gr");
      gr2=getRadCorrW2("h",2,spec.c_str());  
      gr2->SetName("gr2");
      gr3=getRadCorrW2("h",3,spec.c_str());  
      gr3->SetName("gr3");
      tgtDenLen=0.0723*9.9659*.996;
      tgtMass=1.0079;
    }

  Double_t thetacrad=thetac*TMath::Pi()/180;
  Double_t lumdata=tgtDenLen*6.022137e-10/tgtMass*charge/1.602177e-13; // targets/nb
  cout << "lumdata: "<<lumdata<<"targets/cm^2"<<endl;
  outFile <<kin<<"\t"<< ebeam <<"\t"<< thetac<<"\t"<<hsec <<"\t" << lumdata <<"\t";
  Float_t xfoc, yfoc, dxdz, dydz, ztarini, ytarini, delini, xptarini, yptarini;
  Float_t zrec, ytarrec, delrec, yptarrec, xptarrec, xtarini, xstop, ystop, fail_id;
  //      TString fmc = "mc/casey/shms_"+kin+".root";
  //  TString fmc = "mc/"+spec+"_"+kin+"_deb.root";
  TString fmc = "mc/"+spec+"_"+kin+".root";
  //    TString fmc = "mc/deb/shms_39deg_m1p3_h.root";
  //  TString fmc = "mc/abel/h39_1.3.root";
  TFile *fm=new TFile(fmc);
  fm->Print();


  TTree *trm;
  string arm;
  if(spec=="shms")
    {
      trm=(TTree*)fm->Get("h1411");
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

  TString fOut=Form("mcWtOut/pass22/%s_mcWt%s.root",spec.c_str(),kin.c_str());
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
  Float_t hse, hsev, thetaini, sin2, nu, wt, xb, dt, phasespcor, phasespcorCos; 
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
  //  TH1F *delWt=new TH1F("delWt","Monte Carlo Weighted delta",15,-6.,9.);
  TH1F *xpWt=new TH1F("xpWt","Monte Carlo Weighted xp_tar",100,-100.,100.);
  TH1F *ypWt=new TH1F("ypWt","Monte Carlo Weighted yp_tar",100,-100.,100.);
  TH1F *yWt=new TH1F("yWt","Monte Carlo Weighted y_tar",334,-10,10);
  TH1F *w2Wt=new TH1F("w2Wt","Monte Carlo Weighted W2",375,-10,20);
  TH1F *xbWt=new TH1F("xbWt","Monte Carlo Weighted X_{B}",120,0,3.);
  //
  TH1F *centralBorn=new TH1F("avgBorn","#sigma_{born}(#delta,#theta_{central) ",32,-10.,22.);
  TH1F *htemp=new TH1F("htemp","Born for event with 9<#delta<10",100,0,10);
  TH2F *thetaBorn=new TH2F("thetaBorn","Born vs thetaini;#theta (rad);born (nb/GeV/sr)",50,.25,.50,250,0,100);
  TH2F *epBorn=new TH2F("epBorn","Born vs Eprime;E' (GeV);born (nb/GeV/sr)",50,0,8,250,0,100);
  TH2F *thetaRad=new TH2F("thetaRad","Rad vs thetaini;#theta (rad);born (nb/GeV/sr)",50,.25,.50,250,0,100);
  TH2F *epRad=new TH2F("epRad","Rad vs Eprime;E' (GeV);born (nb/GeV/sr)",50,0.,8.0,250,0,100);
  sin2 = sin(thetacrad/2.)*sin(thetacrad/2.);
  nu = ebeam - hsec;
  q2 = 4.*hsec*ebeam*sin2;
  w2= mp2 + 2.*mp*nu-q2;
  sigmac = gr->Interpolate(w2,thetacrad*180./TMath::Pi());// born
  TProfile* hp=new TProfile("hp","Profile of xB Born",100,0,1); 
  TProfile* deltaBornProf=new TProfile("deltaBornProf","Profile of delta Born",32,-10,22); 
  Int_t nEvents=trm->GetEntries();
  //        nEvents=100000;
  Int_t wtf=0;
  cout << "About to loop"<<endl;

  for (Int_t i=0; i<nEvents; i++)
    {

      wt=0;
      if(i%250000==0)cout<<i<<endl;
      trm->GetEntry(i);
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
      w2= mp2 + 2.*mp*nu-q2;
      xb=q2/2./mp/nu;
      //            born = gr->Interpolate(hsev,thetaini*180./TMath::Pi());// born
      //    rci = gr2->Interpolate(hsev,thetaini*180./TMath::Pi());// born/rad
      //    rad = gr3->Interpolate(hsev,thetaini*180./TMath::Pi());// rad
      born = gr->Interpolate(w2,thetaini*180./TMath::Pi());// born
      rci = gr2->Interpolate(w2,thetaini*180./TMath::Pi());// born/rad
      rad = gr3->Interpolate(w2,thetaini*180./TMath::Pi());// rad
      thetaBorn->Fill(thetaini,born);
      thetaRad->Fill(thetaini,rad);
      epBorn->Fill(hsev,born);
      epRad->Fill(hsev,rad);
      hp->Fill(xb,born);
      phasespcor=pow(1+pow(xptarini,2)+pow(yptarini,2),-1.5);
      born_corr=born*phasespcor;
      //Add CSB
      Double_t p0=-2.09 * thetaini*180./TMath::Pi() +12.47;
      Double_t p1=0.2 * thetaini*180./TMath::Pi() -0.6338;
      csb_cx=exp(p0)*(exp(p1*(ebeam-hsev))-1.);
      if(tgt=="d")csb_cx=2*csb_cx;
      if(tgt=="c")csb_cx=0; //need to add. (use rad length)
      csb_cx=0;
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
	 wt = (born/rci+csb_cx)*phasespcor;
	 sigave+=born;
	 ngen++;
	  //	  wt=(rad)/phasespcor;
       }
     //     if(w2>1.2&&fail_id==0 && delrec<delCutHi && delrec >delCutLo)// && coll && fid)
          if(fail_id==0 && delrec<delCutHi && delrec >delCutLo)// && coll && fid)
       {
	 if(abs(xptarrec)<xpCut && abs(yptarrec)<ypCut && abs(ytarrec)<yCut)
	 {
	   if(wt==0)
	     {
	       wtf++;
	       cout<<"WARNING    "<<w2<<"\t"<<thetaini*180./TMath::Pi()<<"\t"<<born<<"\t"<<rci<<"\t"<<rad<<"\t"<<wtf<<endl;
	       //	       wt=.1;
	       //		 wt=100;
	     }
	   deltaBornProf->Fill(delrec,born);
	   delWt->Fill(delrec,wt);
	   xpWt->Fill(xptarrec*1000.,wt);//mr
	   ypWt->Fill(yptarrec*1000.,wt);//mr
	   yWt->Fill(ytarrec,wt);
	   w2Wt->Fill(w2,wt);
	   xbWt->Fill(xb,wt);
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
  cout << "4*dxp*dyp: "<<4.0*dxp*dyp/1000./1000.<<endl;
  cout << "dep" << "\t"<< "phase_space" << "\t"<< "fract" << endl;
  cout << dep << "\t"<< phase_space << "\t"<< fract << endl;
  cout << "Norm: "<< fract << endl;
  cout << "Charge: "<< charge << endl;

  outFile << lummc << "\t"<<charge<<"\t"<<ngen2 <<"\t"<< nacc <<"\t"<< sigave/ngen <<"\t"<< fract/charge << "\t"<<totalWt/nacc<<"\t";

  delWt->Scale(fract/charge);
  ypWt->Scale(fract/charge);
  xpWt->Scale(fract/charge);
  yWt->Scale(fract/charge);
  w2Wt->Scale(fract/charge);
  xbWt->Scale(fract/charge);

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
  xpWt->Write();
  ypWt->Write();
  yWt->Write();
  w2Wt->Write();
  xbWt->Write();
  deltaBornProf->Write();
  thetaBorn->Write();
  epBorn->Write();
  thetaRad->Write();
  epRad->Write();
  hp->Write();
  tree->Write();
  h100->Write();
  h101->Write();
  out->Close();
  return;
}

