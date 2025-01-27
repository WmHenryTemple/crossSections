#include <iostream>
#include <iomanip>
#include "TGraph2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TSystem.h"
#include "src/readReport.cpp"
#include "src/getRadCorr.cpp"
#include "src/getCerEffDelta.cpp"
#include "src/getCSBerr.cpp"
#include "src/getCalEff.cpp"
#include "src/getLivetime.cpp"
#include "src/getKineUnc.cpp"
#include "src/getPionContamination.cpp"
#include "src/fidCut.cpp"
#include "src/collCut.cpp"
#include "src/getRadCorrW2.cpp"
using namespace std;

void dataYield(Int_t run=2868, Double_t ngcCut=2, Double_t betaMin =0.5, Double_t betaMax=1.5, 
	       Double_t deltaMin=-10., Double_t deltaMax=22., Double_t minEdep=0.7, Double_t curCut=5., TString scaleDummy="h",TString fname="pass1_run2868.root"){
  // flags
  bool positron=false;
  bool use_saturation_correction=true;
  bool use_delta_correction=true;
  Double_t target=readReport(run,"target");
  bool use_w2_cut = (target==1.01) || (target>25. && scaleDummy=="h") ;
  use_w2_cut=false;                  

  //cuts
  string arm;
  double xpCut, ypCut, yCut;
  if(spec=="shms"){
    arm="P";
    xpCut= 0.10;
    ypCut= 0.10;
    yCut= 10.;
  }    
  if(spec=="hms"){
    arm="H";
    xpCut= 0.1;
    ypCut= 0.05;
    yCut= 10.;
  }
  
  // ELOG 336
  Double_t entr_fact=.23213;
  Double_t exit_fact=.29571;
  if(scaleDummy=="d"){
      entr_fact=.20118;
      exit_fact=.29107;
    }

  //determine spec, from run number
  string spec="shms";
  if(run<2200)spec="hms";

  //variables
  Double_t beta, delta, etracknorm, ngc, curr, phd, thd, xfp, yfp, xpfp, ypfp, xCer, yCer, xb, el_lo;
  Double_t  q2, w2,cerEff, calEff, mom, xd, yd, goode=0, goode_corr=0, boilCorr, errBoil, wt=0, sime=0,terr_pt2pt=0, terr_glob=0, piC=0;
  Double_t dipole=0;
  TString froot, report, fmc;

  //get corrections 
  TF1* fcer=getCerEffDelta(target);
  TF1 *pionC=getPionContamination(run);
  Double_t charge=readReport(run,"BCM4C charge");
  Double_t livetime=getLivetime(run,"tlt");
  Double_t trackEff=readReport(run,"tr eff");
  Double_t trigEff=readReport(run,"trig eff");
  if (trigEff<.1)trigEff=1;
  Double_t psFact=readReport(run,"Ps2 fact");
  Double_t currentAvg=readReport(run,"BCM4C cut current");

  //positron corrections
  if(psFact<0)positron=true;
  if(positron)livetime=readReport(run,"ps3 clt")/100.;
  if(positron)psFact=readReport(run,"Ps3 fact");

  //kinematic quantities
  Double_t sin2, nu, q2_calc, w2_calc, hse, hstheta, offset;
  Double_t ebeam=10.602;
  Double_t mp = .9382723;
  Double_t mp2 = mp*mp;
  Double_t hsec=readReport(run,"mom");
  cout << "The central momentum is "<<hsec<<endl;

  //hms saturation correction
  if(spec=="hms"&&hsec<5.5){  
    offset = -0.000276*pow(hsec,3) + 0.002585*pow(hsec,2) - 0.008697*hsec+1.0064;
    if(use_saturation_correction)hsec=hsec*offset;
  }
  cout << "The corrected central momentum is "<<hsec<<endl;

  // adjust spec. theta for beam angle
  Double_t thetac=readReport(run,"Spec Angle");
  cout << "The central angle is "<<thetac<<endl;
  Double_t beamTheta=0.00045; //shooting beam right .45mr
  thetac+=beamTheta*180./TMath::Pi();
  cout << "The central angle is "<<thetac<<endl;
  Double_t thetacrad=thetac*TMath::Pi()/180;

  // calculate boiling corrections
  double avgboilCorr=0;
  double avgerrBoil=0;
  double p2perrBoil=0;
  Double_t h_boil =0.0255;
  Double_t h_boil_err =0.0074;
  Double_t d_boil =0.0309;
  Double_t d_boil_err =0.0084;
  Double_t wt_corr = 1.000;
  if(target==1.01){//hydrogen
    boilCorr=1.-currentAvg/100.*h_boil; // ~0.98 +/-
    errBoil=currentAvg/100.*h_boil_err/boilCorr; // ~ 0.0032 / 0.98 +/-
    avgboilCorr=1- 46.94/100.*h_boil; // ~0.98
    avgerrBoil=46.94/100.*h_boil_err/avgboilCorr;  
    p2perrBoil=abs(46.94-currentAvg)/100.*h_boil_err/boilCorr; 
  }
  if(target==2.01){//deuterium
    boilCorr=1.-currentAvg/100.*d_boil; 
    errBoil=currentAvg/100.*d_boil_err/boilCorr; 
    avgboilCorr=1- 49.84/100.*d_boil;
    avgerrBoil=49.84/100.*d_boil_err/avgboilCorr;  
    p2perrBoil=abs(49.84-currentAvg)/100.*d_boil_err/boilCorr; 
  } 
  if(target>2.01){//solid targets
    boilCorr=1; 
    errBoil=0.; 
    avgerrBoil=0;
    p2perrBoil=0;
  }
  // additional correction for deuterium (22.4 K vs 22.0 K)
  if( (spec=="hms" && run >= 1879 && target==2.01) || (spec=="shms" && run >= 2808 && target==2.01) ){
    wt_corr=1/1.006;
    cout << "Correcting density because Deut temp = 22.4 K "<<endl;
    cout << "Yields will go up to match MC where rho = 22.0 K "<<endl;
  }
  cout << "The boiling before the "<<wt_corr<<" density correction is" << boilCorr<<endl;
  boilCorr=boilCorr*wt_corr;
  cout << "The boiling after density correction is" << boilCorr<<endl;

  // Flat corrections
  Double_t scale = (Double_t)1/(livetime)/trackEff/trigEff*psFact;

  // These are relative errors 
  Double_t errCer=.003; //temp.  will get from lookup table
  Double_t errCal=.002; // 
  Double_t errTgL=.015; // from Meekins 
  Double_t errCharge=.005; // 
  Double_t errTrack=.002; //from Deb's residual plot see Winter User Meetinf talk
  Double_t errTrig=.0003;  
  Double_t errPion=.001;
  if(abs(thetac-38.975)<.2)errPion=.002;
  if(abs(thetac+58.98)<.2)errPion=.003;
  Double_t errLive=getLivetime(run,"tlte")/livetime;

  //Output Rootfiles
  TFile *oFile=new TFile("dataYieldOut/pass500/"+fname,"RECREATE");
  TTree *tree=new TTree("tree","Data");
  TTree *tree2=new TTree("tree2","Run Eff.");
  tree2->Branch("boilCorr",&boilCorr);
  tree2->Branch("livetime",&livetime);
  tree2->Branch("trackEff",&trackEff);
  tree2->Branch("trigEff",&trigEff);
  tree2->Branch("psFact",&psFact);
  tree2->Branch("scale",&scale);
  tree2->Fill();
  tree->Branch("calEff",&calEff);
  tree->Branch("piC",&piC);
  tree->Branch("cerEff",&cerEff);
  tree->Branch(Form("%s.gtr.dp",arm.c_str()), &delta);
  tree->Branch(Form("%s.gtr.ph",arm.c_str()), &phd);
  tree->Branch(Form("%s.gtr.thd",arm.c_str()), &thd);
  tree->Branch(Form("%s.gtr.y",arm.c_str()), &yd);
  tree->Branch(Form("%s.kin.W2",arm.c_str()), &w2);
  tree->Branch(Form("%s.kin.Q2",arm.c_str()), &q2);
  tree->Branch("w2_calc", &w2_calc);
  tree->Branch("wt",&wt);

  //get data rootfile
  if(spec=="shms"){
    froot = Form("/lustre/expphy/cache/hallc/E12-10-002/abishek/realpass-3e-shms-data/shms_replay_production_%d_-1.root",run);
  if(spec=="hms"){  
    if(run==1608)froot = Form("/lustre/expphy/cache/hallc/E12-10-002/abishek/realpass-3b-hms-data/hms_replay_production_%d_-1.root",run);
    else froot = Form("/lustre/expphy/cache/hallc/E12-10-002/abishek/realpass-3d-hms-data/hms_replay_production_%d_-1.root",run);
  }
  if (gSystem->AccessPathName(froot)==0)//this if goes to the end of code
    {
      TFile *f=new TFile(froot);
      f->Print();
      TTree *tr=(TTree*)f->Get("T");
      cout << "Setting Branch Addresses"<<endl;
      tr->SetBranchAddress(Form("%s.gtr.beta",arm.c_str()), &beta);
      tr->SetBranchAddress(Form("%s.gtr.dp",arm.c_str()), &delta);
      tr->SetBranchAddress(Form("%s.gtr.p",arm.c_str()), &mom);
      tr->SetBranchAddress(Form("%s.gtr.ph",arm.c_str()), &phd);
      tr->SetBranchAddress(Form("%s.gtr.th",arm.c_str()), &thd);
      tr->SetBranchAddress(Form("%s.gtr.x",arm.c_str()), &xd);
      tr->SetBranchAddress(Form("%s.gtr.y",arm.c_str()), &yd);
      tr->SetBranchAddress(Form("%s.dc.x_fp",arm.c_str()), &xfp);
      tr->SetBranchAddress(Form("%s.dc.y_fp",arm.c_str()), &yfp);
      tr->SetBranchAddress(Form("%s.dc.xp_fp",arm.c_str()), &xpfp);
      tr->SetBranchAddress(Form("%s.dc.yp_fp",arm.c_str()), &ypfp);
      if(spec=="shms")tr->SetBranchAddress(Form("%s.dc.InsideDipoleExit",arm.c_str()), &dipole);
      tr->SetBranchAddress(Form("%s.kin.W2",arm.c_str()), &w2);
      tr->SetBranchAddress(Form("%s.kin.Q2",arm.c_str()), &q2);
      tr->SetBranchAddress(Form("%s.kin.x_bj",arm.c_str()), &xb);
      tr->SetBranchAddress(Form("%s.cal.etracknorm",arm.c_str()), &etracknorm);
      tr->SetBranchAddress(Form("%s.bcm.bcm4c.AvgCurrent",arm.c_str()), &curr);
      if(spec=="shms")tr->SetBranchAddress("P.ngcer.npeSum", &ngc);
      if(spec=="hms"){
	tr->SetBranchAddress("H.cer.npeSum", &ngc);
	tr->SetBranchAddress("T.hms.hEL_LO_tdcTime", &el_lo);
      }
      cout << "Done setting Branch Addresses"<<endl;

      // Define Histograms
      Int_t nbins=60;
      Double_t minBin=-30.;
      Double_t maxBin=30.;
      TH1D *hyld=new TH1D("hyld","Data Raw",nbins,minBin,maxBin);
      TH1D *hdd=new TH1D("hdd","Data Weighted",nbins,minBin,maxBin);
      TH1D *hdd2=new TH1D("hdd2","Statistical + pt2pt",nbins,minBin,maxBin);
      TH1D *hdd3=new TH1D("hdd3","nothing ",nbins,minBin,maxBin);
      TH1D *hdd4=new TH1D("hdd4","Data lt, ps, trk, boil",nbins,minBin,maxBin);
      TH1D *hdd5=new TH1D("hdd5","Data lt, ps, trk, boil, trg",nbins,minBin,maxBin);
      TH1D *hdd6=new TH1D("hdd6","Data lt, ps, trk, boil, trg, cal",nbins,minBin,maxBin);
      TH1D *hdd7=new TH1D("hdd7","Data lt, ps, trk, boil, trg, cal, cer",nbins,minBin,maxBin);
      TH1D *hAvgTheta=new TH1D("hAvgTheta","Average Theta per Bin",nbins,minBin,maxBin);
      TH1D *hAvgDelta=new TH1D("hAvgDelta","Average Delta per Bin",nbins,minBin,maxBin);
      TH1D *hBoilCorr=new TH1D("hBoilCorr","Average Boiling Correction",nbins,minBin,maxBin);
      TH1D *heff=new TH1D("heff","Efficiency",nbins,minBin,maxBin);
      TH1D *heffcal=new TH1D("heffcal","Calo. Efficiency",nbins,minBin,maxBin);
      TH1D *heffcer=new TH1D("heffcer","Cer. Efficiency",nbins,minBin,maxBin);
      TH1D *herrcer=new TH1D("herrcer","Cer. Error",nbins,minBin,maxBin);
      TH1D *heffpion=new TH1D("heffpion","Pion Contamination",nbins,minBin,maxBin);
      TH1D *herr_pt2pt=new TH1D("herr_pt2pt","Point to point error",nbins,minBin,maxBin);
      TH1D *herr_boil=new TH1D("herr_boil","p2p boiling error",nbins,minBin,maxBin);
      TH1D *herr_live=new TH1D("herr_live","p2p livetime error",nbins,minBin,maxBin);
      TH1D *herr_track=new TH1D("herr_track","p2p track eff error",nbins,minBin,maxBin);
      TH1D *herr_trig=new TH1D("herr_trig","p2p trig. eff error",nbins,minBin,maxBin);
      TH1D *herrCSB=new TH1D("herrCSB","CSB error",nbins,minBin,maxBin);
      TH1D *herrKin=new TH1D("herrKin","Kinematic (th,e',Eb) error",nbins,minBin,maxBin);
      TH1D *herrKinRatio=new TH1D("herrKinRatio","Kinematic (th,e',Eb) error on D/H",nbins,minBin,maxBin);
      TH1D *herr_global=new TH1D("herr_global","Total Band error",nbins,minBin,maxBin);
      TH1D *herr_globalR=new TH1D("herr_globalR","Total Band error D/H",nbins,minBin,maxBin);
      TH1D *herrTot=new TH1D("herrTot","Total sys error",nbins,minBin,maxBin);
      TH1D *hxpd=new TH1D("hxpd","Data xptar",100,-100,100);
      TH1D *hypd=new TH1D("hypd","Data yptar",100,-100,100);
      TH1D *hxd=new TH1D("hxd","Data x tar",100,-1.,1.);
      TH1D *hyd=new TH1D("hyd","Data y tar",334,-10,10);
      TH1D *hw2d_calc=new TH1D("hw2d_calc","Data W2 Calc",375,-10,20);//375
      TH1D *hw2d=new TH1D("hw2d","Data W2",1440,-10,26);      
      TH1D *hw2d_calc2=new TH1D("hw2d_calc2","Data W2 Calc",1440,-10,26);//375
      TH1D *hq2d=new TH1D("hq2d","Data Q2",500,-10,50);
      TH1D *hq2d_calc=new TH1D("hq2d_calc","Data Q2 Calc",500,-10,50);
      TH1D *hcerr=new TH1D("hcerr","Cer Eff",100,.9,1.0);      
      TH1D *hpion=new TH1D("hpion","Pion Contamination",200,0,.2);      
      TH1D *hcal=new TH1D("hcal","Cal Eff",100,.995,1.);      
      TH1D *hxb=new TH1D("hxb","xb Good Events",120,0,3);     
      TH1D *hmom=new TH1D("hmom","hmom",80,0.5,3.2);     
      // Pion Contaimination
      TH1D *hmom1=new TH1D("hmom1","hmom1",500,0,7);     
      TH1D *hmom2=new TH1D("hmom2","hmom2",500,0,7);     
      TH1D *hmom3=new TH1D("hmom3","hmom3",500,0,7);     
      TH1D *hcal_e1=new TH1D("hcal_e1","hcal_e1",300,0,2);     
      TH1D *hcal_e2=new TH1D("hcal_e2","hcal_e2",300,0,2);     
      TH1D *hcal_e3=new TH1D("hcal_e3","hcal_e3",300,0,2);     
      TH1D *hcal_pi1=new TH1D("hcal_pi1","hcal_pi1",300,0,2);     
      TH1D *hcal_pi2=new TH1D("hcal_pi2","hcal_pi2",300,0,2);     
      TH1D *hcal_pi3=new TH1D("hcal_pi3","hcal_pi3",300,0,2);     
      // Dummy
      TH2D *hdumFact=new TH2D("hdumFact","Dummy Scale factor vs ytar",50,-10,10,50,.1,.3);
      // Focal Plane Plots
      TH2F *xVy=new TH2F("xVy","x_fp vs y_fp; y_fp (cm); x_fp (cm)",100,-40.,40.0,100,-40.,40.);
      TH2F *xpVyp=new TH2F("xpVyp","xp_fp vs yp_fp; yp_fp (rad); xp_fp (rad)",100,-0.06,0.06,100,-0.1,0.1);
      TH2F *xVxp=new TH2F("xVxp","x_fp vs x_fp; xp_fp (rad); x_fp (cm)",100,-0.1,0.1,100,-40.,40.);
      TH2F *ypVy=new TH2F("ypVy","yp_fp vs y_fp; y_fp (cm); yp_fp (rad)",100,-40.,40.0,100,-0.06,0.06);
      TH2F *yptarVytar=new TH2F("yptarVytar","yp_tar vs y_tar; y_tar (cm); yp_tar (rad)",100,-6,6,100,-0.05,0.05);
      TH2F *yield4acc=new TH2F("yield4acc","yield; theta; delta",30,-65,65,60,-30,30);
      heff->Sumw2();
      hdd->Sumw2();
      //HMS delta correction
      double p1 = 0.0001307595;
      double p2 =-0.0005277879;
      double p3 = 0.0000598111;
      double p4 = 0.0000086922;
      double p5 =-0.0000001957;
      //SHMS acceptance correction
      double p00 = 1.00156;
      double p11 = -0.002473; 
      double p22 = -1.54588e-05;
      double p33 = 6.63986e-06;

      //  MAIN EVENT LOOP
      Int_t nEvents = tr->GetEntries();
      cout << "There are "<<nEvents<<" events"<<endl;
      //      nEvents=1000;
      for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) 
	{
	  if(iEvent%100000==0)cout<<iEvent<<endl;
	  tr->GetEntry(iEvent);
	  //hms delta correction
	  if(spec=="hms" && use_delta_correction)delta=delta-(p1*delta+p2*pow(delta,2)+p3*pow(delta,3)+p4*pow(delta,4)+p5*pow(delta,5));
	  //shms acceptance correction
	  double shms_acc_corr=p00+p11*delta+p22*pow(delta,2)+p33*pow(delta,3);
	  //	  shms_acc_corr=1.;

	  // Calculate E', theta, q2, w2
	  hse=hsec*(1. + delta/100.); // E'
	  if(spec=="shms"){
	  hstheta = acos(cos(thetacrad + phd)*cos(thd));
	  }
	  if(spec=="hms"){
	    // the central angle from the report file is negative
	    hstheta = acos(cos(thetacrad + phd)*cos(thd));
	    shms_acc_corr=1.;
	  }
	  sin2 = sin(hstheta/2.)*sin(hstheta/2.);
	  nu = ebeam - hse;
	  q2_calc = 4.*hse*ebeam*sin2;
	  w2_calc= mp2 + 2.*mp*nu-q2_calc;
	  
	  // fiducial and collimater cut not used
	  bool fid=fidCut(xfp, yfp, xpfp, ypfp);//shms only
	  bool coll=collCut(thd, phd, delta, yd);//shms only
	  // no dipole cut for HMS
	  if(spec=="hms")dipole=1;
	  //   only apply w2 cut for hydrogen analysis
	  bool w2_cut=true;
	  if(use_w2_cut)w2_cut = w2_calc > 1.2;

	  //Histos to calculate Pion Contamination
	  if(abs(thd)<xpCut && abs(phd)<ypCut && abs(yd) < yCut && w2_cut && dipole==1){
	    if(delta > -6 && delta <= -2){
	      hmom1->Fill(mom);
	      if(ngc>1.5)hcal_e1->Fill(etracknorm);
	      if(ngc<1.5&&el_lo>0) hcal_pi1->Fill(etracknorm);
	    }
	    if(delta > -2 && delta <= 3){
	      hmom2->Fill(mom);
	      if(ngc>1.5)hcal_e2->Fill(etracknorm);
	      if(ngc<1.5&&el_lo>0) hcal_pi2->Fill(etracknorm);
	    }
	    if(delta > 3 && delta <= 10){
	      hmom3->Fill(mom);
	      if(ngc>1.5)hcal_e3->Fill(etracknorm);
	      if(ngc<1.5&&el_lo>0) hcal_pi3->Fill(etracknorm);
	    }
	  }
	  // Main Data Cuts
	  if(ngc > ngcCut && delta > deltaMin && delta < deltaMax && etracknorm > minEdep){
	    //	    if(abs(thd)<xpCut && abs(phd)<ypCut && abs(yd) < yCut && w2_cut && dipole==1&&yd>0) //Foil studies
	    if(abs(thd)<xpCut && abs(phd)<ypCut && abs(yd) < yCut && w2_cut && dipole==1)
	      {
		if(curr>curCut){// && fid && coll){
		  //*************************************************************************************
		  //Get event by event corrections
		  //*************************************************************************************
		  // Pion Contamination
		  if(spec=="shms")piC=pionC->Eval(hse);
		  if(spec=="hms"){
		    if(thetac > -40.)piC=pionC->Eval(hse);
		    else piC=0.0;
		  }
		  hpion->Fill(piC);
		  // Cerenkov Effiency
		  cerEff=fcer->Eval(delta);
		  if(spec=="hms")cerEff=0.98;
		  hcerr->Fill(cerEff);
		  //Calorimeter efficiency
		  if(spec=="shms")calEff=getCalEff(hse);
		  hcal->Fill(calEff);
		  if(spec=="shms")calEff=1.0;
		  if(spec=="hms")calEff=0.998;
		  //*************************************************************************************
		  //  Scale Dummy Yields for window thickeness  ELOG 336
		  //*************************************************************************************
		  double dumscale=1;
		  if(target>25.)
		    {
		      if(spec=="shms")
			{
			  if(yd>0.5)dumscale=entr_fact;
			  if(yd<-0.5)dumscale=exit_fact;                 // y=0.5 0    y=-0.5->1
			  if(abs(yd)<=0.5)dumscale = entr_fact + (exit_fact-entr_fact)*(0.5 - yd);
			}
		      if(spec=="hms")
			{
			  if(yd<-0.5)dumscale=entr_fact;
			  if(yd>0.5)dumscale=exit_fact;                 // y=-0.5 0    y=0.5->1
			  if(abs(yd)<=0.5)dumscale = entr_fact + (exit_fact-entr_fact)*(yd + 0.5);
			}
		    }
		  if(scaleDummy=="no")dumscale=1;
		  hdumFact->Fill(yd,dumscale);
		  //*************************************************************************************
		  // ********   Calculate Weight factor *************************************************
		  //*************************************************************************************		  
		  wt=(1.0-piC)/calEff/cerEff/shms_acc_corr*scale*dumscale;
		  if(iEvent%100000==0){
		    cout <<"Event # "<<iEvent<<"\t"; 
		    cout << "Pion Cont "<<piC<<"\t";
		    cout << "P.gtr.p "<<mom<<"\t";
		    cout << "delta "<<delta<<"\t";
		    cout << "phd "<<phd<<"\t";
		    cout << "thd "<<thd<<"\t";
		    cout <<endl;
		  }
		  //*************************************************************************************
		  //********* Fill Histograms   *********************************************************
		  //*************************************************************************************		  
		  // For acceptance method
		  if(spec=="shms")yield4acc->Fill(1000*(hstheta-thetacrad), delta, wt);
		  if(spec=="hms")yield4acc->Fill(1000*(hstheta+thetacrad), delta, wt);
		  hmom->Fill(mom,scale*(1.0-piC));
		  // delta histograms
		  hdd->Fill(delta,wt);
		  hdd2->Fill(delta,wt);
		  hdd3->Fill(delta,wt);
		  hdd4->Fill(delta,wt*trigEff/(1.0-piC));
		  hdd5->Fill(delta,wt/calEff);
		  hdd6->Fill(delta,wt/calEff/(1.0-piC));
		  hdd7->Fill(delta,wt*trigEff/calEff/(1.0-piC));
		  hyld->Fill(delta);
		  heff->Fill(delta,wt);
		  heffcal->Fill(delta,calEff);
		  heffcer->Fill(delta,cerEff);
		  herrcer->Fill(delta,errCer);
		  heffpion->Fill(delta,1-piC);
		  // W2 histograms
		  hw2d->Fill(w2,wt);
		  hw2d_calc->Fill(w2_calc,wt);
		  hw2d_calc2->Fill(w2_calc,wt);
		  // Target variables		  
		  hxpd->Fill(thd*1000,wt);
		  hypd->Fill(phd*1000,wt);
		  hxd->Fill(xd,wt);
		  hyd->Fill(yd,wt);
		  //
		  hAvgDelta->Fill(delta,delta*wt);
		  hAvgTheta->Fill(delta,hstheta*180./TMath::Pi()*wt);
		  hq2d->Fill(q2,wt);
		  hq2d_calc->Fill(q2_calc,wt);
		  hxb->Fill(xb,wt);
		  //2d focal planes
		  xVy->Fill(yfp,xfp,wt);
		  xpVyp->Fill(ypfp,xpfp,wt);
		  xVxp->Fill(xpfp,xfp,wt);
		  ypVy->Fill(yfp,ypfp,wt);
		  yptarVytar->Fill(yd,phd,wt);
		  //*************************************************************************************
		  //********* Pt to Pt and correlated error  handling   *********************************
		  //*************************************************************************************		  
		  // errors should be fractional (%)
		  // Calculate global error
		  if(spec=="hms"&&thetac>-40)errPion=pionC->Eval(hse);
		  terr_glob=0;
		  terr_glob+=pow(errCer,2.);
		  //		    terr+=pow(errCal,2);
		  terr_glob+=pow(avgerrBoil,2.);
		  terr_glob+=pow(errTgL,2.);
		  terr_glob+=pow(errCharge,2.);
		  terr_glob+=pow(errPion,2.);
		  terr_glob=sqrt(terr_glob);
		  // Calculate point to point errors
		  terr_pt2pt=0;
		  terr_pt2pt+=pow(p2perrBoil,2.);
		  terr_pt2pt+=pow(errTrack,2.);
		  terr_pt2pt+=pow(errTrig,2.);
		  terr_pt2pt=sqrt(terr_pt2pt);
		  // Fill histos
		  herr_global->Fill(delta,terr_glob*wt); //after hadd multiple runs, divide by hdd 		  
		  herr_pt2pt->Fill(delta,terr_pt2pt*wt); //after hadd multiple runs, divide by hdd 
		  herr_boil->Fill(delta, p2perrBoil*wt);
		  herr_live->Fill(delta, errLive*wt);
		  herr_track->Fill(delta, errTrack*wt);
		  herr_trig->Fill(delta, errTrig*wt);
		  hBoilCorr->Fill(delta,boilCorr*wt);

		  // Good event Counter		 
		  goode++;
		  goode_corr+=wt;
		  tree->Fill();
		}
	      }
	  }
	}
      //*************************************************************************************
      //********************************** Summary Output   *********************************
      //*************************************************************************************		  
      
      cout << "***************************************************************************************"<<endl;
      cout << "***************************************************************************************"<<endl;
      cout << "ROOTfile: "<<f->GetName()<<endl;
      cout << "Run Number                             "<<run<<endl;
      cout << "File name                              "<<fname<<endl;
      cout << "Delta Cut                              " <<deltaMin<<" to "<<deltaMax <<endl;
      cout << "E/p Cut                                " <<minEdep <<endl;
      cout << "Cerenekov Cut                          " <<ngcCut <<endl;
      cout << "Current Cut (BCM4C)                    " <<curCut <<endl;
      cout << "Prescale Factor                        " <<psFact<<endl;
      cout<<fixed<<setprecision(4);
      cout << "Average Current (BCM4C)                " <<currentAvg<<" uA"<<endl;
      cout << "Boiling Correction                     " <<boilCorr*100.   <<"+/-" << boilCorr*errBoil*100.       <<"%"<< " (de/e="<<errBoil<<")"<<endl;
      cout << "Average boiling correction             " <<avgboilCorr*100.<<"+/-" << avgboilCorr*avgerrBoil*100. <<"%"<< " (de/e="<<avgerrBoil<<")"<<endl;
      cout << "Point to Point boiling correction      " <<boilCorr*100.   <<"+/-" << boilCorr*p2perrBoil*100.    <<"%"<< " (de/e="<<p2perrBoil<<")"<<endl;
      cout << "Livetime                               " <<livetime*100.<<"+/-" <<livetime*errLive*100. <<"%"<< " (de/e="<<errLive<<")"<<endl;
      cout << "Tracking Efficiency                    " <<trackEff*100.<<" %" << " (de/e="<<errTrack<<")"<<endl;
      cout << "Trigger Efficiency                     " <<trigEff*100.<<" %" << " (de/e="<<errTrig<<")"<<endl;
      cout << "Pion Contamintaion (mean +/- rms)      " <<hpion->GetMean()*100<<" +/- "<<hpion->GetRMS()*100.<<" %" << " (de/e="<<errPion<<")"<<endl;
      cout << "Cerenkov Efficiency (mean +/- rms)     " <<hcerr->GetMean()*100.<<" +/- "<<hcerr->GetRMS()*100.<<" %" << " (de/e="<<errCer<<")"<<endl;
      cout << "Calorimeter Efficiency (mean +/- rms)  " <<hcal->GetMean()*100.<<" +/- "<<hcal->GetRMS()*100.<<" %" << " (de/e="<<errCal<<")"<<endl;
      cout << "Good electrons (raw)                   " <<goode<<endl;
      cout << "Good electrons (corr)                  " <<goode_corr<<endl;
      cout << "The integral of hdd                    "<<hdd->Integral()<<endl;
      cout << "Charge (BCM4C)                         "<<charge <<" uC."<<endl;
      cout << "QNY raw with prescale                  " << goode/charge*psFact<<  " e-/mC" << endl;
      cout << "QNY corrected (flat corrections)       "<< goode*scale/charge << " e-/mC" << endl;
      cout << "QNY corrected (all corrections)        "<< goode_corr/charge << " e-/mC" << endl;
      cout << "***************************************************************************************"<<endl;
      cout << "***************************************************************************************"<<endl;

      if(positron){
	ofstream outpos;
	outpos.open("positronQNY.txt",ios::app | ios::out );
	//	outpos << run <<"\t";
	outpos << currentAvg <<"\t";
	outpos << goode_corr/boilCorr/charge*1000 << "\t";
	outpos << sqrt(goode_corr)/boilCorr/charge*1000 << endl;
	outpos.close();
      }

      //output files
      ofstream outFile;
      outFile.open("dataYield_pass500.txt",ios::app | ios::out );

      outFile << "***************************************************************************************"<<endl;
      outFile << "***************************************************************************************"<<endl;
      outFile << "ROOTfile: "<<f->GetName()<<endl;
      outFile << "Run Number                             "<<run<<endl;
      outFile << "File name                              "<<fname<<endl;
      outFile << "Delta Cut                              " <<deltaMin<<" to "<<deltaMax <<endl;
      outFile << "E/p Cut                                " <<minEdep <<endl;
      outFile << "Central Theta                          " <<thetac <<endl;
      outFile << "Central Momentum                       " <<hsec <<endl;
      outFile << "Cerenekov Cut                          " <<ngcCut <<endl;
      outFile << "Current Cut (BCM4C)                    " <<curCut <<endl;
      outFile << "Prescale Factor                        " <<psFact<<endl;
      outFile<<fixed<<setprecision(4);
      outFile << "Average Current (BCM4C)                " <<currentAvg<<" uA"<<endl;
      outFile << "Boiling Correction                     " <<boilCorr*100.<<"+/-" <<boilCorr*errBoil*100. <<"%"<< " (de/e="<<errBoil<<")"<<endl;
      outFile << "Livetime                               " <<livetime*100.<<"+/-" <<livetime*errLive*100. <<"%"<< " (de/e="<<errLive<<")"<<endl;
      outFile << "Tracking Efficiency                    " <<trackEff*100.<<" %" << " (de/e="<<errTrack<<")"<<endl;
      outFile << "Trigger Efficiency                     " <<trigEff*100.<<" %" << " (de/e="<<errTrig<<")"<<endl;
      outFile << "Pion Contamintaion (mean +/- rms)      " <<hpion->GetMean()*100<<" +/- "<<hpion->GetRMS()*100.<<" %" << " (de/e="<<errPion<<")"<<endl;
      outFile << "Cerenkov Efficiency (mean +/- rms)     " <<hcerr->GetMean()*100.<<" +/- "<<hcerr->GetRMS()*100.<<" %" << " (de/e="<<errCer<<")"<<endl;
      outFile << "Calorimeter Efficiency (mean +/- rms)  " <<hcal->GetMean()*100.<<" +/- "<<hcal->GetRMS()*100.<<" %" << " (de/e="<<errCal<<")"<<endl;
      outFile << "Good electrons (raw)                   " <<goode<<endl;
      outFile << "Good electrons (corr)                  " <<goode_corr<<endl;
      outFile << "The integral of hdd                    "<<hdd->Integral()<<endl;
      outFile << "Charge (BCM4C)                         "<<charge <<" uC."<<endl;
      outFile << "QNY raw with prescale                  " << goode/charge*psFact<<  " e-/mC" << endl;
      outFile << "QNY corrected (flat corrections)       "<< goode*scale/charge << " e-/mC" << endl;
      outFile << "QNY corrected (all corrections)        "<< goode_corr/charge << " e-/mC" << endl;
      outFile << "***************************************************************************************"<<endl;
      outFile << "***************************************************************************************"<<endl;
      outFile.close();

      f->Close();
      delete f;
      oFile->cd();
      tree->Write();
      tree2->Write();
      hdumFact->Write();
      yield4acc->Write();
      hdd->Write();
      hmom->Write();
      hmom1->Write();
      hmom2->Write();
      hmom3->Write();
      hdd2->Write();
      hdd3->Write();
      hdd4->Write();
      hdd5->Write();
      hdd6->Write();
      hdd7->Write();
      hAvgTheta->Write();
      hAvgDelta->Write();
      hBoilCorr->Write();
      hyld->Write();
      heff->Write();
      heffcal->Write();
      heffcer->Write();
      herrcer->Write();
      heffpion->Write();
      herr_pt2pt->Write();
      herr_boil->Write();
      herr_live->Write();
      herr_track->Write();
      herr_trig->Write();
      herrCSB->Write();
      herrKin->Write();
      herrKinRatio->Write();
      herrTot->Write();
      herr_global->Write();
      herr_globalR->Write();
      hxpd->Write();
      hypd->Write();
      hxd->Write();
      hyd->Write();
      hw2d->Write();
      hq2d->Write();
      hw2d_calc->Write();
      hw2d_calc2->Write();
      hq2d_calc->Write();
      hxb->Write();
      hcerr->Write();
      hpion->Write();
      hcal->Write();
      xVy->Write();
      xpVyp->Write();
      xVxp->Write();
      ypVy->Write();
      yptarVytar->Write();
      hcal_e1->Write();
      hcal_e2->Write();
      hcal_e3->Write();
      hcal_pi1->Write();
      hcal_pi2->Write();
      hcal_pi3->Write();
      oFile->Close();
      delete oFile;
      
    }
  
  else {cout << "Couldn't find "<<froot<<endl;}
  
  return;
}


