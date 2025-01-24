#include "src/getRadCorrW2.cpp"
#include "src/getChargeError.cpp"
#include "src/getRadError.cpp"
#include "src/getGlobalError.cpp"
#include "src/getKinErrorFromMc.cpp"
#include "src/getAngle.cpp"
#include "src/getMom.cpp"
#include "src/w2cut.cpp"
#include "src/doIrebin.cpp"

// returns cross section w stat error only in a TGraph for a given spectrometer/kinematic if cs==3  (w2==8)
// returns error band in a TGraph for a given spectrometer/kinematic if cs==2 (w2==7)
// returns cross section in a TGraph for a given spectrometer/kinematic if cs==1 (w2==6)
// returns ratio xsec/model in a TGraph for a given spectrometer/kinematic if cs==0  (w2==5)

TGraphErrors* extractCS(string spec="shms", string target="r", string angle="21",string mom="2p7", int cs=1, string pass="pass400", string xaxis="xb", double w2rebin=0){

  cs=cs+5;//W2 binning
  cout << "*****************  extracting .... **************************"<<endl;
  cout <<"spec="<<spec<<endl;
  cout <<"target="<<target<<endl;
  cout <<"angle="<<angle<<endl;
  cout <<"mom= "<<mom<<endl;
  cout <<"cs= "<<cs<<endl;
  cout <<"pass="<<pass<<endl;
  cout <<"xaxis"<<xaxis<<endl;
  w2rebin=doIrebin(spec,angle,mom);
  cout <<"w2rebin: "<<w2rebin<<endl;
  //  w2rebin=0;

  ////////////  need rebin for getKinErrorFromMc ///////////////////////////
  bool rebin=true;
  if(spec=="hms" && angle=="21" && (mom=="5p1" || mom=="5p7"))rebin=false;
  if(spec=="shms" && angle=="21" && (mom=="4p0" || mom=="5p1"))rebin=false;
  if(spec=="shms" && angle=="25" && (mom=="3p5" || mom=="4p4"))rebin=false;
  if(spec=="shms" && angle=="29" && (mom=="3p0" || mom=="3p7"))rebin=false;
  if(spec=="shms" && angle=="33" && (mom=="2p6" || mom=="3p2"))rebin=false;
  if(spec=="shms" && angle=="39" && (mom=="2p0" || mom=="2p5"))rebin=false;
  rebin=false;//to make livetime error work,

  ///////////////////////?/////////////////////////////////////////////////
  //  ofstream ofile3;
  //  ofile3.open("trash.txt",ios::app | ios::out );
  //  ofile3 << spec <<"\t"<< angle <<"\t"<< mom <<"\t"<< rebin<<endl;
  //  ofile3.close();
  ///////////////////////?/////////////////////////////////////////////////
  //  int cnt=0;
  // cout << "Hello"<<cnt<<endl;cnt++;
  //  ofstream ofile;
  //   ofile.open("q2rangeMthn.txt",ios::app | ios::out );
  ///////////////////////?/////////////////////////////////////////////////

  TH1F *hkinErr=getKinErrorFromMc(target, angle, mom, spec,cs);
  //  cout << "Get the kin error hostogram"<<endl;
  if(rebin && cs<5){
    hkinErr->Rebin(3);
    hkinErr->Scale(1/3.);
  }

  double ang=21;
  double spec_flag=0.;
  if(spec=="shms")spec_flag=1.;
  if(angle=="25")ang=25;
  if(angle=="29")ang=29;
  if(angle=="33")ang=33;
  if(angle=="39")ang=39;

  //////////////////////////////
  // Charge Error
  //////////////////////////////
  ////   Just for printing charge error
  //  ofstream ofile2;
  //  ofile2.open("qError.txt",ios::app | ios::out );
  double qh_err=getChargeError(angle.c_str(), "h", mom.c_str(), spec.c_str());
  double qd_err=getChargeError(angle.c_str(), "d", mom.c_str(), spec.c_str());
  double qr_err=sqrt(qh_err*qh_err+qd_err*qd_err);
  double charge_err;
  if(target=="h")charge_err=qh_err;
  if(target=="d")charge_err=qd_err;
  if(target=="r")charge_err=qr_err;
  // add charge error to pt2pt
  //  ofile2 << spec_flag << "\t"<< ang << "\t" << qh_err << "\t"<< qd_err <<endl; 

  //////////////////////////////
  // Radiative Correction Error
  //////////////////////////////
  TGraph *g_rad;
  if(angle=="21")g_rad=getRadError(21, target);
  if(angle=="25")g_rad=getRadError(25, target);
  if(angle=="29")g_rad=getRadError(29, target);
  if(angle=="33")g_rad=getRadError(33, target);
  if(angle=="39")g_rad=getRadError(39, target);
  if(angle=="59")g_rad=getRadError(59, target);

  Float_t deltad, ratiod, errd, ep, modeld, binWidth, w2d;
  Float_t deltah, ratioh, errh, modelh,csbh, csbd, w2h;
  vector <float> cx;
  vector <float> cxe;
  vector <float> eprime;
  
  string kin;
  if(target=="r")kin="h"+angle+"deg"+mom;
  else kin=target+angle+"deg"+mom;
  cout << "The Kinematic is " << kin <<endl;

  //  TString fmc = "../mc/casey/"+spec+"_mc_"+kin+".root";
  //  TFile *fm=new TFile(fmc);

  Double_t pc=getMom(kin,spec);//shms corrected; hms is not
  kin=target+angle+"deg"+mom;
  string skin=angle+"deg"+mom;
  Double_t pOffset=0.;
  Double_t ebeam=10.602;
  Double_t hsec=pc*(1+pOffset);
  //  Double_t thetac=getAngle(angle,spec);
  // extracting cs at nominal shms angle (no beamtheta)
  Double_t thetac=getAngle(angle,"shms");// so we can plot HMS/SHMS 21 deg


  string version="v996t2";
  if(pass=="pass310")version="v0.995";
  if(pass=="pass311")version="v0.990";
  cout << "Going to use version "<<version<<" for "<<pass<<endl;
  TGraph2D *grh=getRadCorrW2("h",1,spec,version);  
  grh->SetName("grh");
  TGraph2D *grd=getRadCorrW2("d",1,spec,version);  
  grd->SetName("grd");
  cout << "The central momentum is "<<hsec<<endl;
  if(spec=="hms"&&hsec<5.5){
    double offset = -0.000276*pow(hsec,3) + 0.002585*pow(hsec,2) - 0.008697*hsec+1.0064;
    hsec=hsec*offset;
  }
    cout << "The corrected central momentum is "<<hsec<<endl;

  // GET HISTOGRAMS FROM RATIOS.CPP
  //*******************************************************************************************
  // Deuterium
  TFile *frd=new TFile(Form("ratiosOut/%s/%s_ratiosd%sdeg%s.root",pass.c_str(),spec.c_str(),angle.c_str(),mom.c_str()));
  TH1F *hrdd, *hmdd;
  hrdd=(TH1F*)frd->Get("hrd");                  // with pt2pt error
  if(cs==3)hrdd=(TH1F*)frd->Get("hrd_stat");    // stat error only  
  if(cs>=5)hrdd=(TH1F*)frd->Get("hsw2");    // W2
  if(cs==8)hrdd=(TH1F*)frd->Get("hsw2_stat");    // W2
  if(cs>=5){
    hmdd=(TH1F*)frd->Get("hmw2");    // W2
    if(w2rebin!=0){
      hmdd->Rebin(w2rebin);
      hrdd->Rebin(w2rebin);
      cout << "In the MC here are "<<hmdd->GetNbinsX()<<endl;
      cout << "In the data here are "<<hrdd->GetNbinsX()<<endl;
      hkinErr->Rebin(w2rebin);
      hkinErr->Scale(1./w2rebin);
    }
      hrdd->Divide(hmdd);    // W2
  }
  // livetime and boiling errors
  TH1F *hlte_d=(TH1F*)frd->Get("herr_live");  
  TH1F *hboil_d=(TH1F*)frd->Get("herr_boil");  
  hrdd->SetDirectory(0);
  hlte_d->SetDirectory(0);
  hboil_d->SetDirectory(0);
  // correctrun run 1640 for temp. fluctuation ELOG xxx
  if( (target=="r" || target=="d") && spec=="hms" && angle=="21" && mom == "3p3"){
    hrdd->Scale(1.0229);
    cout << endl << endl <<"Correcting for run 1640" << endl << endl;
  }
  frd->Close();
  //*******************************************************************************************  
  // Hydrogen
  TFile *frh=new TFile(Form("ratiosOut/%s/%s_ratiosh%sdeg%s.root",pass.c_str(),spec.c_str(),angle.c_str(),mom.c_str()));
  TH1F *hrdh, *hmdh;
  hrdh=(TH1F*)frh->Get("hrd");           // with pt2pt error
  if(cs==3)hrdh=(TH1F*)frh->Get("hrd_stat");  // stat error only
  if(cs>=5)hrdh=(TH1F*)frh->Get("hsw2");    // W2
  if(cs==8)hrdh=(TH1F*)frh->Get("hsw2_stat");    // W2
  if(cs>=5){
    hmdh=(TH1F*)frh->Get("hmw2");    // W2
    if(w2rebin!=0){
      hmdh->Rebin(w2rebin);
      hrdh->Rebin(w2rebin);
    }
    hrdh->Divide(hmdh);    // W2
  }
  // livetime and boiling errors
  TH1F *hlte_h=(TH1F*)frh->Get("herr_live");  
  TH1F *hboil_h=(TH1F*)frh->Get("herr_boil");  
  hrdh->SetDirectory(0);
  hlte_h->SetDirectory(0);
  hboil_h->SetDirectory(0);
  frh->Close();
  //*******************************************************************************************  
  cout << "Got ratio files, about to loop..."<<endl;
  const Int_t nbins=hrdd->GetNbinsX();
  cout << "There are "<<nbins<<" bins in hrdd"<<endl;

  bool otherSide=false;
  double wmin=0;
  double wmax=0;
  bool first=true;
  int leftPts=0;
  int rightPts=0;
  double modErr=0;

  //  Error for model dependence
  TFile *fff=new TFile("modelDepError.root");
  string gname=spec+"_"+target+angle+"deg";
  cout << "gname: "<< gname<< endl;
  TGraphErrors *gmodDep;
  if(angle!="59")gmodDep=(TGraphErrors*)fff->Get(gname.c_str())->Clone();

  //*******************************************************************************************  
  //                                 MAIN LOOP
  //*******************************************************************************************  
  //go bin by bin and wt ratios
  for (Int_t i=1; i<=nbins; i++)
    {
      // deltax: Bin Center
      // ratiox: Data/MC
      // errx:   Errors includes pt2pt, MC, stat

      deltad=hrdd->GetBinCenter(i);
      w2d=hrdd->GetBinCenter(i);
      ratiod=hrdd->GetBinContent(i);
      errd=hrdd->GetBinError(i); 

      deltah=hrdh->GetBinCenter(i);
      w2h=hrdd->GetBinCenter(i);
      ratioh=hrdh->GetBinContent(i);
      errh=hrdh->GetBinError(i); 

      binWidth=hrdh->GetBinWidth(i);//not used atm

      cout <<endl<< "start bin:" << i << "  w2: "<<deltad<<endl;      

      // Calc E' and delta if W2 binned 
      ep=(1+deltah/100)*hsec;
      Float_t mp = .9382723;
      Double_t sin2 = sin(thetac/2/180*TMath::Pi())*sin(thetac/2/180*TMath::Pi());      
      if(cs>=5){
	ep=(2*mp*ebeam+mp*mp-deltah)/(2*mp+4*ebeam*sin2);
	//	cout <<"Eprime: "<<ep<<endl;
	deltah=(ep/hsec-1.)*100;
	deltad=(ep/hsec-1.)*100;
      }
      //      cout << "b4 delta cut, bin:" << i << "  W2: "<<deltad<<endl;      

      // No idea what I was doing here.  Why not a delta check on HMS?
      if( (deltad>-10&&deltad<22&&spec=="shms") || spec=="hms") //This is how I handle rebinned data
	{
	  cout << "after delta cut, bin:" << i <<endl;

	  //Debs delta correction for SHMS  as of data pass 62 its in datayield.cpp 
	  //	  if(spec=="shms")
	  //	    {
	  //	      double p0 = 1.00156;
	  //	      double p1 = -0.002473; 
	  //	      double p2 = -1.54588e-05;
	  //	      double p3 = 6.63986e-06;
	  //	      double shms_delta_corr=p0+p1*deltad+p2*pow(deltad,2)+p3*pow(deltad,3);
	  //	      ratiod=ratiod/shms_delta_corr;
	  //	      ratioh=ratioh/shms_delta_corr;
	  //	}
	  
	  // Add the charge error and radiative dummy area (ELOG 522)
	  // error is .15% for d/h and d. .3% for h 
	  // Will get added to errd for d/h
	  // NEed to think if this is correct for d, ok for d/h
	  if(cs==1 || cs== 6){   
	    if(target=="r"){
	      errh=sqrt(errh*errh+qh_err*qh_err);
	      errd=sqrt(errd*errd+qd_err*qd_err+0.0015*0.0015);
	    }
	    else{
	      errh=sqrt(errh*errh+qh_err*qh_err+0.0030*0.0030);
	      errd=sqrt(errd*errd+qd_err*qd_err+0.0030*0.0030);	      
	    }
	    
	  }

	  Double_t nu= ebeam-ep;
	  Double_t q2 = 4.*ep*ebeam*sin2;
	  Double_t xb=q2/2./mp/nu;
	  Double_t w2= mp*mp + 2.*mp*nu-q2;
	  if(abs(w2-w2d)>.001&&cs>=5){	  
	    cout << "**** ERROR ***, w2 calc check!! ";
	    cout <<"deltah, ep, w2, w2d, :"<<deltah<<"\t";
	    cout<<ep << "\t" << w2 <<"\t"<< w2d <<endl;
	  }

	  bool goodW = false;
	  if(cs>=5&&ratioh!=0)goodW = w2cut(spec,angle,mom,w2);//deltah=w2

	  //  See how many points we are losing with w2 cut
	  if(!otherSide && !goodW && ratioh!=0)leftPts++;
	  if(goodW && ratioh!=0 )otherSide=true;
	  if(otherSide && !goodW && ratioh!=0 )rightPts++;

	  //  get model born xsec
	  modeld=grd->Interpolate(w2,thetac);  //<<"\t"<<
	  modelh=grh->Interpolate(w2,thetac);  //<<"\t"<<

	  if(ratioh!=0 && ratiod!=0 && modelh!=0 && modeld!=0)
	    {
	      float cxh=(ratioh*modelh);
	      float cxd=(ratiod*modeld);
	      double sys_y=0.58;
	      if(angle=="21"){ 
		sys_y=.56;
		//		sys_y=.93;
		//	   	   if(spec=="shms")sys_y=.525;
		//	   if(spec=="hms")sys_y=.56;
	      }
	      if((cs==2 || (cs==7&&goodW)))//error band
		{
		  if(target=="h")cx.push_back(0.);
		  if(target=="d")cx.push_back(0.);
		  if(target=="r")cx.push_back(sys_y);
		  //		  cout << "good 3:  "<<goodW <<endl;	  
		  double lte, boil_err, val;
		  int getBin=30;
		  if(rebin)getBin=10;
		  if(target=="h"){
		    lte=abs(hlte_h->GetBinContent(getBin));
		    boil_err=hboil_h->GetBinContent(getBin);
		  }
		  if(target=="d"){
		    lte=abs(hlte_d->GetBinContent(getBin));
		    boil_err=hboil_d->GetBinContent(getBin);
		  }
		  if(target=="r"){
		    lte=abs(hlte_d->GetBinContent(getBin)-hlte_h->GetBinContent(getBin));
		    boil_err=sqrt(pow(hboil_d->GetBinContent(getBin),2)+pow(hboil_h->GetBinContent(getBin),2) );
		  }		  
		  if(angle!="59")modErr=abs(gmodDep->Eval(w2)/100.);

		  
		  int bin=i;
		  //if(cs>=5)bin=hkinErr->FindBin(w2);
		  //else hkinErr->FindBin(deltah);

		  val=getGlobalError(grd, grh, ep, w2, thetac, hsec, deltah, spec, angle, target, mom, xb, g_rad, hkinErr, bin, lte, charge_err, boil_err, modErr);  
		  if(target=="h")cxe.push_back(val*cxh);
		  if(target=="d")cxe.push_back(val*cxd);
		  if(target=="r")cxe.push_back(val*cxd/cxh/2);		  

		  if(xaxis=="q2")eprime.push_back(q2);	  
		  if(xaxis=="xb")eprime.push_back(xb);
		  if(xaxis=="w2")eprime.push_back(w2);
		  if(xaxis=="ep")eprime.push_back(ep);
		}

	      if(cs==1 || cs==3 || (cs==8 && goodW) || (cs==6&&goodW))//cross section
		{

		  if(target=="h")cx.push_back(cxh);
		  if(target=="d")cx.push_back(cxd);
		  if(target=="r")cx.push_back(cxd/cxh/2);

		  if(target=="h")cxe.push_back(errh*modelh);
		  if(target=="d")cxe.push_back(errd*modeld);
		  if(target=="r")cxe.push_back(sqrt(pow(errd*modeld/cxd,2)+pow(errh*modelh/cxh,2))*cxd/cxh/2.);
		  if(xaxis=="q2")eprime.push_back(q2);	  
		  if(xaxis=="xb")eprime.push_back(xb);
		  if(xaxis=="w2")eprime.push_back(w2);
		  if(xaxis=="ep")eprime.push_back(ep);
		  //		  cout << ".............................xb, d/h: "<<xb<<"\t"<<cxd/cxh/2<<endl;
		}
	      if(cs==0 || (cs==5&& goodW))//data/model
		{
		  //		  cout << "good 3:  "<<goodW <<endl;	  
		  if(first)wmin=w2;first=false;
		  wmax=w2;

		  if(target=="h")cx.push_back(ratioh);
		  if(target=="d")cx.push_back(ratiod);
		  if(target=="r")cx.push_back(ratiod/ratioh);
		  if(target=="h")cxe.push_back(errh);
		  if(target=="d")cxe.push_back(errd);
		  if(target=="r")cxe.push_back(sqrt(pow(errd/ratiod,2)+pow(errh/ratioh,2))*ratiod/ratioh);
		  if(xaxis=="q2")eprime.push_back(q2);	  
		  if(xaxis=="xb")eprime.push_back(xb);
		  if(xaxis=="w2")eprime.push_back(w2);
		  if(xaxis=="ep")eprime.push_back(ep);
		}
	    }
	}
    }
  /*
    ofstream ofile3;
    ofile3.open("w2Cut.txt",ios::app | ios::out );
    ofile3 << kin <<"\t"<< leftPts <<"\t"<< rightPts <<"\t"<<endl;

    ofile3<<"if(kin=="<<skin<<" && w2<"<<wmax<<" && w2>"<<wmin<<")good=true;"<<endl;

    ofile3.close();
  */
  int pts=eprime.size();
  for(int i=0;i<pts;i++)
    {
      //      cout<<eprime.at(i)<<"\t";
      //      cout<<cx.at(i)<<"\t +/-";
      //      cout<<cxe.at(i)<<endl;
    }
  //    ofile.close();
  //  ofile2.close();
  TGraphErrors *gcx=new TGraphErrors(pts,&eprime[0],&cx[0],0,&cxe[0]);
  TGraphAsymmErrors* gcx_glob = new TGraphAsymmErrors(pts,&eprime[0],&cx[0],0,0,0,&cxe[0]);
  //TGraphAsymmErrors* gcx;

  //  if(cs==2||cs==7){
  if(cs==2||cs==7){return (TGraphErrors*)gcx_glob;}
  else{return gcx;}
}

