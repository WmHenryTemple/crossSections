#include "src/getRadCorrW2.cpp"

double marathonCorrection(double m_x=.375, double  m_q2=5.25, double m_w=3.1){
  double m_w2=m_w*m_w;
  string version="v996t2";
  string spec="shms";
  cout << "Going to use version "<<version;
  TGraph2D *grh=getRadCorrW2("h",1,spec,version);  
  grh->SetName("grh");
  TGraph2D *grd=getRadCorrW2("d",1,spec,version);  
  grd->SetName("grd");

  Double_t mp = .9382723;
  Double_t mp2 = mp*mp;
  double eb=10.59;
  double ep=eb-m_q2/(2*mp*m_x);
  double sin2=m_q2/(4*eb*ep);
  double m_th=2*asin(sqrt(sin2));

  m_th=m_th*180/TMath::Pi();
  //Need what ep corresponds to marathons x for SHMS=21.035deg
  double f_eb=10.602;
  double f_sin=sin((21.035*TMath::Pi()/180)/2);
  double f_sin2=f_sin*f_sin;
  double f_ep=(f_eb*m_x*mp)/(2*f_eb*f_sin2+m_x*mp);
  double f_q2=4*f_eb*f_ep*f_sin2;
  double f_w2= mp2 + 2.*mp*(f_eb-f_ep)-f_q2;
  double f_th=21.035;

  double m_d=grd->Interpolate(m_w2,m_th);
  double m_h=grh->Interpolate(m_w2,m_th);
  double m_ratio=m_d/m_h;

  double f_d=grd->Interpolate(f_w2,f_th);
  double f_h=grh->Interpolate(f_w2,f_th);
  double f_ratio=f_d/f_h;
  double correction=f_ratio/m_ratio;

  cout << "Marathon theta/w2:"<<m_th<<"   "<< m_w2 <<endl;
  cout << "For F2 that corresponds to E'="<<f_ep;
  cout << "and W2="<<f_w2<<endl;
  cout<< "The correction is"<<correction<<endl;
  return correction;
}
