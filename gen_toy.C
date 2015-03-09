void gen_toy(){

  //read in the field response
  TFile *file1 = new TFile("convolute.root");
  TH1F *hfield = (TH1F*)file1->Get("hv23");
  Double_t x[5000],y[5000];
  for (Int_t i=0;i!=5000;i++){
    x[i] = hfield->GetBinCenter(i+1);
    y[i] = hfield->GetBinContent(i+1)*2000.;
  }
  TGraph *gfield = new TGraph(5000,x,y);
  // gfield->Draw("AL");

  //generate a toy MC 
  const Int_t nbin = 200;
  TH1F *htrue = new TH1F("htrue","htrue",nbin,0,nbin/2.);

  Double_t sigma = 0.5; //5 us
  for (Int_t i=0;i!=3;i++){
    Double_t time = gRandom->Uniform(30.,70.);
    for (Int_t j=0;j!=nbin;j++){
      Double_t time1 = htrue->GetBinCenter(j+1);
      Double_t content = htrue->GetBinContent(j+1);
      content += 1./sigma*exp(-0.5/sigma/sigma*pow(time-time1,2));
      htrue->SetBinContent(j+1,content);
    }
  }

  TH1F *hfield = new TH1F("hfield","hfield",nbin,0,nbin/2.);
  for (Int_t i=0;i!=nbin;i++){
    Double_t binc = hfield->GetBinCenter(i+1);
    if (binc < nbin/4.){
      hfield->SetBinContent(i+1,gfield->Eval(binc));
    }else{
      hfield->SetBinContent(i+1,gfield->Eval(binc-nbin/2.));
    }
    
  }

  // hfield->Draw();
  // htrue->Draw();
  

  //define a filter
  TF1 *f1  = new TF1("f1","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  Double_t par[10];
  par[0] = 1.73/0.961508;
  par[1] = 1.69;
  par[2] = 1.55;
  par[3] = 0.19;
  par[4] = 3.75;
  f1->SetParameters(par);
  TH1F *hfilter = new TH1F("hfilter","hfilter",nbin,0,nbin);
  
  for (Int_t i=0;i!=nbin;i++){
    Double_t frequency = (nbin/2.-fabs(i-nbin/2.))*2./nbin;
    hfilter->SetBinContent(i+1,f1->Eval(frequency));
  }
  
  //convolute with a filter
  // do a FFT on the h1
  TH1 *h1m = 0;
  TH1 *h1p = 0;
  h1m = htrue->FFT(h1m,"MAG");
  h1p = htrue->FFT(h1p,"PH");

  TH1 *hfm = 0;
  TH1 *hfp = 0;
  
  hfm = hfield->FFT(hfm,"MAG");
  hfp = hfield->FFT(hfp,"PH");
  
  Double_t value_re[nbin],value_im[nbin];
  for (Int_t i=0;i!=nbin;i++){
    Double_t rho = h1m->GetBinContent(i+1);
    Double_t phi = h1p->GetBinContent(i+1);
    
    Double_t rho1 = hfm->GetBinContent(i+1);
    Double_t phi1 = hfp->GetBinContent(i+1);

    rho = rho * rho1;
    phi = phi + phi1;

    //  Double_t filter = hfilter->GetBinContent(i+1);
    value_re[i] = rho*cos(phi)/nbin;
    value_im[i] = rho*sin(phi)/nbin;

  }
  Int_t n = nbin;
  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();
  TH1 *fb;
  fb = TH1::TransformHisto(ifft,fb,"Re");

  //  fb->Draw();
  
  Double_t noise_level = 0.05;

  TH1F *hconv = new TH1F("hconv","hconv",nbin,0,nbin/2.);
  
  Double_t baseline = 0.;//-0.150729;//(fb->GetBinContent(1) + fb->GetBinContent(nbin))/2.;
  for (Int_t i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -= baseline;
    Double_t noise = gRandom->Gaus(0,noise_level);
    content += noise;
    hconv->SetBinContent(i+1,content);
  }
  // hconv->Draw();


  TH1F *hdeconv = new TH1F("hdeconv","hdeconv",nbin,0,nbin/2.);
  
  
  TH1 *hcm = 0;
  TH1 *hcp = 0;
  hcm = hconv->FFT(hcm,"MAG");
  hcp = hconv->FFT(hcp,"PH");

  for (Int_t i=0;i!=nbin;i++){
    Double_t rho = hcm->GetBinContent(i+1);
    Double_t phi = hcp->GetBinContent(i+1);
    
    Double_t rho1 = hfm->GetBinContent(i+1);
    Double_t phi1 = hfp->GetBinContent(i+1);

    if (rho1!=0){
      rho = rho / rho1;
    }else{
      rho = 0.;
    }
    phi = phi - phi1;

    Double_t filter = hfilter->GetBinContent(i+1);
    
    // filter = 1.0;
    // if (i==0) filter = 0.0;
    value_re[i] = rho*cos(phi)*filter/nbin;
    value_im[i] = rho*sin(phi)*filter/nbin;
  }
  
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();
  fb=0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  //baseline = -0.075;
  baseline = -0.075;
  for (Int_t i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -= baseline;
    
    hdeconv->SetBinContent(i+1,content);
  }
  

  //get filter in the time domain
  TH1F *hfilter_time = new TH1F("hfilter_time","hfilter_time",nbin*2,-nbin/4.,nbin/4.);
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  //fb->Draw();

  baseline = -0.0025*2;//fb->GetBinContent(Int_t(nbin/2.));
  //cout << baseline << endl;
  for (Int_t i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    content -=baseline;
    Int_t newbin = i + nbin/2;
    if (newbin > nbin) newbin -=nbin;
    cout << newbin << endl;
    hfilter_time->SetBinContent(2*newbin,content/2.);
    hfilter_time->SetBinContent(2*newbin+1,content/2.);
  }
  //hfilter_time->SetBinContent(1,fb->GetBinContent(nbin/2)/2.);
  hfilter_time->Draw();
  //filter in double binning ...
  //cout << hfilter_time->GetSum() << endl;
  


  TFile *file =  new TFile("toy.root","RECREATE");
  htrue->SetDirectory(file);
  hconv->SetDirectory(file);
  hdeconv->SetDirectory(file);
  hfilter->SetDirectory(file);
  hfilter_time->SetDirectory(file);
  file->Write();
  file->Close();
  
}
