void invert(){
  TFile *file = new TFile("toy.root");
  TH1F *htrue = (TH1F*)file->Get("htrue");
  TH1F *hconv = (TH1F*)file->Get("hdeconv");
  TH1F *hfilter_time = (TH1F*)file->Get("hfilter_time");
  
  
  Int_t nbin = htrue->GetNbinsX();
  Double_t x[10000],y[10000];
  for (Int_t i=0;i!=nbin*2;i++){
    x[i] = hfilter_time->GetBinCenter(i+1);
    y[i] = hfilter_time->GetBinContent(i+1);
  }
  TGraph *gfilter = new TGraph(nbin*2,x,y);
  gfilter->Draw("AL");

  Double_t bin_width = hfilter_time->GetBinWidth(1);
  Double_t time_bin = 2.0/2.; // binning for 2n, binning for n is twice long

  Int_t binnum = time_bin/0.5;

  Int_t n = hconv->GetNbinsX() * hconv->GetBinWidth(1)/time_bin/2.;
  cout << n << " " << 2*n << endl;
  cout << binnum << endl;
 

  //construct A, and AT, and Vy
  TMatrixD *A = new TMatrixD(2*n,n);
  TMatrixD *AT = new TMatrixD(n,2*n);
  TMatrixD *Vy_invert = new TMatrixD(2*n,2*n);
  TMatrixD *MB = new TMatrixD(n,n);
  TMatrixD *FM = new TMatrixD(n,2*n);

  TVectorD rx(n);
  TVectorD ry(2*n);

  for (Int_t i=0;i!=2*n;i++){
    ry[i] = 0.;
    for (Int_t j=0;j!=binnum;j++){
      ry[i] += hconv->GetBinContent(i*binnum+j);
    }
    //cout << ry[i] << endl;
    (*Vy_invert)(i,i) = pow(1./(0.05*sqrt(binnum)),2);
  }
  
  //Vy_invert->Draw("COLZ");



  Double_t sum = 0;
  //rebin filter, assum i=5 is the middle 0
  for (Int_t i=0;i!=20;i++){
    x[i] = (i-10+0.5)*time_bin;
    y[i] = 0.;
    for (Int_t j=0;j!=Int_t(time_bin/bin_width/2.);j++){
      y[i] += gfilter->Eval(x[i]+(2*j+1)*bin_width/2.);
      y[i] += gfilter->Eval(x[i]-(2*j+1)*bin_width/2.);
    }
    //y[i] /= (time_bin/bin_width);
    sum += y[i];
  }
  cout << sum << endl;

  TGraph *g1 = new TGraph(20,x,y);
  g1->Draw("*same");

  for (Int_t i=0;i!=n;i++){
    for (Int_t j=0;j!=20;j++){
      Int_t row = i;
      Int_t column = 2*i-10+j+1;
      if (column >=0 && column < 2*n){
	(*A)(column,row)= y[j];
      }
    }
  }

  AT->Transpose(*A);
  
  //A->Draw("COLZ");
  
  
  *MB = (*AT) * (*Vy_invert) * (*A);
  MB->Invert();
  MB->Draw("COLZ");

  *FM = (*MB) * (*AT) * (*Vy_invert);
  // FM->Draw("COLZ");
  rx = (*FM) * (ry);
  
  TH1F *hmatrix = new TH1F("hmatrix","hmatrix",n,0,n);
  TH1F *htrebin = new TH1F("htrebin","htrebin",n,0,n);
  TH1F *hcrebin = new TH1F("hcrebin","hcrebin",n,0,n);
  for (Int_t i=0;i!=n;i++){
    hmatrix->SetBinContent(i+1,rx[i]);
    Double_t temp = 0;
    for (Int_t j=0;j!=binnum*2;j++){
      temp += htrue->GetBinContent(i*binnum*2+j);
    }
    htrebin->SetBinContent(i+1,temp);
    
    temp = 0;
    for (Int_t j=0;j!=binnum*2;j++){
      temp += hconv->GetBinContent(i*binnum*2+j);
    }
    hcrebin->SetBinContent(i+1,temp);


    hmatrix->SetBinError(i+1,sqrt((*MB)(i,i)));
  }

  hmatrix->Draw();
  htrebin->SetLineColor(2);
  htrebin->Draw("same");

  hcrebin->SetLineColor(4);
  hcrebin->Draw("same");

  hmatrix->SetLineWidth(2.5);
  hcrebin->SetLineWidth(2.5);
  htrebin->SetLineWidth(2.5);

  htrebin->SetLineStyle(2);
  hcrebin->SetLineStyle(3);

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(hmatrix,"Results","l");
  le1->AddEntry(htrebin,"Rebin of True","l");
  le1->AddEntry(hcrebin,"Rebin of Conv.","l");
  le1->Draw();
  

}
