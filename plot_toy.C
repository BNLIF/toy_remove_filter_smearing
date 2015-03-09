void plot_toy(){
  TFile *file = new TFile("toy.root");
  TH1F *htrue = (TH1F*)file->Get("htrue");
  TH1F *hconv = (TH1F*)file->Get("hconv");
  TH1F *hdeconv = (TH1F*)file->Get("hdeconv");
  TH1F *hfilter = (TH1F*)file->Get("hfilter");
  TH1F *hfilter_time = (TH1F*)file->Get("hfilter_time");


  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2,2);
  c1->cd(1);
  htrue->Draw();
  htrue->SetXTitle("Time (us)");
  htrue->SetTitle("True Signal with 0.05 random noise");
  
  c1->cd(2);
  hconv->Draw();
  hconv->SetTitle("Convolution with field response");
  hconv->SetXTitle("Time (us)");

  c1->cd(3);
  // hfilter->Draw();
  // hfilter->SetTitle("Filter in Frequency Domain");
  hdeconv->Draw();
  hdeconv->SetTitle("Deonvolution with Filter");
  hdeconv->SetXTitle("Time (us)");


  c1->cd(4);
  hfilter_time->Draw();
  hfilter_time->SetTitle("Filter in Time Domain");
  hfilter_time->SetXTitle("Time (us)");
 
  
}
