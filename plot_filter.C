void plot_filter(){
  TF1 *f1  = new TF1("f1","(x>0.0)*gaus*exp(-0.5*pow(x/[3],[4]))");
  Double_t par[10];
  par[0] = 1.73;
  par[1] = 1.69;
  par[2] = 1.55;
  par[3] = 0.19;
  par[4] = 3.75;
  f1->SetParameters(par);
  f1->Draw();

}
