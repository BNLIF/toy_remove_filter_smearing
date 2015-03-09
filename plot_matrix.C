void plot_matrix(){
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  TMatrixD *m1 = new TMatrixD(4,2);
  *m1(1,1) =1;
  *m1(0,1) =2;
  *m1(1,0) =3;
  *m1(0,0) =1;

  *m1(2,1) =1;
  *m1(3,1) =2;
  *m1(3,0) =3;
  *m1(2,0) =1;
  m1->Draw("COLZ");
  
  c1->cd(2);
  TMatrixD *m2 = new TMatrixD(2,4);
  m2->Transpose(*m1);
  m2->Draw("COLZ");

  c1->cd(3);
  TMatrixD *m3 =  new TMatrixD(2,2);
  m3->Mult(*m2,*m1);
  m3->Draw("COLZ");

  m3->Print();

  c1->cd(4);
  TVectorD va(2);
  va(0)=1;
  va(1) = 2;

  TVectorD dd = (*m3) * va;
  cout << dd(0)<< " " << dd(1) << endl;
}
