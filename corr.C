#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TCanvas.h>

#define PI 3.14159265

bool gate1D(double value, double min, double max) {
  if ((value < max) && (value > min)) {
    return true;
  }
  return false;
}

void corr(){
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!============== L3T Correlations  =============!"<<endl;
  cout<<"!============== A. Stott - UoYork =============!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;


  TString fname;
  cout<<"Input File Number"<<endl;
  cin>>fname; fname = "../root/R"+fname + "_events.root";
  time_t now = time(0);

  Long64_t *evt=0;
  vector<int> *det=0;
  vector<int> *side=0;

  std::vector<int> *mult0=0;
  std::vector<int> *mult1=0;
  //int mult0[3], mult1[3];

  vector<int> *adc=0;
  vector<Long64_t> *t=0;
  vector<int> *ch=0;
  vector<int> *asic=0;

  int ext_flag;
  Long64_t ext_time;

  TFile *tree = new TFile(fname,"READ");
  TTree *events = (TTree  *) tree->Get("R3B_evt");

  events->SetBranchAddress("evt",&evt);
  events->SetBranchAddress("det",&det);
  events->SetBranchAddress("side",&side);
  events->SetBranchAddress("mult0",&mult0);
  events->SetBranchAddress("mult1",&mult1);
  events->SetBranchAddress("adc",&adc);
  events->SetBranchAddress("timestp_",&t);
  events->SetBranchAddress("ch",&ch);
  events->SetBranchAddress("asic",&asic);
  events->SetBranchAddress("ext_flag",&ext_flag);
  events->SetBranchAddress("tmp_ext",&ext_time);

  double dz = 2330000;
  double x_pos, y_pos;
  int strip0, strip1;
  int strip[2][2];

  int maxstrip = 1536;
  int maxadc = 2000;

  // Prepare hists
  TList *list = new TList;
  TH2D *ss[2], *pp[2], *adc_m[2];
  TH2D *pp_ext[2];
  TH1D *ext_t;
  TH1D *mult_h[2][3], *mult_c[2];
  TH1D *mult_t, *mult_d;
  TH1D *dt_ext;
  TH1D *theta_xz, *theta_yz;
  TH2D *xz_yz;
  TH2D *ext_xz_yz;

  ext_t = new TH1D("ext_t", "external trigger events", 3, -1, 2); list->Add(ext_t); // external trigg
  dt_ext = new TH1D("ext_dt", "ext_dt", 200, -500, 500); list->Add(dt_ext);

  mult_t = new TH1D("mult_t", "Total Strip multiplicity of Event", 20, 0, 20); list->Add(mult_t);
  mult_d = new TH1D("mult_d", "Number of Detectors in Event", 2, 0, 2); list->Add(mult_d);

  theta_xz = new TH1D("theta_xz", "theta_xz", 100, -5, 5); list->Add(theta_xz);
  theta_yz = new TH1D("theta_yz", "theta_yz", 200, -10, 10); list->Add(theta_yz);

  xz_yz = new TH2D("xz_yz", "Angular reconstruction", 100, -5, 5, 200, -10, 10); list->Add(xz_yz);
  ext_xz_yz = new TH2D("ext_xz_yz", "Angular reconstruction (external trigger)", 100, -5, 5, 200, -10, 10); list->Add(ext_xz_yz);

  for (int i=0; i<=1; i++){
    // mult
    mult_c[i] = new TH1D(Form("mult_c_%i", i),Form("side coincidence on det %i", i),2, 0, 2 ); list->Add(mult_c[i]);
    for(int j=0; j<=2; j++){
      mult_h[i][j] = new TH1D(Form("mult_%i_%i", i, j), Form("hit multiplicity side%i det%i", i, j), 10, 0, 10); list->Add(mult_h[i][j]);
    }
    // strip - strip
    ss[i] = new TH2D(Form("ss_%i", i),Form("front-rear strip correlation det%i",i), maxstrip, 0, maxstrip, maxstrip, 0, maxstrip); list->Add(ss[i]);
    // x-y
    pp[i] = new TH2D(Form("pp_%i", i),Form("position reconstruction det%i",i), 270, -20e3, 250e3, 2000, 0, 76000); list->Add(pp[i]);
    // adc_multiplicity
    pp_ext[i] = new TH2D(Form("pp_ext_%i", i), Form("pp_ext_%i", i), 270, -20e3, 250e3, 2000, 0, 76000); list->Add(pp_ext[i]);
    //adc_m[i] = new TH2D(Form("adc_m_%i", i), Form("adc_m_%i", i), 1000, 0, 1000, 50, 0, 50); list->Add(adc_m[i]);

  }

  for(int i = 0; i < events->GetEntries()-1; i++) {
      events->GetEntry(i);
      if(mult0->at(2)+mult1->at(2)>0){
///// Fill multiplicities
      mult_t->Fill(mult0->at(2) + mult1->at(2)); //total strip mult

      if(mult0->at(2)>0 && mult1->at(2)>0){mult_d->Fill(1);} /// Check detector coincidences at least 1 hit on each
      else {mult_d->Fill(0);}

      if(mult0->at(0)>0 && mult0->at(1)>0){mult_c[0]->Fill(1);}  // front-rear det0
      else if (mult0->at(0)>0 || mult0->at(1)>0){mult_c[0]->Fill(0);}

      if(mult1->at(0)>0 && mult1->at(1)>0){mult_c[1]->Fill(1);}  // front-rear det1
      else if (mult1->at(0)>0 || mult1->at(1)>0){mult_c[1]->Fill(0);}

      for (int t=0; t<=2; t++){
          mult_h[0][t]->Fill(mult0->at(t));
          mult_h[1][t]->Fill(mult1->at(t));
        }


      ///Start scanning for correlations
      for(int j=0; j < det->size(); j++){
        strip[det->at(j)][side->at(j)] = ch->at(j)+128*asic->at(j); //load strip
        for(int k=j+1; k< det->size(); k++){


          // Same detector Correlation
          if(((side->at(j)==0 && side->at(k)==1) || ((side->at(k)==0 && side->at(j)==1))) && det->at(j)==det->at(k)){
            strip0 = ch->at(j)+128*asic->at(j);
            strip1 = ch->at(k)+128*asic->at(k);

            ss[det->at(j)]->Fill(strip0, strip1);

            y_pos = (-178.55*real(strip0 - (1536-strip1)));
            x_pos = (25*(strip0 + (1536-strip1)));

            pp[det->at(j)]->Fill(y_pos, x_pos);
            if(ext_flag==1){pp_ext[det->at(j)]->Fill(y_pos, x_pos);}
          }
        }
          // Angular reconstruction
          //if(mult0->at(0)==1 && mult0->at(1)==1 && mult1->at(0)==1 && mult1->at(1)==1){
            if(mult0->at(0)>=1 && mult0->at(1)>=1 && mult1->at(0)>=1 && mult1->at(1)>=1){

            double x0 = (25*(strip[0][0] + (1536-strip[0][1])));
            double x1 = (25*(strip[1][1] + (1536-strip[1][0])));

            double y0 = (-178.55*real(strip[0][0] - (1536-strip[0][1])));
            double y1 = (-178.55*real(strip[1][1] - (1536-strip[1][0])));

            theta_xz->Fill(atan((x1-x0)/dz)* 180 / PI);
            theta_yz->Fill(atan((y1-y0)/dz)* 180 / PI);

            xz_yz->Fill(atan((x1-x0)/dz)* 180 / PI, atan((y1-y0)/dz)* 180 / PI);
            if(ext_flag==1){ext_xz_yz->Fill(atan((x1-x0)/dz)* 180 / PI, atan((y1-y0)/dz)* 180 / PI);}
          }
        }
      

      if (i % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << i << " of " << events->GetEntries() << ", " << 100 * i / events->GetEntries() << "% complete" << "\r" << flush;
    }
    if(ext_flag==1 && adc->size()>0){ext_t->Fill(1); dt_ext->Fill(ext_time-t->at(0));}
    if(ext_flag==1 && adc->size()==0){ext_t->Fill(-1);}
    if(ext_flag==0 && adc->size()>0){ext_t->Fill(0);}
}

  TFile *outfile = new TFile("output.root","recreate");

  outfile->cd();
  list->Write();


  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000); c1->Draw();
  c1->Divide(4, 2);
  c1->cd(3); ss[0]->Draw("colz");
  c1->cd(4); ss[1]->Draw("colz");
  c1->cd(7); pp[0]->Draw("colz");
  c1->cd(8); pp[1]->Draw("colz");

  c1->cd(1); theta_xz->Draw();
  c1->cd(2); theta_yz->Draw();
  c1->cd(5); xz_yz->Draw("colz");
  c1->cd(6); ext_t->Draw("hist text");
  //outfile->Close();

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000); c2->Draw();
  c2->Divide(4, 2);

  c2->cd(1); mult_t->Draw("hist text");
  c2->cd(2); mult_d->Draw("hist text");
  c2->cd(5); mult_h[0][0]->Draw("hist text");
  c2->cd(6); mult_h[1][0]->Draw("hist text");
  c2->cd(3); mult_c[0]->Draw("hist text");
  c2->cd(4); mult_c[1]->Draw("hist text");
  c2->cd(7); mult_h[0][1]->Draw("hist text");
  c2->cd(8); mult_h[1][1]->Draw("hist text");

  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000); c1->Draw();
  c3->Divide(1, 2);
  c3->cd(1); pp[0]->Draw("colz");
  c3->cd(2); pp[1]->Draw("colz");

  time_t now2 = time(0);
  time_t const timediff = now2 - now;
  int minute = timediff / 60;
  int second = timediff % 60;
  cout << "\nTime Taken = " << minute << " Mins " << second << " Secs" << endl;

  }
