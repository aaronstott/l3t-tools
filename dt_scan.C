#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TFile.h>
#include <TCanvas.h>


const int nWorkers = 1U;
TH1D *dt[3][5];
TH1D *dt_det;
TH1D *dt_ext[2];
int t_range = 200, t_range_ext = 2e7;

void read_tree(UInt_t n, TFile *file, TTree* tree){


  TLeaf * t = tree->FindLeaf("timestp");
  TLeaf * te = tree->FindLeaf("tm_stp_ext");
  TLeaf * ch = tree->FindLeaf("ch_id");
  TLeaf * det = tree->FindLeaf("det_id");
  TLeaf * side = tree->FindLeaf("side_id");
  TLeaf * type = tree->FindLeaf("type");
  TLeaf * adc  = tree->FindLeaf("adc_data");
  TLeaf * asic = tree->FindLeaf("asic_id");
  TLeaf * ext = tree->FindLeaf("ext_flag");

  long long int time1, time2, time_diff;
  int det1, det2, asic1, asic2, side1, side2, ch1, ch2, adc1, adc2, type1;
  int hits, ext1, ext2;

  long long int nEntries = tree->GetEntries(); cout << "Thread "<< n << " entries " << nEntries<<endl;
  for(long long int i=0; i<nEntries/10; i++) {
   tree->GetEntry(i);time1=t->GetValue();

     for(long int j = i-100; j < i+100; j++) {
       if(j!=i){
       tree->GetEntry(j); time2=t->GetValue(); type1=type->GetValue();
       if(abs(time2-time1)<t_range_ext){
       det2 = det->GetValue(); asic2=asic->GetValue(); side2=side->GetValue(); ch2=ch->GetValue(); adc2=adc->GetValue();ext2=ext->GetValue();
       if(ext1==1 && ext2!=1){dt_ext[det2]->Fill(time2-time1);}
       if(ext1!=1 && ext2!=1){
       if(det2==det1 && type1==3){dt[side1+side2][det1]->Fill((time2-time1));}
       if(det1==0 && det2==1){dt_det->Fill(time2-time1);}
        }
      }

      }
    }

 if(n==0 && i % 10000 ==0){cout<<100.0*i/nEntries<<"%\r"<<flush;}
  }
  cout<< "Thread "<< n << " Completed"<<endl;
}


void dt_scan() {
 gROOT->SetBatch();
 TString fname, oname;
 cout<<"Input File Number"<<endl;
 cin>>fname; oname = "../root/R"+fname+"_dt.root";
 fname = "../root/R"+fname + "_sorted.root";
 time_t now = time(0);

 for(int det = 0; det<2;det++){
   dt[0][det] = new TH1D(Form("rrt%d",det),Form("rrt%d",det),400,-t_range,t_range);
   dt[0][det]->GetXaxis()->SetTitle("d-time [ns]");
   dt[1][det] = new TH1D(Form("frt%d",det),Form("frt%d",det),400,-t_range,t_range);
   dt[1][det]->GetXaxis()->SetTitle("d-time [ns]");
   dt[2][det] = new TH1D(Form("fft%d",det),Form("fft%d",det),400,-t_range,t_range);
   dt[2][det]->GetXaxis()->SetTitle("d-time [ns]");
   dt_ext[det] = new TH1D(Form("dt_ext_%d", det), Form("dt - ext - det%d", det), 4e3, -t_range_ext, t_range_ext);
   dt_ext[det]->GetXaxis()->SetTitle("d-time [ns]");
   }
   dt_det = new TH1D("d0_d1","dt(det==1) - dt(det==0)", 4e7, -2e9, 2e9);
   dt_det->GetXaxis()->SetTitle("d-time [ns]");

 ROOT::EnableThreadSafety();

   TFile *ofile =  new TFile(fname, "READ");
   TTree *tree = (TTree *) ofile->Get("R3B_sort");
   read_tree(0, ofile, tree);

 TFile *outfile = new TFile(oname,"RECREATE");
 outfile->cd();
 for(int det = 0; det<2;det++){dt[0][det]->Write(); dt[1][det]->Write(); dt[2][det]->Write(); dt_ext[det]->Write();}
 dt_det->Write();
 outfile->Close();

 time_t now2 = time(0);
 time_t const timediff = now2 - now;
 int minute = timediff / 60;
 int second = timediff % 60;
 cout << "\nTime Taken = " << minute << " Mins " << second << " Secs" << endl;
}
