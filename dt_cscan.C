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
TH1D *dt[2][2];
int t_range = 200, t_range_ext = 2e7;

void read_tree(UInt_t n, TFile *file, TTree* tree){


  TLeaf * t = tree->FindLeaf("timestp");
  TLeaf * te = tree->FindLeaf("tm_stp_ext");
  TLeaf * det = tree->FindLeaf("det_id");
  TLeaf * side = tree->FindLeaf("side_id");
  TLeaf * type = tree->FindLeaf("type");
  TLeaf * ext = tree->FindLeaf("ext_flag");

  long long int time1;

  long long int nEntries = tree->GetEntries(); cout << "Thread "<< n << " entries " << nEntries<<endl;
  for(long long int i=1000; i<nEntries; i++) {
    tree->GetEntry(i); time1 = te->GetValue();
    if(ext->GetValue()==1){ //Look for external first
      for(long int j = i-1000; j < i+1000; j++) {
        tree->GetEntry(j);
        if(type->GetValue()==3){ //L3T event
          dt[int(det->GetValue())][int(side->GetValue())]->Fill(abs(t->GetValue() - time1));
        }
      }
    }
  if(n==0 && i % 10000 ==0){cout<<100.0*i/nEntries<<"%\r"<<flush;}
  }
  cout<< "Thread "<< n << " Completed"<<endl;
}

void dt_cscan() {
 gROOT->SetBatch();
 TString fname, oname;
 cout<<"Input File Number"<<endl;
 cin>>fname; oname = "../root/R"+fname+"_dt.root";
 fname = "../root/R"+fname + "_sorted.root";
 time_t now = time(0);

 for(int det = 0; det<2; det++){
   for(int side = 0; side<2; side++){
     dt[det][side] = new TH1D(Form("ext%d_%d",det, side),Form("ext%d_%d",det,side),1e6,-1e6,1e6);
     dt[det][side]->GetXaxis()->SetTitle("d-time [ns]");
   }
 }

 ROOT::EnableThreadSafety();
   TFile *ofile =  new TFile(fname, "READ");
   TTree *tree = (TTree *) ofile->Get("R3B_sort");
   read_tree(0, ofile, tree);

 TFile *outfile = new TFile(oname,"RECREATE");
 outfile->cd();
 for(int det = 0; det<2;det++){
   for(int side = 0; side<2; side++){
     dt[det][side]->Write();
   }
 }

 dt[0][0]->Draw("hist");
 //outfile->Close();

 time_t now2 = time(0);
 time_t const timediff = now2 - now;
 int minute = timediff / 60;
 int second = timediff % 60;
 cout << "\nTime Taken = " << minute << " Mins " << second << " Secs" << endl;
}
