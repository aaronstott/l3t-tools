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

void mult_scan(){
  time_t now = time(0);
  TString infile1, ofile;
  cout<<"input file number"<<endl;
  cin>>infile1; ofile=infile1;
  infile1="../root/R"+infile1+"_sorted.root";

  TFile *file = new TFile(infile1,"READ");
  TTree *tree = (TTree  *) file->Get("R3B_sort");

  TLeaf * t = tree->FindLeaf("timestp");
  TLeaf * type = tree->FindLeaf("type");
  TLeaf * det = tree->FindLeaf("det_id");

  int window = 1e6;
  int min_mult = 1, max_mult = 500;

  TH1D *mult[3];
  mult[0] = new TH1D("mult_ext", "mult_ext", max_mult, 0, max_mult);
  mult[1] = new TH1D("mult_l3t0", "mult_l3t0", max_mult, 0, max_mult);
  mult[2] = new TH1D("mult_l3t1", "mult_l3t1", max_mult, 0, max_mult);

  TH1D *dt[3];
  dt[0] = new TH1D("dt_l3ta_ext","dt_l3ta_ext", 1e4, 0, 1e6);
  dt[1] = new TH1D("dt_l3t0_ext","dt_l3t0_ext", 1e4, 0, 1e6);
  dt[2] = new TH1D("dt_l3t1_ext","dt_l3t1_ext", 1e4, 0, 1e6);

  TH1D *dt_all, *dt_all_ext;
  dt_all = new TH1D("dt_all_test","dt_all_test", 2e4, -1e6, 1e6);
  dt_all_ext = new TH1D("dt_all_ext","dt_all_ext", 2e4, 0, 1e6);


  TH2D *mult_dt;
  mult_dt = new TH2D("mult_dt", "mult_dt", 2e4, 0, 2e6, 10, 0, 10);

  int cmult[3];
  int nEntries = tree->GetEntries();
  int evt, prev;
  long long tl3t[3], text;
  long long te0, t0;

  for(int i=0; i<nEntries; i++){
    cmult[0]=0; cmult[1]=0; cmult[2]=0; //reset 0-ext, 1-det0, 2-det1
    evt=-1;
    tl3t[0]=0; tl3t[1]=0; tl3t[2]=0;

    tree->GetEntry(i);
    cmult[int(type->GetValue()-2)]++;
    evt = int(type->GetValue()-2 + det->GetValue());
    t0 = t->GetValue();

    if(type->GetValue()==2){evt=0; te0  = t->GetValue();}

    for(int j=i+1; j<i+max_mult; j++){
        tree->GetEntry(j);
        if(int(type->GetValue()-2)!= evt) {i=j; break;}
        if(int(type->GetValue()==2)){cmult[0]++;}
        else{cmult[int(type->GetValue()-2 + det->GetValue())]++;}
        if(j==i+max_mult){i=j;}
      }
      if(cmult[2]>=min_mult){mult[2]->Fill(cmult[2]); tl3t[0]=t->GetValue(); tl3t[2]=t->GetValue();}
      if(cmult[1]>=min_mult){mult[1]->Fill(cmult[1]); tl3t[0]=t->GetValue(); tl3t[1]=t->GetValue();}
      if(cmult[0]>=min_mult){mult[0]->Fill(cmult[0]); text=t->GetValue();}

      for(int k=0;k<=2;k++){
        if(cmult[k]>=min_mult && tl3t[k]-text <= 1e8 && evt!=prev){dt[k]->Fill(tl3t[k]-te0); dt[0]->Fill(tl3t[0]-te0);}
      }
	    dt_all->Fill(t->GetValue() - t0);
      dt_all_ext->Fill(t->GetValue() - te0);
      mult_dt->Fill(t->GetValue() - te0,cmult[1]+cmult[2]+cmult[0]);
      prev = evt;
      if(i % 10000 ==0){cout<<100.0*i/nEntries<<"%\r"<<flush;}
    }

  ofile = "../root/R"+ofile+"_extmult.root";
  TFile *outfile = new TFile(ofile,"RECREATE");

  mult[0]->Write();
  mult[1]->Write();
  mult[2]->Write();

  dt[0]->Write();
  dt[1]->Write();
  dt[2]->Write();

  dt_all->Write();
  dt_all_ext->Write();
  mult_dt->Write();
  outfile->Close();



}
