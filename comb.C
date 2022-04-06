#include <vector>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TFile.h>

void comb() {
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!============== L3T tree-merger  ==============!"<<endl;
  cout<<"!============== A. Stott - UoYork =============!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;

  int nfiles;
  cout<<"Input number of files:";
  cin>>nfiles;

  TTree *tree[nfiles];
  TFile *file[nfiles];

  TList *list = new TList;

  for(int i=0; i<nfiles; i++ ){
    TString infile1;
    cout<<"input file number"<<endl;
    cin>>infile1;
    infile1="../root/R"+infile1+"_sorted.root";

    file[i] = new TFile(infile1,"READ");
    tree[i] = (TTree *) file[i]->Get("R3B_sort");

    list->Add(tree[i]);
  }

  TString ofile;
  cout<<"Output File Name:";
  cin>>ofile;
  ofile="../root/R"+ofile+"_sorted.root";
  TFile *output = new TFile(ofile,"RECREATE");

  TTree *newtree = TTree::MergeTrees(list);
  newtree->SetName("R3B_sort");
  newtree->Write();
  output->Close();

}
