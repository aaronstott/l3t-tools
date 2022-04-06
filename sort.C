#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void sort()
{
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!============== L3T time Sorting ==============!"<<endl;
  cout<<"!============== A. Stott - UoYork =============!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;


  time_t now = time(0);
  TString fname, ofile;
  cout<<"Input File Number"<<endl;
  cin>>fname; ofile = fname;
  fname="../root/R"+fname+"_0_sort.root";
  ofile = "../root/R"+ofile+"_sorted.root";

// Collect files and old tree for singles ///
  TFile *file = new TFile(fname,"READ");
  TTree *tree = (TTree  *) file->Get("R3B_sort");

  TLeaf * t_0  = tree->FindLeaf("tm_stp");
  TLeaf * te_0 = tree->FindLeaf("tm_stp_ext");
  TLeaf * type = tree->FindLeaf("type");
  TLeaf * flag_ext = tree->FindLeaf("ext_flag");
  TLeaf * det = tree->FindLeaf("det_id");
	cout<< tree->GetEntries() << endl;

  long int time1;
  int evt = 0;
  vector<pair<long int, int>> sorting;

  TString tmpn; tmpn = "tempf"+to_string(rand()%100);
  TFile *tmpf = new TFile(tmpn,"RECREATE");
  TTree *temp = tree->CloneTree(0);

  temp->Branch("timestp",&time1,"timestp/L");

  int ext_off;
  cout<<"External trigger offset (ns)"<<endl;
  cin>>ext_off;



////////
////    Filtering tree, create pairs
////////
  for(int i=0; i< tree->GetEntries()/1; i++){
	  cout<<100*i/tree->GetEntries()<<"\r"<<flush;
  	tree->GetEntry(i);

	if(type->GetValue()==3){
    time1 = t_0->GetValue();
    sorting.push_back(make_pair(t_0->GetValue(), evt));
    temp->Fill();
    evt++;
    }

	if(type->GetValue()==2 && flag_ext->GetValue()==1){
		time1 = te_0->GetValue() + ext_off; // + 14600; //ext-time-offset
    sorting.push_back(make_pair(te_0->GetValue() + ext_off, evt));
    temp->Fill();
    evt++;
	}
}

tree->DropBuffers(1000000000);
tree->DropBaskets();
file->Close();

  ///////
  /// Sorting
  ///////
cout<<"sorting tree "<<sorting.size()<<endl;
sort(sorting.begin(), sorting.end());
cout<<"sorted data"<<endl;

TFile *outfile = new TFile(ofile, "RECREATE");
outfile->cd();
TTree *newtree;
newtree= temp->CloneTree(0);
newtree->SetDirectory(outfile);
newtree->SetAutoFlush( 30*1024*1024 );	// 30 MB
newtree->SetAutoSave( 100*1024*1024 );	// 100 MB

cout<<"loading baskets"<<endl;
temp->LoadBaskets(1000000);
cout<<"finished baskets"<<endl;

for (int i=0;i<sorting.size();i++){
  temp->GetEntry(sorting[i].second);

  if(i % 100000 ==0 ) cout << 100.0*i/sorting.size()<<"%"<<"\r" << flush; //print the (hopefully sorted) time
  time1 = sorting[i].first;
  newtree->Fill();
}

  outfile->Write();
}
