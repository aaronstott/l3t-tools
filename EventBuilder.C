#include <vector>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TFile.h>

void EventBuilder() {
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!============== L3T EVENTBUILDER ==============!"<<endl;
  cout<<"!============== A. Stott - UoYork =============!"<<endl;
  cout<<"!==============================================!"<<endl;
  cout<<"!==============================================!"<<endl;


  time_t now = time(0);
  TString infile1, ofile;
  cout<<"input file number"<<endl;
  cin>>infile1; ofile=infile1;
  infile1="../root/R"+infile1+"_sorted.root";

  TFile *file = new TFile(infile1,"READ");
  TTree *tree = (TTree  *) file->Get("R3B_sort");

  TLeaf * t = tree->FindLeaf("timestp");
  TLeaf * ti = tree->FindLeaf("tm_stp");
  TLeaf * ch = tree->FindLeaf("ch_id");
  TLeaf * det = tree->FindLeaf("det_id");
  TLeaf * side = tree->FindLeaf("side_id");
  TLeaf * type = tree->FindLeaf("type");
  TLeaf * adc  = tree->FindLeaf("adc_data");
  TLeaf * asic = tree->FindLeaf("asic_id");
  TLeaf * ext = tree->FindLeaf("ext_flag");
  TLeaf * te  = tree->FindLeaf("tm_stp_ext");

  ofile = "../root/R"+ofile+"_events.root";
  TFile *outfile = new TFile(ofile,"RECREATE");
  TTree *newtree = new TTree("R3B_evt", "R3B_evt");

  int evt, det1;
  std::vector<int> det0;
  newtree->Branch("evt",&evt,"evt/L"); //Event Number
  newtree->Branch("det",&det0); //Detector Number

  std::vector<int> side0;
  newtree->Branch("side",&side0);

  std::vector<int> mult0 = {0, 0, 0}; //[side0, side1, total]
  std::vector<int> mult1 = {0, 0, 0}; //[side0, side1, total]
  int mult[2][3];

  // Det0
  newtree->Branch("mult0", &mult0);//multiplicity, [side0, side1, total]
  // Det1
  newtree->Branch("mult1", &mult1);//multiplicity, [side0, side1, total]

  std::vector<int> adc0; //
  newtree->Branch("adc",&adc0); //adc_id

  std::vector<long long int> timestp_0;
  newtree->Branch("timestp_",&timestp_0); //timestamps side0

  std::vector<int> ch0;
  newtree->Branch("ch",&ch0); //ch_id

  std::vector<int> asic0;
  newtree->Branch("asic",&asic0); //asic_id

  int fe;
  newtree->Branch("ext_flag", &fe); //external flag

  Long64_t ext_time;
  newtree->Branch("tmp_ext", &ext_time); // External timestamp

  int em;
  newtree->Branch("ext_mult", &em);

//=================================================================///
  std::vector<int> strip_id; //strip
  newtree->Branch("strip_id",&strip_id); //converted strip_id

  std::vector<int> th0 = {0, 0, 0}; //[side0, side1, total]
  std::vector<int> th1 = {0, 0, 0}; //[side0, side1, total]
  int th[2][3], s0[2][3];
  int tolerence = 5; //strip difference for tolerence.

  newtree->Branch("th0", &th0);//True hits, [side0, side1, total]
  newtree->Branch("th1", &th1);//True hits, [side0, side1, total]
///////////////////////////////////
  //initalise values
  evt = 0;
  int i=0;
  int j=1;
  long long int t_gate=2000; // in ns.
  long long int time_ref0, time_ref1;
  int util, detutil;

  timestp_0.clear();
  det0.clear();
  asic0.clear();
  side0.clear();
  ch0.clear();
  adc0.clear();
  strip_id.clear();
  fe=0; //assume no external
  ext_time=0;

  // Rate Checking
  int ext_rate=0;
  int nds[2][2]; nds[0][0]=0; nds[0][1]=0; nds[1][0]=0; nds[1][1]=0; //Total averaged rates
  int ndss[2][2], ext_srate=0; //dt rates ndss=l3t ext_srate=ext
  long long int t_0, t_f;
  long long int t_0e, t_fe, ite;
  long long int dt, dtn, itn; dt=pow(10,9); dtn=dt/pow(10,9);//dt=time window, dtn=normaliser [nanosecond], itn=Counter for timewindow
  double tl_sec;
  tree->GetEntry(1); t_0=t->GetValue(); t_0e=te->GetValue();
  tree->GetEntry(tree->GetEntries() - 1); t_f=t->GetValue(); t_fe=te->GetValue();
  ite=0; itn=0; //External/Internal rate counters

//////////////////////////////////////////////
//
//Rates Histograms
//
///////////////////////////////////////////////
  TH1D *l3trates[2][2], *ext_rates;
  ext_rates = new TH1D("ext_rate", "external trigger rate [hz]", int((t_fe-t_0e)/dt), 0, (t_fe-t_0e)/pow(10,9));
  for(int d=0; d<2; d++){
    for(int s=0; s<2; s++){
      l3trates[d][s] = new TH1D(Form("rate_d%i_s%i", d, s), Form("rate_d%i_s%i", d, s), int((t_f-t_0)/dt), 0, (t_f-t_0)/pow(10,9));
    }
  }
////////////////////////////////////////////


  tl_sec  =  (t_f - t_0)/(pow(10, 9));
  cout<<"L3T Timestamps   t_0, t_f "<<t_0 <<" , "<<t_f <<endl;
  cout<<"External trigger t_0, t_f "<<t_0e<<" , "<<t_fe<<endl;
  cout<<"Total time of run: "<< (t_f - t_0)/(pow(10, 9)) << " s"<<endl;
  cout<<"Average hit rate: "<< tree->GetEntries()/tl_sec << " hz"  << endl;

  while(i < tree->GetEntries()){ //Read Tree
    tree->GetEntry(i);
    if(ext->GetValue()==1 && fe==0){fe=1; ext_time=t->GetValue(); em=1; ext_rate++;} //// Give external flag
    if(ext->GetValue()==0){
      time_ref0=t->GetValue();
      timestp_0.push_back(t->GetValue());
      det0.push_back(det->GetValue()); nds[int(det->GetValue())][int(side->GetValue())]++;
      asic0.push_back(asic->GetValue());
      side0.push_back(side->GetValue());
      ch0.push_back(ch->GetValue());
      strip_id.push_back(ch->GetValue() + 128*asic->GetValue());
      adc0.push_back(adc->GetValue());
      mult[det0[0]][side0[0]]=1;
      mult[det0[0]][2]=1;

      ///Realistic hit checking
      s0[int(det->GetValue())][int(side->GetValue())] = ch->GetValue() + 128*asic->GetValue(); //fills first strip as an "anchor"
      th[int(det->GetValue())][int(side->GetValue())] = 1; //true hit 1
      th[int(det->GetValue())][2]=1;
    }

    for(j=i+1;j<i+100;j++){
      tree->GetEntry(j);

      //
      // Rate checking every dt [ns]
      //
      //L3T Events
      if(ext->GetValue()==0){
        if((t->GetValue()- (t_0+dt*itn))<dt && mult[int(det->GetValue())][int(side->GetValue())]==0){ //check window + only consider 1st event in window
          ndss[int(det->GetValue())][int(side->GetValue())]++;
        } //L3T event
        else {
          for(int d=0; d<2; d++){
            for(int s=0; s<2; s++){
                l3trates[d][s]->Fill(dt*(itn)/pow(10,9),ndss[d][s]*dtn); //Fill Hist
                ndss[d][s]=0; //Reset
            }
          }
          itn++; //increase time-bin
        }
      }

      if(ext->GetValue()==1){
        if((t->GetValue()- (t_0+dt*ite))<dt){ext_srate++;} //L3T event
        else {
                ext_rates->Fill(dt*(ite)/pow(10,9),ext_srate*dtn); //fill hist
                ext_srate=0; //reset external rate-count
                ite++; //increase time-bin
              }
      }



      time_ref1=t->GetValue();
        if(abs(time_ref0-time_ref1)>t_gate){ //check coincidence

          for(int k=0; k<3; k++){mult0[k]=mult[0][k]; mult1[k]=mult[1][k];}
          for(int k=0; k<3; k++){th0[k]=th[0][k]; th1[k]=th[1][k];}

          newtree->Fill();
          i=j;
          evt++;
          break; // Jump timestamp to next point
        }
          if(ext->GetValue()==1){fe=1; ext_time=t->GetValue(); em++; ext_rate++;} //// Give external flag
          if(ext->GetValue()==0){
            if(abs(ti->GetValue()-(te->GetValue())) <t_gate){fe=1; ext_time=te->GetValue();}

            timestp_0.push_back(t->GetValue());
            det0.push_back(det->GetValue()); nds[int(det->GetValue())][int(side->GetValue())]++;
            asic0.push_back(asic->GetValue());
            side0.push_back(side->GetValue()); util=side->GetValue();
            ch0.push_back(ch->GetValue());
            strip_id.push_back(ch->GetValue() + 128*asic->GetValue());
            adc0.push_back(adc->GetValue());
            mult[int(det->GetValue())][util]++;
            mult[int(det->GetValue())][2]++;

            // Realistic hit checking
            if(s0[int(det->GetValue())][int(side->GetValue())] == 0){
              s0[int(det->GetValue())][int(side->GetValue())] = ch->GetValue() + 128*asic->GetValue();
              th[int(det->GetValue())][int(side->GetValue())]++;
              th[int(det->GetValue())][2]++;
            }
            else if(abs(s0[int(det->GetValue())][int(side->GetValue())]-(ch->GetValue()+128*asic->GetValue()))<=tolerence){
              th[int(det->GetValue())][int(side->GetValue())]++;
              th[int(det->GetValue())][2]++;
            }
        }
    }

	cout<<100.0*i/tree->GetEntries()<<"\r"<<flush;
    //prevent breaks
    if(j==i+100){i=j;}
    timestp_0.clear();
    det0.clear();
    asic0.clear();
    side0.clear();
    ch0.clear();
    strip_id.clear();
    adc0.clear();
    fe=0; //assume no external
    ext_time=0;
    em=0;

    for(int l=0; l<2; l++){
      for(int k=0; k<3; k++){
        mult[l][k] = 0;
        th[l][k]   = 0;
        s0[l][k]   = 0;
      }
    }
  }

  outfile->cd();
  newtree->Write();
  cout<<"Event Rate: "<<newtree->GetEntries()/tl_sec<< " hz"<<endl;

  for(int i=0; i<2; i++){
        cout<<"\nL3T"<<i<<" Rate = "<<(nds[i][0]+nds[i][1])/tl_sec<<" hz"<<endl;
    for(int j=0; j<2; j++){
        cout<<"     side"<<j<<"= "<<nds[i][j]/tl_sec<<" hz"<<endl;
        l3trates[i][j]->Write();
    }
    cout<<endl;
  }

  ext_rates->Write();
  cout<<"\n\nExternal Rate: "<<ext_rate/tl_sec<<" hz"<<endl;


  outfile->Close();

  time_t now2 = time(0);
  time_t const timediff = now2 - now;
  int minute = timediff / 60;
  int second = timediff % 60;
  cout << "\nTime Taken = " << minute << " Mins " << second << " Secs" << endl;
}
