# L3T-tools
Codes for R3B-Si-Tracker (L3T) analysis

Analysis must be performed in a specific order. This is because L3T outputs raw MIDAS files so first we convert to ROOT trees to perform analysis

MarcMIDAS2ROOT is used with midas file to create file:

R#_0_sort.root

This file is raw hits data but is not fully time ordered due to the buffers in each asic emptying at different points.
Any file may be more or less time order than another mostly dependent on noise/rates of l3t but we assume no time ordering for this file.

# sort.C
This takes in the previous file and time-orders it producing R#_sorted.root. It asks for the run number input, you may need to change the input path directory. 
By default it is ../root/
This code also attempts to filter the external events to a degree by removing any duplicate events that can occur. (This code can take a long time for larger files)

# dt_scan.C
This scans time differences between events in several combinations of L3T-L3T hits as well as external hits-L3T to do time synchronisations.

# EventBuilder.C
Main code. This takes input R#_sorted.C and constructs multi-hit events. A time window of 1000 ns is applied in which successive hits in this window are combined into a total "Event". A vector element for each hit in the event is filled with each hit's data.
<br />
This formats an event into a tree with structure and meaning <br />
evt (int) --> Event number from midas <br />
det (vector<int>) --> det id of each hit <br />
side(vector<int>) --> side id of each hit <br />
  <br />
mult0[3] --> Strip multiplicity of det0; [0] side0; [1]side1; [2] total <br />
mult1[3] --> "==================== det1; ==============================" <br />
th0[3] --> "True hits" Strip multiplicity with condition that coincidenct strips must be adjacent +-5 Det0  <br />
th1[3] --> "========================================================================================" Det1 <br />
  <br />
adc --> adc data (energy) of each hit <br />
timestp_ --> Timestamp of each hit <br />
ch  --> channel id from each hit <br />
asic --> ASIC id from each hit <br />
ext_flag --> == 1 if coincident with external event otherwise =0 <br />
tmp_ext  --> external timestamp if ext_flag==1 <br /> 
ext_mult --> Multiplicity of external events <br />
  <br />
strip_id --> Reconstructed strip number for each hit (ch + 128*asic) <br />
  
# mult_scan.C
Custom script for strip multiplicity analysis of Events
  
# corr.C
Constructs correlation plots for online analysis for L3T,
  
  front-rear correlations <br />
  multiple det correlations <br />
  positional reconstruction <br />
  multiplicity plots <br />
  
