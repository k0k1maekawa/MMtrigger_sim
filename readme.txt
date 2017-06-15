Stand Alone MM Trigger Code  by Koki Maekawa  email:koki.maekawa@cern.ch

Please read following procedures
1. and 2.


1. Setting Variables
export BG="1"
#use BackGround(1) or not(0)
export CV="0"
#use Cavern(1) or not(0)
export MU="80"
#pile up events of minBias
export ARTDEADTIME="40"
#ART dead time --ns
export STRIPDEADTIME="300"
#strip dead time ---ns
export THRESHOLD="1250"
#VMM charge threshold in Digitization
export NAMECORE=MC12.NSW.107206.atlasG4.01.0000.Pt.100.generator.HITS_StripDeadtime-300_ARTDeadTime-40_electronicsThreshold-1250_WB
#the name of root file before ".root"
export DRCUT="10"
#using this value of localx (=Rcosdphi) distance from truth hit position like 10mm, recognize signals as coming from muon or secondary particles
export ETASLOPEWIDTH="0.0015"
#slope coincidence width for eta
export STEREOSLOPEWIDTH="0.004"
#slope coincidence width for stereo
export BCstart="0"
#in trackfinding, we start to find tracks with which BC 
export NDATA="1"
#the number of root data file for the situation that there are a number of NAMECORE_*.root files.



2. Using macros

The code is separated into some parts.
readalpha.cc/memo.cc/timetrig.cc/fastalpha.cc/difalpha.cc/resid.cc

Overview

root -q -b -l readalpha.cc+ & root -q	-b -l memo.cc+
root -q -b -l timetrig.cc+
root -q -b -l fastalpha.cc+
root -q -b -l difalpha.cc+
root -q -b -l resid.cc+

readalpha.cc : to read data, translate signals to xyz,localx and calculate timing subtracting TOF from global time
NAMECORE.root -> NAMECORExyz.root 
memo.cc : to memo truth particle eta or vertex_n and so on of each event
NAMECORE.root -> NAMECOREmemo.root
timetrig.cc :  to find tracks and calculate track finding efficiencies of many parameter conditions
NAMECORExyz.root & NAMECOREmemo.root -> NAMECOREtag.root(main) & NAMECOREsel.root(for debugging)
fastalpha.cc : to change format of track segments
NAMECOREtag.root -> NAMECOREfast.root
difalpha.cc : to calculate segment angles (eta,phi,dtheta or elses)
NAMECOREfast.root & NAMECORE.root & NAMECOREmemo.root -> NAMECOREd3.root
resid.cc : to calculate residuals and resolutions
NAMECOREd3.root -> show resolutions and make residual pdfs

Details

there's some switches in macro
Adding more notice to coments written in the macro.

timetrig.cc:
//Overall
using "const int option", you can basically choose the purpose of triggering.
"option = 1" is default and other options don't make an output root file and they are for drawing graphs. 
//selection
there are options to select one from VMM signals (if there is several VMM signals in the same slope road)  
"biason = 1" is default.

difalpha.cc:
"uvfit = 1" was used by Koki by default,
but according to the atlas internal note, we see the "uvfit = 2" option is more actual.

