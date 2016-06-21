# LoneGamma
Standalones for MonoPhoton 13TeV Analysis

export SCRAM_ARCH=slc6_amd64_gcc530
scram pro -n LoneG_slc6_530_CMSSW_8_0_8 CMSSW CMSSW_8_0_8
cd LoneG_slc6_530_CMSSW_8_0_8/src/
cmsenv

git clone https://github.com/tmrhombus/LoneGamma.git

cd LoneGamma

git checkout CMSSW_8_0_8

cd commontools/
bash mklists.sh     --  make lists of input files
bash countlists.sh  --  count total number of events in each of the lists and match with cross section

# QCD Study
cd ../qcdStudy

source setup.sh # set version number here

-- test running - uses list in ../test/ directory, pointed to in callQCDAnalyzer.cc
nohup root -l -b -q callQCDAnalyzer.cc > out.out 2>&1& 

-- after satisfied that everything is running fine, submit with

nohup bash fajSubmit_postAnalyzer_QCD.sh > gitignore/${version}/submit.out 2>&1&
-- there are flags to create the submit directory structure without actually submitting inside fajSubmit_postAnalyzer_QCD.sh

-- purity is basically the same way, but things haven't yet been updated from the 2015 setup


-- Z(nn)G Study
cd ${CMSSW_BASE}/src/LoneGamma/znunugStudy
source setup.sh

-- same thing with callZnunuGAnalyzer.cc to test and fajSubmit_postAnalyzer_ZnunuG.sh to submit
