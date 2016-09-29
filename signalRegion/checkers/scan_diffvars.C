
/hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ggtree_mc_1.root


/hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/ggtree_data_1.root

ggNtuplizer->cd();

TString columns = "event:phoSigmaIEtaIEtaFull5x5:phoSigmaIEtaIPhiFull5x5:phoSigmaIPhiIPhiFull5x5:nPho:phoE:phoSCRawE:phoseedTimeFull5x5";

TString cuts = "event==8545921|| event==13254717 || event==13956965 || event==1266913243 || event==1266666690 || event==1270132990"

((TTreePlayer*)(EventTree->GetPlayer()))->SetScanRedirect(true);
((TTreePlayer*)(EventTree->GetPlayer()))->SetScanFileName("scan_Data.txt");

EventTree->Scan(columns,cuts,"colsize=10")





TString cuts = "event==65464584 || event==5289 || event==5298 || event==5416 || event==124407 || event==124498"
TString cuts = "event==8545921|| event==13254717|| event==13956965|| event==1266913243|| event==1266666690|| event==1270132990|| event==1273281371|| event==1285416440|| event==69113540|| event==72461657|| event==77315330 "



TFile* file_ZnnG = new TFile("/hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ggtree_mc_1.root");
TTree* tree_ZnnG = file_ZnnG.Get("ggNtuplizer/EventTree");
TTree* tree_ZnnG = _file0.Get("ggNtuplizer/EventTree");

TString columns = "event:phoSigmaIEtaIEtaFull5x5:phoSigmaIEtaIPhiFull5x5:phoSigmaIPhiIPhiFull5x5:nPho:phoE:phoSCRawE:phoseedTimeFull5x5";

TString cuts = "event==65464584 || event==5289 || event==5298 || event==5416 || event==124407 || event==124498"

((TTreePlayer*)(New_Tree->GetPlayer()))->SetScanRedirect(true);
((TTreePlayer*)(New_Tree->GetPlayer()))->SetScanFileName("scan_Data.txt");

tree_ZnnG->Scan(columns,cuts,"colsize=10")

