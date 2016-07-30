

#for cut in \
# "n_trig" \
# "n_spike" \
# "n_met" \
# "n_metfilt" \
# "n_dphipho" \
# "n_lepveto" \
# "n_dphijet" 
#  
# do
#sed -n -e 's/^.*'"${cut}"' //p'  /nfs_scratch/tperry/Madrid_met04-SinglePhoton_callpostAnalyzer_Signal/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_*/SinglePhoton_callpostAnalyzer_Signal-ggtree_data_*.out > ${cut}.py
# 
# done


for cut in \
 "n_trig" \
 "n_spike" \
 "n_met" \
 "n_metfilt" \
 "n_dphipho" \
 "n_lepveto" \
 "n_dphijet" 
  
 do
  echo $cut
  python ${cut}.py
 done
