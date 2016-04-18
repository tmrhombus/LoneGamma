echo
echo "submitbase: ${submitbase}"
echo "version: ${version}"

mkdir -p ${submitbase}/gitignore/${version}/eventcomp
printf "${submitbase}/gitignore/${version}/eventcomp \n"

thedir="${submitbase}/gitignore/${version}/eventcomp"

## "ZLLG" \
## "ZJets100To200" \
## "ZJets200To400" \
## "ZJets400To600" \
## "ZJets600ToInf" \
 ##"SinglePhoton"
for samplename in \
 "GJets" \
 "QCD"

do
   grep -h 'Matched    a  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_a.txt
   grep -h 'Matched    b  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_b.txt
   grep -h 'Matched    c  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_c.txt
   grep -h 'Matched    d  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_d.txt
   grep -h 'Matched    e  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_e.txt
   grep -h 'Matched    f  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_f.txt
   grep -h 'Matched    g  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_g.txt
   grep -h 'Matched    h  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_h.txt
   grep -h 'Matched    i  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_M_i.txt

   grep -h 'Unmatched  a  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_a.txt
   grep -h 'Unmatched  b  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_b.txt
   grep -h 'Unmatched  c  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_c.txt
   grep -h 'Unmatched  d  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_d.txt
   grep -h 'Unmatched  e  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_e.txt
   grep -h 'Unmatched  f  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_f.txt
   grep -h 'Unmatched  g  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_g.txt
   grep -h 'Unmatched  h  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_h.txt
   grep -h 'Unmatched  i  --- *[0-9]:*[0-9]*:[0-9]' /nfs_scratch/tperry/${version}*${samplename}*_callpostAnalyzerMC_purity/*/*out > ${thedir}/events_${samplename}_U_i.txt


done # for samplename in "GJets_*" "SingleP"

#sed -i 's/\s*$//' ${thedir}/events_${samplename}.txt
#
#grep -v -F -x -f EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_menothey.txt
#grep -F -x -f EventList.txt ${thedir}/events_SinglePhoton.txt > ${thedir}/events_common.txt
#grep -v -F -x -f ${thedir}/events_SinglePhoton.txt EventList.txt > ${thedir}/events_theynome.txt


