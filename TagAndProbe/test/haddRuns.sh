DATE=Sept20

mkdir -p /data/truggles/TAP_${DATE}_hadd/

for RUN in RunB RunC RunD RunE RunF; do

    hadd -f /data/truggles/TAP_${DATE}_hadd/${RUN}.root /hdfs/store/user/truggles/tagAndProbe${DATE}${RUN}-ConfFile_Data_cfg/*.root

done

# Combined Data Groups
hadd -f /data/truggles/TAP_${DATE}_hadd/ICHEPRuns.root /data/truggles/TAP_${DATE}_hadd/Run[B,C,D].root
hadd -f /data/truggles/TAP_${DATE}_hadd/AllRuns.root /data/truggles/TAP_${DATE}_hadd/Run*.root

# MC Samples
hadd -f /data/truggles/TAP_${DATE}_hadd/DYJets.root /hdfs/store/user/truggles/tagAndProbe${DATE}DYJets-ConfFile_MC_reHLT_cfg/*.root 
