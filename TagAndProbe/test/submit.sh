
DATE=Sept21v3

echo "Do the Data samples"
# Different b/c of Prompt RECO
farmoutAnalysisJobs \
    --assume-input-files-exist \
    --input-dbs-path=/SingleMuon/Run2016F-PromptReco-v1/MINIAOD \
    --input-files-per-job=5 \
    --lumi-mask=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-278808_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt \
    tagAndProbe${DATE}Run2016F $CMSSW_BASE ConfFile_Data_cfg.py

for RUN in Run2016B Run2016C Run2016D Run2016E; do
    farmoutAnalysisJobs \
        --input-dbs-path=/SingleMuon/${RUN}-PromptReco-v2/MINIAOD \
        --assume-input-files-exist \
        --input-files-per-job=5 \
        --lumi-mask=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-278808_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt \
        tagAndProbe${DATE}${RUN} $CMSSW_BASE ConfFile_Data_cfg.py
done



echo "Do the MC samples"

farmoutAnalysisJobs \
    --input-file-list=targetDYJets.txt \
    --assume-input-files-exist \
    --input-files-per-job=5 \
    tagAndProbe${DATE}DYJets $CMSSW_BASE ConfFile_MC_reHLT_cfg.py

farmoutAnalysisJobs \
    --input-file-list=targetAZH220.txt \
    --assume-input-files-exist \
    --input-files-per-job=5 \
    tagAndProbe${DATE}AZH220 $CMSSW_BASE ConfFile_MC_reHLT_cfg.py

farmoutAnalysisJobs \
    --input-file-list=targetGluGluHiggs125.txt \
    --assume-input-files-exist \
    --input-files-per-job=5 \
    tagAndProbe${DATE}ggH125 $CMSSW_BASE ConfFile_MC_reHLT_cfg.py



