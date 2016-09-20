
DATE=Sept20

echo "Do the Data samples"

for RUN in RunB RunC RunD RunE RunF; do
    farmoutAnalysisJobs \
        --input-file-list=targetRunSingleMuon${RUN}.txt \
        --assume-input-files-exist \
        --input-files-per-job=10 \
        --lumi-mask=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-278808_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt \
        tagAndProbe${DATE}${RUN} $CMSSW_BASE ConfFile_Data_cfg.py
done



echo "Do the MC samples"

farmoutAnalysisJobs \
    --input-file-list=targetDYJets.txt \
    --assume-input-files-exist \
    --input-files-per-job=10 \
    tagAndProbe${DATE}DYJets $CMSSW_BASE ConfFile_MC_reHLT_cfg.py

farmoutAnalysisJobs \
    --input-file-list=targetAZH220.txt \
    --assume-input-files-exist \
    --input-files-per-job=2 \
    tagAndProbe${DATE}AZH220 $CMSSW_BASE ConfFile_MC_reHLT_cfg.py

farmoutAnalysisJobs \
    --input-file-list=targetGluGluHiggs125.txt \
    --assume-input-files-exist \
    --input-files-per-job=1 \
    tagAndProbe${DATE}ggH125 $CMSSW_BASE ConfFile_MC_reHLT_cfg.py



