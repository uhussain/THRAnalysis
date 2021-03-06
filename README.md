# Tau Tag and Probe
Get previous HLT_doubleTau35 work

```
cd TagAndProbe
source setup.sh
```


# GenTauAnalysis

`cd GenAnalysis`

Package used to calculate acceptance values for Z->TauTau analysis

To operate:
1. run the analyzer over the DYJets sample file provided

`cmsRun python/ConfFile_cfg.py`

2. add pileup reweighting distribution

`cd puWeights`
`python addPUReweight.py`
`cd ..`

3. print out initial and fiducal cut yields for each channel

`python printWeights.py`

4. for detailed efficiency calculations print out file with all events for each channel passing fiducal cuts + another list with passing genMass window.  This list will be compared to the events in the final signal region

`python genEventInfo.py`

5. copy over your selected final cut based analysis tree to the main directory here (AcceptanceAnalyzer)

`cp .......root .`

6. compare analysis selected events to the gen selected events and dump results to output file

`python recoEventInfo.py > dump.txt`

7. print results of the above dump.  Not sure why I chose this route.

`python finalResults.py`
