import ROOT

chans = {
    'tt' : 'TauTau',
    'tt4030' : 'TauTau4030',
    'et' : 'ETau',
    'mt' : 'MuTau',
    'em' : 'EMu',
    'mm' : 'MuMu',
    }


def eventMap( chan ) :
    #f = ROOT.TFile('DYJetsEff_%s2.root' % chan,'r')
    #f = ROOT.TFile('DYJetsEff_%s.root' % chan,'r')
    f = ROOT.TFile('DYJets_%s.root' % chan,'r')
    #f = ROOT.TFile('dy3_%s.root' % chan,'r')
    t = f.Get('Ntuple')
    ofileName = chan+'RecoEvents.txt'
    ofile = open( ofileName, 'w')

    # List of events passing gen window
    genFile = open( '%sEvents.txt' % chans[chan], 'r' )
    inWindow = set()
    inReco = set()
    notInReco = set()
    for line in genFile :
        info = line.strip().split(' ')
        genEvtInfo = (int(info[0]), int(info[1]), int(info[2]))
        inReco.add( genEvtInfo )
        if float(info[3]) > 60 and float(info[3]) < 120 :
            inWindow.add( genEvtInfo )

    
    totWeight = 0.0
    totWeightWindow = 0.0
    totWeightMatch = 0.0
    cnt = 0
    cntWindow = 0
    cntMatch = 0
    for row in t :
        if chan == 'tt' :
            if row.t1ByVTightIsolationMVArun2v1DBoldDMwLT < 0.5 : continue
            if row.t2ByVTightIsolationMVArun2v1DBoldDMwLT < 0.5 : continue
        run = int(row.run)
        lumi = int(row.lumi)
        evt = int(row.evt)
        evtInfo = (run, lumi, evt)
        #print evtInfo

        # Calc event weight
        trigW = row.trigweight_1
        idIso1 = row.idisoweight_1
        idIso2 = row.idisoweight_2
        puweight = row.puweight
        genW = row.weight
        evtW = trigW * idIso1 * idIso2 * puweight * genW
        #print "evtW: %f" % evtW
        totWeight += evtW

        #if evtInfo in inReco :
        if row.gen_match_1 == 5 and row.gen_match_1 == 5 :
            totWeightMatch += evtW
            cntMatch += 1
            ofile.write('%s %s %s %s %s %s Gen_Match5\n' % (run, lumi, evt, evtW, row.pt_1, row.pt_2 ) )
        else :
            ofile.write('%s %s %s %s %s %s Gen_Match_Not5\n' % (run, lumi, evt, evtW, row.pt_1, row.pt_2 ) )

        if evtInfo in inWindow :
            totWeightWindow += evtW
            cntWindow += 1

        cnt += 1

    print "Number of events: %i    Total weight: %f" % (cnt, totWeight)
    print "Number of events match gen->reco: %i    Total weight match: %f" % (cntMatch, totWeightMatch)
    print "Number of events in window: %i    Total weight in window: %f" % (cntWindow, totWeightWindow)
    print "\n"
    ofile.close()    




# Take our list of non-gen-matched events
# and see if there is a strong low pt issue
def printNonMatchPts( chan ) :

    allAnalysisFile = open( chan+'RecoEvents.txt', 'r' )
    allPassingSet = set()
    weightMap = {}
    for line in allAnalysisFile :
        info = line.strip('\n').strip(' ').split(' ')
        evtInfo = (int(info[0]), int(info[1]), int(info[2]))
        allPassingSet.add( evtInfo )
        weightMap[ evtInfo ] = (float(info[3]), float(info[4]), float(info[5]), info[6] )

   
    f = ROOT.TFile('ttree.root', 'READ')
    t = f.Get('demo/events/Ntuple')
    cnt = 0
    f_out = 0
    f_outW = 0
    not_f_out = 0
    not_f_outW = 0
    f_out_strs = []
    not_f_out_strs = []
    numEvts = 0
    for row in t :
        numEvts += 1
        run = int(row.run)
        lumi = int(row.lumi)
        evt = int(row.eventD)
        evtInfo = (run, lumi, evt)

        if evtInfo in allPassingSet :
            str1 = "AN_Pass_Cuts: GenMass %.2f  Weight %.3f  GenPt1 %.2f  GenPt2 %.2f  GenPt3 %.2f   AnPt1 %.2f  AnPt2 %.2f  %s" % (row.genMass, weightMap[evtInfo][0], row.tauPt1, row.tauPt2, row.tauPt3, weightMap[evtInfo][1], weightMap[evtInfo][2], weightMap[evtInfo][3] )
            if row.genMass < 120 and row.genMass > 60 :
                f_out += 1
                #print "F_out", f_out
                f_out_strs.append( str1 )
                f_outW += weightMap[ evtInfo ][0]
            if row.genMass > 120 or row.genMass < 60 :
                not_f_out += 1
                #print "NOT F_out", not_f_out
                not_f_outW += weightMap[ evtInfo ][0]
                not_f_out_strs.append( str1 )
            allPassingSet.remove( evtInfo )
 

    print "\n F_OUT_STRS"
    for stri in f_out_strs :
        print stri
    print "\n\n NOT F_OUT_STRS"
    for stri in not_f_out_strs :
        print stri
    print "\n\n NO MATCHES between An and Reco"
    for evt in allPassingSet :
        print evt

    print "F_out: ",f_out,"    F_out Weight: ",f_outW
    print "NOT F_out: ",not_f_out,"    NOT F_out Weight: ",not_f_outW
            
    print "TOTAL GEN ANALYSIS EVTS: ",numEvts

        


if __name__ == '__main__' :
    
    channels = ['tt4030',]#'em']
    
    for chan in channels :
        print chan
        eventMap( chan )

        printNonMatchPts( chan )



