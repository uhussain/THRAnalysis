import ROOT

chans = {
    'tt' : 'TauTau',
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
    for line in genFile :
        info = line.strip().split(' ')
        genEvtInfo = (float(info[0]), float(info[1]), float(info[2]))
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

        if evtInfo in inReco :
            totWeightMatch += evtW
            cntMatch += 1

        if evtInfo in inWindow :
            totWeightWindow += evtW
            cntWindow += 1

        # run is redundant for MC, keeping for consistency w/ data style
        ofile.write('%s %s %s\n' % (run, lumi, evt) )


        cnt += 1

    print "Number of events: %i    Total weight: %f" % (cnt, totWeight)
    print "Number of events match gen->reco: %i    Total weight match: %f" % (cntMatch, totWeightMatch)
    print "Number of events in window: %i    Total weight in window: %f" % (cntWindow, totWeightWindow)
    print "\n"
    ofile.close()    


if __name__ == '__main__' :
    
    channels = ['tt','em']
    
    for chan in channels :
        print chan
        eventMap( chan )
