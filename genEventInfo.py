import ROOT


def eventMap( chan ) :
    ofileName = chan+'Events.txt'
    ofile = open( ofileName, 'w')
    
    for row in t :
        passing = getattr( row, chan+'Pass' )
        if passing != 0 :
            run = int(row.run)
            lumi = int(row.lumi)
            evt = int(row.eventD)
            genMass = row.genMass 
            pt1 = row.tauPt1
            pt2 = row.tauPt2
            # run is redundant for MC, keeping for consistency w/ data style
            ofile.write('%i %i %i %f %f %f\n' % (run, lumi, evt, genMass, pt1, pt2) )
    ofile.close()    


if __name__ == '__main__' :
    f = ROOT.TFile('ttree.root','r')
    t = f.Get('demo/events/Ntuple')
    
    channels = ['ETau', 'MuTau', 'EMu', 'TauTau', 'TauTau4030', 'MuMu']
    
    for chan in channels :
        print chan
        eventMap( chan )
