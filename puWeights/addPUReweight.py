import ROOT
from array import array
from pileUpVertexCorrections import PUreweight
import math



puDict = PUreweight()
#print puDict

f = ROOT.TFile('../ttree.root','UPDATE')
d = f.Get('demo/events')

t = d.Get('Ntuple')

puweight = array('f', [0])
puweightB = t.Branch('puweight', puweight, 'puweight/F')

count = 0
for i in range( t.GetEntries() ) :
    t.GetEntry( i )
    if count % 10000 == 0 : print "Event:",count
    nTrPu = ( math.floor(t.nTruePU * 10))/10
    puweight[0] = puDict[ nTrPu ]
    puweightB.Fill()
    count += 1

print "DONE!"

d.cd()
t.Write('', ROOT.TObject.kOverwrite)
f.Close()
