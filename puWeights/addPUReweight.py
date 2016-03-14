import ROOT
from array import array
from pileUpVertexCorrections import PUreweight
import math



puDict = PUreweight()
#print puDict

f = ROOT.TFile('ttree.root','UPDATE')
d = f.Get('demo/events')

t = d.Get('Ntuple')

puweight = array('f', [0])
puweightB = t.Branch('puweight', puweight, 'puweight/F')

count = 0
firstMass = 0
for i in range( t.GetEntries() ) :
#for row in t :
    t.GetEntry( i )
    #if row.genMass == firstMass : 
    #if t.genMass == firstMass : 
    #    print "Looped: count:",count
    #    break
    #if count == 0 : 
        #firstMass = row.genMass
    #    firstMass = t.genMass
    #    print "First Mass: %f" % firstMass
    #if round(row.genMass,9) == round(79.6630020142,9) :
    #if round(t.genMass,9) == round(79.6630020142,9) :
        #print "%f !!!! count: %i " % (row.genMass, count)
    #    print "%f !!!! count: %i " % (t.genMass, count)
    if count % 1000 == 0 : print "Event:",count
    #nTrPu = ( math.floor(row.nTruePU * 10))/10
    nTrPu = ( math.floor(t.nTruePU * 10))/10
    puweight[0] = puDict[ nTrPu ]
    #if count % 1000 == 0 : print "Event:",count,"puweight:",puweight[0]
    t.Fill()
    count += 1
t.Write()
f.Close()
