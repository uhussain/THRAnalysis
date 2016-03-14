import ROOT
from ROOT import gPad
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)

f = ROOT.TFile('ttree.root','r')

t = f.Get('Ntuple')

'''
w/o weights
w/o + f window
w/ weights
w/ + f windo
'''

fWindow = '*(genMass > 60 && genMass < 120)'
puweight = '*(puweight)'


channels = ['ETau', 'MuTau', 'EMu', 'TauTau']
mapperAll = OrderedDict()
mapperPass = OrderedDict()

for chan in channels :
    mapperPass[ chan+' w/o weights' ] = [ chan+'Pass', ]
    mapperPass[ chan+' w/o weights f_out' ] = [ chan+'Pass', fWindow ]
    mapperPass[ chan+' w/ weights' ] = [ chan+'Pass', puweight ]
    mapperPass[ chan+' w/ weights f_out' ] = [ chan+'Pass', fWindow, puweight ]
    mapperAll[ chan+' w/o weights' ] = [ chan+'D', ]
    mapperAll[ chan+' w/o weights f_out' ] = [ chan+'D', fWindow ]
    mapperAll[ chan+' w/ weights' ] = [ chan+'D', puweight ]
    mapperAll[ chan+' w/ weights f_out' ] = [ chan+'D', fWindow, puweight ]


for i,name in enumerate(mapperPass) :
    chan = name.split(' ')[0]
    #print chan
    cutsAndWeightAll = ''
    for cut in mapperAll[ name ] :
        cutsAndWeightAll += cut
    h1a = ROOT.TH1F('h1a%i' % i, 'h1a%i' % i, 100, 0, 3500 )
    t.Draw('genMass>>h1a%i' % i, cutsAndWeightAll)
    h1a = gPad.GetPrimitive( 'h1a%i' % i )
    h1aInt = h1a.Integral()

    cutsAndWeightPass = ''
    for cut in mapperPass[ name ] :
        cutsAndWeightPass += cut
    h1p = ROOT.TH1F('h1p%i' % i, 'h1p%i' % i, 100, 0, 3500 )
    t.Draw('genMass>>h1p%i' % i, cutsAndWeightPass)
    h1p = gPad.GetPrimitive( 'h1p%i' % i )
    h1pInt = h1p.Integral()


    #print "key: %s,  cut all: %s    cut pass: %s" % (name, cutsAndWeightAll, cutsAndWeightPass)
    # divide by 2 because the Ntuple is double written, could be fixed in addPUReweight.py
    #print " ---  all:",h1aInt/2.
    #print " --- pass:",h1pInt/2.
    #print " --- pass/all: %f" % (h1pInt/(h1aInt*2))
    print "key: %s" % name
    print h1aInt/2.,h1pInt/2.,(h1pInt/h1aInt)
    
