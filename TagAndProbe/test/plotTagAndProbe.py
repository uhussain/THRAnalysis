import ROOT
from ROOT import gROOT
import pyplotter.tdrstyle as tdr
from array import array
from crystalBall import CBeff
from array import array
import json
import math


ROOT.gROOT.SetBatch(True)
tdr.setTDRStyle()
DATE='Sept20'

def buildLegend( items, names ) :
    legend = ROOT.TLegend(0.45, 0.75, 0.95, 0.93)
    legend.SetMargin(0.3)
    legend.SetBorderSize(0)
    for item, name in zip(items, names) : #range(0, stack.GetStack().GetLast() + 1) :
        legend.AddEntry( item, name, 'l')
    return legend

def decorate(cmsLumi) :
    logo = ROOT.TText(.2, .88,"CMS Preliminary")
    logo.SetTextSize(0.03)
    logo.DrawTextNDC(.2, .89,"CMS Preliminary")
    
    lumi = ROOT.TText(.7,1.05,"X fb^{-1} (13 TeV)")
    lumi.SetTextSize(0.03)
    lumi.DrawTextNDC(.7,.96,"%.1f / fb (13 TeV)" % cmsLumi )

#def fitArcTan( graph ) :
#    xVals = array('d', [i for i in range(0, 161)])
#    yVals = array('d', [0 for i in range(0, 161)])
#    g = ROOT.TGraph(len(xVals), xVals, yVals)
#    print graph.GetTitle()
#    print graph.GetRMS()
#    #func = ROOT.TF2( 'func', '[0] + (x * [1]) +(y *[2])' )
#    #func = ROOT.TF1( 'func', '[0] + [1] * TMath::ATanH(x)' )
#    #func = ROOT.TF1( 'func', '([0] + [1] * TMath::Log((1+x)/(1-x))/2)', 0, 160 )
#
#    #func = ROOT.TF1( 'func', '([0] + [1] * TMath::ATan([2] - x))', 0, 160 )
#    func = ROOT.TF1( 'func', '([0] + [1] * ROOT::Math::crystalball_function(x,[2],[3],[4],[5]) + TMath::Ceil([6]-x))', 0, 160 )
#    f1 = gROOT.GetFunction( 'func' )
#    f1.SetParName( 0, 'vertical' )
#    f1.SetParName( 1, 'scale' )
#    f1.SetParName( 2, 'horizontal' )
#    f1.SetParName( 3, 'gaus mean' )
#    f1.SetParName( 4, 'gaus sigma' )
#    f1.SetParName( 5, 'gaus sigmax' )
#    f1.SetParName( 6, 'ceil' )
#    f1.SetParameter( 0, 0. )
#    f1.SetParameter( 1, 1. )
#    f1.SetParameter( 2, 40. )
#    f1.SetParameter( 3, 40. )
#    f1.SetParameter( 4, 10. )
#    f1.SetParameter( 5, 10. )
#    f1.SetParameter( 6, 40. )
#    g.Apply( func)
#    graph.Fit('func', 'R' )
#    fit = graph.GetFunction('func')
#    print fit.GetParameter(0)
#    print fit.GetParameter(1)
#    g.Apply( graph.GetFunction('func') )
#    g.SetLineWidth(2)
#    g.SetTitle('All Data - ATan')
#    g.SetName('All Data - ATan')
#    return g


def getCBEffGraph( iso ) :
    with open('../triggerSF/di-tau/real_taus_cumulative.json') as f1 :
        inJson = json.load(f1)
    iso += 'Iso'
    alpha = inJson[iso]['alpha']
    m0 = inJson[iso]['m_{0}']
    sigma = inJson[iso]['sigma']
    norm = inJson[iso]['norm']
    n = inJson[iso]['n']
    xVals = array('d', [i for i in range(0, 161)])
    yVals = array('d', [CBeff(i, m0, sigma, alpha, n, norm) for i in range(0, 161)])
    g = ROOT.TGraph(len(xVals), xVals, yVals)
    g.SetLineWidth(2)
    g.SetTitle(iso+' Tau POG')
    g.SetName(iso+' Tau POG')
    return g
    

def getHist( tree, var, cut, name, run, iso ) :
    binning = array('d', [20,22.5,25,27.5,30,32.5,35,37.5,40,\
        42.5,45,47.5,50,55,60,67.5,80,100,160])
    #h = ROOT.TH1F( name, name, 20, 0, 100)
    h = ROOT.TH1F( name, name, len(binning)-1, binning)
    doCut = '(IsoMu22 == 1 && transMass < 30 && m_vis > 40 && m_vis < 80)'
    doCut += '*(tMVAIso%s==1)*(mTrigMatch>0)*(leptonDR>0.5)*' % iso
    doCut += cut

    # Gen matching additions
    if 'RealTau' in run : doCut += '*(t_gen_match == 5)'
    if 'FakeTau' in run : 
        doCut = '(IsoMu22 == 1)'
        doCut += '*(tMVAIsoTight==1)*(mTrigMatch>0)*(leptonDR>0.5)*'
        doCut += cut
        doCut += '*(t_gen_match != 5)'


    t.Draw( var+' >> '+name, doCut )
    print name, h.Integral()
    h.GetXaxis().SetTitle('#tau p_{T} (GeV)')
    h.GetYaxis().SetTitle('Number of Events')
    h.SetDirectory( 0 )
    return h

def subtractTH1( h1, h2 ) :
    h3 = h1
    h3.Add( -1 * h2 )
    h3.SetTitle( h1.GetTitle()+'_Minus_'+h2.GetTitle() )
    return h3

def saveLoop( c, map ) :
    for name, h in map.iteritems() :
        h.Draw()
        c.SaveAs('/afs/cern.ch/user/t/truggles/www/TAP/'+name+'.png')
        c.Clear()

def divideTH1( h1, h2 ) :
    ### FIXME Check bins to make sure Pass <= All
    for b in range( 1, h1.GetNbinsX()+1 ) :
        b1 =  h1.GetBinContent( b )
        b2 =  h2.GetBinContent( b )
        print b1, b2
        if b1 > b2 :
            h2.SetBinContent( b, b1 )
    g = ROOT.TGraphAsymmErrors( h1, h2 )
    return g

if __name__ == '__main__' :

    #colors = [i for i in range(1, 10)]
    effPlots = {}
    #for run in ['RunB', 'RunC', 'RunD', 'RunE', 'RunF', 'AllRuns', 'ICHEPRuns'] :
    runs = ['RunB', 'RunC', 'RunE', 'RunF', 'AllRuns', 'ICHEPRuns', 'DYJets','DYJetsRealTau','DYJetsFakeTau']# 'ggH125']
    isolations = ['VLoose','Loose','Medium','Tight','VTight']
    #isolations = ['VLoose',]#'Loose','Medium','Tight','VTight']
    for iso in isolations :
        effPlots[iso] = {}
        for i, run in enumerate(runs) :
            openName = run
            if 'DYJets' in run : openName = 'DYJets'
            f = ROOT.TFile('/data/truggles/TAP_'+DATE+'_hadd/'+openName+'.root','r')
            t = f.Get('tagAndProbe/tagAndProbe/Ntuple')

            cuts = {
                'SSPass'+run: '(SS == 1 && IsoMu21MediumIsoTau32 == 1\
                    && tTrigMatch>0.5 && tL1Match>0.5)',
                'OSPass'+run: '(SS == 0 && IsoMu21MediumIsoTau32 == 1\
                    && tTrigMatch>0.5 && tL1Match>0.5)',
                #'SSFail'+run: '(SS == 1 && IsoMu21MediumIsoTau32 == 0)',
                #'OSFail'+run: '(SS == 0 && IsoMu21MediumIsoTau32 == 0)',
                'SSAll'+run: '(SS == 1)',
                'OSAll'+run: '(SS == 0)',
                }

            hists = {}
            for name, cut in cuts.iteritems() :
                print name, cut
                hists[ name ] = getHist( t, 'tPt', cut, name, run, iso )
                
            ### Save all
            c = ROOT.TCanvas('c','c',600,600)
            #saveLoop( c, hists )

            ### Do OS - SS
            #groups = ['Pass'+run,'Fail'+run,'All'+run]
            groups = ['Pass'+run,'All'+run]
            subMap = {}
            for group in groups :
                subMap[ group ] = subtractTH1( hists['OS'+group], hists['SS'+group] )
            #saveLoop( c, subMap )


            ### Make Eff Plot
            g = divideTH1( subMap['Pass'+run], subMap['All'+run] )    
            c.SetGrid()
            g.SetMaximum( 1.2 )
            g.GetXaxis().SetTitle('#tau p_{T} (GeV)')
            g.GetYaxis().SetTitle('L1 + HLT Efficiency')
            #g.SetTitle(run+' HLT MediumIso35Tau Eff. per Tau')
            g.SetTitle(run)
            #g.SetLineColor( colors[i] )
            g.SetLineWidth(2)
            g.Draw()
            c.SaveAs('/afs/cern.ch/user/t/truggles/www/TAP/FinalEff_'+run+'.png')
            c.Clear()
            effPlots[iso][run] = g

        del c
        c = ROOT.TCanvas('c','c',900,900)
        c.SetGrid()
        effPlots[iso]['AllRuns'].SetMaximum(1.5)
        effPlots[iso]['AllRuns'].Draw()
        #finalRuns = ['ggH125', 'AllRuns', 'ICHEPRuns', 'DYJets']
        #finalRuns = ['AllRuns', 'ICHEPRuns',]# 'DYJets', 'DYJetsRealTau']
        finalRuns = ['AllRuns', 'DYJets', 'DYJetsRealTau']
        #colors = [ROOT.kBlack, ROOT.kGray, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+1, ROOT.kYellow-2]
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+1, ROOT.kYellow-2]
        for i, run in enumerate(finalRuns) :
            print i, run
            effPlots[iso][run].SetLineColor( colors[i] )
            if run == 'AllRuns' : continue
            effPlots[iso][run].Draw('SAME')
        oldEff = getCBEffGraph( iso )
        oldEff.SetLineColor(ROOT.kBlue)
        oldEff.Draw('SAME')
        #bestFit = fitArcTan( effPlots[iso]['AllRuns'] )
        #bestFit.SetLineColor(ROOT.kCyan)
        #bestFit.Draw('SAME')
        legItems = [effPlots[iso]['AllRuns'],
            #effPlots[iso]['ICHEPRuns'],
            oldEff,#]
            effPlots[iso]['DYJets'],
            effPlots[iso]['DYJetsRealTau'],]
        legNames = ['20/fb - Wisc.',
            #'ICHEP - Wisc.',
            'ICHEP - Tau POG',#]
            'DYJets - All Taus',
            'DYJets - Real Taus',]
        leg = buildLegend( legItems, legNames )
        leg.Draw()
        decorate(20.0)
        #c.BuildLegend()
        c.SaveAs('/afs/cern.ch/user/t/truggles/www/TAP/_Combined_FinalEff_%s.png' % iso)
        c.SaveAs('/afs/cern.ch/user/t/truggles/www/TAP/_Combined_FinalEff_%s.pdf' % iso)
    








