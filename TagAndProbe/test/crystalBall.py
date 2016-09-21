import ROOT
import math




# Directly from Riccardo
def CBeff(x, m0, sigma, alpha, n, norm):
    sqrtPiOver2 = math.sqrt(ROOT.TMath.PiOver2())
    sqrt2       = math.sqrt(2.)
    sig         = abs(sigma)
    t           = (x - m0)/sig * alpha / abs(alpha)
    absAlpha    = abs(alpha/sig)
    a           = ROOT.TMath.Power(n/absAlpha, n) * ROOT.TMath.Exp(-0.5 * absAlpha * absAlpha)
    b           = absAlpha - n/absAlpha
    arg         = absAlpha / sqrt2;
    if   arg >  5.: ApproxErf =  1.
    elif arg < -5.: ApproxErf = -1.
    else          : ApproxErf = ROOT.TMath.Erf(arg)
    leftArea    = (1. + ApproxErf) * sqrtPiOver2
    rightArea   = ( a * 1./ROOT.TMath.Power(absAlpha-b, n-1) ) / (n-1)
    area        = leftArea + rightArea
    if t <= absAlpha:
        arg = t / sqrt2
        if   arg >  5.: ApproxErf =  1.
        elif arg < -5.: ApproxErf = -1.
        else          : ApproxErf = ROOT.TMath.Erf(arg)
        return norm * (1. + ApproxErf) * sqrtPiOver2 / area
    else:
        return norm * (leftArea + a * (1./ROOT.TMath.Power(t-b,n-1.) - \
1./ROOT.TMath.Power(absAlpha - b,n-1.)) / (1.-n)) / area
