import numpy as np
import pylab as plt
from scipy.optimize import root_scalar

mmin = 0.03
mmax = 120
p1 = 0.3
p2 = 1.3
p3 = 2.3
break1=0.08
break2=0.5

slopes = [0.3, 1.3, 2.3]
breaks = [mmin, break1, break2, mmax]
normfactor = 1

def BrokenPowerLaw(slopes, breaks, p):
    
    
    def PowerLaw(slope, m1, m2):
        """ Power law with slope slope in the interval m1,m2 """
        slope = slope
        m1 = float(m1)
        m2 = float(m2)
        assert (m1 < m2)
        assert (m1 > 0)
        assert (m1 != -1)
        
       	return slope,m1,m2
    
    def calcpows(slopes,breaks):
        if not (len(slopes) == len(breaks) - 1):
            raise ValueError(
                'The length of array of slopes must be equal to length of ' +
                'array of break points minus 1')
        if not ((np.diff(breaks) > 0).all()):
            raise ValueError('Power law break-points must be monotonic')
        nsegm = len(slopes)
        pows =[]
        print(slopes,breaks) 
        for ii in range(nsegm):
            pows.append(
                PowerLaw(slopes[ii], breaks[ii],
                         breaks[ii + 1]))
        return pows, nsegm
        
    pows, nsegm = calcpows(slopes,breaks)
    print(pows)
    
    def pdf(pows_ii, breaks_ii, nsegm, weights, pows, breaks):
        x1 = np.asarray(breaks_ii)
        ret = np.atleast_1d(x1) * 0.
        for ii in range(nsegm):
            xind = (x1 < breaks[ii + 1]) & (x1 >= breaks[ii])
            if xind.sum() > 0:
                ret[xind] = weights[ii] * pdf(pows[ii], x1[xind], nsegm, weights, pows, breaks)
        return ret.reshape(x1.shape)
    
    
    def calcweights(slopes,breaks,pows):
        nsegm = len(slopes)
        weights =[1]
        for ii in range(1, nsegm):
            rat = pdf(pows[ii], breaks[ii], nsegm, weights, pows, breaks) / pdf(
                pows[ii - 1], breaks[ii], nsegm, weights, pows, breaks)
            weights.append(weights[-1] / rat)
        weights = np.array(weights)
        weights = weights / np.sum(weights)  # relative normalizations
        
        return weights
    
    weights = calcweights(slopes, breaks, pows)
    #print(p,weights,nsegm,pows)
    def ppf(x0, weights, nsegm, pows):
        x = np.asarray(x0)
        x1 = np.atleast_1d(x)
        edges = np.r_[[0], np.cumsum(weights)]
        # edges of powerlaw in CDF scale from 0 to 1
        pos = np.digitize(x1, edges)  # bin positions, 1 is the leftmost
        pos = np.clip(pos, 1, nsegm)
        #  we can get zeros here if input is corrupt
        left = edges[pos - 1]
        w = weights[pos - 1]
        x2 = np.clip((x1 - left) / w, 0, 1)  # mapping to 0,1 on the segment

        # must force float b/c int dtypes can result in truncation
        ret = np.zeros_like(x1, dtype='float')
        for ii in range(x.size):
            ret[ii] = pows[pos[ii] - 1].ppf(x2[ii])

        isnan = (x1 < 0) | (x1 > 1)
        if any(isnan):
            ret[isnan] = np.nan
        return ret.reshape(x.shape)
    
    return ppf(p, weights, nsegm, pows)
    
def inverse_imf(p,
                mmin=None,
                mmax=None,
                massfunc='kroupa',
                **kwargs):
    """
    Inverse mass function.  Given a likelihood value in the range [0, 1),
    return the appropriate mass.  This just calls the mass function's ppdf
under the hood.


    Parameters
    ----------
    p: np.array
        An array of floats in the range [0, 1).  These should be uniformly random
        numbers.
    mmin: float
    mmax: float
        Minimum and maximum stellar mass in the distribution
    massfunc: string or function
        massfunc can be 'kroupa', 'chabrier', 'salpeter', 'schechter', or a
        function
    """


    # this should be the entirety of "inverse-imf".  The rest is a hack
    
    return BrokenPowerLaw(slopes, breaks, p)


def get_massfunc_name(massfunc):
    if massfunc in reverse_mf_dict:
        return reverse_mf_dict[massfunc]
    elif type(massfunc) is str:
        return massfunc
    elif hasattr(massfunc, '__name__'):
        return massfunc.__name__
    else:
        raise ValueError("invalid mass function")


mass = inverse_imf(np.random.rand(1000000), mmin, mmax , massfunc= 'kroupa')


plt.hist(mass, bins=np.geomspace(mmin,mmax))
plt.xscale('log')
plt.yscale('log')
plt.show()
