def enoise(data,navg):
    """    [noise,stdv,i_max,index]=enoise3(data,navg)
    Routine to estimate noise level based on Hildebrand and Sekhon
    Modified by Pablo Reyes 15apr2009:
     - count=min(npts,max(int(npts/50.),10))# in case few data points
    Modified by Pablo Reyes 10oct2011:
     - use : range(count, npts) instead of np.arange(count,npts)
     - use : "n" instead of "id" since "id" is a python built-in function
    """
    index = data.argsort()
    data = data[index]
    npts = data.size    
    count=min(npts,max(int(npts/50.),10))
    #count=min(10.,npts/2)
    s = data[0:count].sum()
    s2 = (data[0:count] ** 2.).sum()
    for n in range(count, npts): #last id=npts-1
#    for id in np.arange(count, npts): #last id=npts-1
        s = s + data[n]
        s2 = s2 + data[n] ** 2.
        p = s / (n + 1.)
        dp2 = (s2 - (n + 1.) * p ** 2.) / n
        if (p) ** 2./navg < dp2:
            noise = data[0:n].sum() / n
            stdv = dp2 ** 0.5
            i_max = n
            index = index[0:i_max]
            return noise, stdv, i_max, index
    noise = data.sum() / npts
    stdv = dp2 ** 0.5
    stdv = dp2
    i_max = npts
    return noise, stdv, i_max, index
