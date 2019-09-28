import sys,optparse
from math import log,exp,floor # natural log
from datetime import datetime

def getAd(alpha):
    return(-log(alpha))

''' Read genome-wide substitution:polymorphism ratio
#format: N sub poly; no header
#0 for sub, 1 for poly
'''
def getPolyDen(spectfile):
    g=[]; G=[]
    with open(spectfile,'r') as spect:
        l=spect.next()
        l=l.strip().split('\t')
        N=int(l[0]); s=float(l[1]); p=float(l[2])
        print('Substitutions: %s ; polymorphisms: %s' %(s,p))
    try:
        assert s+p == 1
    except:
        print('s+p= %s + %s = %s != 1' % (s,p,s+p))
        sys.exit()
    g={(0,N):s,(1,N):p}
    G={(0,N):log(s), (1,N):log(p)}
    return(g,G,N)

'''Read unormalized neutral spect
##format: x n freq; no header
'''
def getLogSpect(spectfile,nosub,MAF):
    g={}; G={}; N=[]; checksum=0.
    with open(spectfile,'r') as spect:
        for l in spect:
            l=l.strip().split("\t")
            x=int(l[0]); n=int(l[1])
            f=float(l[2])
            g[(x,n)]=f
            G[(x,n)]=log(f)
            checksum += f
            N.append(n)
    N = list(set(N))
    print('Total prob is %f in the spectrum file.' %(checksum))
    if len(N) >= 2:
        print("Current implementation only supports uniform sample size. Please try again.")
        sys.exit()
    else:
        N=N[0]
    #Make sure there's no missing keys:
    fold=False
    for k in range(N+1):
        if (k,N) not in g:
            g[(k,N)] = 0
        elif MAF and k > N/2:
            fold = True
    #if fold;     
    if MAF and fold:
        fg = {}; fG = {}; checksum=0
        for i in range(N/2):
            if nosub and i==0:
                print('Skip substitutions')
                continue
            fg[(i,N)] = g[(i,N)] + g[(N-i,N)]
            fG[(i,N)] = log(fg[(i,N)])
            checksum += fg[(i,N)]
        if N%2 == 0:
            fg[(N/2,N)] = g[(N/2,N)]
            fG[(N/2,N)] = log(fg[(N/2,N)])
            checksum += fg[(N/2,N)]
        print('After folding, total probability is %s.'%(checksum))
        return(fg,fG,N)
    #if not folded:
    return(g,G,N) #spect, logspect, size

''' Precompute binomial coefficients for sample size N
#log of Binom(n choose k)
# B = N!/k!(N-k)!
# logB = sum_{n-k+1}^n{logi} - sum_1^{n-k}{logi}
'''
def logBinCoeff(N):
    logIs = [log(i) for i in xrange(1,N+1)]
    logBs={}
    logBs[0]=0 ; logBs[N]=0
    for k in range(1,N):
        logBs[k] = sum(logIs[N-k:])-sum(logIs[:k])
    return logBs

''' Precompute binomial probability for sample size N, chance x for B2, B2maf or B1
# Ancestral sites (k=0) not considered. Normalize probability across k=1...N
# For each x, return a dash table of k:prob[k]
'''
def getBinom(N,x,logBs):
    #logIs = [log(i) for i in xrange(1,N+1)]
    logProb = {}; Prob={}
    logP0 = N*log(1-x)
    #logPn = N*log(x)
    norm = log( 1-exp(logP0) )
    #print('log(1-(1-x)^N) = %e'%(norm))
    for k in xrange(N+1):
        logB = logBs[k]
        logP = logB + k*log(x) + (N-k)*log(1-x)
        logProb[k] = logP - norm
        Prob[k] = exp(logProb[k])
    return(Prob)#,logProb

''' Precompute binomial probability distribution for B0.
# Monomorphic sites (k=0 || k=N) not considered. Normalize probability across k=1...N-1 '''
def getNonZeroBinom(N,x,logBs):
    #logIs = [log(i) for i in xrange(1,N+1)]
    logProb = {}; Prob={}
    logP0 = N*log(1-x)
    logPn = N*log(x)
    norm = log(1-exp(logPn)-exp(logP0))
    #print('log(1 - x^N - (1-x)^N) = %e'%(norm))
    for k in xrange(N+1):
        logB = logBs[k]
        logP = logB + k*log(x) + (N-k)*log(1-x)
        logProb[k] = logP - norm
        Prob[k] = exp(logProb[k])
    return(Prob)#,logProb

''' Return folded probability distribution for B0maf or B2maf'''
def foldBinom(Prob,N,nosub):
    foldedProb={}
    for k in range(N/2):
        if nosub and k==0:
            #print('Skip substitutions')
            continue
        if k not in Prob:
            Prob[k] = 0
            #print('Assign f[%s] = 0'%(k))
        if N-k not in Prob:
            Prob[N-k] = 0
            #print('Assign f[%s] = 0'%(N-k))
        foldedProb[k] = Prob[k] + Prob[N-k]
    if N%2 == 0:
        foldedProb[N/2] = Prob[N/2]
    return(foldedProb)

''' Generate config_alt for B1
# Normalized by 1-h(0)'''
def getSubLogConfig(N,x):#,logBs
    #no need to get full prob; just get h(N) and h(0)
    logPn= N*log(x) ; logP0 = N*log(1-x)

    norm = 1 - exp(logP0)
    p=exp(logPn)/norm
    Config=[p,1-p]
    return(Config)#,logConfig

''' Return log likelihood for the given window
# testSite is genetic position (in recomb unit)
# logP_neut = logSpect[(k,N)]; logP_alt = Fx[x][a][k]'''
def calcBaller(pos_list,Ks,Ns, testSite, testSite_i, logSpect, biFx, xGrid, AGrid, AdGrid, aGrid,MAF):
    #note that the test site is included in Ks
    #Optimize over xGrid and AGrid
    L = len(AdGrid) ; numSites = len(pos_list)
    #testSite = pos_list[testSite_i]
    biTmax = [-100,0,0] #LR, x, A
    # Go through all A values
    for A in AGrid:
        # i indexes positions. pos_list ascending
        # c indexes Ad value in the pre-computed descendingly sorted list
        i_index=[] ; c_index=[]
        # going through sites from the center to outward
        i = max(testSite_i-1,0) ; c = len(AdGrid)-1
        # Leftward:
        while i >= 0 and c > 0 and testSite_i != 0:
            pos = pos_list[i]
            dist = testSite - pos
            Ad = A*dist
            while AdGrid[c] < Ad and c >= 0:
                c -= 1
            #now AdGrid[c] >= Ad or c=0
            if c >=0:
                i_index.append(i)
                c_index.append(c)
            i -= 1
        # now i==0. Starting rightward from center
        #i=min(testSite_i+1,pos_list[-1]) 
        i = min(testSite_i+1,numSites) ; c = len(AdGrid)-1
        while i < len(pos_list) and c > 0:
            pos = pos_list[i]
            dist = pos - testSite
            Ad = A*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c >= 0:
                i_index.append(i)
                c_index.append(c)
            i += 1
        # if noCenter, consider testSite_i too
        if testSite != pos_list[testSite_i]:
            pos = pos_list[testSite_i]
            dist = abs(pos - testSite)
            Ad = A*dist
            while AdGrid[c] < Ad and c >=0:
                c -= 1
            if c>=0:
                i_index.append(testSite_i)
                c_index.append(c)
        #go through the grid of x
        #note that this grid only spans 0-0.5 
        for x in xGrid:
            La1=0; La2=0; La_bi=0; L0=0
            #go through all the sites
            for j in range(len(i_index)):
                i = i_index[j] #index in the pos_list
                c = c_index[j]
                alpha = aGrid[c]
                k = Ks[i] ; N = Ns[i]
                La_bi += biFx[x][alpha][k]
                L0 += logSpect[(k,N)]
            #get LR
            Tbi = 2*(La_bi - L0)
            if Tbi >=biTmax[0]:
                biTmax = [Tbi,x,A]
    return(biTmax)

''' Generate the grids of x and A values to optimize over. The grid of alpha is fixed.'''
def getGrids(x,MAF, seqA, listA):
    #define the default list of alpha values:
    aGrid = [0,1e-8,1e-7,1e-6,1e-5]+[m*1e-4 for m in range(1,10001)]
    #generate sorted AdGrid based on the dense grid of alpha:
    AdGrid = [-log(1e-32)]+[ getAd(a) for a in aGrid[1:] ] 
    assert len(aGrid) == len(AdGrid)
    #get xGrid
    if x:
        xGrid = [float(x)]
    else:
        xGrid = [.05*i for i in range(1,11)]
    #print(xGrid)
    if not seqA and not listA: 
        #AGrid = [500*i for i in xrange(2,21)]
        AGrid = [100*i for i in range(1,12)] + [200*i for i in range(6,13)] + [500*i for i in range(5,10)] + [1000*i for i in range(5,11)] + [1e6,1e8]
        #The grid is: 100-1100, 1200-2400, 2500-4000, 5000-9000, 1e4,1e6,1e8
        #total of 29 values
    elif listA:
        AGrid = [float(x) for x in listA.split(',')]
    else:
        Amin,Amax,Astep = [float(x) for x in seqA.split(',')]
        n=(Amax-Amin)/Astep
        AGrid=[Amin+Atep*i for i in range(n+1)]
    return(xGrid,AGrid,AdGrid,aGrid)

''' Generate the look-up table for g(k) and f(k,d; x,A) distributions '''
#By default, the bin width of alpha values is 1e-4
def initialize(spectfile,xGrid,aGrid,MAF,nofreq,nosub):
    Prob1={}; Fx={}; biFx={}
    if nofreq: #B1
        print('Only consider polymorphism density...')
        Spec,logSpec,N = getPolyDen(spectfile)
        print('Sample size: %s'%(N))
        #0 for subs, 1 for poly
        for x in xGrid:
            biFx[x]={}
            Prob1[x] = getSubLogConfig(N,x) #a list of [s,p]
            Prob1[1-x] = getSubLogConfig(N,1-x)
            for a in aGrid:
                biFx[x][a]={}
                for k in [0,1]:
                    gx = Spec[(k,N)]
                    bifx = a*(.5*Prob1[x][k]+.5*Prob1[1-x][k]) + (1-a)*gx
                    biFx[x][a][k] = log(bifx)
    else: 
        Spec,logSpec,N = getLogSpect(spectfile,nosub,MAF)
        print('Sample size: %s'%(N))
        logBs = logBinCoeff(N)
        if nosub: #B0s
            print('Only consider polymorphisms in the input...')
            for x in xGrid:
                biFx[x]={}
                Prob1[x] = getNonZeroBinom(N,x,logBs)
                Prob1[1-x] = getNonZeroBinom(N,1-x,logBs)
                if MAF:
                    Prob1[x] = foldBinom(Prob1[x],N,nosub)
                    Prob1[1-x] = foldBinom(Prob1[1-x],N,nosub)
                for a in aGrid:
                    biFx[x][a]={}
                    if MAF: #B0maf
                        for k in range(1,N/2 +1):
                            gx = Spec[(k,N)]
                            bifx = a*(.5*Prob1[x][k]+.5*Prob1[1-x][k]) + (1-a)*gx
                            biFx[x][a][k] = log(bifx)
                    else: #B0
                        for k in range(1,N):
                            gx = Spec[(k,N)]
                            bifx = a*(.5*Prob1[x][k]+.5*Prob1[1-x][k]) + (1-a)*gx
                            biFx[x][a][k] = log(bifx)
        else: #B2s
            for x in xGrid:
                biFx[x]={}
                Prob1[x] = getBinom(N,x,logBs)
                Prob1[1-x] = getBinom(N,1-x,logBs)
                if MAF:
                    Prob1[x] = foldBinom(Prob1[x],N,nosub)
                    Prob1[1-x] = foldBinom(Prob1[1-x],N,nosub)
                for a in aGrid:
                    biFx[x][a]={}
                    if MAF: #B2maf
                        for k in range(N/2 +1):
                            gx = Spec[(k,N)]
                            bifx = a*(.5*Prob1[x][k]+.5*Prob1[1-x][k]) + (1-a)*gx
                            biFx[x][a][k] = log(bifx)
                    else: #B2
                        for k in range(1,N+1):
                            gx = Spec[(k,N)]
                            bifx = a*(.5*Prob1[x][k]+.5*Prob1[1-x][k]) + (1-a)*gx
                            biFx[x][a][k] = log(bifx)
        if not MAF:
            print('Using derived allele frequency...')
        else:
            print('Using minor allele frequency...')
    print('%s Initialization finished.\n' % (datetime.now()))
    return(Spec,logSpec,biFx,N)


''' Return parsed data
# all sitePos are genetic positions (by cM)
# default value of Rrate is 1e-6 cM/site
# Input format: physPos, genPos, x, n
'''
def readInput(infile,nosub,nofreq,MAF,phys=False,Rrate=1e-6):
    phys_pos = []
    Ks=[]; Ns = []; pos_list=[]; numSites=0 
    postype=1-int(phys) #index; 0 is physical position, 1 is genetic position
    with open(infile,'r') as sites:
        l=sites.next()#skipping header
        if nosub:
            #make sure to only read frequencies
            for l in sites:
                l=l.strip().split('\t')
                if float(l[2])/float(l[3]) not in [0.,1.]:
                    numSites += 1
                    physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                    #genPos = float(l[1])
                    physPos = int(float(l[0]))
                    sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                    k = MAF*min(k,n-k) + (1-MAF)*k
                    Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                    phys_pos.append(physPos)
        elif nofreq: # Take all none N cases as 1
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = (k!=n) #k=1 if not sub, k=0 if sub
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
        else:
            for l in sites:
                l=l.strip().split('\t')
                numSites += 1
                physPos,k,n = [float(l[0]),int(l[2]),int(l[3])]
                #genPos = float(l[1])
                physPos = int(float(l[0]))
                sitepos = float(l[postype])*phys*Rrate + float(l[postype])*(1-int(phys))
                k = MAF*min(k,n-k) + (1-MAF)*k
                Ks.append(k); Ns.append(n); pos_list.append(sitepos)
                phys_pos.append(physPos)
    return(phys_pos,pos_list,Ks,Ns,numSites)


def scan(xGrid,AGrid,AdGrid,aGrid,outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec,biFx,N,size=False,r=0,s=1,phys=False,nofreq=False,MAF=False,noCenter=False,Rrate=1e-6):
    print("writing output to %s" % (outfile))
    with open(outfile,'w') as scores:
        scores.write('physPos\tgenPos\tLR\txhat\tAhat\n')#\tLRpos1\txhat_p1\tAhat_p1\tLRpos2\txhat_p2\tAhat_p2
        #fixed distance
        if size:
            print('You\'ve chosen to fix the scan window size.')
            if r==0:
                print('Please set a window width in bp with \"-w\" or \"--window\" command.')
                sys.exit()
            if not phys:
                print('Please make sure to use physical positions as coordinates if fixed-length windows are chosen. Scan will continue with physical positions.')
                phys = True
                #sys.exit()
            if noCenter:
                w = float(r) 
                print('Computing LR on %.3f kb windows on every %s bp. Using physical positions by default.' % (w/1e3, s))
                start = int( floor(2*float(phys_pos[0])/w)*(w/2) )
                end = start + s ; midpos = start + s/2
                start_i=0; end_i=0; pos_i=0
                while midpos <= phys_pos[-1]:
                    #define the window
                    while phys_pos[start_i] < start:
                        start_i +=1 
                    while phys_pos[pos_i] < pos_i:
                        pos_i += 1
                    while end_i < numSites:
                        if phys_pos[end_i] < end:
                            end_i += 1
                        else:
                            break
                    if start_i >= end_i:
                        scores.write('%g\tNA\t0\tNA\tNA\n' % (midpos))
                        start+=s; midpos+=s; end+=s
                    else:
                        Window = Ks[start_i:end_i+1]
                        #calcBaller args: Window, N, testSite, logSpect, biFx, xGrid, AGrid, AdGrid, aGrid, all=False, phys=False, nofreq=False, MAF=False)
                        #pos_list and phys_pos has matching indice
                        Tmax = calcBaller(pos_list[start_i:end_i+1], Ks[start_i:end_i+1], Ns[start_i:end_i+1], midpos*Rrate, pos_i, logSpec, biFx, xGrid, AGrid, AdGrid, aGrid, MAF)
                        scores.write('%g\tNA\t%s\t%.2f\t%g\n' % (midpos, Tmax[0], Tmax[1], Tmax[2]))#
                        #read in, take next step
                        start+=s; midpos+=s; end+=s 
            else: #site-centered
                w = float(r) 
                print('Computing LR on %.2f kb windows on every %g informative sites. Using physical positions by default.' % (w/1e3,s))
                i=0; #testSite=phys_pos[i]
                #start =  max(0, testSite-r/2); end = min(testSite+r/2,phys_pos[-1])
                start_i=0; end_i=0
                while i < numSites:
                    testSite=phys_pos[i]
                    start =  max(0, testSite-r/2); end = min(testSite+r/2,phys_pos[-1])
                    while phys_pos[start_i] < start:
                        start_i += 1
                    while end_i < numSites:
                        if phys_pos[end_i] < end:
                            end_i += 1
                        else:
                            break
                    assert end_i >= start_i
                    end_i = min(end_i, numSites-1)
                    #Window = Ks[start_i:end_i+1]
                    Tmax = calcBaller(pos_list[start_i:end_i+1], Ks[start_i:end_i+1], Ns[start_i:end_i+1], pos_list[i], i, logSpec, biFx, xGrid, AGrid, AdGrid, aGrid, MAF)
                    scores.write('%s\t%s\t%s\t%.2f\t%g\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2]))#
                    i+=int(s)
        # When fixSize == False, and radius (-r) provided. Scan with fixed number of sites
        elif r != 0: 
            print('Computing LR on every %s site/s, with %s sites on either side.' % (s, r))
            i=0
            while i < numSites:
                testSite = pos_list[i]
                start_i = max(0, i-r) ; end_i = min(numSites,i+r+1)
                Tmax = calcBaller(pos_list[start_i:end_i], Ks[start_i:end_i], Ns[start_i:end_i], testSite, i, logSpec, biFx, xGrid, AGrid, AdGrid, aGrid, MAF)
                scores.write('%s\t%s\t%s\t%.2f\t%g\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2]))
                i+=s
        #window size not given 
        #then use all data (but the test site)
        else:
            print('Computing LR on every %s site/s, using all the data with alpha >= 1e-8.' % (s))
            i=0
            while i < numSites:
                testSite = pos_list[int(i)]
                Tmax = calcBaller(pos_list, Ks, Ns, testSite, int(i), logSpec, biFx, xGrid, AGrid, AdGrid, aGrid, MAF)
                scores.write('%s\t%s\t%s\t%.2f\t%g\n' % (phys_pos[i], pos_list[i], Tmax[0], Tmax[1], Tmax[2]))
                i+=int(s)
        scores.close()

#give the config file given the concatenated input
def getConfig(infile,configfile):
    Config={}; numSites=0# N: [s,p]
    with open(infile,'r') as sites:
        l=next(sites)#skip the header by default
        for l in sites:
            x,n = [int(x) for x in l.strip().split('\t')[2:] ]
            if x==0:
                print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) will be ignored.\n')
                continue
            if n not in Config:
                Config[n] = [0,0]
            Config[n][0] += int(x==n)
            Config[n][1] += 1-int(x==n)
            numSites+=1
    sizes = sorted(Config.keys())
    with open(configfile,'w') as config:
        for N in sizes:
            config.write('%s\t%s\t%s\n' % ( N, Config[N][0]/float(numSites) , Config[N][1]/float(numSites) ))
    sites.close(); config.close()
    print('Done')

#give spectrum file given the concatenated input
def getSpect(infile,spectfile,nosub=False,MAF=False):
    Spect={}; numSites=0
    if nosub:
        print('Generating spectrum for all polymorphic sites. Substitutions (x=n or x=0) won\'t be considered.')
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                if x == 0 or x == n:
                    continue
                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    else:
        with open(infile,'r') as sites:
            l=next(sites)
            for l in sites:
                (x,n)=[ int(i) for i in l.strip().split('\t')[2:] ]
                if MAF:
                    x = min(x, n-x)
                elif x==0:
                    print('Please make sure the input has derived allele frequency. Sites with 0 observed allele count (k=0) should not be included.\n')
                    sys.exit()

                if (x,n) in Spect:
                    Spect[(x,n)] += 1
                else:
                    Spect[(x,n)] = 1
                numSites+=1
    #write out
    pairs = sorted(Spect.keys())
    with open(spectfile,'w') as spec:
        for x,n in pairs:
            spec.write('%s\t%s\t%s\n' % (x,n,float(Spect[(x,n)])/float(numSites) ))
    sites.close(); spec.close()
    print('Done.')

#main function to scan through the input
#only work with minor allele frequency
def main(): 
    #parsing arguments
    parser = optparse.OptionParser(usage='python %prog -i <input file> -o <output file> --spect <spect/config file> [--help] [--nofreq] [--nosub] [--MAF] [--getSpect] [--getConfig] [--fixSize] [--physPos] [--rec <recomb rate>] [-w <window size>] [--noCenter] [-s <step size>] [--rangeA <min,max,step>] [--listA <A1,A2,..,Ak>] [--fixX <x>]')
    parser.add_option('-i','--input', dest='infile', help = 'Path and name of your input file.\n')
    parser.add_option('-o','--output', dest='outfile', help = 'Path and name of your output file.\n')
    parser.add_option('--spect',dest='spectfile', help = 'Path and name of the allele frequency spectrum file or configuration file.\n')
    parser.add_option('--getSpect', dest='getSpec', action='store_true', default=False, help='Option to generate frequency spectrum file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively. Indicate the input type with \"--MAF\" and/or \"--nosub\".\n')
    parser.add_option('--getConfig', dest='getConfig', action='store_true', default=False, help='Option to generate configuration file from the concatenated input file. Use \"-i\" and \"--spect\" commands to provide names and paths to input and output files, respectively.\n\n')

    parser.add_option('--nofreq', dest='nofreq', action='store_true', default=False, help = 'Option to ignore allele frequency information. All polymorphic sites will be considered as equivalent.')
    parser.add_option('--nosub', dest='nosub', action='store_true', default=False, help = 'Option to not include substitution in input data.')
    parser.add_option('--MAF', dest='MAF', action='store_true', default=False, help = 'Option to use minor allele frequency, instead of polarized allele frequency. The latter is default.')
    parser.add_option('--physPos', action='store_true', dest = 'phys', default = False, help = 'Option to use physical positions instead of genetic positions (in cM). Default is using genetic positions.\n')
    parser.add_option('--rec', dest='Rrate', default = 1e-6 , help='The uniform recombination rate in cM/nt. Default value is 1e-6 cM/nt. Only useful when choose to use physical positions as coordinates.\n\n')

    parser.add_option('--fixSize', action='store_true', dest = 'size', default = False, help = 'Option to fix the size of scanning windows. When true, provide the length of window in neucleotide (nt) with \"-w\" or \"--window\" command.\n')
    parser.add_option('-w','--window', dest='r', type = 'int', default=0, help='Number of sites flanking the test locus on either side. When choose to fix window size (\"--fixSize\"), input the length of window in bp.\n')
    parser.add_option('--noCenter', action='store_true', dest='noCenter', default=False, help = 'Option to have the scanning windows not centered on informative sites. Require that the window size (\"-w\") in physical positions (\"--physPos\") is provided. Default is True.\n')
    parser.add_option('-s','--step', dest='step', type = 'float', default=1, help='Step size in bp (when using \"--noCenter\") or the number of informative sites. Default value is one site or one nucleotide.\n\n')

    parser.add_option('--fixX', dest='x', help='Option to fix the presumed equilibrium frequency.\n')
    parser.add_option('--rangeA', dest='seqA', help='Range of the values of the parameter A to optimize over. Format should follow <Amin>,<Amax>,<Astep> with no space around commas.\n')
    parser.add_option('--listA', dest='listA', help='Manually provide a list of A values to optimize over. Please separate the values with comma, no space.\n')

    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    opt,v = parser.parse_args(sys.argv[1:])

    if opt.getSpec:
        print('You\'ve chosen to generate site frequency spectrum...')
        print('Concatenated input: %s \nSpectrum file: %s' % (opt.infile, opt.spectfile ))
        getSpect(opt.infile, opt.spectfile, opt.nosub, opt.MAF)
        sys.exit()
    elif opt.getConfig:
        print('You\'ve chosen to generate the substitution-polymorphism configuration...')
        print('Concatenated input: %s \nConfiguration file: %s' % (opt.infile, opt.spectfile ))
        getConfig(opt.infile, opt.spectfile)
        sys.exit()

    #generate the grids to optimize over/with
    xGrid,AGrid,AdGrid,aGrid = getGrids(opt.x, opt.MAF, opt.seqA, opt.listA)
    #print xGrid
    print('\nOptimizing over x= '+', '.join([str(x) for x in xGrid]))
    print('Optimizing over A= '+', '.join([str(x) for x in AGrid]))

    #initialization for a grid of x
    print('\n%s. Initializing...'%(datetime.now()))
    (Spec, logSpec, biFx,N) = initialize(opt.spectfile,xGrid,aGrid,opt.MAF,opt.nofreq,opt.nosub)#initialize(spectfile,xGrid,aGrid,MAF,nofreq,nosub)

    #start reading data
    print("\n%s. Reading input file: %s" % (datetime.now(),opt.infile) )
    (phys_pos,pos_list,Ks,Ns,numSites) = readInput(opt.infile,opt.nosub,opt.nofreq,opt.MAF,opt.phys,opt.Rrate)

    #finished reading file, start writing output
    print("\n%s. Start computing likelihood raito..." %(datetime.now()))
    scan(xGrid,AGrid,AdGrid,aGrid,opt.outfile,phys_pos,pos_list,Ks,Ns,numSites,Spec, logSpec, biFx,N, opt.size,opt.r,opt.step,opt.phys,opt.nofreq,opt.MAF,opt.noCenter,opt.Rrate)#, ProbsStar

    #pipeline finished
    print('\n%s. Pipeline finished.'%(datetime.now()))



if __name__ == '__main__':
    main()
