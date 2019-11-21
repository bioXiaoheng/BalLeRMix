#Nov 21, 2019: Modified to generalize the path names.
#Jan 9, 2018: Modified to parse all MS-formatted slim output.
#Nov 10, 2017: Modified script to run more simulations on cluster

import sys
from datetime import datetime

###parameters list as below:
#"chromosome" length, 50kb
L = 50000 
# 1e-6 centimorgan/nt
recDist = 1e-6 


#get allsites into each pop dictionary with sorted keys
def even(dictA, allsites): 
    for k in allsites:
        if k not in dictA:
            dictA[k] = 0
    return dictA

#polarize the freq based on Gorilla's sequence
def polarize(p1,p,N,drv=True):
    for k in p1:
        if p1[k] == int(drv):
            p[k] = N - p[k]
    return p

#read the slim output and return dictionaries of all populations
#N is the list of sample sizes, by default it's of length 3
def readslim(L, N, Hfile):
    slim = open(Hfile,'r')
    Pops = []; allpos=[]
    l=slim.readline()
    #Start recording counts to each population
    for samples in N:
        while l != '//\n':
            l=slim.readline()
        l=slim.readline()
        nsites = int(l.strip().split(' ')[1])
        pos = [ float(x)*L for x in slim.readline().strip().split(' ')[1:] ]
        pos = [ int(round(x)) for x in pos]
        allpos += pos
        dup=[];  Pop={}; counts = [0]*nsites
        #read first lineage
        l = slim.readline().strip()
        #record sampled lineages in this population
        #for j in range(samples):
        j=0
        while l[0] in ['0','1']:
            i=0
            for s in l:
                counts[i] += int(s)
                #Pop[pos[i]] += int(s)
                i += 1
            try:
                assert i == nsites
            except(AssertionError):
                print i, nsites
                sys.exit()
            l = slim.readline().strip(); j+=1
            if l == '': break#record dups
        try:
            assert j == samples
        except(AssertionError):
            print j, samples,
            sys.exit()
        #get dups
        for p in range(nsites):
            if pos[p] not in Pop and counts[p] != 0:
                Pop[pos[p]] = counts[p]
            else:
                dup.append(pos[p])
        #Count the number of sites for this pop
        print 'Within-pop segsites #:', (nsites-counts.count(0)), 'nsites=', nsites,''
        Pops.append(Pop)
    allpos = list(set(allpos))
    #remove dups 
    dup = list(set(dup))
    for d in dup:
        allpos = [p for p in allpos if p != d]
        for i in range(len(N)):
            if d in Pops[i]:
                del Pops[i][d]
    allpos.sort()
    return Pops,allpos

#Integrated parsing function: 
#Assuming uniform recombination rate
def parse(slmout, path, condition, rep, N, MAF=False, drv=True, anc=2):
    if MAF:
        polar='MAF'
    elif drv:
        polar='DAF'
    else:
        polar='AAF'
    if anc==1:
        outfile = '%sinput/HC-%s_%s_%s.txt' % (path, condition, rep, polar)
        N = [N[0],1]
    elif anc==2:
        #N by default is Nh, Hc, 1
        outfile = '%sinput/HG-%s_%s_%s.txt' % (path, condition, rep, polar)
    #slmout=prefix
    Pops,allpos = readslim(L, N, slmout)

    #prep the dict
    if not MAF:
        for num in range(anc):
            Pops[num] =  even(Pops[num],allpos)
            Pops[anc] = even(Pops[anc], allpos)
            Pops[num] = polarize( Pops[anc], Pops[num], N[num], drv)
    else:
        for num in range(anc):
            Pops[num] =  even(Pops[num],allpos)
            Pops[anc] = even(Pops[anc], allpos)

    #order the sites
    allpos.sort()

    #write to output
    #only considers single-species input for now
    with open(outfile,'w') as parsed:
        parsed.write('physPos\tgenPos\tx\tn\n')
        for p in allpos:
            if Pops[0][p] != (1-drv)*50:
                ol = '%s\t%.6f\t%s\t%s\n' % (p, p*recDist, Pops[0][p], N[0])
                parsed.write(ol)
    print('End of pipeline.')

#==============================================

def main():
    infile = sys.argv[1]
    #assume $condition contains "HCG_"
    condition = sys.argv[2]
    rep = int(sys.argv[3])
    env_path = sys.argv[4]
    #N=[50,50,1]
    N = [int(x) for x in sys.argv[4].split(',')]
    path = env_path+'50kb_'+condition+'/'
    parse(infile,path,condition,rep,N, MAF=False, drv=True, anc=2)

if __name__ == '__main__':
    main()
