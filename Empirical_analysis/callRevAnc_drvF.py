#Jan 30, 2018: Modified for parsing bonobo data. Use hg38 as ancestral and panPan2 as derived
#May 24, 2018: Modified s.t. compatible with other population sizes
#May 9, 2018: Modified script to make sure to output derived freq
#go through alignment, vcf, and HWE pvalues at the same time and generate derived frequency

import gzip, sys
from datetime import datetime

Bases = ['A','T','C','G']

def findMatch(matchfile, vcffile, pvalfile, outfile):
	#sample size of %s
	#N = popSize[pop]
	N=26
	#read and store pvalues
	pval={}
	with open(pvalfile,'r') as P:
		for l in P:
			l=l.strip().split('\t')
			pos = int(l[0])
			pval[pos] = float(l[1])

	vcf = open(vcffile,'r')
	lv=vcf.readline()#header: 0pos 1id 2ref 3alt 4x 5n 6AA 7Aa 8aa
	lv=vcf.readline().strip().split('\t')
	vpos = int(lv[0]); popSize=int(lv[5])

	out = open(outfile,'w')
	out.write('posH\trefH\trefB\tpval\tx\tn\n')
	#out.write('posH\tdrv\tanc\tpval\tx\tn\n')

	with gzip.open(matchfile,'rb') as match:
		l=next(match) #header: 0chr 1posH 2refH 3refB
		for l in match:
			l=l.strip().split('\t')
			Mpos = int(l[1]); refH=l[2].upper(); refB=l[3].upper() #
			if Mpos < vpos:
				if refH != refB: #monomorphic derived
					out.write('%s\t-\t%s\tNA\t%s\t%s\n'%(Mpos,refB,N,N))
				else: #monomorphic ancestral
					out.write('%s\t%s\t%s\tNA\t0\t%s\n'%(Mpos,refH,refB,N))
			elif Mpos == vpos:
				ref = lv[1].upper(); altH=lv[2].upper()
				if (ref in Bases) and (altH in Bases) and len(set([altH,ref,refB,refH]))==2: #bi-allelic
					x=int(lv[3]); n=int(lv[4])
					if Mpos in pval:
						p = pval[Mpos]
					else:
						p='NA'
					if refH == altH: #alt=refH is ancestral
						out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(Mpos,altH,refB,p,n-x,n))
					elif refH == ref: #ref=refH, the ancestral
						out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(Mpos,ref,refB,p,x,n))
					else:
						print([altH,ref,refB,refH])
						sys.exit()
				#else:
					#print(lv)
				#	print refH, ref, altH, refB, set([refH, ref, altH, refB]), len(set([refH, ref, altH, refB])), (ref in Bases), (altH in Bases)
			else: #Mpos > vpos:
				while Mpos > vpos:
					lv=vcf.readline().strip().split('\t')
					if lv == ['']:
						break
					vpos = int(lv[0])
				#if Mpos > vpos:
				#	print l,'\n',lv
	out.close()


def main():
	ch = sys.argv[1]
	path = sys.argv[2]
	vcffile='%sparsed_vcf/Panpan_hg38_Chr%s_counts-from-vcf.txt' % (path,ch)
	pvalfile = '%shg38hwe/pval_onetail_n26_chr%s.out'%(path,ch)
	matchfile = '%shg38matchPP2/chr%s.hg38.pP2.txt.gz' % (path,ch)
	outfile = '%sinput/Chr%s.hg38.panPan2_drvFreq.txt'%(ch)
	print datetime.now()
	findMatch(matchfile, vcffile, pvalfile, outfile)
	print datetime.now()

if __name__ == '__main__':
	main()
