#Jan 30, 2018: Modified for parsing pilot bonobo data. Use panTro2 as ancestral call. Ignore hg38 ref

#May 24, 2018: Modified s.t. compatible with other population sizes
#May 9, 2018: Modified script to make sure to output derived freq

#popSize = {'YRI':216, 'CEU':198, 'CHB':206, 'CHS':210, 'CDX':186}

#go through alignment and vcf at the same time and generate derived frequency
#use hg19.pT5 as match template
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
					#if ref == refB: #alt=refH is ancestral
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
	vcffile='/gpfs/group/mxd60/default/Xiaoheng_temp/BallerMix/Empirical/parsed_vcf/Panpan_hg38_Chr%s_counts-from-vcf.txt' % (ch)
	pvalfile = '/gpfs/group/mxd60/default/Xiaoheng_temp/BallerMix/Empirical/hg38hwe/pval_onetail_n26_chr%s.out'%(ch)
	matchfile = '/gpfs/group/mxd60/default/Xiaoheng_temp/BallerMix/Empirical/hg38matchPP2/chr%s.hg38.pP2.txt.gz' % (ch)
	outfile = '/gpfs/group/mxd60/default/Xiaoheng_temp/BallerMix/Empirical/input/redo-Chr%s.hg38.panPan2_drvFreq.txt'%(ch)
	print datetime.now()
	findMatch(matchfile, vcffile, pvalfile, outfile)
	print datetime.now()

if __name__ == '__main__':
	main()
