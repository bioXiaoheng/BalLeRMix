#Jan 30, 2019: Adopted to use panPan2 to call substitution to hg18 alignment
#python script to call ancestral with pantro5 on hg19 coordinates.

import gzip, os, sys
from datetime import datetime

ptChrs = ['chr1','chr2a','chr2b'] + ['chr%s'%(i) for i in xrange(3,23)]

pair = {'A':'T','T':'A','G':'C','C':'G','-':'-'}

Bases = ['A','T','C','G']

def match(ch,axtfile, hfile):
	Hanc = gzip.open(hfile,'wb')
	Hanc.write('Chr\tposition\tref_hg\tref_pt\n')
	axt = gzip.open(axtfile, 'rb')
	l = axt.readline()
	while l[:1] == '#':
		l = axt.readline()
	l=l.strip().split(' ') #format: 0# 1chr_hg 2start_hg 3end_hg 4chr_pt 5start_pt 6end_pt 7strand 8blastscore
	#start reading through the file
	while l != ['']:
		while l[4] not in ptChrs:
			#read on
			for meh in xrange(4):
				l=axt.readline().strip().split(' ')
			if l == ['']: break
		#here they're a matched line
		if l == ['']: break
		start = int(l[2]); end = int(l[3])
		#size = end-start+1
		hgseq = axt.readline().strip().upper()
		ptseq = axt.readline().strip().upper()
		pos=start-1
		#with gaps, size < len(hgseq)
		for i in xrange(len(hgseq)):
			if hgseq[i] in Bases:
				pos+=1
				if ptseq[i] in Bases:
					ptRef = ptseq[i]
					hgRef = hgseq[i]
					Hanc.write('%s\t%s\t%s\t%s\n'%(ch,pos,hgRef,ptRef))

			''' #Gaps in human and in chimp cannot be processed equally 
			if hgseq[i] in Bases and ptseq[i] in Bases: #ignore gaps
				anc = ptseq[i]
				ref = hgseq[i]
				pos += 1
				Hanc.write('%s\t%s\t%s\t%s\n'%(ch,pos,ref,anc))'''
		#read the next alignment
		l=axt.readline()
		l=axt.readline().strip().split(' ')
	axt.close()
	Hanc.close()




def main():
	ch = sys.argv[1]
  path = sys.argv[2]
	axtfile = '%s/hg18_to_panTro2/chr%s.hg18.panTro2.net.axt.gz' % (path,ch)
	hfile = '%s/hg18matchPT2/chr%s.hg18.pT2.txt.gz' %(path,ch)
	match(ch, axtfile, hfile)
	print 'Chrmosome', ch, 'done at', datetime.now()

if __name__ == '__main__':
	main()
