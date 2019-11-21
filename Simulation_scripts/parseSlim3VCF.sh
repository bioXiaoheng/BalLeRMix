#shell script to parse VCFs outputs by slim3.3+

cd $( pwd )

#Path and name of input and output files
infile=$1
outfile=$2

#Index No. of the individual to use as ``ancestral'' sequence
outseq=$3

# Read & parse the vcf
cat $infile | awk -v O=$outseq 'BEGIN{
	OFS="\t"; 
	print "physPos\tgenPos\tx\tn" ;
}{ 
	if(NF > 9 && $9 =="GT" && $5 ~ /^[A|T|C|G]$/ && $4 != $5 ){  
		split($O,g,"|"); anc = g[1]
		drv=0 ; total=0
		for(x=10;x<=34;x++){
			split($x,h,"|")
			if(h[1] != anc) drv++ ;
			if(h[2] != anc) drv++ ;
			total = total + 2
		}
		if(drv>0) print $2, 1e-8*$2,drv,total 
	}
}' > $outfile
