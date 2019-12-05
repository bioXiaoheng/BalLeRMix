# script to get average coverage for each population, by each chromosome
# Request 1 processors on 1 node
#PBS -l nodes=1:ppn=1

# Request 4 hours of wall clock time
#PBS -l walltime=200:00:00

# Request 1 gigabyte of memory per process
#PBS -l pmem=500gb

# Request regular output go to the same file
#PBS -j oe

module load gcc/5.3.1
#module load bamtools/2.4.1
module load bedtools/2.26.0

#pop='KHV'
ch=$chr

date
cd /gpfs/scratch/xzc5154/coverage_${pop}/

ids=$( cat ${pop}_IDs.txt | awk '{printf $1" "}' )
num=$( cat ${pop}_IDs.txt | awk 'END{print NR}')
#grab chromosomes
echo 'Making sure the files are tab delimited'
for i in $( seq 1 ${num} );
do
	id=$( head -$i ${pop}_IDs.txt | tail -1 | cut -f 1 )
	echo extract chromosome $ch from individual $id
	cat beds/${id}.bg | awk -v c=$ch 'BEGIN{OFS="\t"}{if($1=c) print $0}' > tempChrs/${id}_chr${ch}.bg
	#make sure they're sorted
	bedtools sort -i tempChrs/${id}_chr${ch}.bg > tempChrs/${id}_chr${ch}.sorted.bg
done
files=$( cat ${pop}_IDs.txt | awk -v c=$ch '{printf "tempChrs/"$1"_chr"c"sorted.bg "}' )
echo 'files:' $files
echo 'ids:'$ids
#make union file -header 
bedtools unionbedg -i ${files} -names ${ids} > tempChrs/${pop}_chr${ch}.bg
bedtools sort -i tempChrs/${pop}_chr${ch}.bg > tempChrs/${pop}_chr${ch}.sorted.bg
bedtools merge -i tempChrs/${pop}_chr${ch}.sorted.bg > tempChrs/${pop}_chr${ch}.sorted.merged.bg

#get average
echo -e 'chrom\tstart\tend\tavgDepth\tnumZeros\tSD\tnonZeroSD' > avgDepth/hg19.chr${ch}.all${pop}_avgDepth.begGraph

tail -n+2 tempChrs/${pop}_chr${ch}.sorted.merged.bg | awk -v N=${num} 'BEGIN{OFS="\t"}{
	sum=0
	for(i=4;i<=NF;i++){ sum+=$i }
	avg=sum/(NF-3)
	var=0; NZvar=0; NZ=0
	for(i=4;i<=NF;i++){
		var+=($i - avg)*($i - avg)
		if($i != 0){
			NZvar+=($i - avg)*($i - avg)
			NZ+=1
		}
	}
	sd=sqrt(var/(NF-3))
	if (NZ>0) {
		NZsd=sqrt(NZvar/NZ)
	}else{NZsd = 0}
	
	print $1,$2,$3,avg,NF-3-NZ,sd,NZsd
}' >>  avgDepth/hg19.chr${ch}.all${pop}_avgDepth.begGraph

echo 'done'
date

