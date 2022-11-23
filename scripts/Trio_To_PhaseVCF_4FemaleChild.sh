#!/bin/bash
if [[ -n "$1" && -n "$2" && -n "$3" && -n "$4" ]]; then
	echo "Job started: `date`"

	awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } !($1,$2,$5) in seen { print $0 }' $1 $2 | awk -F'\t' '$1 !~ /X/ && $1 !~ /Y/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /1\/1/ || $10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /1\|1/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' > $4"_PternalSNVs.vcf"

	awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]=$10; next } ($1,$2,$5) in seen { print $0, seen[$1,$2,$5] }' $1 $2 | awk -F'\t' '$1 !~ /X/ && $1 !~ /Y/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1 && $11 !~ /1\/1/ && $11 !~ /1\|1/ && $11 !~ /1\/\./ && $11 !~ /1\|\./ && $11 !~ /\.\/1/ && $11 !~ /\.\|1/' | awk -F'\t' '$10 ~ /1\/1/ || $10 ~ /1\|1/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' >> $4"_PternalSNVs.vcf"

        awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } !($1,$2,$5) in seen { print $0, seen[$1,$2,$5] }' $1 $2 | awk -F'\t' '$1 ~ /X/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /\/1/ || $10 ~ /\|1/ || $10 ~ /^1/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' >> $4"_PternalSNVs.vcf"

	awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]=$10; next } ($1,$2,$5) in seen { print $0, seen[$1,$2,$5] }' $1 $2 | awk -F'\t' '$1 ~ /X/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1 && $11 !~ /1\/1/ && $11 !~ /1\|1/ && $11 !~ /1\/\./ && $11 !~ /1\|\./ && $11 !~ /\.\/1/ && $11 !~ /\.\|1/' | awk -F'\t' '$10 ~ /\/1/ || $10 ~ /\|1/ || $10 ~ /^1/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' >> $4"_PternalSNVs.vcf" 



	awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } !($1,$2,$5) in seen { print $0 }' $2 $1 | awk -F'\t' '$1 !~ /Y/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /1\/1/ || $10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /1\|1/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' > $4"_MaternalSNVs.vcf"

	awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]=$10; next } ($1,$2,$5) in seen { print $0, seen[$1,$2,$5] }' $2 $1 | awk -F'\t' '$1 !~ /X/ && $1 !~ /Y/ && $4 != "." && $5 != "." && length($4)==1 && length($5)==1 && $11 !~ /1\/1/ && $11 !~ /1\|1/ && $11 !~ /1\/\./ && $11 !~ /1\|\./ && $11 !~ /\.\/1/ && $11 !~ /\.\|1/' | awk -F'\t' '$10 ~ /1\/1/ || $10 ~ /1\|1/' | awk 'BEGIN{ FS=OFS="\t" } NR==FNR { seen[$1,$2,$5]; next } ($1,$2,$5) in seen { print $0 }' - $3 | awk -F'\t' '$4 != "." && $5 !=  "." && length($4)==1 && length($5)==1' | awk -F'\t' '$10 ~ /0\/1/ || $10 ~ /1\/0/ || $10 ~ /0\|1/ || $10 ~ /1\|0/' >> $4"_MaternalSNVs.vcf" 


        awk -F'\t' '{sub(/0\/1/,"0|1",$10); sub(/1\/0/,"0|1",$10);sub(/1\|0/,"0|1",$10); print}' OFS='\t' $4"_PternalSNVs.vcf" > $4"_Ready_4_NanoMethPhase_Maternal-As-HP1_Paternal-As-HP2.vcf"  

	awk -F'\t' '{sub(/0\/1/,"1|0",$10); sub(/1\/0/,"1|0",$10);sub(/0\|1/,"1|0",$10); print}' OFS='\t' $4"_MaternalSNVs.vcf" >> $4"_Ready_4_NanoMethPhase_Maternal-As-HP1_Paternal-As-HP2.vcf" 
	rm $4"_MaternalSNVs.vcf"
        rm $4"_PternalSNVs.vcf" 
	sort -k1,1 -k2,2n $4"_Ready_4_NanoMethPhase_Maternal-As-HP1_Paternal-As-HP2.vcf" -o $4"_Ready_4_NanoMethPhase_Maternal-As-HP1_Paternal-As-HP2.vcf"


	echo "Job finished: `date`" 

else
	echo -e "Please insert arguments:\n
	First argument is maternal's vcf Second is paternal's vcf, third is child's vcf forth is output prefixt.
	Example: Trio_To_PhaseVCF_4FemaleChild.sh Maternal_SNVs.vcf Paternal_SNVs.vcf Child_SNVs.vcf Output
	Note: Make sure your input vcf files contain only SNVs. We also recommend using high-quality PASS SNVs only 
	because there will be no filtering step to keep high quality Passed variants."
fi
