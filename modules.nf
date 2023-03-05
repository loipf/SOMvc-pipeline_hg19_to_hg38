
//TODO
// dont public index_reference

// need params.data_dir declared in main.nf, not implemented in DSL2 yet
params.data_dir	= "$launchDir/data"



process INDEX_REFERENCE { 

	input:
		path reference_genome
		path bed_file

	output:
		tuple path("*.fa"), path("*.fa.fai"), path("*.dict"), emit: reference_genome
		tuple path("*.bed_sorted.gz"), path("*.bed_sorted.gz.tbi"), emit: bed_file

	shell:
	'''
	reference_name=!{reference_genome}
	reference_name=${reference_name%.*}  # removes last file extension .gz

	gunzip -c !{reference_genome} > $reference_name 
	samtools faidx $reference_name -o $reference_name.fai
	samtools dict $reference_name -o ${reference_name%.*}.dict  # remove .fa so name is only .dict

	bedtools sort -i !{bed_file} | bgzip -c > !{bed_file}_sorted.gz
	tabix --zero-based -b 2 -e 3 !{bed_file}_sorted.gz
	'''
}



process BAM_TO_FASTQ_PREPROCESS { 
	tag "$sample_id"
	cache false
	
	input:
		tuple val(sample_id), path(bam_file)
		val num_threads

	output:
		tuple val(sample_id), path("${sample_id}_prepro_1.fastq.gz"), path("${sample_id}_prepro_2.fastq.gz"), emit: reads_prepro
		path "${sample_id}_cutadapt_output.txt", emit: cutadapt

	shell:
	'''
	
	### extract reads from bam files
	samtools sort -@ !{num_threads} -n -o sorted_n.bam !{bam_file} 
	bedtools bamtofastq -i sorted_n.bam -fq !{sample_id}_reads_1.fastq -fq2 !{sample_id}_reads_2.fastq

	cutadapt --max-n 0.1 --discard-trimmed --pair-filter=any --minimum-length 10 -o !{sample_id}_prepro_1.fastq.gz -p !{sample_id}_prepro_2.fastq.gz !{sample_id}_reads_1.fastq !{sample_id}_reads_2.fastq > !{sample_id}_cutadapt_output.txt

	'''
}



process CREATE_BWA_INDEX { 
	publishDir "$params.data_dir/bwa_index", mode: "copy"

	input:
		path reference_genome

	output:
		path "*.{amb,ann,bwt,pac,sa}", emit: bwa_index

	shell:
	'''
	bwa index !{reference_genome} 
	'''
}




process MAPPING_BWA { 
	tag "$sample_id"
	publishDir "$params.data_dir/reads_mapped", mode: 'copy', pattern:"*_stats.txt", saveAs: { filename -> "${sample_id}/$filename" }
	cache false

	input:
		tuple val(sample_id), path(reads_prepro_1), path(reads_prepro_2) 
		val num_threads
		path reference_genome
		path bwa_index  // to ensure index is created
		val file_name_extension


	output:
		tuple val(sample_id), path("${sample_id}${file_name_extension}.bam"), path("${sample_id}${file_name_extension}.bam.bai"), emit: reads_mapped
		path "*", emit: all


	shell:
	'''
	sample_name=!{sample_id}!{file_name_extension}
	
	bwa mem -Y -R "@RG\\tID:$sample_name\\tSM:$sample_name" -t !{num_threads} -K 100000000 !{reference_genome} !{reads_prepro_1} !{reads_prepro_2}\
	| samtools view -@ !{num_threads} -h -b -u - \
	| samtools sort -n -@ !{num_threads} - \
	| samtools fixmate -m -@ !{num_threads} - - \
	| samtools sort -@ !{num_threads} - \
	| samtools markdup -@ !{num_threads} -f  $sample_name"_markdup_stats.txt" - $sample_name".bam"

	samtools index -b -@ !{num_threads}  $sample_name".bam"
	samtools stats -@ !{num_threads}  $sample_name".bam" >  $sample_name"_stats.txt"

	'''
}







process SOMVC_LOFREQ { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/lofreq", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(normal_file), path(normal_file_index), path(tumor_file), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		tuple val("$sample_id"), path("lofreq_somatic_final.snvs.vcf.gz"), path("lofreq_somatic_final.snvs.vcf.gz.tbi"), path("lofreq_somatic_final.indels_vt.vcf.gz"), path("lofreq_somatic_final.indels_vt.vcf.gz.tbi"), emit: lofreq_output
		tuple path("lofreq_somatic_raw.snvs.vcf.gz"), path("lofreq_somatic_raw.snvs.vcf.gz.tbi")
		tuple path("lofreq_somatic_raw.indels.vcf.gz"), path("lofreq_somatic_raw.indels.vcf.gz.tbi")

	shell:
	'''
	### https://csb5.github.io/lofreq/commands/#somatic
	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # vardict need unzipped

	lofreq viterbi -f !{reference_genome[0]} !{normal_file} | samtools sort -@ !{num_threads} -o normal_file_viterbi.bam -
	lofreq indelqual --dindel -f !{reference_genome[0]} -o normal_file_indelqual.bam normal_file_viterbi.bam
	samtools index -b -@ !{num_threads} normal_file_viterbi.bam

	lofreq viterbi -f !{reference_genome[0]} !{tumor_file} | samtools sort -@ !{num_threads} -o tumor_file_viterbi.bam -
	lofreq indelqual --dindel -f !{reference_genome[0]} -o tumor_file_indelqual.bam tumor_file_viterbi.bam
	samtools index -b -@ !{num_threads} tumor_file_viterbi.bam

	lofreq somatic -n normal_file_viterbi.bam -t tumor_file_viterbi.bam -f !{reference_genome[0]} --threads !{num_threads} -o lofreq_ -l bed_file_unzipped.bed --call-indels
	
	### vt normalization for indels
	vt normalize lofreq_somatic_final.indels.vcf.gz -r !{reference_genome[0]} -o lofreq_somatic_final.indels_vt.vcf.gz
	tabix -p vcf lofreq_somatic_final.indels_vt.vcf.gz

	'''
}


process SOMVC_MUTECT2 { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/mutect2", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(normal_file), path(normal_file_index), path(tumor_file), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		tuple val("$sample_id"), path("mutect2_filtered_vt.vcf.gz"), path("mutect2_filtered_vt.vcf.gz.tbi"), emit: mutect2_output
		path "mutect2_filtered.vcf.filteringStats.tsv"

	shell:
	'''
	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # mutect2 need unzipped
	
	gatk GetSampleName -I !{normal_file} -O sample_normal.txt
	normal_sample_name=$(cat sample_normal.txt)
	gatk GetSampleName -I !{tumor_file} -O sample_tumor.txt
	tumor_sample_name=$(cat sample_tumor.txt)
	
	gatk Mutect2 -R !{reference_genome[0]} -I !{normal_file} -normal $normal_sample_name -I !{tumor_file} -tumor $tumor_sample_name --native-pair-hmm-threads !{num_threads} --germline-resource /usr/src/mutect2_genome/af-only-gnomad_ensembl.hg38.vcf.gz --panel-of-normals /usr/src/mutect2_genome/1000g_pon_ensembl.hg38.vcf.gz --intervals bed_file_unzipped.bed -O mutect2_unfiltered.vcf
	gatk FilterMutectCalls -R !{reference_genome[0]} -V mutect2_unfiltered.vcf -O mutect2_filtered.vcf

	### rename for somatic-combiner
	printf '%s\n' "$normal_sample_name !{sample_id}_normal" "$tumor_sample_name !{sample_id}_tumor" > sample_names.txt 
	bcftools reheader --samples sample_names.txt -o mutect2_filtered_name.vcf mutect2_filtered.vcf

	bgzip -c mutect2_filtered_name.vcf > mutect2_filtered.vcf.gz
	vt normalize mutect2_filtered.vcf.gz -r !{reference_genome[0]} -o mutect2_filtered_vt.vcf.gz
	tabix -p vcf mutect2_filtered_vt.vcf.gz
	'''
}



process SOMVC_STRELKA { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/strelka", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(normal_file), path(normal_file_index), path(tumor_file), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		tuple val("$sample_id"), path("somatic.snvs.vcf.gz"), path("somatic.snvs.vcf.gz.tbi"), path("somatic.indels_vt.vcf.gz"), path("somatic.indels_vt.vcf.gz.tbi"), emit: strelka_output
		path manta_sv

	shell:
	'''
	configManta.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome[0]} --runDir manta_dir
	manta_dir/runWorkflow.py -m local -j !{num_threads}	

	configureStrelkaSomaticWorkflow.py --tumorBam !{tumor_file} --normalBam !{normal_file} --referenceFasta !{reference_genome[0]} --exome --runDir strelka_dir --indelCandidates manta_dir/results/variants/candidateSmallIndels.vcf.gz --callRegions !{bed_file[0]}
	strelka_dir/runWorkflow.py -m local -j !{num_threads}

	printf '%s\n' "NORMAL !{sample_id}_normal" "TUMOR !{sample_id}_tumor" > sample_names.txt 
	bcftools reheader --samples sample_names.txt -o somatic.snvs.vcf.gz strelka_dir/results/variants/somatic.snvs.vcf.gz
	bcftools reheader --samples sample_names.txt -o somatic.indels.vcf.gz strelka_dir/results/variants/somatic.indels.vcf.gz
	tabix -p vcf somatic.snvs.vcf.gz

	### vt normalization for indels
	vt normalize somatic.indels.vcf.gz -r !{reference_genome[0]} -o somatic.indels_vt.vcf.gz
	tabix -p vcf somatic.indels_vt.vcf.gz
	
	mv manta_dir/results/variants manta_sv
	'''
}


process SOMVC_VARDICT { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/vardict", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(normal_file), path(normal_file_index), path(tumor_file), path(tumor_file_index) 
		path reference_genome
		path bed_file
		val num_threads

	output:
		tuple val("$sample_id"), path("vardict_output_filtered.vcf.gz"), path("vardict_output_filtered.vcf.gz.tbi"), emit: vardict_output

	shell: 
	'''
	gunzip -c !{bed_file[0]} > bed_file_unzipped.bed   # vardict need unzipped

	# /usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome[0]} -k 1 -b "!{tumor_file}|!{normal_file}" -Q 5 -z 1 -c 1 -S 2 -E 3 -g 4 -th !{num_threads} bed_file_unzipped.bed | /usr/src/VarDict-1.8.2/bin/testsomatic.R | /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M -N "!{sample_id}_tumor|!{sample_id}_normal" > vardict_output.vcf
	
	### bugfix https://github.com/aryarm/varCA/issues/42 - but should be undone - removes <dup-10> lines
	/usr/src/VarDict-1.8.2/bin/VarDict -G !{reference_genome[0]} -k 1 -b "!{tumor_file}|!{normal_file}" -Q 5 -z 1 -c 1 -S 2 -E 3 -g 4 -th !{num_threads} bed_file_unzipped.bed | /usr/src/VarDict-1.8.2/bin/testsomatic.R | /usr/src/VarDict-1.8.2/bin/var2vcf_paired.pl -P 0.9 -m 4.25 -f 0.01 -M -N "!{sample_id}_tumor|!{sample_id}_normal" | awk -F $"\t" -v 'OFS=\t' '/^#/ || $5 !~ /<dup/' > vardict_output.vcf

	### https://github.com/bcbio/bcbio-nextgen/blob/5cfc02b5974d19908702fa21e6d2f7a50455b44c/bcbio/variation/vardict.py#L248
	gatk VariantFiltration -R !{reference_genome[0]} -V vardict_output.vcf -O vardict_output_filtered.vcf --filter-name bcbio_advised --filter-expression "((AF*DP<6)&&((MQ<55.0&&NM>1.0)||(MQ<60.0&&NM>2.0)||(DP<10)||(QUAL<45)))" 

	bgzip -c vardict_output_filtered.vcf > vardict_output_filtered.vcf.gz
	tabix -p vcf vardict_output_filtered.vcf.gz
	'''
}



process SOMATIC_COMBINER { 
	tag "$sample_id"
	publishDir "$params.data_dir/vc_caller/somatic_combiner", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(lofreq_snv_vcf), path(lofreq_snv_vcf_index), path(lofreq_indel_vcf), path(lofreq_indel_vcf_index), path(mutect2_vcf), path(mutect2_vcf_index), path(strelka_snv_vcf), path(strelka_snv_vcf_index), path(strelka_indel_vcf), path(strelka_indel_vcf_index), path(vardict_vcf), path(vardict_vcf_index)

	output:
		tuple val("${sample_id}"), path("*_somatic_combiner_all.vcf.gz"), path("*_somatic_combiner_all.vcf.gz.tbi"), emit: somatic_combiner_vcf


	shell:
	'''
	java -jar /usr/src/somaticCombiner.jar -L !{lofreq_indel_vcf} -l !{lofreq_snv_vcf} -M !{mutect2_vcf} -s !{strelka_snv_vcf} -S !{strelka_indel_vcf} -D !{vardict_vcf} -o somatic_combiner_raw.vcf

	printf '%s\n' "TUMOR !{sample_id}_tumor" "NORMAL !{sample_id}_normal" > sample_names.txt
	bcftools reheader --samples sample_names.txt -o somatic_combiner_sample.vcf somatic_combiner_raw.vcf
	
	### somatic combiner defined but not added to vcf
	sed -i '7 i ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed, according to somatic-combiner">' somatic_combiner_sample.vcf
	
	bgzip -c somatic_combiner_sample.vcf > !{sample_id}_somatic_combiner_all.vcf.gz
	tabix -p vcf !{sample_id}_somatic_combiner_all.vcf.gz
	'''
}




process CONPAIR_CONTAMINATION { 
	tag "$sample_id"
	publishDir "$params.data_dir/conpair", mode: 'copy', saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(normal_file), path(normal_file_index), path(tumor_file), path(tumor_file_index) 
		path reference_genome

	output:
		path "*_stats.txt", emit: conpair_info

	shell:
	'''
	MARKER_FILE_BED="/usr/src/conpair/conpair/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
	sed -e 's/chr//g' $MARKER_FILE_BED > markers_GRCh38_snv_formatted.bed  ### rename chromosomes

	MARKER_FILE_TXT="/usr/src/conpair/conpair/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
	sed -e 's/chr//g' $MARKER_FILE_TXT > markers_GRCh38_snv_formatted.txt  ### rename chromosomes

	run_gatk_pileup_for_sample.py -B !{tumor_file} -O tumor_pileup --reference !{reference_genome[0]} --conpair_dir /usr/src/conpair/ --markers markers_GRCh38_snv_formatted.bed
	run_gatk_pileup_for_sample.py -B !{normal_file} -O normal_pileup --reference !{reference_genome[0]} --conpair_dir /usr/src/conpair/ --markers markers_GRCh38_snv_formatted.bed

	verify_concordance.py -T tumor_pileup -N normal_pileup --markers markers_GRCh38_snv_formatted.txt --outfile !{sample_id}_concordance_stats.txt

	estimate_tumor_normal_contamination.py -T tumor_pileup -N normal_pileup --markers markers_GRCh38_snv_formatted.txt --outfile !{sample_id}_contamination_stats.txt
	'''
}



process VARIANT_CALLING_STATS { 
	container "dnavc-pipeline:latest"
	tag "$sample_id"
	publishDir "$params.data_dir/variants_vcf", mode: "copy", overwrite: false, saveAs: { filename -> "${sample_id}/$filename" }

	input:
		tuple val(sample_id), path(vcf_file), path(vcf_file_index)
		val num_threads

	output:
		path "*_ADJ_PASS_vcfstats.txt", emit: vcf_stats

	shell:
	'''
	bcftools stats -f ADJ_PASS -s !{sample_id}_normal --threads !{num_threads} !{vcf_file} > !{sample_id}_normal_ADJ_PASS_vcfstats.txt
	bcftools stats -f ADJ_PASS -s !{sample_id}_tumor --threads !{num_threads} !{vcf_file} > !{sample_id}_tumor_ADJ_PASS_vcfstats.txt
	'''
}


process MERGE_VCF { 
	container "dnavc-pipeline:latest"
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false

	input:
		path vcf_files
		val num_threads

	output:
		path "all_samples_vcf_merged.vcf.gz", emit: vcf_all

	shell:
	'''
	bcftools merge -Oz -o all_samples_vcf_merged.vcf.gz --threads !{num_threads} *_somatic_combiner_all.vcf.gz
	'''
}




process VARIANT_ANNOTATION { 
	container "dnavc-pipeline:latest"
	publishDir "$params.data_dir/variants_vcf/_all", mode: "copy", overwrite: false

	input:
		path vcf_all
		val num_threads

	output:
		path "*"

	shell:
	'''
	oc module ls -t annotator > oc_databases_version.txt  ### output OpenCRAVAT database versions
	
	oc run -l hg38 -t csv -x --mp !{num_threads} !{vcf_all}
	'''
}


process MULTIQC_VCF { 
	container "dnavc-pipeline:latest"
	publishDir "$params.data_dir/quality_reports", mode: "copy"

	input:
		path stat_files
		path conpair_files

	output:
		path "*"

	shell:
	'''
	multiqc -f -o variants_vcf .
	'''
}

















