
/* 
 * SOMATIC-VC PIPELINE 
 * for paired normal-tumor mapped reads
 */


/* 
 * import modules 
 */

nextflow.enable.dsl=2

include { 
	INDEX_REFERENCE;
	CREATE_BWA_INDEX;
	BAM_TO_FASTQ_PREPROCESS as BAM_TO_FASTQ_PREPROCESS_NORMAL;
	BAM_TO_FASTQ_PREPROCESS as BAM_TO_FASTQ_PREPROCESS_TUMOR;
	MAPPING_BWA as MAPPING_BWA_NORMAL;
	MAPPING_BWA as MAPPING_BWA_TUMOR;
	SOMVC_LOFREQ;
	SOMVC_MUTECT2;
	SOMVC_STRELKA;
	SOMVC_VARDICT;
	SOMATIC_COMBINER;
	CONPAIR_CONTAMINATION;
	VARIANT_CALLING_STATS;
	MERGE_VCF;
	MULTIQC_VCF;
	VARIANT_ANNOTATION;
} from './modules.nf' 




/*
 * default parameters
 */ 

params.dev_samples = -1

params.project_dir		= "$projectDir"
params.sample_match_file	= "$params.project_dir/normal_tumor_pairs.csv"
params.data_dir		= "$params.project_dir/data"
params.scripts_dir		= "$params.project_dir/scripts"


/*
 * other parameters
 */

params.num_threads		= 3
params.reference_genome	= "$params.project_dir/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
params.bed_file		= "$params.project_dir/data/Homo_sapiens.GRCh38.cds.all.bed"



log.info """\
SOMATIC-VC PIPELINE
===================================================
sample_match_file	: $params.sample_match_file
data_dir		: $params.data_dir
reference_genome	: $params.reference_genome
bed_file		: $params.bed_file

===================================================

"""



/* 
 * main pipeline logic
 */
workflow {

	channel_sample_match = Channel
			.fromPath(params.sample_match_file)
			.splitCsv(header:true)
			.map{ row -> tuple(row.sample_id, row.normal_file, row.tumor_file) }
			.ifEmpty { error "cannot read in sample_match_file correctly: ${params.sample_match_file}" }
			.take( params.dev_samples )  // only consider a few files for debugging
	//channel_sample_match.view()

	// preprocessing
	INDEX_REFERENCE(params.reference_genome, params.bed_file)
	CREATE_BWA_INDEX(params.reference_genome)


	// mapping
	BAM_TO_FASTQ_PREPROCESS_NORMAL(channel_sample_match.map{ it -> tuple(it[0], it[1])}, params.num_threads)
	MAPPING_BWA_NORMAL(BAM_TO_FASTQ_PREPROCESS_NORMAL.out.reads_prepro, params.num_threads, params.reference_genome, CREATE_BWA_INDEX.out.bwa_index.collect(), "_normal")
	
	BAM_TO_FASTQ_PREPROCESS_TUMOR(channel_sample_match.map{ it -> tuple(it[0], it[2])}, params.num_threads)
	MAPPING_BWA_TUMOR(BAM_TO_FASTQ_PREPROCESS_TUMOR.out.reads_prepro, params.num_threads, params.reference_genome, CREATE_BWA_INDEX.out.bwa_index.collect(), "_tumor")
	
	channel_sample_match_mapped = MAPPING_BWA_NORMAL.out.reads_mapped.join(MAPPING_BWA_TUMOR.out.reads_mapped)
	//channel_sample_match_mapped.view()
	

	// variant calling
	SOMVC_LOFREQ(channel_sample_match_mapped, INDEX_REFERENCE.out.reference_genome, INDEX_REFERENCE.out.bed_file, params.num_threads)
	SOMVC_MUTECT2(channel_sample_match_mapped, INDEX_REFERENCE.out.reference_genome, INDEX_REFERENCE.out.bed_file, params.num_threads)
	SOMVC_STRELKA(channel_sample_match_mapped, INDEX_REFERENCE.out.reference_genome, INDEX_REFERENCE.out.bed_file, params.num_threads)
	SOMVC_VARDICT(channel_sample_match_mapped, INDEX_REFERENCE.out.reference_genome, INDEX_REFERENCE.out.bed_file, params.num_threads)

	channel_all_vc = SOMVC_LOFREQ.out.lofreq_output
			.join(SOMVC_MUTECT2.out.mutect2_output, by: 0)
			.join(SOMVC_STRELKA.out.strelka_output, by: 0)
			.join(SOMVC_VARDICT.out.vardict_output, by: 0)
	//SOMATIC_COMBINER(channel_all_vc)
	
	//CONPAIR_CONTAMINATION(channel_sample_match_mapped, INDEX_REFERENCE.out.reference_genome)
	
	
	//VARIANT_CALLING_STATS(SOMATIC_COMBINER.out.somatic_combiner_vcf, params.num_threads); 
	//MERGE_VCF(SOMATIC_COMBINER.out.somatic_combiner_vcf.collect{[it[1],it[2]]}, params.num_threads) 
	//MULTIQC_VCF(VARIANT_CALLING_STATS.out.vcf_stats.collect(), CONPAIR_CONTAMINATION.out.conpair_info.collect())

}



workflow.onComplete { 
	println ( workflow.success ? "\ndone! check the quality reports in --> $params.data_dir/quality_reports\n" : "oops .. something went wrong" ) } 




