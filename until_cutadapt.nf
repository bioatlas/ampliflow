#!/usr/bin/env nextflow
/*
========================================================================================
						 nf-core/ampliseq
========================================================================================
 nf-core/ampliseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/ampliseq
----------------------------------------------------------------------------------------
*/


def helpMessage() {
	log.info nfcoreHeader()
	log.info"""
	Usage:

	The minimal command for running the pipeline is as follows:
	nextflow run nf-core/ampliseq -profile singularity --reads "data" --FW_primer GTGYCAGCMGCCGCGGTAA --RV_primer GGACTACNVGGGTWTCTAAT


	Main arguments:
	  -profile [strings]            Use this parameter to choose a configuration profile. If not specified, runs locally and expects all software
	                                to be installed and available on the `PATH`. Otherwise specify a container engine, "docker" or "singularity" 
	                                and a specialized profile such as "binac".
	  --reads [path/to/folder]      Folder containing paired-end demultiplexed fastq files
	                                Note: All samples have to be sequenced in one run, otherwise also specifiy "--multipleSequencingRuns"
	  --FW_primer [str]             Forward primer sequence
	  --RV_primer [str]             Reverse primer sequence
	  --metadata [path/to/file]     Path to metadata sheet, when missing most downstream analysis are skipped (barplots, PCoA plots, ...)

	Other input options:
	  --extension [str]             Naming of sequencing files (default: "/*_R{1,2}_001.fastq.gz"). 
	                                The prepended "/" is required, also one "*" is required for sample names and "{1,2}" indicates read orientation
	  --multipleSequencingRuns      If samples were sequenced in multiple sequencing runs. Expects one subfolder per sequencing run
	                                in the folder specified by "--reads" containing sequencing data of the specific run. These folders 
	                                may not contain underscores. Also, fastQC is skipped because multiple sequencing runs might 
	                                create overlapping file names that crash MultiQC.
	  --split [str]                 A string that will be used between the prepended run/folder name and the sample name. (default: "-")
	                                May not be present in run/folder names and no underscore(s) allowed. Only used with "--multipleSequencingRuns"
	  --phred64                     If the sequencing data has PHRED 64 encoded quality scores (default: PHRED 33)

	Filters:
	  --exclude_taxa [str]          Comma separated list of unwanted taxa (default: "mitochondria,chloroplast")
	                                To skip taxa filtering use "none"
	  --min_frequency [int]         Remove entries from the feature table below an absolute abundance threshold (default: 1)
	  --min_samples [int]           Filtering low prevalent features from the feature table (default: 1)                   

	Cutoffs:
	  --retain_untrimmed            Cutadapt will retain untrimmed reads
	  --max_error_rate              Cutadapt level of error tolerance. Allowed errors are mismatches, insertion and deletions. Default 0.2= 20%
	  --trunclenf [int]             DADA2 read truncation value for forward strand
	  --trunclenr [int]             DADA2 read truncation value for reverse strand
	  --trunc_qmin [int]            If --trunclenf and --trunclenr are not set, 
	                                these values will be automatically determined using 
	                                this mean quality score (not preferred) (default: 25)

	Options for R-DADA2:
	  --trunclenF [int]    		Value for truncation position for forward reads. No default, mandatory.
	  --trunclenR [int]             Value for truncation position for reverse reads. No default, mandatory.
	  --fwdmark [str]    		Name pattern to identify forward reads in file names, default "_R1".
	  --revmark [str]    		Name pattern to identify reverse reads in file names, default "_R2".
	  --nsamples [int]              Number of samples to select from DADA2 learning algorithm (process dada2errmodels).
	                                Should be high enough to make sure at least 1 milion reads are included. If 0, an eight of the number of samples will be used,
	                                assuming you have at least 10 milion reads in total, fairly evenly distributed over samples.
	                                Default 24, testing proved that the learning algorithm cannot keep up with more than 24 samples at once.
	  --maxconsist [int]		Maximum number of iterations to achieve self consistency (process dada2errmodels), default "10".
	  --concatenate [str]           Concatenate sequences instead of try to merge (use when overlap is too short), default FALSE. Setting this overrides --maxmismatch and --minoverlap.
	  --minoverlap [int]            Minimum number of overlapping bases in merge, default 20.
	  --maxmismatch [int]           Maximum number of non matching bases in merge, default 0.
	  --method [str] 		Method used for bimera detection, default "consensus". See R doc for removeBimeraDenovo. 
	  --minab [int]                 Minimum parent abundance, default 8. See R doc for isBimeraDenovo.
	  --overab [int]                Parent overabundance multiplier, default 1. See R doc for isBimeraDenovo.
	  --bimeraoff [str]             Handles the off bimeras. If '--oneoff', allows one off bimeras; if '--nooneoff' does not allow one off bimeras.
	                                Default '--oneoff'. See R doc for isBimeraDenovo.

	References:                     If you have trained a compatible classifier before
	  --classifier [path/to/file]   Path to QIIME2 classifier file (typically *-classifier.qza)
	  --classifier_removeHash       Remove all hash signs from taxonomy strings, resolves a rare ValueError during classification (process classifier)

	Statistics:
	  --metadata_category [str]     Comma separated list of metadata column headers for statistics (default: false)
	                                If not specified, all suitable columns in the metadata sheet will be used.
	                                Suitable are columns which are categorical (not numerical) and have multiple  
	                                different values that are not all unique.

	Other options:
	  --untilQ2import               Skip all steps after importing into QIIME2, used for visually choosing DADA2 parameter
	  --Q2imported [path/to/file]   Path to imported reads (e.g. "demux.qza"), used after visually choosing DADA2 parameter
	  --onlyDenoising               Skip all steps after denoising, produce only sequences and abundance tables on ASV level
	  --keepIntermediates           Keep additional intermediate files, such as trimmed reads or various QIIME2 archives
	  --outdir [path/to/folder]     The output directory where the results will be saved
	  --email [email]               Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	  --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, 
	                                it will not be attached (Default: 25MB)
	  -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

	Skipping steps:
	  --skip_fastqc                 Skip FastQC
	  --skip_alpha_rarefaction      Skip alpha rarefaction
	  --skip_taxonomy               Skip taxonomic classification
	  --skip_barplot                Skip producing barplot
	  --skip_abundance_tables       Skip producing any relative abundance tables
	  --skip_diversity_indices      Skip alpha and beta diversity analysis
	  --skip_ancom                  Skip differential abundance testing     
	""".stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
	helpMessage()
	exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")
Channel.fromPath("$baseDir/assets/matplotlibrc")
	.into { ch_mpl_for_make_classifier; ch_mpl_for_qiime_import; ch_mpl_for_ancom_asv; ch_mpl_for_ancom_tax; ch_mpl_for_ancom; ch_mpl_for_beta_diversity_ord; ch_mpl_for_beta_diversity; ch_mpl_for_alpha_diversity; ch_mpl_for_metadata_pair; ch_mpl_for_metadata_cat; ch_mpl_for_diversity_core; ch_mpl_for_alpha_rare; ch_mpl_for_tree; ch_mpl_for_barcode; ch_mpl_for_relreducetaxa; ch_mpl_for_relasv; ch_mpl_for_export_dada_output; ch_mpl_filter_taxa; ch_mpl_classifier; ch_mpl_dada; ch_mpl_dada_merge; ch_mpl_for_demux_visualize; ch_mpl_for_classifier }

// Defines all parameters that are independent of a test run
params.trunc_qmin = 25 //to calculate params.trunclenf and params.trunclenr automatically
params.trunclenf = false
params.trunclenr = false
params.metadata_category = false
params.retain_untrimmed = false
params.exclude_taxa = "mitochondria,chloroplast"
params.keepIntermediates = false
params.classifier_removeHash = false
params.min_frequency = false
params.min_samples = false
params.multipleSequencingRuns = false
params.phred64 = false
params.split = "-"

//Database specific parameters
//currently only this is compatible with process make_SILVA_132_16S_classifier
params.reference_database = "https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip"
params.dereplication = 99


/*
 * Define pipeline steps
 */
params.untilQ2import = false

params.Q2imported = false
if (params.Q2imported) {
	params.skip_fastqc = true
	params.skip_multiqc = true
} else {
	params.skip_multiqc = false
}

//Currently, fastqc doesnt work for multiple runs when sample names are identical. These names are encoded in the sequencing file itself.
if (params.multipleSequencingRuns) {
	params.skip_fastqc = true
} else {
	params.skip_fastqc = false
}

params.onlyDenoising = false
if (params.onlyDenoising || params.untilQ2import) {
	params.skip_abundance_tables = true
	params.skip_barplot = true
	params.skip_taxonomy = true
	params.skip_alpha_rarefaction = true
	params.skip_diversity_indices = true
	params.skip_ancom = true
} else {
	params.skip_abundance_tables = false
	params.skip_barplot = false
	params.skip_taxonomy = false
	params.skip_alpha_rarefaction = false
	params.skip_diversity_indices = false
	params.skip_ancom = false
}

params.idtaxa_db = false
params.species_db = false
params.rdp_db = false

/*
 * Import input files
 */
if (params.metadata) {
	Channel.fromPath("${params.metadata}", checkIfExists: true)
		.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
} else {
	Channel.from()
		.into { ch_metadata_for_barplot; ch_metadata_for_alphararefaction; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_metadata_category_all; ch_metadata_for_metadata_category_pairwise; ch_metadata_for_beta_diversity; ch_metadata_for_beta_diversity_ordination; ch_metadata_for_ancom; ch_metadata_for_ancom_tax; ch_metadata_for_ancom_asv }
}

if (params.Q2imported) {
	Channel.fromPath("${params.Q2imported}", checkIfExists: true)
		   .into { ch_qiime_demux_import; ch_qiime_demux_vis; ch_qiime_demux_dada }
}

if (params.classifier) {
	Channel.fromPath("${params.classifier}", checkIfExists: true)
		   .set { ch_qiime_classifier }
}

/*
 * Sanity check input values
 */
if (!params.Q2imported) { 
	if (!params.FW_primer) { exit 1, "Option --FW_primer missing" }
	if (!params.RV_primer) { exit 1, "Option --RV_primer missing" }
	if (!params.reads) { exit 1, "Option --reads missing" }
	if (!params.trunclenF || !params.trunclenR) { exit 1, "Values for mandatory options --trunclenF and --trunclenR missing"}
}

if (params.Q2imported && params.untilQ2import) {
	exit 1, "Choose either to import data into a QIIME2 artefact and quit with --untilQ2import or use an already existing QIIME2 data artefact with --Q2imported."
}

if ("${params.split}".indexOf("_") > -1 ) {
	exit 1, "Underscore is not allowed in --split, please review your input."
}

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
	if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
	if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")


// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/ampliseq'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

if( !params.trunclenF || !params.trunclenR ){
	if ( !params.untilQ2import ) log.info "\n######## WARNING: No DADA2 cutoffs were specified, therefore reads will not be truncated.Truncating reads where median quality drops below ${params.trunc_qmin} is deprecated, and therefore disabled.\n This is because cutoffs based on median quality do not account for required overlap for merging, therefore DADA2 might have poor merging efficiency or even fail.\n The pieplien will now stop.\n"
}

// add if condition to warn users to select a DADA2 taxonomy DB

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
	def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
	yaml_file.text  = """
	id: 'nf-core-ampliseq-summary'
	description: " - this information is collected when the pipeline is started."
	section_name: 'nf-core/ampliseq Workflow Summary'
	section_href: 'https://github.com/nf-core/ampliseq'
	plot_type: 'html'
	data: |
		<dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
		</dl>
	""".stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
	publishDir "${params.outdir}/pipeline_info", mode: 'copy',
	saveAs: {filename ->
		if (filename.indexOf(".csv") > 0) filename
		else null
	}

	output:
	file 'software_versions_mqc.yaml' into ch_software_versions_yaml
	file "software_versions.csv"

	script:
	// TODO nf-core: Get all tools to print their version number here
	"""
	echo $workflow.manifest.version > v_pipeline.txt
	echo $workflow.nextflow.version > v_nextflow.txt
	fastqc --version > v_fastqc.txt
	multiqc --version > v_multiqc.txt
	cutadapt --version > v_cutadapt.txt
	qiime --version > v_qiime.txt
	scrape_software_versions.py &> software_versions_mqc.yaml
	"""
}


if (!params.Q2imported){

	/*
	* Create a channel for input read files
	*/
	if(params.readPaths && params.reads == "data${params.extension}" && !params.multipleSequencingRuns){
		//Test input for single sequencing runs, profile = test

		Channel
			.from(params.readPaths)
			.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
			.ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
			.into { ch_read_pairs; ch_read_pairs_fastqc; ch_read_pairs_name_check }

	} else if ( !params.readPaths && params.multipleSequencingRuns ) {
		//Standard input for multiple sequencing runs

		//Get files
		Channel
			.fromFilePairs( params.reads + "/*" + params.extension, size: 2 )
			.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}/*${params.extension}\nNB: Path needs to be enclosed in quotes!" }
			.into { ch_extract_folders; ch_rename_key }

		//Get folder information
		ch_extract_folders
			.flatMap { key, files -> [files[0]] }
			.map { it.take(it.findLastIndexOf{"/"})[-1] }
			.unique()
			.into { ch_count_folders; ch_check_folders; ch_report_folders }

		//Report folders with sequencing files
		ch_report_folders
			.collect()
			.subscribe {
				String folders = it.toString().replace("[", "").replace("]","") 
				log.info "\nFound the folder(s) \"$folders\" containing sequencing read files matching \"${params.extension}\" in \"${params.reads}\".\n" }

		//Stop if folder count is 1
		ch_count_folders
			.count()
			.subscribe { if ( it == 1 ) exit 1, "Found only one folder with read data but \"--multipleSequencingRuns\" was specified. Please review data input." }
		
		//Stop if folder names contain "_" or "${params.split}"
		ch_check_folders
			.subscribe { 
				if ( it.toString().indexOf("${params.split}") > -1 ) exit 1, "Folder name \"$it\" contains \"${params.split}\", but may not. Please review data input or choose another string using \"--split [str]\" (no underscore allowed!)."
				if ( it.toString().indexOf("_") > -1 ) exit 1, "Folder name \"$it\" contains \"_\", but may not. Please review data input." 
			}

		//Add folder information to sequence files
		ch_rename_key
			.map { key, files -> [ key, files, (files[0].take(files[0].findLastIndexOf{"/"})[-1]) ] }
			.into { ch_read_pairs; ch_read_pairs_fastqc }

	} else if ( params.readPaths && params.multipleSequencingRuns ) {
		//Test input for multiple sequencing runs, profile = test_multi

		Channel
			.from(params.readPaths)
			.map { row -> [ row[0], [file(row[1][0]), file(row[1][1])], row[2] ] }
			.ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
			.into { ch_read_pairs; ch_read_pairs_fastqc }
			
	} else {
		//Standard input

		Channel
			.fromFilePairs( params.reads + params.extension, size: 2 )
			.ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}${params.extension}\nNB: Path needs to be enclosed in quotes!" }
			.into { ch_read_pairs; ch_read_pairs_fastqc }
	}

	/*
	 * fastQC
	 */
	if (!params.multipleSequencingRuns){
		process fastqc {
			tag "${pair_id}"
			publishDir "${params.outdir}/fastQC", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

			input:
			set val(pair_id), file(reads) from ch_read_pairs_fastqc

			output:
			file "*_fastqc.{zip,html}" into ch_fastqc_results

			when:
			!params.skip_fastqc

			script: 
			"""
			fastqc -q ${reads}
			"""
		}
	} else {
		process fastqc_multi {
			tag "${folder}${params.split}${pair_id}"
			publishDir "${params.outdir}/fastQC", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

			input:
			set val(pair_id), file(reads), val(folder) from ch_read_pairs_fastqc

			output:
			file "*_fastqc.{zip,html}" into ch_fastqc_results

			when:
			!params.skip_fastqc

			script: 
			"""
			fastqc -q ${reads}
			"""
		}
	}

	/*
	 * Trim each read-pair with cutadapt
	 */
	if (!params.multipleSequencingRuns){
		process trimming {
			tag "${pair_id}"  
			publishDir "${params.outdir}/trimmed", mode: 'copy',
				saveAs: {filename -> 
				if (filename.indexOf(".gz") == -1) "logs/$filename"
				else if(params.keepIntermediates) filename 
				else null}
			
			input:
			set val(pair_id), file(reads) from ch_read_pairs
		
			output:
			file "trimmed/*.*" into ch_fastq_trimmed, ch_fastq_trimmed_manifest
			file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

			script:
			discard_untrimmed = params.retain_untrimmed ? '' : '--discard-untrimmed'
			"""
			mkdir -p trimmed
			cutadapt -g ${params.FW_primer} -G ${params.RV_primer} ${discard_untrimmed} \
				-o trimmed/${reads[0]} -p trimmed/${reads[1]} -e ${max_error_rate}\
				${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
			"""
		}
	} else {
		process trimming_multi {
			tag "$folder${params.split}$pair_id"  
			publishDir "${params.outdir}/trimmed", mode: 'copy',
				saveAs: {filename -> 
				if (filename.indexOf(".gz") == -1) "logs/$filename"
				else if(params.keepIntermediates) filename 
				else null}
		
			input:
			set val(pair_id), file(reads), val(folder) from ch_read_pairs
		
			output:
			set val(pair_id), file ("trimmed/*.*") into ch_fastq_trimmed
			file "trimmed/*.*" into ch_fastq_trimmed_manifest
			file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

			script:
			discard_untrimmed = params.retain_untrimmed ? '' : '--discard-untrimmed'
			"""
			mkdir -p trimmed
			cutadapt -g ${params.FW_primer} -G ${params.RV_primer} ${discard_untrimmed} \
				-o trimmed/$folder${params.split}${reads[0]} -p trimmed/$folder${params.split}${reads[1]} \
				${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
			"""
		}
	}

	/*
	 * multiQC
	 */
	process multiqc {
		publishDir "${params.outdir}/MultiQC", mode: 'copy'

		input:
		file ('fastqc/*') from ch_fastqc_results.collect()
		file ('cutadapt/logs/*') from ch_fastq_cutadapt_log.collect()

		output:
		file "*multiqc_report.html" into multiqc_report
		file "*_data"

		when:
		!params.skip_multiqc

		script:
		"""
		multiqc --force --interactive .
		"""
	}

	/*
	 *(eemis-dada2) [https://github.com/erikrikarddaniel/eemisdada2]
	 *dada2filter: quality and length trimming 
	 */
	process dada2filter {

		publishDir "${params.outdir}/dada2filter", mode: 'copy',
			saveAs: {filename -> 
			if (filename.indexOf(".gz") == -1) "logs/$filename"
			else if(params.keepIntermediates) filename 
			else null}
		input:
		file(reads) from ch_fastq_trimmed.collect()
	
		output:
		file "dada2filter/*fastq.gz" into ch_fastq_dada2filter, ch_fastq_dada2filter2
		file "dada2filter_log" into ch_fastq_dada2filter_log

		script:
		/*--filterdir: Directory for quality truncated reads, dada2filter will create it if it does not exist.
		*See more on its settings on https://github.com/erikrikarddaniel/eemsdada2/blob/master/src/R/dada2filter#L20
		*/
		//correcting --trunclen by subtracting the length of primers from the truncation position
		trunclenR1 = params.trunclenF - params.FW_primer.size()
		trunclenR2 = params.trunclenR - params.RV_primer.size()
		"""
		dada2filter.R --verbose --trimleft=0,0 --trunclen=${trunclenR1},${trunclenR2} \
		--filterdir=dada2filter --fwdmark=${params.fwdmark} --revmark=${params.revmark} \
		>dada2filter_log 2>&1 
		"""
	}
	//dada2errmodels: calculation of error models from the data
        
	process dada2errmodels  {
	
	 	publishDir "${params.outdir}/dada2errmodels", mode: 'copy',
	 		saveAs: {filename -> 
	 		if (filename.indexOf(".rds") == -1) "logs/$filename"
	 		else if(params.keepIntermediates) filename 
	 		else null}
	  
	 	input:
	 	file(dada2filtered) from ch_fastq_dada2filter.collect()
	 
	 	output:
		file ("*.rds") into ch_fastq_dada2errmodels
	 	file ("dada2errmodels_log") into ch_fastq_dada2errmodels_log
	 	script:
	 	"""
	 	dada2errmodels.R --verbose --nsamples=${params.nsamples} --maxconsist=${params.maxconsist} \
	 	--filterdir=. --fwdmark=${params.fwdmark} --revmark=${params.revmark} \
	 	>dada2errmodels_log 2>&1
	 	"""
	 }

	//dada2cleanNmerge: sequence correction and merge
	process dada2cleanNmerge  {
	
	 	publishDir "${params.outdir}/dada2cleanNmerge", mode: 'copy',
	 		saveAs: {filename -> 
	 		if (filename.indexOf(".merged.rds") == -1) "logs/$filename"
	 		else if(params.keepIntermediates) filename 
	 		else null}
	  
	 	input:
	 	file(dada2filtered) from ch_fastq_dada2filter2.collect()
		file(dada2errmodels) from ch_fastq_dada2errmodels.collect()
	 
	 	output:
		file ("*merged.rds") into ch_fastq_dada2cleanNmerge
	 	file ("dada2cleanNmerge_log") into ch_fastq_dada2cleanNmerge_log
	 	script:
	 	"""
	 	dada2cleanNmerge.R --verbose ${params.concatenate} --filterdir=. \
		--minoverlap=${params.minoverlap} --maxmismatch=${params.maxmismatch} --fwderrmodel=seq.dada2errmodels.fwd.errorates.rds \
		--reverrmodel=seq.dada2errmodels.rev.errorates.rds --fwdmark=${params.fwdmark} --revmark=${params.revmark} >dada2cleanNmerge_log 2>&1
	 	"""
	 }
	//dada2bimeras: filter bimeras, create the final table
	process dada2bimeras  {
	
	 	publishDir "${params.outdir}/dada2bimeras", mode: 'copy',
	 		saveAs: {filename -> 
	 		if (filename.indexOf(".bimeras.rds") == -1) "logs/$filename"
	 		else if(params.keepIntermediates) filename 
	 		else null}
	  
	 	input:
		file(dada2cleanNmerged) from ch_fastq_dada2cleanNmerge.collect()
	 
	 	output:
		file ("*bimeras.tsv.gz") into ch_fastq_dada2bimeras
	 	file ("dada2bimeras_log") into ch_fastq_dada2bimeras_log
	 	script:
	 	"""
	 	dada2bimeras.R --verbose --method=${params.method} --minab=${params.minab} --overab=${params.overab} \
		${params.bimeraoff} --prefix=dada2.cleaned.merged.bimeras --seqtabfile=dada2.cleaned.merged.rds \
		>dada2bimeras_log 2>&1
		"""
	 }

	//dada2idseq:
	process dada2idseq  {
	
	 	publishDir "${params.outdir}/dada2idseq", mode: 'copy',
	 		saveAs: {filename -> 
	 		if (filename.indexOf(".tsv.gz") == -1) "logs/$filename"
	 		else if(params.keepIntermediates) filename 
	 		else null}
	  
	 	input:
		file(dada2idseq) from ch_fastq_dada2bimeras.collect()
	 
	 	output:
		file ("*.fna") into ch_fastq_dada2taxonomy
	 	file ("dada2idseq_log") into ch_fastq_dada2taxonomy_log
	 	script:
	 	"""
	 	dada2idseq.R --verbose --fnafile=unique_seqs.fna --outtable asv_table.tsv.gz \
		dada2.cleaned.merged.tsv.gz > dada2idseq_log 2>&1
		"""
	 }

	//dada2taxonomy:
	process dada2taxonomy  {
	
	 	publishDir "${params.outdir}/dada2taxonomy", mode: 'copy',
	 		saveAs: {filename -> 
	 		if (filename.indexOf(".") == -1) "logs/$filename"
	 		else if(params.keepIntermediates) filename 
	 		else null}
	  
	 	input:
		file(dada2seqid) from ch_fastq_dada2idseq.collect()
	 
	 	output:
		file ("*.tsv") into ch_fastq_dada2taxonomy
	 	file ("dada2taxonomy_log") into ch_fastq_dada2taxonomy_log
	 	
		when (idtaxa_db == SILVA) {
		
		script:
	 	"""
	 	dada2taxonomy.R --verbose idtaxa_rdata /data/SILVA/SILVA_SSU_r132_March2018.RData  > dada2taxonomy_log 2>&1
		"""

	 } else if (idtaxa_db == GTDB) {

		script:
	 	"""
	 	dada2taxonomy.R --verbose idtaxa_rdata /data/GTDB/GTDB_r89-mod_June2019.RData  > dada2taxonomy_log 2>&1
		"""

	 } else if (idtaxa_db == UNITE) {
		
		script:
	 	"""
	 	dada2taxonomy.R --verbose idtaxa_rdata /data/UNITE/UNITE_v2019_July2019.RData  > dada2taxonomy_log 2>&1
		"""

	} else 
		error "Invalid taxonomy database name: ${idtaxa_db}"



workflow.onComplete {
	// from line 1781 of main.nf: On success try attach the multiqc report

	def mqc_report = null
	try {
		if (workflow.success) {
			mqc_report = multiqc_report.getVal()
			if (mqc_report.getClass() == ArrayList){
				log.warn "[nf-core/ampliseq] Found multiple reports from process 'multiqc', will use only one"
				mqc_report = mqc_report[0]
			}
		}
	} catch (all) {
		log.warn "[nf-core/ampliseq] Could not attach MultiQC report to summary email"
	}

	// Render the TXT template
	def engine = new groovy.text.GStringTemplateEngine()
	def tf = new File("$baseDir/assets/email_template.txt")
	def txt_template = engine.createTemplate(tf).make(email_fields)

	// Render the HTML template
	def hf = new File("$baseDir/assets/email_template.html")
	def html_template = engine.createTemplate(hf).make(email_fields)


	if(workflow.success){
		log.info "${c_purple}[nf-core/ampliseq]${c_green} Pipeline completed successfully${c_reset}"
	} else {
		checkHostname()
		log.info "${c_purple}[nf-core/ampliseq]${c_red} Pipeline completed with errors${c_reset}"
	}

}

}

def nfcoreHeader(){
	// Log colors ANSI codes
	c_reset = params.monochrome_logs ? '' : "\033[0m";
	c_dim = params.monochrome_logs ? '' : "\033[2m";
	c_black = params.monochrome_logs ? '' : "\033[0;30m";
	c_green = params.monochrome_logs ? '' : "\033[0;32m";
	c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
	c_blue = params.monochrome_logs ? '' : "\033[0;34m";
	c_purple = params.monochrome_logs ? '' : "\033[0;35m";
	c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
	c_white = params.monochrome_logs ? '' : "\033[0;37m";

return """${c_dim}----------------------------------------------------${c_reset}
	                                ${c_green},--.${c_black}/${c_green},-.${c_reset}
${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
	                                ${c_green}`._,._,\'${c_reset}
${c_purple}  nf-core/ampliseq v${workflow.manifest.version}${c_reset}
${c_dim}----------------------------------------------------${c_reset}
""".stripIndent()
}

def checkHostname(){
	def c_reset = params.monochrome_logs ? '' : "\033[0m"
	def c_white = params.monochrome_logs ? '' : "\033[0;37m"
	def c_red = params.monochrome_logs ? '' : "\033[1;91m"
	def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
	if(params.hostnames){
		def hostname = "hostname".execute().text.trim()
		params.hostnames.each { prof, hnames ->
			hnames.each { hname ->
				if(hostname.contains(hname) && !workflow.profile.contains(prof)){
					log.error "====================================================\n" +
							"  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
							"  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
							"  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
							"============================================================"
				}
			}
		}
	}
}

