#!/usr/bin/env nextflow
/*
========================================================================================
    HiTE
========================================================================================
    Github : https://github.com/CSU-KangHu/HiTE
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

def echo(text) {
  println text
}

def helpMessage() {
    log.info"""
    HiTE - Nextflow PIPELINE (v$workflow.manifest.version)
    =================================
    Usage:
    The typical command is as follows:
    nextflow run main.nf --genome ../demo/genome.fa --thread 48 --outdir ../demo/test --plant 0

    Mandatory arguments:
      --genome      Genome assembly path (format: fasta, fa, and fna)
      --outdir      Output directory; It is recommended to use a new directory to avoid automatic deletion of important files.
    General options:
      --chunk_size                      The chunk size of large genome, default = [ 400 MB ]
      --plant                           Is it a plant genome, 1: true, 0: false. default = [ 1 ]
      --recover                         Whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false. default = [ 0 ]
      --miu                             The neutral mutation rate (per bp per ya). default = [ 1.3e-8 ]
      --classified                      Whether to classify TE models, HiTE uses RepeatClassifier from RepeatModeler to classify TEs, 1: true, 0: false. default = [ 1 ]
      --debug                           Open debug mode, and temporary files will be kept, 1: true, 0: false. default = [ 0 ]
      --flanking_len                    The flanking length of candidates to find the true boundaries, default = [ 50 ]
      --fixed_extend_base_threshold     The length of variation can be tolerated during pairwise alignment, default = [ 1000 ]
      --tandem_region_cutoff            Cutoff of the candidates regarded as tandem region, default = [ 0.5 ]
      --max_repeat_len                  The maximum length of a single repeat, default = [ 30000 ]
      --chrom_seg_length                The length of genome segments, default = [ 500000 ]
    """.stripIndent()
}

def printSetting() {
    log.info"""
    ====================================Parameter settings========================================
      [Setting] Reference sequences / assemblies path = [ $params.genome ]
      [Setting] The chunk size of large genome = [ $params.chunk_size ] MB
      [Setting] Is plant genome = [ $params.plant ]
      [Setting] recover = [ $params.recover ]
      [Setting] Is classified = [ $params.classified ]
      [Setting] The neutral mutation rate (per bp per ya) = = [ $params.miu ]
      [Setting] debug = [ $params.debug ]
      [Setting] Output Directory = [ $params.outdir ]

      [Setting] Fixed extend bases threshold = [ $params.fixed_extend_base_threshold ]
      [Setting] Flanking length of TE = [ $params.flanking_len ]
      [Setting] Cutoff of the repeat regarded as tandem sequence = [ $params.tandem_region_cutoff ]
      [Setting] Maximum length of TE = [ $params.max_repeat_len ]
      [Setting] The length of genome segments = [ $params.chrom_seg_length ]
    """.stripIndent()
}

echo "Welcome to HiTE nextflow!"

// check input =====================================================================
// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (!params.genome){
    exit 1, "--genome option not specified!"
}

if (!params.outdir){
    exit 1, "--outdir option not specified!"
}

printSetting()
projectDir = workflow.projectDir
ch_module = "${projectDir}/module"
ch_classification = "${projectDir}/classification"
tools_module = "${projectDir}/tools"
HSDIR = "${projectDir}/bin/HelitronScanner/TrainingSet"
HSJAR = "${projectDir}/bin/HelitronScanner/HelitronScanner.jar"
sh_dir = "${projectDir}/module"
member_script_path = "${projectDir}/tools/make_fasta_from_blast.sh"
subset_script_path = "${projectDir}/tools/ready_for_MSA.sh"
lib_module = "${projectDir}/library"
ch_EAHelitron = "${projectDir}/bin/EAHelitron-master"
ch_ltrfinder = "${projectDir}/bin/LTR_FINDER_parallel-master"
ch_protein = "${projectDir}/library/RepeatPeps.lib"
genome_name = file(params.genome).getName()
out_genome = "${params.outdir}/${genome_name}"

filePrefix = genome_name.substring(0, genome_name.lastIndexOf('.'))
out_genome_rename = "${params.outdir}/${filePrefix}.rename.fa"

//parameters of HiTE
tmp_output_dir = "${params.outdir}"
chrom_seg_length = "${params.chrom_seg_length}"
chunk_size = "${params.chunk_size}"
fixed_extend_base_threshold = "${params.fixed_extend_base_threshold}"
max_repeat_len = "${params.max_repeat_len}"
flanking_len = "${params.flanking_len}"
tandem_region_cutoff = "${params.tandem_region_cutoff}"
recover = "${params.recover}"
plant = "${params.plant}"
classified = "${params.classified}"
domain = "${params.domain}"
debug = "${params.debug}"
miu = "${params.miu}"
//parameters of Evaluation


// Check all tools work well
process EnvCheck {
    tag "envcheck"
    errorStrategy 'terminate'

    label 'process_low'

    input:
    val dependency

    output:
    stdout

    script:
    """
    echo "---Check if ${dependency} is installed."
    command -v ${dependency} >/dev/null 2>&1 || { echo >&2 "${dependency} is required but not installed.  Aborting."; exit 1; }
    """
}

// Check all tools work well
process splitGenome {
    tag "${ref}"

    label 'process_low'

    input:
    path ref

    output:
    path "genome.cut*.fa"
    

    script:
    cores = task.cpus
    ref_name = ref.getName()
    """
    cp ${ref} ${tmp_output_dir}
    python3 ${ch_module}/split_genome_chunks.py \
     -g ${tmp_output_dir}/${ref_name} --tmp_output_dir ${tmp_output_dir} --chrom_seg_length ${chrom_seg_length} \
     --chunk_size ${chunk_size}
    cp ${tmp_output_dir}/genome.cut*.fa ./
    """
}

process coarseBoundary {
    tag "${cut_ref}"

    label 'process_high'

    input:
    path cut_ref

    output:
    path "longest_repeats_*.flanked.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]
    """
    python3 ${ch_module}/coarse_boundary.py \
     -g ${cut_ref} --tmp_output_dir ${tmp_output_dir} \
     --fixed_extend_base_threshold ${fixed_extend_base_threshold} \
     --max_repeat_len ${max_repeat_len} --thread ${cores} \
     --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
     --ref_index ${ref_index} -r ${out_genome} --recover ${recover} --debug ${debug}

    ## Since nextflow will look for output files in the work directory, we need to copy the script output files to the work directory.
    cp ${tmp_output_dir}/longest_repeats_${ref_index}.flanked.fa ./
    """
}

def groupTuple = {
    input_ch.collectFileGroups { file -> file.baseName.replaceAll("[^\\d]", "") }
}

process TIR {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(lrf)
    path ltrs

    output:
    path "confident_tir_*.fa"


    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]

    script:
    """
    python3 ${ch_module}/judge_TIR_transposons.py \
    -g ${cut_ref} --seqs ${lrf} --confident_ltr_cut_path ${ltrs} \
    -t ${cores} --TRsearch_dir ${tools_module}  \
    --tmp_output_dir ${tmp_output_dir} \
    --flanking_len ${flanking_len} --tandem_region_cutoff ${tandem_region_cutoff} \
    --ref_index ${ref_index} --member_script_path ${member_script_path} \
    --subset_script_path ${subset_script_path} \
    --plant ${plant} --recover ${recover} --debug ${debug}

    cp ${tmp_output_dir}/confident_tir_${ref_index}.fa ./
    """
    
}

process Helitron {
    tag "${cut_ref}"

    label 'process_high'

    input:
    tuple path(cut_ref), path(lrf)


    output:
    path "confident_helitron_*.fa"

    script:
    cores = task.cpus
    (full, ref_index) = (cut_ref =~ /genome.cut(\d+)\.fa/)[0]

    script:
    """
    python3 ${ch_module}/judge_Helitron_transposons.py \
    -g ${cut_ref} --seqs ${lrf} \
    -t ${cores} --tmp_output_dir ${tmp_output_dir} \
    --HSDIR ${HSDIR} --HSJAR ${HSJAR} --sh_dir ${sh_dir} \
    --member_script_path ${member_script_path} --subset_script_path ${subset_script_path} \
    --flanking_len ${flanking_len} --debug ${debug} \
    --ref_index ${ref_index} --recover ${recover} \

    cp ${tmp_output_dir}/confident_helitron_${ref_index}.fa ./
    """
}

process OtherTE {
    label 'process_high'

    input:
    tuple path(cut_ref), path(lrf)


    output:
    path "confident_other_*.fa"

    script:
    cores = task.cpus
    (full, lrf_index) = (lrf =~ /longest_repeats_(\d+)\.flanked\.fa/)[0]
    """
    python3 ${ch_module}/judge_Other_transposons.py \
     -g ${cut_ref} --member_script_path ${member_script_path} --subset_script_path ${subset_script_path} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --ref_index ${lrf_index} --library_dir ${lib_module} --recover ${recover}

    cp ${tmp_output_dir}/confident_other_${lrf_index}.fa ./
    """
}

process LTR {
    tag "${ref}"

    label 'process_high'

    input:
    path ref

    output:
    path "confident_ltr_cut.fa"

    script:
    cores = task.cpus
    """
    python3 ${ch_module}/judge_LTR_transposons.py \
     -g ${ref} --ltrfinder_home ${ch_ltrfinder} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --recover ${recover} --miu ${miu}

    cp ${tmp_output_dir}/confident_ltr_cut.fa ./
    """
}

process UnwrapNested {
    tag "${ref}"

    label 'process_high'

    input:
    path ref
    path ltr
    path tir
    path helitron
    path other


    output:
    path "confident_TE.fa"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/remove_nested.py \
     -g ${ref} --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_other ${other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --global_flanking_filter ${global_flanking_filter} \
     --remove_nested ${remove_nested} --test_home ${ch_module}

    cp ${tmp_output_dir}/confident_TE.fa ./
    """
}

process BuildLib {
    label 'process_high'

    input:
    path ltr
    path tir
    path helitron
    path other

    output:
    path "confident_TE.cons.fa"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/get_nonRedundant_lib.py \
     --confident_ltr_cut ${ltr} \
     --confident_tir ${tir} \
     --confident_helitron ${helitron} \
     --confident_other ${other} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --test_home ${ch_module}

    cp ${tmp_output_dir}/confident_TE.cons.fa ./
    """
}

process ClassifyLib {
    label 'process_high'

    input:
    path lib

    output:
    path "${lib}.classified"


    script:
    cores = task.cpus
    """
    python3 ${ch_module}/get_classified_lib.py \
     --confident_TE_consensus ${lib} \
     -t ${cores} --tmp_output_dir ${tmp_output_dir} \
     --classified ${classified} --domain ${domain} --TEClass_home ${ch_classification} \
     --protein_path ${ch_protein} \
     --debug ${debug}
    """
}

process CleanLib {
    label 'process_low'

    input:
    path lib

    output:
    stdout


    script:
    """
    python3 ${ch_module}/clean_lib.py \
     --tmp_output_dir ${tmp_output_dir} \
     --debug ${debug}
    """
}

process BM_RM2 {
    tag "${TE}"

    label 'process_high'

    input:
    path TE
    path curatedLib
    path rm2_script

    output:
    path "${TE}.out"
    path "res.log"


    script:
    cores = task.cpus
    """
    RepeatMasker -lib ${curatedLib} -nolow -pa ${cores} ${TE}
    mkdir rm2_test
    cd rm2_test && rm -rf * && sh ../${rm2_script} ../${TE}.out > ../res.log
    """
}

process BM_EDTA {
    tag "${TE}"

    label 'process_high'

    input:
    path TE
    path curatedLib
    val reference
    val EDTA_home

    output:
    path "repbase.out"
    path "HiTE.out"
    path "HiTE.out.*"


    script:
    cores = task.cpus
    """
    RepeatMasker -e ncbi -pa ${cores} -q -no_is -norna -nolow -div 40 -lib ${curatedLib} -cutoff 225 ${tmp_output_dir}/genome.rename.fa
    mv ${tmp_output_dir}/genome.rename.fa.out repbase.out
    cp repbase.out ${tmp_output_dir}/

    RepeatMasker -e ncbi -pa ${cores} -q -no_is -norna -nolow -div 40 -lib ${TE} -cutoff 225 ${tmp_output_dir}/genome.rename.fa
    mv ${tmp_output_dir}/genome.rename.fa.out HiTE.out
    cp HiTE.out ${tmp_output_dir}/

    perl ${EDTA_home}/lib-test.pl -genome ${tmp_output_dir}/genome.rename.fa -std repbase.out -tst HiTE.out -cat Total
    """
}

// allow multi-thread
process test {
    tag "${cut_reference}"

    label 'process_low'

    input:
    val cut_reference

    output:
    stdout

    script:
    cores = task.cpus


    """
    # genome.fa.cut0.fa
    echo "${cut_reference}"
    """

}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

workflow {
    if (!params.skip_HiTE) {
        if ( !file(params.genome).exists() )
                exit 1, "genome reference path does not exist, check params: --genome ${params.genome}"


            // dependency check
            dependencies = Channel.from(params.dependencies)
            EnvCheck(dependencies) | view { "$it" }

            // split genome into chunks
            Channel.fromPath(params.genome, type: 'any', checkIfExists: true).set{ ch_g }

            //LTR identification
            ch_ltrs = LTR(ch_g)
            //test(ch_ltrs) | view { "$it" }

            cut_genomes = splitGenome(ch_g)
            ch_cut_g = cut_genomes.flatten()

            // coarse-grained Boundary identification
            longest_repeats = coarseBoundary(ch_cut_g)

            // merge files with the same index 
            ch_cut_g_map = ch_cut_g.map {file -> tuple(file.baseName.replaceAll("[^\\d]", ""), file)}
            longest_repeats_map = longest_repeats.map {file -> tuple(file.baseName.replaceAll("[^\\d]", ""), file)}
            merged_channel = ch_cut_g_map.combine(longest_repeats_map, by: 0).map {file -> tuple(file[1], file[2])}
            //test(merged_channel) | view { "$it" }

            //TIR identification
            ch_tirs = TIR(merged_channel, ch_ltrs).collectFile(name: "${params.outdir}/confident_tir.fa")

            //Helitron identification
            ch_h = Helitron(merged_channel).collectFile(name: "${params.outdir}/confident_helitron.fa")

            //Other identification
            ch_o = OtherTE(merged_channel).collectFile(name: "${params.outdir}/confident_other.fa")
            //test(ch_o) | view { "$it" }

            //Unwrap nested TE
            //ch_TE = UnwrapNested(params.genome, ch_ltrs, ch_tirs, ch_h, ch_o)
            //test(ch_TE) | view { "$it" }

            //Build TE library
            ch_lib = BuildLib(ch_ltrs, ch_tirs, ch_h, ch_o)
            //test(ch_lib) | view { "$it" }

            //Classify TE library
            ch_lib.splitFasta(by: params.classify_chunk_size, file:true).set { ch_fasta }
            ch_classified_lib = ClassifyLib(ch_fasta)
            ch_final = ch_classified_lib.collectFile(name: "${params.outdir}/confident_TE.cons.fa.classified")
            //test(ch_lib) | view { "$it" }

            //Clean TE library
            CleanLib(ch_final)
    }
        
    
    if (params.BM_RM2){
        if (!params.species)
            exit 1, "--BM_RM2 is set as true, but there is no --species specified! Choose from test, dmel, rice, cb, zebrafish, and maize."
        
        if (params.skip_HiTE)
            Channel.fromPath("${params.outdir}/confident_TE.cons.fa.classified", type: 'any', checkIfExists: true).set{ ch_final }

        if (params.species == "dmel"){
            lib_path = "${lib_module}/drorep.ref"
        } else if (params.species == "rice"){
            lib_path = "${lib_module}/oryrep.ref"
        } else if (params.species == "cb"){
            lib_path = "${lib_module}/cbrrep.ref"
        } else if (params.species == "zebrafish"){
            lib_path = "${lib_module}/zebrep.ref"
        } else if (params.species == "maize"){
            lib_path = "${lib_module}/maize.ref"
        } else{
            lib_path = "${lib_module}/test.ref"
        }

        Channel.fromPath("${lib_path}", type: 'any', checkIfExists: true).set{ curatedLib }
        Channel.fromPath("${projectDir}/bin/get_family_summary_paper.sh", type: 'any', checkIfExists: true).set{ rm2_script }
        (ch_out, ch_log) = BM_RM2(ch_final, curatedLib, rm2_script)
        ch_log.collectFile(name: "${params.outdir}/BM_RM2.log")
    }

    if (params.BM_EDTA){
        if (!params.species)
            exit 1, "--BM_EDTA is set as true, but there is no --species specified! Choose from test, dmel, rice, cb, zebrafish, and maize."
        
        if (!params.EDTA_home)
            exit 1, "--BM_EDTA is set as true, but there is no --EDTA_home specified!"

        if (params.skip_HiTE)
            Channel.fromPath("${params.outdir}/confident_TE.cons.fa.classified", type: 'any', checkIfExists: true).set{ ch_final }

        if (params.species == "dmel"){
            lib_path = "${lib_module}/drorep.ref"
        } else if (params.species == "rice"){
            lib_path = "${lib_module}/oryrep.ref"
        } else if (params.species == "cb"){
            lib_path = "${lib_module}/cbrrep.ref"
        } else if (params.species == "zebrafish"){
            lib_path = "${lib_module}/zebrep.ref"
        } else if (params.species == "maize"){
            lib_path = "${lib_module}/maize.ref"
        } else if (params.species == "test"){
            lib_path = "${lib_module}/test.ref"
        }

        Channel.fromPath("${lib_path}", type: 'any', checkIfExists: true).set{ curatedLib }
        (ch_rep_out,ch_hi_out,ch_report) = BM_EDTA(ch_final, curatedLib, out_genome_rename, params.EDTA_home)
        ch_report.collectFile(name: "${params.outdir}/BM_EDTA.log")
    }
       

}


/*
========================================================================================
    THE END
========================================================================================
*/
