#!/usr/bin/env nextflow

params.fqs = "/projects/b1059/projects/20180119-RNA-Seq-Analysis/Ye/subsample_data/results/**.gz"
params.transcriptome = "/projects/b1059/projects/20180119-RNA-Seq-Analysis/Ye/SEmRNA-seq-nf/test_data/c.elegans.cdna.ncrna.fa"
params.fragment_len = '250'
params.fragment_sd = '50'
params.bootstrap = '100'


log.info """\
         R N A S E Q - N F   P I P E L I N E  (Kallisto plus QC)
         ===================================
         transcriptome: ${params.transcriptome}
         fqs          : ${params.fqs}
         fragment_len : ${params.fragment_len}
         fragment_sd  : ${params.fragment_sd}
         bootstrap    : ${params.bootstrap}

         """
         .stripIndent()

transcriptome_file = file(params.transcriptome)

if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"


Channel
    .fromFilePairs( params.fqs, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.fqs}" }
    .set { reads }

process kal_index {

    input:
        file transcriptome_file

    output:
        file "transcriptome.index" into transcriptome_index

    script:
    //
    // Kallisto mapper index
    //
    """
    kallisto index -i transcriptome.index ${transcriptome_file}
    """
}

process kal_mapping {

    tag "reads: $name"

    input:
        file index from transcriptome_index
        set val(name), file(fq) from reads

    output:
        file "kallisto_${name}" into kallisto_out_dirs

    script:
    //
    // Kallisto tools mapper
    //
    def single = fq instanceof Path
    if( !single ){
        """
        mkdir kallisto_${name}
        kallisto quant --bootstrap ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${fq}
        """
    }
    else {
        """
        mkdir kallisto_${name}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} --bootstrap ${params.bootstrap} -i ${index} -t ${task.cpus} -o kallisto_${name} ${fq}
        """
    }
}

