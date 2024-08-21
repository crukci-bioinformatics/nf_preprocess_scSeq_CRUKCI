#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

// include functions
include { displayParameters } from './functions/configuration'

// Download the reference data
process downloadReferences {
    publishDir "${launchDir}/references", mode: 'link'
    errorStrategy 'retry'
    maxRetries 5

    when:
    params.species != null && params.referenceDir == null

    input:
        val species

    output:
        path "refdata-gex*", emit: referenceDir

    script:
        """
        if [[ "${species}" == "homo_sapiens" ]]; then
            curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
            tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz
            rm refdata-gex-GRCh38-2024-A.tar.gz
        elif [[ "${species}" == "mus_musculus" ]]; then
            curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
            tar -xzvf refdata-gex-GRCm39-2024-A.tar.gz
            rm refdata-gex-GRCm39-2024-A.tar.gz
        else
            echo "Species ${species} not supported"
            exit 1
        fi
        """
}

process downloadData {
    publishDir "${launchDir}/${slxid}", mode: 'link'
    input:
        val slxid

    output:
        path outdir 

    script:
        outdir = "fastq"
        """
        curl -o clarity-tools.jar http://internal-bioinformatics.cruk.cam.ac.uk/software/clarity-tools.jar 
        java -jar clarity-tools.jar --library ${slxid}
        mv ${slxid} ${outdir}
        """
}

process renameFastq {
    input:
        path rawdir

    output:
        path 'fastq_renamed/**', emit: fastqFiles
        path 'fastq_renamed', emit: fastqDir

    shell:
        srchTerm = 'SLX-[0-9]*.\\([A-Z0-9]*\\).[A-Z0-9]*.s_\\([1-4]\\).\\([ri]\\)_\\([12]\\).fq.gz'
        replTerm = '\\1_S\\2_L001_\\u\\3\\4_001.fastq.gz'
        ''' 
        mkdir -p fastq_renamed
        for fqFile in !{rawdir}/*[12].fq.gz; do
            newName=`basename \${fqFile} | sed "s/!{srchTerm}/!{replTerm}/"`
            ln -rsT ${fqFile} fastq_renamed/${newName}
        done
        '''
}

process cellRangerCount {
    publishDir "${launchDir}/${params.slxid}", mode: 'link'
    cpus 16
    memory 32.GB
    time 24.hour

    input:
        tuple val(barcode), path(fastqDir), path(referenceDir)

    output:
        path outputDir, emit: crDir
        tuple val(barcode), path(metricsFile), emit: metricsFile
        tuple val(barcode), path(webSummary), emit: webSummary

    script:
        outputDir = barcode
        metricsFile = "${outputDir}/outs/metrics_summary.csv"
        webSummary = "${outputDir}/outs/web_summary.html"
        """
        cellranger count \
            --id=${barcode} \
            --fastqs=${fastqDir} \
            --sample=${barcode} \
            --transcriptome=${referenceDir} \
            --create-bam ${params.createBam} \
            --localcores=16 \
            --localmem=32
        """
}

process collectWebSummaries {
    publishDir "${launchDir}/${params.slxid}/reports", mode: 'link'

    executor 'local'

    input:
        tuple val(barcode), path(inHtml)

    output:
        path outHtml

    script:
        outHtml = "${barcode}.web_summary.html"
        """
        cp ${inHtml} ${outHtml}
        """
}

process modifyMetrics {
    executor 'local'

    input:
        tuple val(barcode), path(metricsFile)

    output:
        path modifiedMetrics 

    script:
        modifiedMetrics = "${barcode}.appended.csv"
        """
        sed -e '1s/^/Barcode,/' -e '2s/^/${barcode},/' ${metricsFile} > ${modifiedMetrics}
        """
}

process plotMetrics {
    publishDir "${launchDir}/${params.slxid}/reports", mode: 'link'

    executor 'local'

    input:
        path metricsFile

    output:
        path plotFile

    script:
        plotFile = "summary_metrics.pdf"
        template 'PlotMetrics.R'

}

displayParameters(params)

workflow {
    speciesChannel = params.species
    slxChannel = params.slxid

    downloadReferences(speciesChannel)
    downloadData(slxChannel)
    renameFastq(downloadData.out)
    
    firstPartChannel = renameFastq.out.fastqFiles
        .flatten() 
        .map { file ->
            def basename = file.getName()
            def firstPart = basename.split("_")[0]
            return firstPart
        }
        .unique()

    // if referencDir is not set in the params the refChannel comes from
    // downloadReferences, otherwise it is set to the referenceDir
    refChannel = params.referenceDir == null ? downloadReferences.out.referenceDir : channel.fromPath(params.referenceDir)

    combChannel = firstPartChannel
        .combine(renameFastq.out.fastqDir)
        .combine(refChannel)

    cellRangerCount(combChannel)
    collectWebSummaries(cellRangerCount.out.webSummary)
    modifiedMetricsChannel = modifyMetrics(cellRangerCount.out.metricsFile)
    collatedChannel = modifiedMetricsChannel
        .collectFile(name: 'collated_metrics_file.csv', skip: 1, keepHeader: true, storeDir: "${params.slxid}/reports")
    plotMetrics(collatedChannel)
}
