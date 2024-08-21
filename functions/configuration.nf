/*
 * Write a log message summarising how the pipeline is configured and the
 * locations of reference files that will be used.
 */

def displayParameters(params)
{
    params.with
    {
        log.info "SLX ID: ${slxid}"
        log.info "species: ${species}"
        log.info "reference directory: ${referenceDir}"
        log.info "cellranger directory: ${cellrangerDir}"
    }
}