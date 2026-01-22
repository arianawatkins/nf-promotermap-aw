#!/usr/bin/env nextflow

/*
========================================================================================
   Variant calling Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-promotermap
   Contact  : Eachan Johnson <user@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/

def pipeline_name = """\
         S C B I R   P R O M O T E R   M A P P I N G   P I P E L I N E
         =============================================================
         """.stripIndent()

if ( params.help ) {
   println """${pipeline_name}
         Nextflow pipeline to map Illumina sequences to bacterial genomes and call peaks.
         Usage:
            nextflow run scbirlab/nf-promotermap --sample_sheet <csv> --fastq_dir <dir> --control_label <text> [--inputs <dir> --mapper (bowtie2|minimap2)]
            nextflow run scbirlab/nf-promotermap -c <config-file>
         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 to be processed
            fastq_dir         Path to where FASTQ files are stored
            control_label     The bin ID of background controls
         Optional parameters (with defaults):  
            inputs            Directory containing inputs. Default: "./inputs".
            outputs           Directory to contain outputs. Default: "./outputs".
            trim_qual         Minimum base-call quality for trimming. Default: 5.
            min_length        Discard reads shorter than this number of bases after trimming. Default: 9.
            mapper            Alignment tool. Default: "bowtie2"
         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.control_label ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a control label (--control_label)")
}
if ( !params.from_sra ) {
   if ( !params.fastq_dir ) {
      throw new Exception("!!! PARAMETER MISSING: Please provide a path to fastq_dir")
   }
}

working_dir = params.outputs

log.info """${pipeline_name}
         inputs
            input dir.     : ${params.inputs}
            sample sheet   : ${params.sample_sheet}
            FASTQ dir.     : ${params.fastq_dir}
            input sample   : ${params.control_label}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         Aligner           : ${params.mapper}
         output            : ${params.outputs}
         """
         .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/
include { 
   annotate_nearest_gene;
   extract_peak_sequences;
   get_peak_coverage;
   get_per_base_coverage_within_peaks;
   gff2bed;
   make_peak_summary_table;
} from './modules/bedtools.nf'
include { 
   bowtie2_index; 
   bowtie2_align;
} from './modules/bowtie.nf'
include { 
   bam2wig;
   plot_peaks;
} from './modules/deeptools.nf'
include { 
   macs3 as MACS3_all_peaks;
} from './modules/macs3.nf'
include { 
   minimap_index; 
   minimap_align;
} from './modules/minimap.nf'
include { multiQC } from './modules/multiqc.nf'
include { 
   fetch_genome_from_NCBI; 
   fetch_FASTQ_from_SRA;
} from './modules/ncbi.nf'
include { fastQC } from './modules/qc.nf'
include { 
   sort_and_index_bam;
   plot_bamstats;
   SAMtools_stats;
   SAMtools_coverage;
   SAMtools_flagstat;
} from './modules/samtools.nf'
include { 
   trim_using_cutadapt; 
} from './modules/trimming.nf'
include { 
   get_peak_variance; 
} from './modules/variance.nf'

workflow {

   Channel.fromPath( 
      params.sample_sheet, 
      checkIfExists: true 
   )
      .splitCsv( header: true )
      .set { csv_ch }

   csv_ch
      .map { tuple( 
         it.sample_id,
         tuple( it.adapter_read1_5prime, it.adapter_read2_5prime ),
         tuple( it.adapter_read1_3prime, it.adapter_read2_3prime ),
      ) }
      .set { adapter_ch }  // sample_id, [adapt5], [adapt3]

   if ( params.from_sra ) {
      Channel.of( params.ncbi_api_key ).set { ncbi_api_key }
      csv_ch
         .map { tuple( it.sample_id, it.Run ) }
         .unique()
         | fetch_FASTQ_from_SRA
         | set { reads_ch }  // sample_id, reads
   }

   else {
      csv_ch
         .map { tuple( 
            it.sample_id,
            file( 
               "${params.fastq_dir}/*${it.fastq_pattern}*",
               checkIfExists: true
            ).sort()
         ) }
         .unique()
         .transpose()
         .groupTuple( by:0 )
         .set { reads_ch }  // sample_id, [reads]
   }

   csv_ch
      .map { tuple( it.sample_id, it.genome_accession ) }
      .unique()
      .set { genome_ch }  // sample_id, genome_acc

   csv_ch
      .map { tuple( it.sample_id, it.expt_id ) }
      .unique()
      .set { expt_ch }  // sample_id, expt_id

   csv_ch
      .map { tuple( it.sample_id, it.expt_id, it.bin_id ) }
      .filter { it[2] == "${params.control_label}" }
      .map { it[0..1] }
      .unique()
      .set { ctrl_ch }  // sample_id, expt_id

   csv_ch
      .map { tuple( it.sample_id, it.expt_id, it.bin_id ) }
      .filter { it[2] != "${params.control_label}" }
      .map { it[0..1] }
      .unique()
      .set { treat_ch }  // sample_id, expt_id
   
   // Purpose: Add an orientation channel containing promoter direction
   csv_ch
      .map { tuple( it.expt_id, it.promoter_direction) }
      .unique()
      .set { orientation_ch } // expt_id, promoter_direction

   /*
   ========================================================================================
      Processing
   ========================================================================================
   */

   reads_ch | fastQC 

   genome_ch
      .map { it[1] }  // genome_acc
      .unique()
      | fetch_genome_from_NCBI   // genome_acc, genome, gff
   
   genome_ch
      .map { it[1..0] }  // genome_acc, sample_id
      .combine( fetch_genome_from_NCBI.out, by: 0 )  // genome_acc, sample_id, genome, gff
      .map { tuple( it[1], it[-1] ) }  // sample_id, gff
      .set { genome_gff }

   trim_using_cutadapt(
      reads_ch.combine( adapter_ch, by: 0 ),  // sample_id, [reads], [adapt5], [adapt3]
      Channel.value( params.trim_qual ),
      Channel.value( params.min_length ),
   )  // sample_id, [reads]
   trim_using_cutadapt.out.main.set { trimmed }
   trim_using_cutadapt.out.logs.set { trim_logs }

   if ( params.mapper == "bowtie2" ) {

      fetch_genome_from_NCBI.out
         .map { it[0..1] } 
         | bowtie2_index
         | set { genome_idx0 }

   } 

   else if ( params.mapper == "minimap2" ) {

      minimap_index(
         fetch_genome_from_NCBI.out
            .map { it[0..1] },
         Channel.value( params.nanopore ),

      )
         | set { genome_idx0 }

   }

   else {
      error "Unsupported mapper: ${params.mapper}. Choose from: bowtie2, minimap2, star."
   }

   genome_idx0
      .combine( 
         genome_ch.map { it[1..0] }, 
         by: 0,
      )  // genome_acc, [genome_idx], sample_id
      .map { it[-1..0] }  // sample_id, [genome_idx], genome_acc
      .set { genome_idx }

   trimmed
      .combine( genome_idx, by: 0 )  // sample_id, [reads], [genome_idx], genome_acc
      .set { pre_mapper }

   if ( params.mapper == "bowtie2" ) {

      pre_mapper
         | bowtie2_align
         | set { mapped_reads }

   } 
   
   else if ( params.mapper == "minimap2" ) {

      minimap_align(
         pre_mapper,
         Channel.value( params.nanopore ),

      )
         | set { mapped_reads }

   }

   mapped_reads.main
      | sort_and_index_bam
      | (
         SAMtools_stats
         & SAMtools_coverage
         & SAMtools_flagstat
         & bam2wig
      )
   SAMtools_stats.out | plot_bamstats

   fetch_genome_from_NCBI.out  // genome_acc, genome, gff
      .map { tuple( it[0], it[2] ) }  // genome_acc, gff
      | gff2bed

   gff2bed.out
      .combine( 
         genome_ch.map { it[1..0] }, 
         by: 0,
      )  // genome_acc, bed, sample_id
      .map { it[-1..1] }  // sample_id, bed
      .combine( treat_ch, by: 0 )  // sample_id, bed, expt_id
      .map { it[-1..1] }  // expt_id, bed
      .unique()
      .set { expt2gffbed }

   sort_and_index_bam.out
      .combine( bam2wig.out, by: 0 )
      .combine( expt_ch, by: 0 )     // sample_id, bam_bai, bw, expt_id
      .map { tuple( it[0], it[-1], it[1], it[2] ) }  // sample_id, expt_id, bam_bai, bw
      .set { alignments }

   alignments
      .map { it[0..-2] }
      .join( treat_ch, by: [0, 1] )  // sample_id, expt_id, bam_bai_treat
      .set { treat_alignments }

   alignments
      .map { it[0..-2] }
      .join( ctrl_ch, by: [0, 1] )  // sample_id, expt_id, bam_bai_ctrl
      .set { ctrl_alignments }

   treat_alignments
      .map { tuple( it[1], it[2] ) }  // expt_id, bam_bai
      .unique()
      .groupTuple( by: 0 )  // expt_id, [bam_bai, ..]
      .combine(
         ctrl_alignments.map { tuple( it[1], it[-1] ) }.unique(),  // expt_id, bam_bai_ctrl
         by: 0,
      )  // expt_id, [bam_bai, ..], bam_bai_ctrl
      .set { expt2bams }
   
   // Purpose: Combine with orientation channel to attach promoter direction 
   expt2bams
      .combine( orientation_ch, by: 0 ) // expt_id, promoter_orientation
      | MACS3_all_peaks

   alignments
      .map { tuple( it[1], it[3] ) } // expt_id, bw
      .unique()
      .groupTuple( by: 0 )  // expt_id, [bw, ..]
      .set { expt2bigwig }
   expt2bigwig
      .combine( expt2gffbed, by: 0 )
      | plot_peaks

   MACS3_all_peaks.out.peaks  // expt_id, bed
      .combine(
         expt2bams.map { it[0..1] }.unique(),
         by: 0,
      )
      | (
         get_peak_coverage
         & get_per_base_coverage_within_peaks
      )

   get_peak_coverage.out.tsv
      | get_peak_variance
   
   MACS3_all_peaks.out.summits  // expt_id, bed
      .combine( 
         expt2gffbed,
         by: 0 
      )  // expt_id, bed, gff
      | annotate_nearest_gene

   MACS3_all_peaks.out.peaks  // expt_id, bed
      .combine( 
         fetch_genome_from_NCBI.out  // genome_acc, genome, gff
            .map { tuple( it[0], it[1] ) }  // genome_acc, genome
            .unique()
            .combine( 
               genome_ch.map { it[1..0] }, 
               by: 0,
            )  // genome_acc, genome, sample_id
            .map { it[-1..1] }  // sample_id, genome
            .unique()
            .combine( treat_ch, by: 0 )  // sample_id, genome, expt_id
            .map { it[-1..1] }  // expt_id, genome
            .unique(),
         by: 0,
      )  // expt_id, bed, genome
      | extract_peak_sequences

   annotate_nearest_gene.out
      .combine(
         extract_peak_sequences.out.tsv,
         by: 0,
      )
      .combine(
         get_peak_coverage.out.tsv,
         by: 0,
      )
      | make_peak_summary_table

   trim_logs
      .map { it[1] }
      .concat(
         fastQC.out.multiqc_logs,
         mapped_reads.logs,
         SAMtools_coverage.out.map { it[1] },
         SAMtools_flagstat.out.map { it[1] },
         SAMtools_stats.out.map { it[1] },
      )
      .flatten()
      .unique()
      .collect()
      | multiQC

}


/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/