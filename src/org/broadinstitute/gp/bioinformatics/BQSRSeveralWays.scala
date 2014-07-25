
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.{CommandLineGATK, BaseRecalibrator, PrintReads}


class BQSRSeveralWays extends QScript {
  qscript =>

  @Input(shortName = "i", required = true, doc = "Input Bam to be recalibrated") var inputFile: List[File] = _

  @Argument(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  @Argument(required = true, doc = "Intervals to use for BQSR") var Intervals: File = _

  @Argument(shortName = "v", required = true, doc = "Standard Variants to avoid") var Variants: List[File] = _

  @Argument(required = true, doc = "Standard Variants to avoid") var NA12878Variants: List[File] = _

  @Argument(shortName = "sc", doc = "scatter count") var scatterCount: Int = 1

  def script() {
    for (file <- inputFile)
      for (intervals <- Map("Intervals" -> Intervals, "NoIntervals" -> null))
        for (variants <- Map("GenericVariants" -> Variants, "NA12878_variants" -> NA12878Variants))
          BQSRFile(file, intervals, variants)
  }


  def BQSRFile(file: File, intervals: (String, File), variants: (String, List[File])) {

    trait CommonArguments extends CommandLineGATK {
      this.reference_sequence = referenceFile
      this.intervals = intervals
      this.input_file :+= file
      this.useOriginalQualities = true
    }

    val bqsr = new BaseRecalibrator with CommonArguments
    bqsr.knownSites = variants._2
    bqsr.out = swapExt(file, ".bam", "." + intervals._1 + "." + variants._1 + ".tbl")
    bqsr.scatterCount

    add(bqsr)

    val pr = new PrintReads with CommonArguments
    pr.out = swapExt(bqsr.out, "tbl", "bam")
    pr.BQSR = bqsr.out

    add(pr)
  }

}