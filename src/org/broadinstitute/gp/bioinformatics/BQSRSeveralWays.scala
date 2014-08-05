
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk.{CommandLineGATK, BaseRecalibrator, PrintReads}
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction


class BQSRSeveralWays extends QScript {
  qscript =>

  @Input(shortName = "i", required = true, doc = "Input Bam to be recalibrated") var inputFile: List[File] = _

  @Input(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")

  @Input(required = true, doc = "Intervals to use for BQSR") var Intervals: List[File] = _

  @Input(shortName = "v", required = true, doc = "Standard Variants to avoid") var Variants: List[File] = _

  @Input(required = true, doc = "Standard Variants to avoid") var NA12878Variants: List[File] = _

  @Argument(shortName = "sc", required = false, doc = "scatter count") var scatterCount: Int = 1
  @Argument(shortName = "P", required = false, doc="probability of using base in quality check") var probability:Double = _

  def script() {
    for (file <- inputFile)
      for (intervals <- Map("Intervals" -> Intervals, "NoIntervals" -> Nil))
        for (variants <- Map("GenericVariants" -> Variants, "NA12878_variants" -> NA12878Variants))
          BQSRFile(file, intervals, variants)
  }

  trait CommonArguments extends CommandLineGATK {
    reference_sequence = referenceFile
    useOriginalQualities = true
  }

  def BQSRFile(file: File, intervals: (String, List[File]), variants: (String, List[File])) {

    val bqsr = new BaseRecalibrator with CommonArguments
    bqsr.knownSites = variants._2
    bqsr.out = swapExt(file, ".bam", "." + intervals._1 + "." + variants._1 + ".tbl")
    bqsr.intervals ++= intervals._2
    bqsr.scatterCount = qscript.scatterCount
    bqsr.input_file :+= file

    add(bqsr)

    val pr = new PrintReads with CommonArguments
    pr.out = swapExt(bqsr.out, "tbl", "bam")
    pr.BQSR = bqsr.out
    pr.input_file :+= file

    add(pr)

    val cbem2 = new CollectBamErrorMetrics
    cbem2.bam = pr.out
    cbem2.out = swapExt(pr.out, ".bam", "")
    cbem2.variants = NA12878Variants(0)
    cbem2.reference = referenceFile

    add(cbem2)
  }

  class CollectBamErrorMetrics  extends JavaCommandLineFunction {

    @Input var bam: File = _
    @Input var variants: File = _
    @Output var out: File = _
    @Argument var reference: File = _

    val nist_interval = "/seq/tng/giab/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.interval_list"
    this.jarFile = "/seq/tng/farjoun/temp/CollectBamErrorMetrics.jar"
    var tempDir: List[File] = List(new File("/local/scratch/"), new File("/seq/picardtemp3"))
    this.memoryLimit = Option(8)

    this.jobResourceRequests :+= "virtual_free=%dM".format(9000)

    override def commandLine = super.commandLine +
      required("I=", bam, spaceSeparated = false) +
      required("O=", out, spaceSeparated = false) +
      required("V=", variants, spaceSeparated = false) +
      required("R=", reference, spaceSeparated = false) +
      optional("P=", probability, spaceSeparated = false)+
      repeat("L=", List(nist_interval,"20"), spaceSeparated = false)
  }

}