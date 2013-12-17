
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import scala.annotation.tailrec
import org.broadinstitute.gp.bioinformatics.utils.RScriptFunction

class ExtractSubPopulationMAF extends QScript {
  @Input(shortName = "vcfSrcDir", doc = "The source directory of VCFs to process", required = false)
  var vcfSrcDir: File = "/humgen/1kg/DCC/ftp/release/20110521/"

  @Output(shortName = "out", doc = "The output VCF", required = true)
  var out: File = _

  @Input(shortName = "panel", doc = "The file containing the population description of the different samples", required = false)
  var SamplePanel: File = "/humgen/1kg/DCC/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel"

  @Argument(shortName = "pop", doc = "The list of target populations", required = false)
  var TargetPopulations: List[String] = _

  @Argument(shortName = "scatterCount", doc = "The number of ways to scatter this job", required = false)
  val scatterCount: Int = 25

  @Argument(shortName = "L", required = false)
  var interval: String = _

  @Argument(shortName = "R", required = false)
  var reference: File = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.reference_sequence = reference
    this.intervalsString :+= interval
  }

  class GetPopulationSamples extends RScriptFunction {
    @Argument var pop: String = _
    @Input var SamplePanel: File = _
    @Output var SampleList: File = _


    lazy val rscript: String = s"""
        |read.table($SamplePanel ,col.names=c("SAMPLE","SUB.POP","POP","SEQUENCING.CENTER","SEQUENCING.CENTER.2"),fill=TRUE)->data
        |write.table(file=$SampleList,subset(data,POP="AFR")$$SAMPLE)
        """.stripMargin

  }

  @tailrec
  final def findallVCF(f: Seq[File], acc:Seq[File]=Nil): Seq[File] = {

    val dirs=f.filter(_.isDirectory).filter(x=> !x.getName.startsWith("."))
    val newFiles=f.filter(x=> x.isFile & (x.getName.endsWith(".vcf")| x.getName.endsWith(".vcf.gz")))

    if(!dirs.isEmpty) findallVCF(dirs.flatMap(_.listFiles),acc=acc++newFiles)
    else acc++newFiles
  }

  def script() {


    val cv = new CombineVariants with UNIVERSAL_GATK_ARGS
    cv.variant = findallVCF(List(vcfSrcDir))
    cv.out=swapExt(out,"vcf","subsetted.vcf")

    add(cv)

    for (pop:String <- TargetPopulations) {

      val gps = new GetPopulationSamples
      gps.pop = pop
      gps.SampleList = "sample." + pop + ".list"
      gps.SamplePanel = SamplePanel
      add(gps)

      val sv = new SelectVariants with UNIVERSAL_GATK_ARGS

      sv.variant=cv.out
      sv.sample_file :+= gps.pop
      sv.excludeNonVariants = false
      sv.keepOriginalAC = false
      sv.out = swapExt(out, "vcf", pop + ".vcf")

      add(sv)

    }


  }
}
