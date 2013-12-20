
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

  @Argument(shortName = "pop", doc = "The list of target populations", required = true)
  var TargetPopulations: List[String] = _

  @Argument(shortName = "scatterCount", doc = "The number of ways to scatter this job", required = false)
  val scatterCount: Int = 25

  @Argument(shortName = "L", required = false)
  var intervals: List[String] = _

  @Argument(shortName = "R", required = false)
  var reference: File = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"


  class GetPopulationSamples extends RScriptFunction {
    @Argument var pop: String = _
    @Input var SamplePanel: File = _
    @Output var SampleList: File = _


    lazy val rscript: String = s"""
        |read.table("$SamplePanel" ,col.names=c("SAMPLE","SUB.POP","POP","SEQUENCING.CENTER","SEQUENCING.CENTER.2"),fill=TRUE)->data
        |show(unique(data$$POP))
        |write.table(file="$SampleList",subset(data,POP=="$pop")$$SAMPLE,row.names=FALSE,quote=FALSE,col.names=FALSE)
        """.stripMargin

  }

  class MakeSampleFreeVCF extends CommandLineFunction {
    @Input var vcf: File = _
    @Output var out: File = _

    override def commandLine =
      required("grep") +
        required("^##") +
        required(vcf) +
        required(">", escape = false) +
        required(out) +
        required(";", escape = false) +
        required("grep") +
        required("-v", "^##") +
        required(vcf) +
        required("|", escape = false) +
        required("cut") +
        required("-f", "1-8") +
        required(">>", escape = false) +
        required(out)
  }

  final def findallVCF(f: Seq[File]): Seq[File] = {

    val dirs=f.filter(_.isDirectory).filter(x=> !x.getName.startsWith("."))
    println("dir:"+dirs)
    (f++dirs.flatMap(_.listFiles())).filter(x=> x.isFile & (x.getName.endsWith(".vcf")| x.getName.endsWith(".vcf.gz")))
  }

  def script() {


    val cv = new CombineVariants
    cv.intervalsString++=intervals
    cv.reference_sequence=reference
    cv.variant = findallVCF(List(vcfSrcDir))
    cv.out=swapExt(out,"vcf","subsetted.vcf")
    cv.scatterCount=scatterCount
    println("intervals:"+cv.intervalsString)

    add(cv)

    for (pop:String <- TargetPopulations) {

      val gps = new GetPopulationSamples
      gps.pop = pop
      gps.SampleList = "sample." + pop + ".list"
      gps.SamplePanel = SamplePanel
      add(gps)

      val sv = new SelectVariants
      sv.reference_sequence=reference
      sv.intervalsString++=intervals
      sv.variant=cv.out
      sv.sample_file :+= gps.SampleList
      sv.excludeNonVariants = false
      sv.keepOriginalAC = false
      sv.out = swapExt(out, "vcf", pop + ".vcf")
      add(sv)

      val msfvcf = new MakeSampleFreeVCF
      msfvcf.vcf=sv.out
      msfvcf.out=swapExt(msfvcf.vcf,"vcf","sites.only.vcf")
      add(msfvcf)

    }


  }
}
