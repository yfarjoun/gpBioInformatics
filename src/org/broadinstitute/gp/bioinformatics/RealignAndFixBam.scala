
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.commandline.ClassType
import org.broadinstitute.gp.bioinformatics.utils.PicardCommandLineFunction


class RealignAndFixBam extends QScript {
  qscript =>


  @Input(shortName = "i", required = true, doc = "Input Bam to be fixed") var inputFile:File=_
 // @Output(required = false) var output:File = _

  @Argument(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  @ClassType(classOf[Int])
  @Argument(shortName = "t", required = false, doc = "Thread count for bwa") var threads: Option[Int] = _

  @Argument(shortName = "bq", required = false, doc = "Base quality to reset all the qualities to.") var baseQuality: Int = 40


//  class PicardCommandLineFunction extends JavaCommandLineFunction{
//    var tempDir:List[File]=List("/local/scratch/","/seq/picardtemp3")
//    var jarPath:File="/seq/software/picard/current/bin"
//    var jarName:String=null
//
//    this.memoryLimit=2
//
//    override def freezeFieldValues(): Unit = {
//      super.freezeFieldValues()
//      jarFile=new File(jarPath,jarName)
//    }
//
//    override def commandLine: String = super.commandLine + repeat("TMP_DIR=",tempDir)
//  }

  /*class SamToFastQ extends PicardCommandLineFunction{
    @Input
    var bam:File=_
    @Output
    var fasta:File=_

    jarName="SamToFastq.jar"
    override def commandLine: String = super.commandLine +
      required("I=",bam)+
      required("F=",fasta)
  } **/

  class SamToFastQAndBWAMem extends PicardCommandLineFunction {
    @Input
    var in: File = _
    @Output
    var out: File = _

    @Argument(shortName = "r", required = false, doc = "Reference sequence")
    var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

    if (hasValue(threads)) {
      this.memoryLimit = 2 * threads+1
      this.jobNativeArgs :+= "-pe smp_pe " + (threads+1)
    } else
      this.memoryLimit = 2 * 2



    jarName = "SamToFastq.jar"

    override def commandLine: String = super.commandLine +
      required("I=", in) +
      required("F=", "/dev/stdout") + required("|", escape = false) +
      required("/seq/software/picard/current/3rd_party/bwa_mem/bwa", "mem") +
      required("-p") +
      optional("-t", threads) +
      required(referenceFile) +
      required("/dev/stdin") +
      required(">", escape = false) +
      required(out)
  }


  class ChangeBQ extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_

    jarPath="/seq/tng/farjoun/temp/"
    jarName="ChangeSAMReadQuality.jar"

    override def commandLine: String = super.commandLine + required("I=",in)+required("O=",out) + required("BQ=",baseQuality)
  }

  class AddSyntheticRGAndSort extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_
    @Argument
    var rgName:String="Synthetic"
    jarName="AddOrReplaceReadGroups.jar"

    override def commandLine: String = super.commandLine +
      required("I=",in)+
      required("O=",out)+
      required("PL=",rgName)+
      required("LB=",rgName)+
      required("SM=",rgName)+
      required("PU=",rgName)+
      required("SO=","coordinate")
  }



  class IndexSam(@Input var in:File) extends CommandLineFunction{

    @Output(doc="BAM file index to output", required=false)
    var index: File =  new File(in.getPath + ".bai")

    def commandLine = required("samtools") +
      required("index") +
      required(in) +
      required(index)
  }

  def script() {


    var stfabm=new SamToFastQAndBWAMem()
    stfabm.in=inputFile
    stfabm.out=swapExt(stfabm.in,".bam",".aligned.bam")
    stfabm.referenceFile= referenceFile
    add(stfabm)

    var cbq=new ChangeBQ()
    cbq.in=stfabm.out
    cbq.out=swapExt(cbq.in,".bam",".modBQ.bam")
    add(cbq)

    var asrg=new AddSyntheticRGAndSort()
    asrg.in=cbq.out
    asrg.out=swapExt(asrg.in,".bam",".modRG.sorted.bam")
    add(asrg)

    var ib=new IndexSam(asrg.out)
    add(ib)


  }
}