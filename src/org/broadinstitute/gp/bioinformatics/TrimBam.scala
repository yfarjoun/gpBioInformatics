
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import java.io.InputStream
import java.util.Scanner
import org.broadinstitute.sting.commandline.Argument

class TrimBam extends QScript {
  qscript =>


  @Input(shortName = "I", doc = "A bam to aggregate into one.", required = true)
  var Input: File =_

  @Output(shortName="O",doc = "the name of the final output file")
  var Output:File = _

  @Argument(shortName="otherBam",doc="Bam to which we would like to match the read-length and coverage of the input bam.",required=true)
  var OtherBAM:File=_

  @Argument
  var bamReadLength:Integer=_
  @Argument
  var otherBamReadLength:Integer=_

  @Argument
  var bamReads:Integer=_
  @Argument
  var otherReads:Integer=_


  def script(){


    val sampleRatio=Math.min(1,otherReads.toDouble / bamReads)

    val trimLeft=Math.max(0,bamReadLength-otherBamReadLength)


   /* val msf=new MergeSamFiles
    msf.Input=Input
    msf.Output=swapExt(Output,".bam",".merged.bam")
    add(msf)
     */

    val rsf = new RevertSamFile
    rsf.Input=Input
    rsf.Output=swapExt(rsf.Input,".bam",".reverted.bam")
    add(rsf)

    val tr=new TrimReads
    tr.Input=rsf.Output
    tr.Trim=trimLeft
    tr.Output=swapExt(tr.Input,".bam",".trimmed.bam")
    add(tr)

    val dss= new DownSampleSam
    dss.Input=swapExt(if(trimLeft!=0){tr.Output}else {rsf.Output},".bam",".downsampled")
    dss.Fraction= sampleRatio
    dss.Output=Output
    add(dss)

  }

  class MergeSamFiles extends PicardCommandLineFunction{

    @Input var Input:List[File]=_
    @Output var Output:File=_

    jarName="MergeSamFiles.jar"

    override def commandLine: String =
      repeat("I=",Input) +
      required("O=",Output)+
      required("MERGE_SEQUENCE_DICTIONARIES=")
  }
  class RevertSamFile extends PicardCommandLineFunction{
    @Input var Input:File=_
    @Output var Output:File=_

    jarName="RevertSam.jar"

    override def commandLine: String =
      required("I=",Input) +
      required("O=", Output)

  }

  class TrimReads extends PicardCommandLineFunction{

    @Input var Input:File=_
    @Output var Output:File=_
    @Argument var Trim:Int=0

    jarPath="/home/unix/tfennell/bin/"
    jarName="TrimSamFile.jar"

    override def commandLine: String =
      required("I=",Input) +
      required("O=",Output)+
      required("TRIM_RIGHT=",Trim)
  }


  class DownSampleSam extends  PicardCommandLineFunction{

    @Input var Input:File=_
    @Output var Output:File=_

    @Argument var Fraction:Double=1

    jarName="DownsampleSam.jar"


    override def commandLine: String =
      required("I=",Input) +
      required("O=",Output)+
      required(s"P=$Fraction)")
  }







}



