
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import java.io.InputStream
import java.util.Scanner

class TrimBam extends QScript {
  qscript =>


  @Input(shortName = "I", doc = "A bam to aggregate into one.", required = true)
  var Input: List[File] =_

  @Output(shortName="O",doc = "the name of the final output file")
  var Output:File = _

  @Argument(shortName="otherBam",doc="Bam to which we would like to match the read-length and coverage of the input bam.",required=true)
  var OtherBAM:File=_




  def script(){


    val inputReadLengths=Input.map(x=>{
      val findInputReadLength=new FindReadLength
      findInputReadLength.Input=x
      findInputReadLength.getOutput.toInt
    })

    if(inputReadLengths.reduce(Math.min)!=inputReadLengths.reduce(Math.max)){
      throw new RuntimeException("Not all input bam reads have same length. Not currently implemented")
    }
    val inputReadLength=inputReadLengths(0)

    val findOtherReadLength=new FindReadLength
    findOtherReadLength.Input=OtherBAM
    val otherReadLength=findOtherReadLength.getOutput.toInt
    val trimLeft=Math.max(0,inputReadLength-otherReadLength)



    val inputNumberOfRecords=Input.map(x=>{
      val temp=new FindNumberOfRecords
      temp.Input=x
      temp.getOutput.toInt
    }).reduce((x,y)=> x + y)


    val findOtherNumberOfRecords=new FindNumberOfRecords
    findOtherNumberOfRecords.Input=OtherBAM
    val otherNumberOfRecords=findOtherNumberOfRecords.getOutput.toInt
    val sampleRatio=Math.min(1,otherNumberOfRecords.toDouble / inputNumberOfRecords)



    val msf=new MergeSamFiles
    msf.Input=Input
    msf.Output=swapExt(Output,".bam",".merged.bam")
    add(msf)

    val rsf = new RevertSamFile
    rsf.Input=msf.Output
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


  abstract class GetRuntimeCommand extends CommandLineFunction{

    val pipe = required("|" ,escape=false)

    def execCmd(cmd:String):String={
      val process=Runtime.getRuntime.exec(cmd)
      if(process.exitValue()==0){
        getFromStream(process.getInputStream)
      }else{
        val temp=getFromStream(process.getErrorStream)
        logger.error(s"command returned an error:\n $temp")
        val retval=process.exitValue()
        throw new RuntimeException(s"command returned an error $retval")
      }
    }


    def getFromStream(is:InputStream):String={
      val s:Scanner = new java.util.Scanner(is) useDelimiter "\\A"
      if (s.hasNext) s.next() else ""
    }


    def getOutput:String={
      logger.debug(s"calling INLINE command:\n $commandLine")
      val retval=execCmd(commandLine)
      logger.debug(s"got INLINE command output:\n $retval")
      retval
    }

  }

  class FindReadLength extends GetRuntimeCommand{
    @Input var Input:File=_


    override def commandLine:String =
      required("samtools", "view",escape=false)+
      required(Input.getAbsolutePath)+ pipe +
      required("head", "-n1",escape=false)+ pipe +
      required("cut", "-f","10",escape=false)+ pipe +
      required("wc", "-c",escape=false) + pipe +
      required("awk", "'{print $1-1}'",escape=false)
  }
  class FindNumberOfRecords extends GetRuntimeCommand{
    @Input var Input:File=_


    override def commandLine:String =
      required("samtools", "idxstats",escape=false)+
      required(Input.getAbsolutePath)+ pipe +
      required("cut", "-f","3,4",escape=false)+ pipe +
      required("tr", "\t", "\n",escape=false ) + pipe +
      required("paste", "-sd","+",escape=false)+ pipe+
      required("bc",escape=false)
  }
  class MergeSamFiles extends PicardCommandLineFunction{

    @Input var Input:List[File]=_
    @Output var Output:File=_

    jarName="MergeSamFiles.jar"

    override def commandLine: String =
      repeat("I=",Input) +
      required("O=",Output)+
      required("MERGE_SEQUENCE_DICTIONARIES=",escape = true)
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



