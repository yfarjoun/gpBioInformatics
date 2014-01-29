/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.utils.interval.IntervalSetRule
import org.broadinstitute.sting.commandline.{ArgumentException, ClassType}



class FindNoKmerCoverageScript extends QScript {
  qscript =>

  @Argument(shortName = "i", required = false, doc = "Intervals") var interval: String = _
  @Argument(shortName = "xi", required = false, doc = "Excluded Intervals") var excludeInterval: String = _

  @Argument(shortName = "t", required = false, doc = "Names of bams to be placed in plots") var inputTag:List[String]=_
  @Argument(shortName = "b", required = true, doc = "List of BAM files") var bamList: List[File] = _
  @Argument(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
  @Argument(shortName = "sc", required = false, doc = "Scatter count") var jobs: Int = 25
  @Argument(shortName = "isr", required = false, doc = "Interval set rule") var intervalSetRule = IntervalSetRule.UNION
  @Argument(doc="A flag informing that the bams are to be processed separately")var KeepBamsSeparate:Boolean=false
  @ClassType(classOf[Int])
  @Argument(shortName = "mkl", required = false, doc = "minimum kmer length(s)") var minKmerLength:List[Int] = List(5)
  @ClassType(classOf[Int])
  @Argument(shortName = "mbq", required = false, doc = "minimum base quality(ies)") var minBaseQuality: List[Int] = List(20)
  @ClassType(classOf[Int])
  @Argument(shortName = "mmq", required = false, doc = "minimum mapping quality(ies)") var minMappingQuality: List[Int] = List(20)

  @Argument(shortName = "o", required = true, doc ="output base name") var outputBaseName:File = _
  @Input(shortName="amountOfN",required=false,doc="a file containing the number of N bases each contig has") var amountOfN:File=null

  @Argument(shortName = "baits", required = false, doc ="coding regions") var baitsfile:String = _
  @Argument(shortName = "targets", required = false, doc ="targets") var targetsfile:String = _
  @Argument(shortName = "gff", required = false, doc ="gene features") var genomeFeature:String = null // "/seq/references/Homo_sapiens_assembly19/v1/annotation/gencode.v12.annotation.patched_contigs.gff"
  @Argument(shortName="rpath",required=false,doc="the name of the path to where the R scripts for this qscript are stored") var rFilesPath:File="./private/R/scripts/org/broadinstitute/sting/techdev/NoCoverageScripts"
  @Argument(shortName="centLoc",doc="A file documenting the location of various elements in the genome (from http://genome.ucsc.edu/cgi-bin/hgTables, choose table:gap)",required = false)
  var centromereLocation:File=null

  val tempIntervalFile="targetInterval.interval_list"

  trait CommonArguments extends CommandLineGATK{
      this.reference_sequence=referenceFile
 
  }


  def script() {


    //write input interval as interval_list
    val mil=new MergeIntervalLists()
    val ciltr2 = new ConvertIntervalListToR(false)
    if(interval!=null){
      mil.reference_sequence=referenceFile
      mil.targetIntervalsString:+=interval
      mil.intervalsString:+=interval
      mil.excludeIntervalsString:+=excludeInterval
      mil.out=tempIntervalFile
      add(mil)

      ciltr2.intervalList=mil.out
      ciltr2.output=swapExt(mil.out,"",".forR")
      add(ciltr2)
    }

    var outputs:List[String]=null

    val minlength=List(minKmerLength.length,minBaseQuality.length,minMappingQuality.length) reduce math.min
    val maxlength=List(minKmerLength.length,minBaseQuality.length,minMappingQuality.length) reduce math.max

    if(minlength!=maxlength & minlength!=1)  throw new ArgumentException("the lists minKmerLength, minBaseQuality, and minMappingQuality should be of equal length or length 1")

    if(minBaseQuality.length==1)    minBaseQuality    = List.fill(4)(minBaseQuality).flatten
    if(minMappingQuality.length==1) minMappingQuality = List.fill(4)(minMappingQuality).flatten
    if(minKmerLength.length==1)     minKmerLength     = List.fill(4)(minKmerLength).flatten

    val listOfArgs=( minKmerLength zip minBaseQuality) zip minMappingQuality map {
      case ((x,y),z) => (x,y,z)
    }

    for( (mkl,mbq,mmq) <- listOfArgs) {
      if (KeepBamsSeparate) outputs = bamList.map(x => process_bams(mkl, mbq,mmq, List(x)))
      else outputs = List(process_bams(mkl ,mbq,mmq, bamList))
    }


    //R scripts to plot stuff

    val mp=new MakePlots()
    mp.rpath=qscript.rFilesPath
    mp.input=outputs
    mp.inputTags=inputTag
    mp.output=swapExt(outputBaseName,"",".plots.pdf")
    if(interval!=null) mp.interval=ciltr2.output
    add(mp)

  }


  //the output is the name of the r-readable file that is the output of ConvertIntervalListToR
  def process_bams(mkl:Int,mbq:Int,mmq:Int, bamList: List[File]):String={


    val fncbl_output = swapExt(outputBaseName, "",bamList.map(x=>x.getName).mkString(".")+ ".kmer." + mkl + ".mbq."+ mbq+ ".mmq." + mmq+ ".interval_list")


    //Find No Coverage
    val fncbl = new FindNoKmerCoverageByLocus() with CommonArguments
    if(interval!=null) fncbl.intervalsString:+=interval

    if(excludeInterval!=null) fncbl.excludeIntervalsString:+=excludeInterval

    fncbl.interval_set_rule=intervalSetRule
    fncbl.out = fncbl_output
    fncbl.mkl = mkl
    fncbl.minimumBaseQuality=mbq
    fncbl.minimumMappingQuality=mmq

    fncbl.input_file = bamList
    fncbl.memoryLimit = 8
    fncbl.scatterCount = qscript.jobs
    add(fncbl)

    //Convert To R-Readable

    val ciltr = new ConvertIntervalListToR(false)
    ciltr.intervalList=fncbl.out
    ciltr.output=swapExt(fncbl.out,"",".forR")
    add(ciltr)


    //Diagnose Missing Intervals
      val qmi=new QualifyMissingIntervals()
      qmi.out = swapExt(fncbl_output,"interval_list","missing.grp")

      qmi.baitsfile=if(baitsfile!=null)baitsfile else fncbl_output
      qmi.targetsfile=if(targetsfile!=null)targetsfile else fncbl_output
      qmi.reference_sequence=qscript.referenceFile
      qmi.intervals:+=fncbl_output
      qmi.input_file=qscript.bamList
      qmi.scatterCount=1
      add(qmi)

    if (genomeFeature != null) {
      // Convert To BED
      val ciltb = new ConvertIntervalListToBed()
      ciltb.intervalList = fncbl.out
      ciltb.output = swapExt(fncbl.out, "interval_list", "bed")
      add(ciltb)

      // bedtools intersect
      val bti = new BEDToolsIntersect()
      bti.bedInput = ciltb.output
      bti.output = swapExt(ciltb.output, "bed", "annotated.bed")
      bti.genomeFeature = genomeFeature
      add(bti)

      // R script to expose the uncovered genes
      val cgrs=new CleanGFFRScript()
      cgrs.input=bti.output
      cgrs.output=swapExt(bti.output,"annotated.bed","genes")
      cgrs.rpath=rFilesPath
      add(cgrs)

    }
    qmi.out
  }


  class CleanGFFRScript() extends CommandLineFunction {
    @Argument var rpath:File=_
    @Input var input:File=_
    @Output var output:File=_
    def commandLine = "Rscript %s/clean.gff.R %s %s".format(rpath, input, output)
  }

  class MakePlots() extends CommandLineFunction {
    @Argument var rpath:File=_
    @Input var input:List[File]=_
    @Argument var inputTags:List[String]=null
    @Output var output:File=_
    @Input(required=false) var interval:File=_
    def commandLine = required("Rscript")+
      required("%s/AnalyseNoCoverageData.R".format(rpath)) +
 //   required("--args")+
    required("--sourceF",input.mkString(","))+
    required("--sourceN",(if(inputTags!=null)inputTags else input.map(x=>x.toString)).mkString(",") ) +
    optional("--centromereLocation", centromereLocation)+
    optional("--intervalList",interval)+
    optional("--amountOfN",amountOfN)+
    required("--dictionary",swapExt(referenceFile,".*",".dict"))
    required("--output",output)
  }



  class BEDToolsIntersect extends CommandLineFunction{
    @Input(shortName="bedin",doc="A bed file to intersect with the genome definition file") var bedInput:File=_
    @Argument(shortName="gff",doc="A General feature format containing the gene definition one would like to use") var genomeFeature:File=_
    @Output(shortName="bedOut",doc="The name of the output file")var output:File=_

    def commandLine=
      required("bedtools")+
      required("intersect")+
      required("-wb")+
      required("-a",bedInput)+
      required("-b",genomeFeature)+
      required(">",escape=false)+
      required(output)
  }


  class ConvertIntervalListToR(withPlus:Boolean=true) extends CommandLineFunction{
    @Input(shortName = "il",doc="An interval list to convert to R readable format") var intervalList:File = null
    @Output(shortName="o",doc="The name of the output file") var output:File=null

    def commandLine =
        required("sed")+
        required("s/:/\t/; s/-/\t/")+
        required(intervalList)  +
    required("|",escape=false)+
      required("awk")+
          conditional(withPlus,"BEGIN{OFS=\"\\t\"}{if($4==\"\") { $4=$3; $3=$2;} ; print $0}" )+
          conditional(!withPlus,"BEGIN{OFS=\"\\t\"}{if($3==\"\") { $3=$2;} ; print $0}" )+
          required(">",escape = false) +
      required(output)
  }

  class ConvertIntervalListToBed extends CommandLineFunction{
    @Input(shortName = "il",doc="An interval list to convert to R readable format") var intervalList:File = null
    @Output(shortName="o",doc="The name of the output file",required = false) var output:File=null


    def commandLine =
      required("sed")+
        required("s/:/\t/; s/-/\t/ ; s/$/\t+\tinterval/g")+
        required(intervalList)+
        required("|",escape=false)+
        required("awk")+
        required("BEGIN{OFS=\"\\t\"}{if($3==\"+\") { $5=$4; $4=$3; $3=$2 } ; $2=$2-1;print $0\"-\"NR}")+
        required(">",escape = false)+
        required(output)
  }



}


