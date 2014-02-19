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
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import org.broadinstitute.sting.commandline
import org.broadinstitute.sting.queue.extensions.picard.AddOrReplaceReadGroups


class RealignAndFixBam extends QScript {
  qscript =>


  @Input(shortName = "i", required = true, doc = "Input Bam to be fixed") var inputFile:File=_
 // @Output(required = false) var output:File = _

  @Argument(shortName = "r", required = false, doc = "Reference sequence") var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")

  @Argument(shortName = "t", required = false, doc = "Thread count for bwa") var threads: Int = 1

  @Argument(shortName = "bq", required = false, doc = "Base quality to reset all the qualties to.") var baseQuality: Int = 40




  class PicardCommandLineFunction extends JavaCommandLineFunction{
    var tempDir:List[File]=List("/local/scratch/","/seq/picardtemp3")
    var jarPath:File="/seq/software/picard/current/bin"
    var jarName:String=null

    override def freezeFieldValues(): Unit = {
      super.freezeFieldValues()
      jarFile=new File(jarPath,jarName)
    }

    protected def lRequired(prefix: String, param: Any, suffix: String = "", spaceSeparated: Boolean = false,
                                    escape: Boolean = true, format: String = "%s"): String =
      super.required(prefix, param, suffix, spaceSeparated, escape, format)


    protected def lOptional(prefix: String, param: Any, suffix: String = "", spaceSeparated: Boolean = false,
                                    escape: Boolean = true, format: String = "%s"): String =
      super.optional(prefix, param, suffix, spaceSeparated, escape, format)


    protected def lRepeat(prefix: String, params: Traversable[_], suffix: String = "", separator: String = " ",
                                  spaceSeparated: Boolean = false, escape: Boolean = true, format: String = "%s",
                                  formatPrefix: (String, Any) => String = (prefix, value) => prefix): String =
      super.repeat(prefix, params, suffix, separator, spaceSeparated, escape, format, formatPrefix)

    override def commandLine: String = super.commandLine + lRepeat("TMP_DIR=",tempDir)
  }

  class SamToFastQ extends PicardCommandLineFunction{
    @Input
    var bam:File=_
    @Output
    var fasta:File=_

    jarName="SamToFastq.jar"
    override def commandLine: String = super.commandLine +
      lRequired("I=",bam)+
      lRequired("F=",fasta)
  }

  class SamToFastQAndBWAMem extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_

    @Argument(shortName = "r", required = false, doc = "Reference sequence")
    var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    @Argument
    var threads:Option[Int]=None

    this.memoryLimit=2048*threads



    jarName="SamToFastq.jar"
    override def commandLine: String = super.commandLine +
      lRequired("I=",in)+
      lRequired("F=","/dev/stdout")  + required("|",escape = false)+
      required("/seq/software/picard/current/3rd_party/bwa_mem/bwa","mem")+
      required("-p")+
      optional("-t",threads)+
      required(referenceFile)+
      required("/dev/stdin")+
      required(">",escape = false)+
      required(out)
  }



  class ChangeBQ extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_
    @Argument
    var bq:Int=0

    jarPath="/seq/tng/farjoun/temp/"
    jarName="ChangeSAMReadQuality.jar"

    override def commandLine: String = super.commandLine + lRequired("I=",in)+lRequired("O=",out) + lRequired("BQ=",bq)
  }

  class AddSyntheticRG extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_
    @Argument
    var rgName:String="Synthetic"
    jarName="AddOrReplaceReadGroups.jar"
    override def commandLine: String = super.commandLine +
      lRequired("I=",in)+
      lRequired("O=",out)+
      lRequired("PL=",rgName)+
      lRequired("LB=",rgName)+
      lRequired("SM=",rgName)+
      lRequired("PU=",rgName)
      lRequired("SO=","unsorted")
  }

  class SortSam extends PicardCommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_
    @Argument
    var sortOrder="coordinate"
    jarName="SortSam.jar"
    override def commandLine: String = super.commandLine + lRequired("I=",in)+lRequired("O=",out) + lRequired("SO=",sortOrder)
  }

  class BWAMem extends CommandLineFunction{
    @Input
    var in:File=_
    @Output
    var out:File=_
    @Argument(shortName = "r", required = false, doc = "Reference sequence")
    var referenceFile: File = new File("/humgen/1kg/reference/human_g1k_v37_decoy.fasta")
    @Argument
    var threads:Option[Int]=None

    this.memoryLimit=2048*threads


    def commandLine: String =  required("/seq/software/picard/current/3rd_party/bwa_mem/bwa")+
    required("mem")+
    required("-p")+
    optional("-t",threads)+
    required(referenceFile)+
    required(in)+
    required(">",escape = false)+
    required(out)
  }

  class IndexSam extends CommandLineFunction{
    @Input
    var in:File=_

    @Output(doc="BAM file index to output", required=false)
    var index: File = _

    override def freezeFieldValues() {
      super.freezeFieldValues()
      if (index == null && in != null)
        index = new File(in.getPath + ".bai")
    }

    def commandLine = required("samtools") +
      required("index") +
      required(in) +
      required(index)
  }



  def script() {


//    var stf=new SamToFastQ()
//    stf.bam=inputFile
//    stf.fasta=swapExt(inputFile,".bam",".fasta")
//    add(stf)
//
//    var  bwa=new BWAMem()
//    bwa.in=stf.fasta
//    bwa.threads=this.threads
//    bwa.out=swapExt(bwa.in,".fasta",".aligned.bam")
//    add(bwa)
//
    var stfabm=new SamToFastQAndBWAMem()
    stfabm.in=inputFile
    stfabm.out=swapExt(stfabm.in,".bam",".aligned.bam")
    stfabm.threads=this.threads
    stfabm.referenceFile= referenceFile
    add(stfabm)

    var cbq=new ChangeBQ()
    cbq.in=stfabm.out
    cbq.out=swapExt(cbq.in,".bam",".modBQ.bam")
    cbq.bq=this.baseQuality
    add(cbq)
    
    var asrg=new AddSyntheticRG()
    asrg.in=cbq.out
    asrg.out=swapExt(asrg.in,".bam",".modRG.bam")
    add(asrg)

    var sb=new SortSam()
    sb.in=asrg.out
    sb.out=swapExt(sb.in,".bam",".sorted.bam")
    add(sb)

    var ib=new IndexSam()
    ib.in=sb.out
    add(ib)


  }
}