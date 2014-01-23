
package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine.OUTPUT_MODE._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model._
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE._

import java.io.File
import org.broadinstitute.sting.queue.function.RetryMemoryLimit

/**
 *
 */
class GenotypeProjects extends QScript {
  qscript =>

  @Input(doc = "The reference file for the bam files.", shortName = "R")
  var referenceFile: File = _

  @Input(doc = "The dbSNP file to use", shortName = "D", required = false)
  var dbSNP: File = _

  @Input(doc = "Bam list files to genotype.", shortName = "bs")
  var bamFiles: List[File] = _

  trait CommonArguments extends CommandLineGATK {
    this.reference_sequence = qscript.referenceFile
    this.intervals :+= qscript.dbSNP
  }

  class myUG extends UnifiedGenotyper with CommonArguments with RetryMemoryLimit {
    this.dbsnp = qscript.dbSNP
    this.memoryLimit = 3.5
    this.gt_mode = GENOTYPE_GIVEN_ALLELES
    this.out_mode = EMIT_ALL_SITES
    this.genotype_likelihoods_model = SNP
    this.alleles = qscript.dbSNP
    this.isIntermediate=false
    this.contamination=0
  }

  class myVT extends VariantsToTable with CommonArguments {
    this.memoryLimit = 1
    this.fields ++= List("CHROM","POS","REF","ALT","ID")
    this.genotypeFields ++= List("GT","GQ","PL")
  }


  def script() {
    bamFiles.map(x => {
      val ug = new myUG()
      ug.input_file :+= x
      ug.out = swapExt(x.getParent,x, "bam.list", "vcf")
      ug.contaminationFile= swapExt(x.getParent,x,"bam.list","contam")
      add(ug)

      val vt = new myVT()
      vt.variant :+= ug.out
      vt.out = swapExt(ug.out.getParent,ug.out, "vcf", "table")
      add(vt)
    })

  }
}
