


package org.broadinstitute.gp.bioinformatics.utils

import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Created by farjoun on 12/17/13.
 */


trait RScriptFunction extends CommandLineFunction{

  val rscript:String

  override def commandLine=required("echo", rscript) +
    required("|",escape=false) +
    required("R")+
    required("--vanilla")+
    required("--slave")

}
