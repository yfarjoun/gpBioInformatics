


package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.function.CommandLineFunction

/**
 * Created by farjoun on 12/17/13.
 */


abstract trait RScriptFunction extends CommandLineFunction{

  abstract val rscript:String

  override def commandLine=required("echo", rscript) +
    required("|",escape=false) +
    required("R")+
    required("--vanilla")+
    required("--slave")

}