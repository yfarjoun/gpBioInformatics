package org.broadinstitute.gp.bioinformatics

import org.broadinstitute.sting.queue.function.JavaCommandLineFunction
import java.io.File

/**
 * Created by farjoun on 3/7/14.
 */

class PicardCommandLineFunction extends JavaCommandLineFunction{
    var tempDir:List[File]=List(new File("/local/scratch/"),new File("/seq/picardtemp3"))
    var jarPath:File=new File("/seq/software/picard/current/bin")
    var jarName:String=null

    this.memoryLimit=Option(2)

    override def freezeFieldValues(): Unit = {
      super.freezeFieldValues()
      jarFile=new File(jarPath,jarName)
    }

    override def commandLine: String = super.commandLine + repeat("TMP_DIR=",tempDir)

}
