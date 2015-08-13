package abeel.hadeda

import java.io.File
import atk.io.PatternFileFilter
import atk.util.Tool
import org.apache.commons.math3.stat.inference.TTest
/**
 * This should probably be redone with Kolmogorov-Smirnov test
 * 
 */
object DifferentMean extends Tool {

  // WARNING: NON-FUNCTIONAL PROGRAM
  def main(args: Array[String]): Unit = {
    val date = null
    val folder = new File(".")

    val files = new File(folder + "/clone_identification/").listFiles(new PatternFileFilter("clusterDvalues.1000.*.txt"))
    
    for (a <- files) {
      val x = tLines(a).map(_.toDouble)
      for (b <- files) {
        val y = tLines(b).map(_.toDouble)
        val p=new TTest().tTest(x.toArray,y.toArray)
        println(a.getName+"\t"+b.getName+"\t"+p+"\t"+math.min(1,p*files.size*files.size))
      }
    }

  }

}