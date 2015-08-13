package abeel.hadeda

import atk.util.Tool
import java.io.PrintWriter
import java.io.File

object DeepDiveMutationsWithinClusters extends Tool {

  override val version = """
    2015/08/13        First version on GitHub
    
    """

  case class Config(val input: File = null, val output: String = "", val matrix:File=null)
  def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar hadeda.jar clones") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Cluster membership file")
      opt[File]("matrix")required () action { (x, c) => c.copy(matrix = x) } text ("Matrix file containing SNP: first header row = positions, second header row = annotations, remainder is AP per sample.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("File prefix for the output files.")

    }
    parser.parse(args, Config()) map { config =>

      run(config)
    }
  }

  private def run(config: Config) = {

    val pairs = tLines(config.input).filterNot(_.startsWith("$$")).map { f => val pair = f.split("\t"); pair(0) -> pair(1) }
    val clusters = pairs.groupBy(_._2)

    val input = tLines(config.matrix, false, true).toList
    val position = input(0).split("\t").toList.drop(1)
    val label = input(1).split("\t").toList.drop(1)

    val matrix = tMap(input.drop(2)).mapValues(_.split("\t").toList)
    val log = new PrintWriter(config.output+"deepdive.txt")
    log.println(generatorInfo)
    for (cluster <- clusters) {
      val (identifier, itemsI) = cluster
      println("Processing: " + identifier)
      println("\t" + cluster._2)
      val items = itemsI.map(_._1).toSet
      val subset = matrix.filterKeys(items.contains(_)).toList

      println(subset.map((_._1)))

      val keys = List("position", "label") ++ subset.map(_._1) //gnumbers
      val ssm = subset.map(_._2)
      val values = List(position) ++ List(label) ++ ssm //AP values

      val transposed = values.transpose
      val filtered = transposed.filterNot(list => {
        val numbers = list.drop(2).map(_(0)).filter(_ != '-').groupBy(identity)
        assume(!numbers.contains('-'))
        numbers.size <= 1

      })
      log.println("#=====================================")
      log.println("#  CLONE " + identifier + " -- " + ssm.size + " strains")
      log.println("#=====================================")
      log.println("pos\t" + subset.map((_._1)).mkString("\t") + "\tlabelling")
      val PGRS = filtered.filter(list => list(1).contains("PGRS"))
      val PPE = filtered.filter(list => list(1).contains("PPE"))
      val transposon = filtered.filter(list => list(1).contains("transposase"))
      val filteredVariants = PGRS.size + PPE.size + transposon.size
      log.println("## Excluded: " + nfP.format(filteredVariants / (filtered.size.toDouble)) + "\n#\tPGRS=" + PGRS.size + "\n#\tPPE=" + PPE.size + "\n#\ttransposon=" + transposon.size)
      val removeDeletionInsertion = filtered.filterNot(list => list(1).contains("PGRS") || list(1).contains("PPE") || list(1).contains("transposase")) //filtered//.filterNot(list=>list(1).contains("DELETION")||list(1).contains("INSERTION"))
      log.println("## Retained: " + nfP.format(removeDeletionInsertion.size / filtered.size.toDouble) + "\t" + removeDeletionInsertion.size)

      log.println(removeDeletionInsertion.map(f => f(0) + "\t" + f.drop(2).mkString("\t") + "\t" + f(1)).mkString("\n"))
      log.println("Unique variants: " + removeDeletionInsertion.size)

    }
    log.close

  }

}
