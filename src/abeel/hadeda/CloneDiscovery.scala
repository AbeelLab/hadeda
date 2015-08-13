package abeel.hadeda

import atk.util.Tool
import be.abeel.util.HashMap2D
import scala.collection.JavaConversions._
import scala.util.Random
import scala.collection.mutable.MutableList
import scala.collection.mutable.HashSet
import java.io.PrintWriter
import be.abeel.util.FrequencyMap
import java.io.File
import scala.collection.mutable.LinkedList
import be.abeel.util.FrequencyMapUtils

object CloneDiscovery extends Tool {

  
  override val version="""
    2014-2015     Initial version that identifies clusters based on DBSCAN
    2015/04/07    Made all detailed files optional verbose output
                  Improved clustering naming scheme  
    2015/08/13    First version on GitHub
    """
  
  /**
   * Specifies the density (the range-query must contain at least minPoints
   * instances)
   */
  private val minPoints = 2;

  class Wrapper(val key: String) {
    /*
     * ClusterIDX
     *
     * 0 non classified
     * -1 noise
     * 1...  cluster
     */
    var clusterIdx = 0
  }

  case class Config(val verbose: Boolean = false, val outputPairs: Boolean = true, val input: String = null, val output: String = null)

  def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar hadeda.jar clones") {
      opt[String]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file that contains the pairwise snp count matrix (*.apr15.snpmatrix.txt).") 

      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("File prefix for the output files.")

      opt[Unit]("no-pairs") action { (x, c) => c.copy(outputPairs = false) } text ("Do not output clusters that only contain two strains.")
      opt[Unit]("verbose") action { (x, c) => c.copy(verbose = true) } text ("Output a lot of additional files with detailed cluster information.")

    }
    parser.parse(args, Config()) map { config =>

      run(config)
    }
  }

  private def obfus(id: Int, range: Int) = {
   
    range + "-" + range.toHexString + "." + id.toHexString
  }

  private def run(config: Config) {
    //    val dbs = new DensityBasedSpatialClustering();
    val map = tMap(tLines(config.input))
    val keys = map.getOrElse("Taxon", null).split("\t")

    val matrix = new HashMap2D[String, String, Int]()
    map.filterNot(f => f._1.equals("Taxon")).map(p => {
      val second = keys.zip(p._2.split("\t").map(_.toInt)).toMap

      matrix.put(p._1, second)

    })

    val rg = new Random()

    val outputFile = new File(config.output + "clusterlog.txt")
    outputFile.getAbsoluteFile().getParentFile().mkdirs()
    val pw = new PrintWriter(outputFile)
    pw.println(generatorInfo)
    pw.println("##SNP distance\tno clusters\tavg. cluster size\tavg. pairwise snpd averages \tmin pairwises snpd averages\tmax...\tactuals (cluster Id, avg. intra pairwise)")
    val eps = List(5, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 100, 200, 500, 1000, 2000, 3000)
    for (epsilon <- eps) {
      val wrapMap = keys.map(f => f -> new Wrapper(f)).toMap
      val outFile = new File(config.output + "clustermembership." + epsilon + ".txt")
      val out = new PrintWriter(outFile)
      out.println(generatorInfo)
      out.println("$$\t" + epsilon)
      var currentIdx = 1
      var clusterCount = 0
      val clusters: MutableList[List[Wrapper]] = new MutableList[List[Wrapper]]

      val fmClusterSize = new FrequencyMap
      for (item <- wrapMap) {
        //      println("Processing clustering: "+currentIdx)
        if (item._2.clusterIdx == 0) {

          expand(item._1, currentIdx)
          val cluster = wrapMap.filter(_._2.clusterIdx == currentIdx).toList

          if (cluster.size > 0) {
            val link: MutableList[Wrapper] = new MutableList[Wrapper]

            clusters += cluster.map(_._2).toList
            clusterCount += 1
            fmClusterSize.count(cluster.size)

            println("Writing: " + currentIdx + "\t" + cluster.size + "\t")
            val ob = obfus(currentIdx, epsilon)
            out.println(cluster.map(_._1).mkString("\t" + ob + "\n") + "\t" + ob)

          }

          currentIdx += 1
        }
      }
      /**
       * Extract pair clones
       */
      if (config.outputPairs) {
        for (item <- wrapMap.filter(_._2.clusterIdx == -1)) {

          val neighbors = matrix.getOrElse(item._1, null).filter(p => p._2 <= epsilon && !p._1.equals(item._1)).map(f => f._1).toList

          println("Unclustered: " + item + "\t" + neighbors.size)
          assume(neighbors.size < minPoints, "Too many neighbors, clustering likely failed: " + neighbors)
          if (neighbors.size > 0) {
            wrapMap.getOrElse(item._1, null).clusterIdx = currentIdx
            for (n <- neighbors) {
              wrapMap.getOrElse(n, null).clusterIdx = currentIdx
            }

            val cluster = wrapMap.filter(_._2.clusterIdx == currentIdx).toList

            if (cluster.size > 0) {
              val link: MutableList[Wrapper] = new MutableList[Wrapper]

              clusters += cluster.map(_._2).toList
              clusterCount += 1
              fmClusterSize.count(cluster.size)

              println("Writing: " + currentIdx + "\t" + cluster.size + "\t")
              val ob = obfus(currentIdx, epsilon)
              out.println(cluster.map(_._1).mkString("\t" + ob + "\n") + "\t" + ob)

            }

            currentIdx += 1
          }

        }
      }

      out.close

      def expand(item: String, currentClusterIdx: Int) {
        val initialSeeds: HashSet[String] = new HashSet()
        initialSeeds.addAll(matrix.getOrElse(item, null).filter(p => p._2 <= epsilon && !p._1.equals(item)).map(f => f._1).toList)
        /* Seeds that have been used in expansion */
        val usedSeeds: HashSet[String] = new HashSet()
        if (initialSeeds.size < minPoints)
          wrapMap.getOrElse(item, null).clusterIdx = -1
        else {
          initialSeeds.map(f => wrapMap.getOrElse(f, null).clusterIdx = currentClusterIdx)

          while (initialSeeds.size > 0) {
            val currentSeed = initialSeeds.head
            usedSeeds.add(currentSeed)
            val neighbors = matrix.getOrElse(currentSeed, null).filter(p => p._2 <= epsilon && !p._1.equals(currentSeed)).map(f => f._1).toList
            if (neighbors.size > minPoints) {
              for (n <- neighbors) {
                if (!usedSeeds.contains(n)) {
                  initialSeeds.add(n)
                }
                wrapMap.getOrElse(n, null).clusterIdx = currentClusterIdx
              }
            } else {
              for (n <- neighbors) {
                wrapMap.getOrElse(n, null).clusterIdx = currentClusterIdx
              }
            }
            initialSeeds.remove(currentSeed)
          }
        }
      }

      val averageDistances = for {
        cluster <- clusters

        val fm = new FrequencyMap
        val distances = for {
          i <- cluster;
          j <- cluster
          if (i != j)
        } yield {
          fm.count(matrix.get(i.key, j.key))
          matrix.get(i.key, j.key)
        }

      } yield {
        val idx = if (cluster.size > 0) cluster.head.clusterIdx else "-"
        if (config.verbose) {
          if (fm.size() > 1)
            FrequencyMapUtils.plot(fm, config.output + "clusterDdistribution." + epsilon + "." + idx + ".png")

          val pw = new PrintWriter(config.output + "clusterDvalues." + epsilon + "." + idx + ".txt")
          pw.println(generatorInfo)
          pw.println("# Pairwise distances for epsilon: " + epsilon + ", and cluster: " + idx)
          pw.println(distances.mkString("\n"))
          pw.close
        }
        (distances.sum.toDouble / distances.size)

      }

      val clusterIdx = clusters.map(cluster => if (cluster.size > 0) cluster.head.clusterIdx else -1)

      if (averageDistances.size == 0) {
        averageDistances += 0
        clusterIdx += -1
      }
      assume(averageDistances.size == clusterIdx.size)
      nf.setMaximumFractionDigits(0)
      pw.println(epsilon + "\t" + clusterCount + "\t" + nf.format(fmClusterSize.average()) + "\t" + nf.format(averageDistances.sum / averageDistances.size) + "\t" + nf.format(averageDistances.min) + "\t" + nf.format(averageDistances.max) + "\t" + (clusterIdx.zip(averageDistances.map(nf.format(_)))).mkString(";"))

    }
    pw.close
   
  }
}



