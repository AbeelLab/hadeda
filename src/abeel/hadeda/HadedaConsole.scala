package abeel.hadeda

object HadedaConsole {
  def main(args: Array[String]): Unit = {

    if (args.length == 0) {

      listInstructions
    } else {
      args(0) match {
        case "list" => listInstructions
        case "help" => listInstructions
        case "clones" => CloneDiscovery.main(args.drop(1))
        case "clone-details" => DeepDiveMutationsWithinClusters.main(args.drop(1))
        case _ => listInstructions
      }
    }

  }

  def listInstructions() {
    println("Usage:java -jar hadeda.jar [instruction] [instruction options...]")
    println("Instructions:")
    println("  clones             Detect clones within data set")
    println("  clone-details      Output detailed information about individual clones")
    println("")
    println("  list               Show list of all instructions")
    println("  help               Show this help message")

  }

}