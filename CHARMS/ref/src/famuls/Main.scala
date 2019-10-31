/*
 * Main.scala
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package famuls

object Main {


  def main(args: Array[String]) {

    var test_node_1 = new FEN(1)
    var test_node_2 = new FEN(2)
    var test_node_3 = new FEN(3)

    var my_mesh = new GMESH("my_mesh")

    println("mesh name: " + my_mesh.name)

    println("number of nodes: " + my_mesh.fen_total)

    println("adding some nodes")

    my_mesh.add_fen(test_node_1)
    my_mesh.add_fen(test_node_2)
    my_mesh.add_fen(test_node_3)

    println("1st node: " + my_mesh.find_fen(3).id)
    println("number of nodes: " + my_mesh.fen_total)
    println("end test program")

  }

}


