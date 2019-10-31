/*
 * GMESH.scala
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package famuls
import scala.collection.mutable.HashMap
import scala.collection.jcl.LinkedList
import scala.Math._

class GMESH(s: String) {

  // mesh name
  private var _name: String = s
  
  // max node id
  private var _max_fen_id: Int = 0

  // list of nodes
  private var _fens = new LinkedList[FEN]()

  // id to node map
  private var _fen_map = new HashMap[ Int, FEN ]

  //list <GSUBMESH *>                      _gsubmeshes;

  // ------------------------------------------------------------
  //  private methods
  // ------------------------------------------------------------
  // test method - should read from file later
  private def read_fens(fen: FEN): Unit = {
      _fen_map += fen.id -> fen
  }

  // ------------------------------------------------------------
  //  public methods
  // ------------------------------------------------------------
  // get mesh name
  def name(): String = _name

  // Get the total of nodes associated with mesh.
  def fen_total(): Int = _fens.length

  // Find a finite element node by id
  def find_fen (id: Int): FEN = _fen_map(id)
      
  // Add a finite element node.
  def add_fen(fen: FEN): Unit = {
      _fens.add(fen)
      read_fens(fen)
      _max_fen_id = max( _max_fen_id, fen.id );
  }

  // return max node number
  def max_fen_id(): Int = _max_fen_id
  
}
