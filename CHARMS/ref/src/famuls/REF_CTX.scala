/*
 * REF_CTX.scala
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package famuls
import scala.collection.mutable.HashSet

class REF_CTX(gmesh: GMESH) {

  private var _gmesh = gmesh
  private var _activated_fens = new HashSet[FEN]

  //REF_FEN_SRC          *_ref_fen_src;
  //set <FEN *>           _activated_fens;
  //set <FEN *>           _deactivated_fens;
  //bool                  _true_hierarchical;

}
