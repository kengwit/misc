/*
 * FEN.scala
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package famuls

class FEN(id: Int) {

    // node id
    private var _id = id

    // coordinates
    private var _ref_loc: Array[Double] = new Array[Double](3)

    // return node id
    def id(): Int = _id

    // return nodal coordinates
    def ref_loc(): Array[Double] = _ref_loc

    // set coordinates
    def set_loc(x: Double, y: Double, z: Double) = {
        _ref_loc(0) = x
        _ref_loc(1) = y
        _ref_loc(2) = z
    }
}
