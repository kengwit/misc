/*
 * GCELL.scala
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package famuls

abstract class GCELL {

    // At which level in the hierarchy is the cell?
    // The coarsest level is zero (0).
    def level(): Int

   // Return the parent of this gcell.  May be null if the
   // cell has no parent (meaning it is at the bottom).
   def parent(): GCELL

   // Get the i-th child of the cell.
   def child(i: Int): GCELL


   /**
     Get the number of children.  If the gcell is not refined, zero should be
     returned; otherwise the number of child gcells is returned
     (depending on the type of the gcell).
  */
  def  nchildren(): Int = 0


  /**
     Has the gcell been divided into children?  Note we're not asking if they
     *could* exist, rather whether they actually exist.
  */
 def divided(): Boolean = false
 
  /**
     Physically divide the gcell by creating its children (and all the nodes
     these children connect on the higher level).  Needs to be defined
     by the implementation of the specialization.
  */
 def divide()



    
    //virtual void divide (class REF_CTX *ref_ctx) = 0;
  /**
     Return the set of detail functions that are needed for the refinement
     of the function associated with the node on input.
     This function may be called *ONLY* if the gcell had been divided
     into children before.
  */
  //virtual void detail_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf) = 0;
  /**
     Return the set of all refinement functions (i.e. detail functions plus
     the one private function) that are needed for the refinement
     of the function associated with the node on input.
     This function may be called *ONLY* if the gcell had been divided
     into children before.
  */
  //virtual void  complete_refinement_set (FEN *fen_of_bfun_to_refine,set <FEN*> *rf) = 0;

}
