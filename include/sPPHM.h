#ifndef sPPHM_H
#define sPPHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 08-03-2011\n\n
 * This class, sPPHM, is a class written for spinsymmetrical, translationally invaraint two-particle-one-hole matrices. 
 * It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give
 * the relationship between the pph (two-particle one hole) and the sp basis. This matrix has M blocks, M/2 for S = 1/2 block with degeneracy 2
 * and M/2 for S = 3/2 block with degeneracy 4.
 */
class sPPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the sPPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const sPPHM &pphm_p);

   public:
      
      //constructor
      sPPHM();

      //copy constructor
      sPPHM(const sPPHM &);

      //destructor
      virtual ~sPPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      void trace_x(const PPHM &);

      static void init();

      static void clear();

   private:

      //!static list that takes in a pph index i and a blockindex B for spin and momentum, and returns three sp indices a,b,c and intermediate spin S_ab:
      static vector< vector<int> > *pph2s;

      //!static list that takes three sp indices a, b and c, a blockindex B, and an intermediate spinindex S_ab, and returns a pph index i
      static int *****s2pph;

      //!list of block characteristics: 
      static int **block_char;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
