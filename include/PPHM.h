#ifndef PPHM_H
#define PPHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 08-03-2011\n\n
 * This class, PPHM, is a class written for spinsymmetrical, translationally invaraint two-particle-one-hole matrices. 
 * It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give
 * the relationship between the pph (two-particle one hole) and the sp basis. This matrix has M blocks, M/2 for S = 1/2 block with degeneracy 2
 * and M/2 for S = 3/2 block with degeneracy 4.
 */
class PPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PPHM &pphm_p);

   public:
      
      //constructor
      PPHM();

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const;

      int get_inco(int B,int S_ab,int a,int b,int c,int &i) const;

      //maak een PPHM van een TPM via de T2 conditie
      void T(const TPM &);

      void convert(double **) const;

      static void convert_st(double **);

      static void convert_st2(double **);

      static int gblock_char(int,int);

      static int gpph2s(int,int,int);

      static int gs2pph(int,int,int,int,int);

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
