#ifndef sDPM_H
#define sDPM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 04-03-2011\n\n
 */
class sDPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpm_p the sDPM you want to print
    */
   friend ostream &operator<<(ostream &output,const sDPM &dpm_p);

   public:
      
      //constructor
      sDPM();

      //copy constructor
      sDPM(const sDPM &);

      //destructor
      virtual ~sDPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();
      
      void trace_x(const DPM &);

      static void init();

      static void clear();

   private:

      //!static list that takes in a dp index i for block B and returns an intermediate spin: S_ab, and three sp indices: a, b and c
      static vector< vector<int> > *dp2s;

      //!static list that takes a block index B, an intermediate spin-index S_ab and three sp indices a,b and c "and" returns a dp index i for block B:
      static int *****s2dp;

      //!static list that takes a blockindex B and returns the dp spin S and the dp momenta K_x and K_y.
      static int **block_char;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
