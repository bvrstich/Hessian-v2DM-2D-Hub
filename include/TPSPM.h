#ifndef TPSPM_H
#define TPSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "RecMat.h"

/**
 * @author Brecht Verstichel
 * @date 14-01-2013\n\n
 * This class TPSPM is a class written for the singly-traced TPTPM matrix object,
 * being a rectangular matrix on TP for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix and adds some special member functions
 */
class TPSPM : public RecMat {

   public:
      
      //constructor
      TPSPM();

      //copy constructor
      TPSPM(const TPSPM &);

      //destructor
      virtual ~TPSPM();

      using RecMat::operator=;

      using RecMat::operator();

      void dpt(double,const TPM &);

      void dpt(double,double **);

      void dpt3(double,double **);

   private:

};

#endif
