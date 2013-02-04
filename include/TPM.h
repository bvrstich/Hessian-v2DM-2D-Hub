#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::ifstream;
using std::vector;

#include "BlockMatrix.h"

class Gradient;
class PHM;
class DPM;
class SUP;

/**
 * @author Brecht Verstichel
 * @date 10-05-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry and translational symemtry included,
 * It inherits alle the function from its mother BlockMatrix.
 * Some special member functions and two lists that give the relationship between the sp and the tp basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and with tp spin and momentum quantumnumbers
      double operator()(int S,int a,int b,int c,int d) const;

      void hubbard(double U);

      void unit();

      double S_2() const;

      double line_search(double,SUP &,const TPM &);

      void convert(const Gradient &);

      void Q(int,double,double,double,const TPM &);

      void Q(int,const TPM &);

      void G(const PHM &);

      void T(const DPM &);

      void bar(double,const DPM &);

      static int gt2s(int,int,int);

      static int gs2t(int,int,int);

      static int gblock_char(int,int);

      static int gdim(int);

      static double gnorm(int,int);

      static void init();

      static void clear();

   private:

      //!static list that takes in a tp index i and a blockindex B, and returns two sp indices: a = t2s[B][i][0] and b = t2s[B][i][1]
      static vector< vector<int> > *t2s;

      //!static list that takes two sp indices a,b and a blockindex B, and returns a tp index i: i = s2t[B][a][b]
      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S and the tp momenta K_x
      static int **block_char;

      static double **norm;

};

#endif
