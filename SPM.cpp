#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::cout;
using std::endl;

#include "include.h"

/**
 * constructor, makes vector of dimension Tools::gL()*Tools::gL(), the SPM is completely diagonal in momentum space, both in k_x and k_y.
 */
SPM::SPM() {

   spm = new double [Tools::gL2()];
   
}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) {

   spm = new double [Tools::gL2()];

   for(int a = 0;a < Tools::gL2();++a)
      spm[a] = spm_copy[a];
   
}

/**
 * destructor
 */
SPM::~SPM(){ 

   delete [] spm;
   
}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int a = 0;a < Tools::gL2();++a)
      output << a << "\t" << spm_p[a] << endl;

   return output;

}

/**
 * overload the [] operator
 */
double SPM::operator[](int a) const{

   return spm[a];

}

/**
 * overload the [] operator
 */
double &SPM::operator[](int a){

   return spm[a];

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int a = 0;a < Tools::gL2();++a){

      (*this)[a] = 0.0;

      //S = 0
      //b < a
      for(int b = 0;b < a;++b)
         (*this)[a] += tpm(0,a,b,a,b);

      //b == a
      (*this)[a] += 2.0*tpm(0,a,a,a,a);

      //b > a
      for(int b = a + 1;b < Tools::gL()*Tools::gL();++b)
         (*this)[a] += tpm(0,a,b,a,b);

      //S = 1
      for(int b = 0;b < Tools::gL()*Tools::gL();++b)
         (*this)[a] += 3.0*tpm(1,a,b,a,b);

      (*this)[a] *= 0.5 * scale;

   }

}