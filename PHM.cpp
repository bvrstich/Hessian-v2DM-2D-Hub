#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::endl;

#include "include.h"

vector< vector<int> > *PHM::ph2s;
int ***PHM::s2ph;

int **PHM::block_char;

/**
 * initializes the statics
 */
void PHM::init(){

   //allocate stuff
   ph2s = new vector< vector<int> > [Tools::gM()];

   s2ph = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2ph[B] = new int * [Tools::gL2()];

      for(int a = 0;a < Tools::gL2();++a)
         s2ph[B][a] = new int [Tools::gL2()];

   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [2];

   vector<int> v(2);

   int block = 0;

   //ph index
   int ph;

   //loop over the K_x K_y blocks
   for(int K = 0;K < Tools::gL2();++K){

         ph = 0;

         //S = 0
         block_char[block][0] = 0;
         block_char[block][1] = K;

         //S = 1
         block_char[Tools::gL2() + block][0] = 1;
         block_char[Tools::gL2() + block][1] = K;

         for(int a = 0;a < Tools::gL2();++a)
            for(int b = 0;b < Tools::gL2();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == Hamiltonian::ga_xy(K,0) )
                  
                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == Hamiltonian::ga_xy(K,1) ) ){

                  v[0] = a;
                  v[1] = b;

                  ph2s[block].push_back(v);//S = 0
                  ph2s[Tools::gL2() + block].push_back(v);//S = 1

                  s2ph[block][a][b] = ph;
                  s2ph[Tools::gL2() + block][a][b] = ph;

                  ++ph;

               }

            }

         ++block;

      }

}

/**
 * static function that deallocates the static variables
 */
void PHM::clear(){

   delete [] ph2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL2();++a)
         delete [] s2ph[B][a];

      delete [] s2ph[B];

      delete [] block_char[B];

   }

   delete [] s2ph;

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 */
PHM::PHM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL2();++B)//S = 0
      setMatrixDim(B,ph2s[B].size(),1);

   for(int B = Tools::gL2();B < Tools::gM();++B)//S = 1
      setMatrixDim(B,ph2s[B].size(),3); 

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks of dimension M*M/4 and copies the content of phm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor
 */
PHM::~PHM(){ }

ostream &operator<<(ostream &output,const PHM &phm_p){

   int S,K_x,K_y;

   for(int B = 0;B < phm_p.gnr();++B){

      S = phm_p.block_char[B][0];

      K_x = Hamiltonian::ga_xy(phm_p.block_char[B][1],0);
      K_y = Hamiltonian::ga_xy(phm_p.block_char[B][1],1);

      output << "S =\t" << S << "\tK_x =\t" << K_x << "\tK_y =\t" << K_y << "\tdimension =\t" << phm_p.gdim(B) << "\tdegeneracy =\t" << phm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(B);++i)
         for(int j = 0;j < phm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[B][i][0] << "\t" << phm_p.ph2s[B][i][1]

               << "\t" << phm_p.ph2s[B][j][0] << "\t" << phm_p.ph2s[B][j][1] << "\t" << phm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access the elements of the matrix in sp mode, 
 * @param S The tp spin quantumnumber
 * @param a first sp index that forms the ph row index i in block B together with b
 * @param b second sp index that forms the ph row index i in block B together with a
 * @param c first sp index that forms the ph column index j in block B together with d
 * @param d second sp index that forms the ph column index j in block B together with c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int K_x =  ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%Tools::gL();
   int K_y =  ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%Tools::gL();

   //check if momentum is conserved
   if( K_x != ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%Tools::gL() )
      return 0;

   if( K_y != ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%Tools::gL() ) 
      return 0;

   int B = Tools::gL2()*S + Hamiltonian::gxy_a(K_x,K_y);

   int i = s2ph[B][a][b];
   int j = s2ph[B][c][d];

   return (*this)(B,i,j);

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = ph2s[B][i][0];

         //transform k_b to tpm sp-momentum:
         b = Hamiltonian::bar(ph2s[B][i][1]);

         for(int j = i;j < gdim(B);++j){

            c = ph2s[B][j][0];

            //transform k_d to tpm sp-momentum:
            d = Hamiltonian::bar(ph2s[B][j][1]);

            (*this)(B,i,j) = - Tools::g6j(0,0,0,S) * tpm(0,a,d,c,b) - 3.0 * Tools::g6j(0,0,1,S) * tpm(1,a,d,c,b);

            if(a == d)
               (*this)(B,i,j) *= std::sqrt(2.0);

            if(b == c)
               (*this)(B,i,j) *= std::sqrt(2.0);

         }

         (*this)(B,i,i) += spm[a];

      }

   }

   this->symmetrize();

}

/**
 * convert a PHM matrix to a double array for faster access to the number
 */
void PHM::convert(double **array) const {

   int K,S;

   int L2 = Tools::gL2();

   for(int B = 0;B < 2*L2;++B){

      S = block_char[B][0];
      K = block_char[B][1];

      for(int a = 0;a < L2;++a)
         for(int c = 0;c < L2;++c)
            array[B][a + c*L2] = (*this)(S,a,Hamiltonian::adjoint(K,a),c,Hamiltonian::adjoint(K,c));

   }

}

/**
 * access to the lists from outside of the class
 */
int PHM::gph2s(int B,int i,int option){

   return ph2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int PHM::gs2ph(int B,int a,int b){

   return ph2s[B][a][b];

}

/**
 * access to the lists from outside of the class
 */
int PHM::gblock_char(int B,int option){

   return block_char[B][option];

}

/* vim: set ts=3 sw=3 expandtab :*/
