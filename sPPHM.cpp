#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > *sPPHM::pph2s;
int *****sPPHM::s2pph;

int **sPPHM::block_char;

/**
 * allocate and initialize the statics
 */
void sPPHM::init(){

   //allocate
   pph2s = new vector< vector<int> > [2*Tools::gL()];

   s2pph = new int **** [2*Tools::gL()];

   //S = 1/2
   for(int B = 0;B < Tools::gL();++B){

      s2pph[B] = new int *** [2];

      for(int Z = 0;Z < 2;++Z){

         s2pph[B][Z] = new int ** [Tools::gL()];

         for(int a = 0;a < Tools::gL();++a){

            s2pph[B][Z][a] = new int * [Tools::gL()];

            for(int b = 0;b < Tools::gL();++b)
               s2pph[B][Z][a][b] = new int [Tools::gL()];

         }
      }
   }

   //S = 3/2
   for(int B = Tools::gL();B < 2*Tools::gL();++B){

      s2pph[B] = new int *** [1];

      s2pph[B][0] = new int ** [Tools::gL()];

      for(int a = 0;a < Tools::gL();++a){

         s2pph[B][0][a] = new int * [Tools::gL()];

         for(int b = 0;b < Tools::gL();++b)
            s2pph[B][0][a][b] = new int [Tools::gL()];

      }
   }


   block_char = new int * [2*Tools::gL()];

   for(int B = 0;B < 2*Tools::gL();++B)
      block_char[B] = new int [2];

   int block = 0;

   int pph;

   vector<int> v(4);

   for(int K = 0;K < Tools::gL();++K){

      //S = 1/2
      block_char[block][0] = 0;//means 1/2
      block_char[block][1] = K;

      pph = 0;

      //S_ab = 0
      for(int a = 0;a < Tools::gL();++a)
         for(int b = a;b < Tools::gL();++b)
            for(int c = 0;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                     v[0] = 0;
                     v[1] = a;
                     v[2] = b;
                     v[3] = c;

                     pph2s[block].push_back(v);

                     s2pph[block][0][a][b][c] = pph;

                     ++pph;

                  }

            }

      //S_ab = 1
      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b)
            for(int c = 0;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                  v[0] = 1;
                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  pph2s[block].push_back(v);

                  s2pph[block][1][a][b][c] = pph;

                  ++pph;

               }

            }

      //S = 3/2
      block_char[Tools::gL() + block][0] = 1;//means 3/2
      block_char[Tools::gL() + block][1] = K;

      pph = 0;

      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b)
            for(int c = 0;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                  v[0] = 1;
                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  pph2s[Tools::gL() + block].push_back(v);

                  s2pph[Tools::gL() + block][0][a][b][c] = pph;

                  ++pph;

               }

            }

      ++block;

   }

}

/**
 * static function that deallocates the static lists.
 */
void sPPHM::clear(){

   delete [] pph2s;

   for(int B = 0;B < Tools::gL();++B){

      for(int Z = 0;Z < 2;++Z){

         for(int a = 0;a < Tools::gL();++a){

            for(int b = 0;b < Tools::gL();++b)
               delete [] s2pph[B][Z][a][b];

            delete [] s2pph[B][Z][a];

         }

         delete [] s2pph[B][Z];

      }

      delete [] s2pph[B];

   }

   for(int B = Tools::gL();B < 2*Tools::gL();++B){

      for(int a = 0;a < Tools::gL();++a){

         for(int b = 0;b < Tools::gL();++b)
            delete [] s2pph[B][0][a][b];

         delete [] s2pph[B][0][a];

      }

      delete [] s2pph[B][0];

      delete [] s2pph[B];

   }

   delete [] s2pph;

   for(int B = 0;B < 2*Tools::gL();++B)
      delete [] block_char[B];

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 */
sPPHM::sPPHM() : BlockMatrix(2*Tools::gL()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL();++B)//S = 1/2
      setMatrixDim(B,pph2s[B].size(),2);

   for(int B = Tools::gL();B < 2*Tools::gL();++B)//S = 3/2
      setMatrixDim(B,pph2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with M blocks, M/2 for S=1/2 and M/2 for S=3/2, and copies the content of the pphm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param pphm_c sPPHM object to be copied into (*this)
 */
sPPHM::sPPHM(const sPPHM &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor, if counter = 1 the lists will be deallocated.
 */
sPPHM::~sPPHM(){ }

void sPPHM::trace_x(const PPHM &pphm){

   int a,b,c,d,e,z;
   int S_ab,S_de;

   for(int B = 0;B < 2*Tools::gL();++B){

      int S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         S_ab = pph2s[B][i][0];

         a = pph2s[B][i][1];
         b = pph2s[B][i][2];
         c = pph2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            S_de = pph2s[B][j][0];

            d = pph2s[B][j][1];
            e = pph2s[B][j][2];
            z = pph2s[B][j][3];

            (*this)(B,i,j) = 0.0;

            for(int k_a = 0;k_a < Tools::gL();++k_a)
               for(int k_b = 0;k_b < Tools::gL();++k_b)
                  for(int k_c = 0;k_c < Tools::gL();++k_c)
                     for(int k_d = 0;k_d < Tools::gL();++k_d)
                        for(int k_e = 0;k_e < Tools::gL();++k_e)
                           for(int k_z = 0;k_z < Tools::gL();++k_z)
                              (*this)(B,i,j) += pphm(S,S_ab,Hamiltonian::gxy_a(k_a,a),Hamiltonian::gxy_a(k_b,b),Hamiltonian::gxy_a(k_c,c),
                                    S_de,Hamiltonian::gxy_a(k_d,d),Hamiltonian::gxy_a(k_e,e),Hamiltonian::gxy_a(k_z,z));

            (*this)(B,i,j) /= (double) Tools::gL() * Tools::gL() * Tools::gL();

         }

      }

   }

   this->symmetrize();

}
