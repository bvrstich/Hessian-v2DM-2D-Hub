#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;
using std::vector;

#include "include.h"

vector< vector<int> > *sDPM::dp2s;
int *****sDPM::s2dp;

int **sDPM::block_char;

/**
 * initialize the statics and allocate and construct all the lists
 */
void sDPM::init(){

   dp2s = new vector< vector<int> > [2*Tools::gL()];

   s2dp = new int **** [2*Tools::gL()];

   //S = 1/2 part
   for(int B = 0;B < Tools::gL();++B){

      s2dp[B] = new int *** [2];

      for(int S = 0;S < 2;++S){

         s2dp[B][S] = new int ** [Tools::gL()];

         for(int a = 0;a < Tools::gL();++a){

            s2dp[B][S][a] = new int * [Tools::gL()];

            for(int b = 0;b < Tools::gL();++b)
               s2dp[B][S][a][b] = new int [Tools::gL()];

         }
      }
   }

   //S = 3/2 part
   for(int B = Tools::gL();B < 2*Tools::gL();++B){

      s2dp[B] = new int *** [1];

      s2dp[B][0] = new int ** [Tools::gL()];

      for(int a = 0;a < Tools::gL();++a){

         s2dp[B][0][a] = new int * [Tools::gL()];

         for(int b = 0;b < Tools::gL();++b)
            s2dp[B][0][a][b] = new int [Tools::gL()];

      }
   }

   block_char = new int * [2*Tools::gL()];

   for(int B = 0;B < 2*Tools::gL();++B)
      block_char[B] = new int [2];

   int block = 0;

   vector<int> v(4);

   int dp;

   //loop over the blocks
   for(int K = 0;K < Tools::gL();++K){

      //S = 1/2
      block_char[block][0] = 0;//0 means spin 1/2
      block_char[block][1] = K;

      dp = 0;

      //first S = 1/2, S_ab = 0, a = b != c
      for(int a = 0;a < Tools::gL();++a){

         for(int b = 0;b < a;++b){

            if( (2*a + b)%Tools::gL() == K ){

               v[0] = 0;//S_ab
               v[1] = a;
               v[2] = a;
               v[3] = b;

               dp2s[block].push_back(v);

               s2dp[block][0][a][a][b] = dp;

               ++dp;

            }

         }

         for(int b = a + 1;b < Tools::gL();++b){

            if( (2*a + b)%Tools::gL() == K ){

               v[0] = 0;//S_ab
               v[1] = a;
               v[2] = a;
               v[3] = b;

               dp2s[block].push_back(v);

               s2dp[block][0][a][a][b] = dp;

               ++dp;

            }

         }

      }

      //then S = 1/2, S_ab = 0, a < b < c
      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b)
            for(int c = b + 1;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                  v[0] = 0;//S_ab
                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[block].push_back(v);

                  s2dp[block][0][a][b][c] = dp;

                  ++dp;

               }

            }

      //then S = 1/2, S_ab = 1, a < b < c
      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b)
            for(int c = b + 1;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                  v[0] = 1;//S_ab
                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[block].push_back(v);

                  s2dp[block][1][a][b][c] = dp;

                  ++dp;

               }

            }

      //S = 3/2
      block_char[Tools::gL() + block][0] = 1;//0 means spin 1/2
      block_char[Tools::gL() + block][1] = K;

      dp = 0;

      //then S = 3/2, S_ab = 1, a < b < c
      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + 1;b < Tools::gL();++b)
            for(int c = b + 1;c < Tools::gL();++c){

               if( (a + b + c)%Tools::gL() == K ){

                  v[0] = 1;//S_ab
                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[Tools::gL() + block].push_back(v);

                  s2dp[Tools::gL() + block][0][a][b][c] = dp;

                  ++dp;

               }

            }

      ++block;

   }

}

/**
 * deallocate the static lists
 */
void sDPM::clear(){

   delete [] dp2s;

   for(int B = 0;B < Tools::gL();++B){

      for(int S = 0;S < 2;++S){

         for(int a = 0;a < Tools::gL();++a){

            for(int b = 0;b < Tools::gL();++b)
               delete [] s2dp[B][S][a][b];

            delete [] s2dp[B][S][a];

         }

         delete [] s2dp[B][S];

      }

      delete [] s2dp[B];

   }

   for(int B = Tools::gL();B < 2*Tools::gL();++B){

      for(int a = 0;a < Tools::gL();++a){

         for(int b = 0;b < Tools::gL();++b)
            delete [] s2dp[B][0][a][b];

         delete [] s2dp[B][0][a];

      }

      delete [] s2dp[B][0];

      delete [] s2dp[B];

   }


   delete [] s2dp;

   for(int B = 0;B < 2*Tools::gL();++B)
      delete [] block_char[B];

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 * (Tools::gL()*Tools::gL()) blocks, (Tools::gL()*Tools::gL()) for S = 1/2 and (M/2) 3/2.
 */
sDPM::sDPM() : BlockMatrix(2*Tools::gL()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL();++B)//S = 1/2
      setMatrixDim(B,dp2s[B].size(),2);

   for(int B = Tools::gL();B < 2*Tools::gL();++B)//S = 3/2
      setMatrixDim(B,dp2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with 2 * (M/2) blocks, on (M/2) for S=1/2 and (M/2) for S=3/2,
 * and copies the content of the dpm_c blocks into it.
 * @param dpm_c sDPM to be copied into (*this)
 */
sDPM::sDPM(const sDPM &dpm_c) : BlockMatrix(dpm_c) { }

/**
 * destructor: 
 */
sDPM::~sDPM(){ }

void sDPM::trace_x(const DPM &dpm){

   int a,b,c,d,e,z;
   int S_ab,S_de;

   for(int B = 0;B < 2*Tools::gL();++B){

      int S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         S_ab = dp2s[B][i][0];

         a = dp2s[B][i][1];
         b = dp2s[B][i][2];
         c = dp2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            S_de = dp2s[B][j][0];

            d = dp2s[B][j][1];
            e = dp2s[B][j][2];
            z = dp2s[B][j][3];

            (*this)(B,i,j) = 0.0;

            for(int k_a = 0;k_a < Tools::gL();++k_a)
               for(int k_b = 0;k_b < Tools::gL();++k_b)
                  for(int k_c = 0;k_c < Tools::gL();++k_c)
                     for(int k_d = 0;k_d < Tools::gL();++k_d)
                        for(int k_e = 0;k_e < Tools::gL();++k_e)
                           for(int k_z = 0;k_z < Tools::gL();++k_z)
                              (*this)(B,i,j) += dpm(S,S_ab,Hamiltonian::gxy_a(k_a,a),Hamiltonian::gxy_a(k_b,b),Hamiltonian::gxy_a(k_c,c),
                                    S_de,Hamiltonian::gxy_a(k_d,d),Hamiltonian::gxy_a(k_e,e),Hamiltonian::gxy_a(k_z,z));

            (*this)(B,i,j) /= (double) Tools::gL() * Tools::gL() * Tools::gL();

         }

      }

   }

   this->symmetrize();

}
