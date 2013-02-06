#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > *PPHM_ns::pph2s;
int ****PPHM_ns::s2pph0;
int ***PPHM_ns::s2pph1;

int **PPHM_ns::block_char;

/**
 * allocate and initialize the statics
 */
void PPHM_ns::init(){

   //allocate
   pph2s = new vector< vector<int> > [Tools::gM()];

   s2pph0 = new int *** [Tools::gL2()];

   //S = 1/2
   for(int B = 0;B < Tools::gL2();++B){

      s2pph0[B] = new int ** [2];

      for(int Z = 0;Z < 2;++Z){

         s2pph0[B][Z] = new int * [Tools::gL2()];

         for(int a = 0;a < Tools::gL2();++a)
            s2pph0[B][Z][a] = new int [Tools::gL2()];

      }
   }

   s2pph1 = new int ** [Tools::gL2()];

   //S = 3/2
   for(int K = 0;K < Tools::gL2();++K){

      s2pph1[K] = new int * [Tools::gL2()];

      for(int a = 0;a < Tools::gL2();++a)
         s2pph1[K][a] = new int [Tools::gL2()];

   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [2];

   int block = 0;

   int pph;

   vector<int> v(4);

   for(int K = 0;K < Tools::gL2();++K){

      //S = 1/2
      block_char[block][0] = 0;//means 1/2
      block_char[block][1] = K;

      pph = 0;

      //S_ab = 0
      for(int a = 0;a < Tools::gL2();++a)
         for(int b = 0;b < Tools::gL2();++b)
            for(int c = 0;c < Tools::gL2();++c){

               if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == Hamiltonian::ga_xy(K,0) )
                  if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == Hamiltonian::ga_xy(K,1) ){

                     v[0] = 0;
                     v[1] = a;
                     v[2] = b;
                     v[3] = c;

                     pph2s[block].push_back(v);

                     s2pph0[block][0][a][b] = pph;

                     ++pph;

                  }

            }

      //S_ab = 1
      for(int a = 0;a < Tools::gL2();++a)
         for(int b = 0;b < Tools::gL2();++b)
            for(int c = 0;c < Tools::gL2();++c){

               if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == Hamiltonian::ga_xy(K,0) )
                  if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == Hamiltonian::ga_xy(K,1) ){

                     v[0] = 1;
                     v[1] = a;
                     v[2] = b;
                     v[3] = c;

                     pph2s[block].push_back(v);

                     s2pph0[block][1][a][b] = pph;

                     ++pph;

                  }

            }

      //S = 3/2
      block_char[Tools::gL2() + block][0] = 1;//means 3/2
      block_char[Tools::gL2() + block][1] = K;

      pph = 0;

      for(int a = 0;a < Tools::gL2();++a)
         for(int b = 0;b < Tools::gL2();++b)
            for(int c = 0;c < Tools::gL2();++c){

               if( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) + Hamiltonian::ga_xy(c,0))%Tools::gL() == Hamiltonian::ga_xy(K,0) )
                  if( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) + Hamiltonian::ga_xy(c,1))%Tools::gL() == Hamiltonian::ga_xy(K,1) ){

                     v[0] = 1;
                     v[1] = a;
                     v[2] = b;
                     v[3] = c;

                     pph2s[Tools::gL2() + block].push_back(v);

                     s2pph1[block][a][b] = pph;

                     ++pph;

                  }

            }

      ++block;

   }

}

/**
 * static function that deallocates the static lists.
 */
void PPHM_ns::clear(){

   delete [] pph2s;

   for(int B = 0;B < Tools::gL2();++B){

      for(int Z = 0;Z < 2;++Z){

         for(int a = 0;a < Tools::gL2();++a)
            delete [] s2pph0[B][Z][a];

         delete [] s2pph0[B][Z];

      }

      delete [] s2pph0[B];

   }

   for(int K = 0;K < Tools::gL2();++K){

      for(int a = 0;a < Tools::gL2();++a)
            delete [] s2pph1[K][a];

      delete [] s2pph1[K];

   }

   delete [] s2pph0;
   delete [] s2pph1;

   for(int B = 0;B < Tools::gM();++B)
      delete [] block_char[B];

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 */
PPHM_ns::PPHM_ns() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL2();++B)//S = 1/2
      setMatrixDim(B,pph2s[B].size(),2);

   for(int B = Tools::gL2();B < Tools::gM();++B)//S = 3/2
      setMatrixDim(B,pph2s[B].size(),4);

}

/**
 * copy constructor: constructs BlockMatrix object with M blocks, M/2 for S=1/2 and M/2 for S=3/2, and copies the content of the pphm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param pphm_c PPHM_ns object to be copied into (*this)
 */
PPHM_ns::PPHM_ns(const PPHM_ns &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor
 */
PPHM_ns::~PPHM_ns(){ }

ostream &operator<<(ostream &output,const PPHM_ns &pphm_p){

   for(int B = 0;B < pphm_p.gnr();++B){

      int S = pphm_p.gblock_char(B,0);
      int K = pphm_p.gblock_char(B,1);

      int K_x = Hamiltonian::ga_xy(K,0);
      int K_y = Hamiltonian::ga_xy(K,1);

      output << S << "\t" << K_x << "\t" << K_y << "\t" << pphm_p.gdim(B) << "\t" << pphm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(B);++i)
         for(int j = 0;j < pphm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               pphm_p.pph2s[B][i][0] << "\t" << pphm_p.pph2s[B][i][1] << "\t" << pphm_p.pph2s[B][i][2] << "\t" << pphm_p.pph2s[B][i][3] << 

               "\t" << pphm_p.pph2s[B][j][0] << "\t" << pphm_p.pph2s[B][j][1] << "\t" << pphm_p.pph2s[B][j][2] << "\t" << pphm_p.pph2s[B][j][3] 

               << "\t" << pphm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access to the lists from outside of the class
 */
int PPHM_ns::gblock_char(int B,int option){

   return block_char[B][option];

}

/**
 * access to the lists from outside of the class
 */
int PPHM_ns::gpph2s(int B,int i,int option){

   return pph2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int PPHM_ns::gs2pph0(int B,int S_ab,int a,int b){

   return s2pph0[B][S_ab][a][b];

}

/**
 * access to the lists from outside of the class
 */
int PPHM_ns::gs2pph1(int B,int a,int b){

   return s2pph1[B][a][b];

}
