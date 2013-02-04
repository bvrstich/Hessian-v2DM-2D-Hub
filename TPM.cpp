#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

vector< vector<int> > *TPM::t2s;

int ***TPM::s2t;

int **TPM::block_char;
int ***TPM::char_block;

/**
 * static function that initializes the static variables
 */
void TPM::init(){

   //allocate stuff
   t2s = new vector< vector<int> > [Tools::gM()];

   s2t = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2t[B] = new int * [Tools::gL2()];

      for(int a = 0;a < Tools::gL2();++a)
         s2t[B][a] = new int [Tools::gL2()];

   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [3];

   char_block = new int ** [2];

   for(int S = 0;S < 2;++S){

      char_block[S] = new int * [Tools::gL()];

      for(int x = 0;x < Tools::gL();++x)
         char_block[S][x] = new int [Tools::gL()];

   }

   vector<int> v(2);

   int block = 0;

   //tp index
   int t;

   //loop over the K_x K_y blocks
   for(int K_x = 0;K_x < Tools::gL();++K_x)
      for(int K_y = 0;K_y < Tools::gL();++K_y){

         t = 0;

         //S = 0
         block_char[block][0] = 0;
         block_char[block][1] = K_x;
         block_char[block][2] = K_y;

         char_block[0][K_x][K_y] = block;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a;b < Tools::gL()*Tools::gL();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == K_x ) 

                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == K_y ) ){

                  v[0] = a;
                  v[1] = b;

                  t2s[block].push_back(v);

                  s2t[block][a][b] = t;
                  s2t[block][b][a] = t;

                  ++t;

               }

            }

         t = 0;

         //S = 1
         block_char[Tools::gL()*Tools::gL() + block][0] = 1;
         block_char[Tools::gL()*Tools::gL() + block][1] = K_x;
         block_char[Tools::gL()*Tools::gL() + block][2] = K_y;

         char_block[1][K_x][K_y] = Tools::gL()*Tools::gL() + block;

         for(int a = 0;a < Tools::gL()*Tools::gL();++a)
            for(int b = a + 1;b < Tools::gL()*Tools::gL();++b){

               if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == K_x ) 

                     && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == K_y ) ){

                  v[0] = a;
                  v[1] = b;

                  t2s[Tools::gL()*Tools::gL() + block].push_back(v);

                  s2t[Tools::gL()*Tools::gL() + block][a][b] = t;
                  s2t[Tools::gL()*Tools::gL() + block][b][a] = t;

                  ++t;

               }

            }

         ++block;

      }

}

/**
 * static function that deallocates the static variables
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL2();++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;

   delete [] block_char;

   for(int S = 0;S < 2;++S){

      for(int x = 0;x < Tools::gL();++x)
         delete [] char_block[S][x];

      delete [] char_block[S];

   }

   delete [] char_block;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice: 
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 */
TPM::TPM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL2();++B)//S = 0
      setMatrixDim(B,t2s[B].size(),1);

   for(int B = Tools::gL2();B < Tools::gM();++B)//S = 1
      setMatrixDim(B,t2s[B].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice:
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K_x,K_y;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K_x = tpm_p.block_char[B][1];
      K_y = tpm_p.block_char[B][2];

      output << "S =\t" << S << "\tK_x =\t" << K_x << "\tK_y =\t" << K_y <<

         "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;

      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(B);++i)
         for(int j = 0;j < tpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[B][i][0] << "\t" << tpm_p.t2s[B][i][1]

               << "\t" << tpm_p.t2s[B][j][0] << "\t" << tpm_p.t2s[B][j][1] << "\t" << tpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The tp spin quantumnumber
 * @param a first sp index that forms the tp row index i of block B(S,K_x,K_y), together with b
 * @param b second sp index that forms the tp row index i of block B(S,K_x,K_y), together with a
 * @param c first sp index that forms the tp column index j of block B(S,K_x,K_y), together with d
 * @param d second sp index that forms the tp column index j of block B(S,K_x,K_y), together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   int K_x = ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%Tools::gL();
   int K_y = ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%Tools::gL();

   //check if momentum is ok
   if( K_x != ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%Tools::gL() )
      return 0;

   if( K_y !=  ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%Tools::gL() )
      return 0;

   int B = char_block[S][K_x][K_y];

   if(S == 0){

      int i = s2t[B][a][b];
      int j = s2t[B][c][d];

      return (*this)(B,i,j);

   }
   else{//S = 1

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[B][a][b];
         int j = s2t[B][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase * (*this)(B,i,j);

      }

   }

}

/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp indices

   double ward = 1.0/(Tools::gN() - 1.0);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j){

               (*this)(B,i,i) = -2.0 * ward * ( cos( 2.0 * Hamiltonian::ga_xy(a,0) * 3.141592653589793238462 / (double) Tools::gL())  

                     + cos( 2.0 * Hamiltonian::ga_xy(a,1) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,0) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,1) * 3.141592653589793238462 / (double) Tools::gL()) );

            }

            //on-site repulsion
            if(B < Tools::gL()*Tools::gL()){//only spin zero

               double hard = 2.0*U / (double) (Tools::gL()*Tools::gL());

               if(a == b)
                  hard /= std::sqrt(2.0);

               if(c == d)
                  hard /= std::sqrt(2.0);

               (*this)(B,i,j) += hard;

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(Tools::gM()*(Tools::gM() - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::S_2() const{

   double ward = 0.0;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}

/**
 * perform a line search what step size in along the Tools::gN()ewton direction is ideal.
 * @param t potential scaling factor
 * @param Z SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &Z,const TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   Z.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(Z,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * access to the lists from outside of the class
 */
int TPM::gt2s(int B,int i,int option){

   return t2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gs2t(int B,int a,int b){

   return s2t[B][a][b];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gblock_char(int B,int option){

   return block_char[B][option];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gdim(int B){

   return t2s[B].size();

}

/**
 * convert a 2DM object from vector to matrix/TPM form
 * @param grad input Gradient object
 */
void TPM::convert(const Gradient &grad){

   int tpmm;
   int S;

   for(int B = 0;B < Tools::gM();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i)
         for(int j = i;j < gdim(B);++j){

            tpmm = TPTPM::gt2tpmm(B,i,j);

            (*this)(B,i,j) = grad[tpmm]/ ( 2.0 * Gradient::gnorm(tpmm) * (2.0 * S + 1.0) );

         }

   }

   this->symmetrize();

}
