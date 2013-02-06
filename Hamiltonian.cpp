#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int **Hamiltonian::xy_a;
int **Hamiltonian::a_xy;

int *Hamiltonian::bar;

int *Hamiltonian::adjoint2;
int *Hamiltonian::adjoint3;

int *Hamiltonian::add2;
int *Hamiltonian::add3;

int *Hamiltonian::adjoint_sum;

/**
 * function that allocates and constructs the lists.
 */
void Hamiltonian::init(){

   //allocate
   xy_a = new int * [Tools::gL()];

   for(int x = 0;x < Tools::gL();++x)
      xy_a[x] = new int [Tools::gL()];

   a_xy = new int * [Tools::gL2()];

   for(int i = 0;i < Tools::gL2();++i)
      a_xy[i] = new int [2];

   //construct:
   int a = 0;

   for(int x = 0;x < Tools::gL();++x)
      for(int y = 0;y < Tools::gL();++y){

         a_xy[a][0] = x;
         a_xy[a][1] = y;

         xy_a[x][y] = a;

         a++;

      }

   int L2 = Tools::gL2();

   //construct the list which finds the sp index which combines with sp-index k to sum up to sp index K.
   bar = new int [L2];

   for(int a = 0;a < L2;++a)
      bar[a] = xy_a[(-a_xy[a][0] + Tools::gL())%Tools::gL()][(-a_xy[a][1] + Tools::gL())%Tools::gL()];

   adjoint2 = new int [L2*L2];
   
   for(int K = 0;K < L2;++K)
      for(int k = 0;k < L2;++k)
         adjoint2[K + L2*k] = xy_a[(a_xy[K][0] - a_xy[k][0] + Tools::gL())%Tools::gL()][(a_xy[K][1] - a_xy[k][1] + Tools::gL())%Tools::gL()];

   adjoint_sum = new int [L2*L2*L2];

   for(int a = 0;a < L2;++a)
      for(int k = 0;k < L2;++k)
         for(int e = 0;e < L2;++e){

            adjoint_sum[a + k*L2 + e*L2*L2] = 
            
               xy_a[(a_xy[a][0] + a_xy[k][0] - a_xy[e][0] + Tools::gL())%Tools::gL()][(a_xy[a][1] + a_xy[k][1] - a_xy[e][1] + Tools::gL())%Tools::gL()];

         }

   adjoint3 = new int [L2*L2*L2];

   for(int K = 0;K < L2;++K)
      for(int a = 0;a < L2;++a)
         for(int b = 0;b < L2;++b){

            adjoint3[K + a*L2 + b*L2*L2] = 
            
               xy_a[(a_xy[K][0] - a_xy[a][0] - a_xy[b][0] + 2*Tools::gL())%Tools::gL()][(a_xy[K][1] - a_xy[a][1] - a_xy[b][1] + 2*Tools::gL())%Tools::gL()];

         }

   add2 = new int [L2*L2];

   for(int a = 0;a < L2;++a)
      for(int b = 0;b < L2;++b)
         add2[a + b*L2] = xy_a[(a_xy[a][0] + a_xy[b][0])%Tools::gL()][(a_xy[a][1] + a_xy[b][1])%Tools::gL()];

   add3 = new int [L2*L2*L2];

   for(int a = 0;a < L2;++a)
      for(int b = 0;b < L2;++b)
         for(int c = 0;c < L2;++c)
            add3[a + b*L2 + c*L2*L2] = xy_a[(a_xy[a][0] + a_xy[b][0] + a_xy[c][0])%Tools::gL()][(a_xy[a][1] + a_xy[b][1] + a_xy[c][1])%Tools::gL()];

}

/**
 * deallocates the lists.
 */
void Hamiltonian::clear(){

   //delete xy_a
   for(int x = 0;x < Tools::gL();++x)
      delete [] xy_a[x];

   delete [] xy_a;

   //delete a_xy
   for(int i = 0;i < Tools::gL2();++i)
      delete [] a_xy[i];

   delete [] a_xy;

   delete [] bar;
   delete [] adjoint2;
   delete [] adjoint_sum;
   delete [] adjoint3;
   delete [] add2;
   delete [] add3;

}

/**
 * print the list
 */
void Hamiltonian::print(){

   for(int i = 0;i < Tools::gL2();++i)
      std::cout << i << "\t" << a_xy[i][0]<< "\t" << a_xy[i][1] << std::endl;

}

/**
 * access to the a_xy list from outside this class
 * @param a the sp-index
 * @param option can be 0 or 1
 * @return k_x if option == 0 , k_y if option == 1
 */
int Hamiltonian::ga_xy(int a,int option){

   return a_xy[a][option];

}

/**
 * access to the list xy_a from outside this class
 * @param k_x the x momentum
 * @param k_y the y momentum
 * @return the sp index corresponding to k_x and k_y
 */
int Hamiltonian::gxy_a(int k_x,int k_y){
   
   return xy_a[k_x][k_y];

}

/**
 * transform the "particle momentum" sp-index to the "hole momentum" sp-index
 * @param a the input sp-index
 */
int Hamiltonian::gbar(int a){

   return bar[a];

}

/**
 * find the sp index which combines with sp-index k to sum up to sp index K.
 */
int Hamiltonian::gadjoint(int K,int k){

  return adjoint2[K + k*Tools::gL2()];

}

/**
 * find the sp index which combines with sp-index e to sum up to sp index formed by (a + k).
 */
int Hamiltonian::gadjoint_sum(int a,int k,int e){

  return adjoint_sum[a + k*Tools::gL2() + e*Tools::gL2()*Tools::gL2()];

}

/**
 * add two sp indices together to form a new sp index (i.e. add the separate x and y momenta and recombine)
 */
int Hamiltonian::gadd(int a,int b){

   return add2[a + b*Tools::gL2()];

}

/**
 * find the sp index which combines with sp-indices a and b to sum up to sp index K.
 */
int Hamiltonian::gadjoint(int K,int a,int b){

  return adjoint3[K + a*Tools::gL2() + b*Tools::gL2()*Tools::gL2()];

}

/**
 * add three sp indices together to form a new sp index (i.e. add the separate x and y momenta and recombine)
 */
int Hamiltonian::gadd(int a,int b,int c){

   return add3[a + b*Tools::gL2() + c*Tools::gL2()*Tools::gL2()];

}
