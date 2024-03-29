#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>
#include <cstdlib>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 28-02-2011\n\n
 * This is a class written for the 2 dimension hubbard model, it translates the physical degrees of freedom of the
 * 2D hubbard model: (x,y,sigma) to the one particle basis of the regular program.
 */

class Hamiltonian{

   public:

      //initializes the lists.
      static void init();

      //clears the lists;
      static void clear();

      //print the list
      static void print();

      //access the lists from outside the class
      static int gxy_a(int,int);

      //access the lists from outside the class
      static int ga_xy(int,int);

      static int gadjoint(int,int);

      static int gadjoint_sum(int,int,int);

      static int gadjoint(int,int,int);

      static int gadd(int,int,int);

      static int gbar(int);

      static int gadd(int,int);

   private:

      static int L2;

      static int L4;
      
      //!static list that translates the two indices of the 2D hubbard model to one sp index.
      static int **xy_a;

      //!static list that translates the sp index alpha to the two physical indices of the 2D hubbard model.
      static int **a_xy;

      //!static list relating sums of vector-momenta (x,y) to sp-indices
      static int *adjoint2;

      static int *adjoint3;

      static int *add2;

      static int *add3;

      static int *adjoint_sum;

      static int *bar;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
