// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Richard Comblen, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_LINEARSYSTEMPETSC_MAD
#define _H_LINEARSYSTEMPETSC_MAD

// -------------------------------------------------------------------
//  Interface to PETSc
// -------------------------------------------------------------------

#ifdef _HAVE_PETSC_

#include "MAdLinearSystem.h"
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include <iostream>
#include <string>

namespace MAd {
  class MAdLinearMatrixPETSc {
  private:
    int * nnz;
    int local_size, global_size;
    bool allocated, assembled;
  public:
    Mat mat;
    MAdLinearMatrixPETSc(int _local_size, int _global_size) {
      local_size  = _local_size;
      global_size = _global_size;
      nnz = new int[local_size];
      mat = NULL;
      allocated = false;
      assembled = false;
    }
    ~MAdLinearMatrixPETSc() {
      if(mat) MatDestroy(mat);
      delete [] nnz;
    }
    void set_nnz(int row, int nz){
      nnz[row] = nz;
    }
    // You must call set_nnz before!
    void allocate(){
      if(!allocated){
        MatCreate(PETSC_COMM_WORLD, &mat);
        MatSetSizes(mat,local_size,local_size,global_size,global_size);
        MatSetFromOptions(mat);
        //    MatMPIAIJSetPreallocationCSR(mat,sparsity->rowstart,sparsity->colind,NULL);
        MatSeqAIJSetPreallocation(mat,0,nnz);
        allocated=true;
      }
    }
    void add(int row, int col, double val){ 
      MatSetValue(mat, row, col, val, ADD_VALUES);
      assembled = false;
    }
    void assemble()
    {
      if(assembled) return;
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd  (mat, MAT_FINAL_ASSEMBLY);
#if PETSC_VERSION_MAJOR==3
      MatSetOption(mat, MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
      MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
#else
      MatSetOption(mat, MAT_NO_NEW_NONZERO_LOCATIONS);
#endif 
      assembled = true;
    }
    void zero(){
      //assemble();
      MatZeroEntries(mat);
    }
    void print(std::string name)
    {
      printf("\nPrinting matrix %s:\n\n",name.c_str());
      if ( !allocated ) {printf("Not allocated!\n"); return;}
      if ( !assembled ) {printf("Not assembled!\n"); return;}

      MatView(mat,PETSC_VIEWER_STDOUT_WORLD);
    }
  };
  class MAdLinearVectorPETSc {
  private:
    int local_size, global_size, nghosts;
    double * entries;
  public:
    Vec vec;
    MAdLinearVectorPETSc(int _local_size, int _global_size){
      global_size = _global_size;
      local_size  = _local_size;
      nghosts = 0;
      entries = new double[local_size+nghosts];
      VecCreateMPIWithArray(MPI_COMM_WORLD,local_size,global_size,entries,&vec);
      VecAssemblyBegin(vec);
      VecAssemblyEnd(vec);
      VecSet(vec, 0);
      get_array();
    }
    ~MAdLinearVectorPETSc(){
      VecDestroy(vec);
      delete []entries;// Check whether Petsc does destroy it or not...
    }
    inline double& operator()(int row){
      return entries[row];
    }
    void zero(){
      for(int i=0;i<local_size+nghosts;i++){
        entries[i]=0;
      }
    }
    void gather(){
      // To be implemented - see slim_vector.cc
    }
    void scatter(){
      // To be implemented - see slim_vector.cc
    }
    void get_array(){
      VecGetArray(vec,&entries);
    }
    void restore_array(){
      VecRestoreArray(vec,&entries);
    }
    void print(std::string name)
    {
      printf("\nPrinting vector %s:\n\n",name.c_str());
      VecView(vec,PETSC_VIEWER_STDOUT_WORLD);
    }
  };

  // -------------------------------------------------------------------
  class MAdLinearSystemPETSc : public MAdLinearSystem {

  private:
    MAdLinearMatrixPETSc * matrix;
    MAdLinearVectorPETSc * rhs;
    MAdLinearVectorPETSc * X;
    KSP ksp;

  public:

    MAdLinearSystemPETSc():
      MAdLinearSystem()
    {
      PetscInitialize(0,NULL, NULL, NULL);
      KSPCreate(PETSC_COMM_WORLD, &ksp);
    }
    ~MAdLinearSystemPETSc ()
    {
      KSPDestroy(ksp);
      delete matrix;
      delete rhs;
      delete X;
    }
    bool isAllocated () const {throw;}
    void allocate (int nbRows){
      int local_size, global_size;
      local_size = global_size = nbRows;
      matrix = new MAdLinearMatrixPETSc(local_size, global_size);
      rhs    = new MAdLinearVectorPETSc(local_size, global_size);
      X      = new MAdLinearVectorPETSc(local_size, global_size);
    }
    void addToMatrix (int row, int col, double val) {
      matrix->add(row,col,val);
    }
    double getFromMatrix (int _row, int _col) const {
      throw;
    }
    void addToRightHandSide (int row, double val) {
      (*rhs)(row) += val;
    }
    double getFromRightHandSide (int row) const {
      return (*rhs)(row);
    }
    double getFromSolution (int row) const {
      return (*X)(row);
    }
    void zeroMatrix () {
      matrix->zero();
    }
    void zeroRightHandSide () {
      rhs->zero();
    }
    void reorder () {
    }
    void set_nnz(int row,int nz){
      matrix->set_nnz(row,nz);
    }
    void allocate_matrix(){
      matrix->allocate();
    }
    void setSolver(SolverType _type) {
    }
    int systemSolve () {
      rhs->restore_array();
      X->restore_array();
      matrix->assemble();

// #warning "debug"
//       printMatrix();
//       printRhs();

      std::string options;
      options="-pc_type bjacobi -sub_pc_type icc -ksp_type cg -ksp_rtol 1e-6";
//      options="-pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_print_statistics -ksp_type gmres -ksp_rtol 1e-7 -ksp_monitor";
      PetscOptionsInsertString(options.c_str());
      KSPSetFromOptions(ksp);
      KSPSetOperators (ksp, matrix->mat, matrix->mat, DIFFERENT_NONZERO_PATTERN);
      KSPSolve(ksp, rhs->vec, X->vec);
      
      X->get_array();
      rhs->get_array();
      KSPConvergedReason reason;
      return KSPGetConvergedReason(ksp,&reason);
    }

    void printMatrix(std::string name="")
    {
      matrix->print("Left hand size");
    }

    void printRhs(std::string name="")
    {
      rhs->print("Right hand size");
    }

  };
}

#endif

#endif
