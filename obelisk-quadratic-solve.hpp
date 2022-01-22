#ifndef _OBELISK_QUADRATIC_SOLVE_
#define _OBELISK_QUADRATIC_SOLVE_

#include <map>

#include <cfloat>

#include <tgmath.h>

#include <assert.h>

#include <coin/ClpSimplex.hpp>

#include <x86_64-linux-gnu/cblas.h>

#include "obelisk-compressed-rrmat.hpp"

#include "obelisk-quadratic-solve.hpp"

using namespace std;

int swap_pointers(double** x, double** y);
int get_rr_inv_pair(double** input_mat, double** inv, compressed_rrmat* lhs, int dim);
int print_nz(int size, double* array);
int column_collapse(compressed_rrmat* lhs, map <int, double*>::iterator lhs_col, double** inv, int inv_cols);
int active_set_alg(compressed_rrmat* lhs, int num_rows, double** inv, double* solution, double* kkt_rhs, int variables, int constraints);
int update_inv(compressed_rrmat* lhs, int* num_rows, double** inv, int inv_cols, int new_active);

class quadratic_programming_problem {
  public:
    pair<int, double>** problem_matrix;
    int* constraint_cols;
    double* rhs;
    double* quadratic_coefficients;
    double* linear_coefficients;
    int variables;
    int constraints;
    
    int solve(double* this_solution) {

      // Generate a linear programming problem to find an interior point for the quadratic
      // problem solve. It will have a single additional variable t and constraints requiring
      // all other problem and slack variables to be greater than t. Maximize t.
  
      /* If uncommented, print the tech matrix for debugging instead.
  
      double** print_matrix = new double*[this->constraints];
      for (int i = 0; i < this->constraints; i++) {
        print_matrix[i] = new double[this->variables];
        for (int j = 0; j < this->variables; j++) {
          print_matrix[i][j] = 0;
        }
      }
  
      for (int i = 0; i < this->constraints; i++) {
        for (int j = 0; j < this->constraint_cols[i]; j++) {
          print_matrix[i][this->problem_matrix[i][j].first] = this->problem_matrix[i][j].second;
        }
      }
  
      for (int i = 0; i < this->constraints; i++) {
        for (int j = 0; j < this->variables; j++) {
          cout << "\t" << print_matrix[i][j] << "\t";
        }
        cout << "\n";
      }
  
      return 0;
  
      */
  
      double* init_solution_rhs_l = new double[this->constraints + this->variables];
      memset(init_solution_rhs_l, 0, sizeof(double) * (this->constraints + this->variables) );
      double* init_solution_rhs_u = new double[this->constraints+this->variables];
      memset(init_solution_rhs_u, 0, sizeof(double) * (this->constraints + this->variables) );
      double* init_solution_objs = new double[this->variables+1];
      memset(init_solution_objs, 0, sizeof(double) * (this->variables + 1) );
 
      CoinPackedMatrix init_solution_mat;
      init_solution_mat.setDimensions(this->constraints+this->variables, this->variables+1);
  
      for (int i = 0; i < this->constraints; i++) {
        for (int j = 0; j < this->constraint_cols[i]; j++) {
         init_solution_mat.modifyCoefficient(i,this->problem_matrix[i][j].first,this->problem_matrix[i][j].second);
        }
      }
    
      for (int i = 0; i < this->variables; i++) {
        init_solution_mat.modifyCoefficient(this->constraints+i,i,1);
        init_solution_mat.modifyCoefficient(this->constraints+i,this->variables,-1);
      }
  
      for (int i = 0; i < this->constraints; i++) {
        init_solution_rhs_l[i] = this->rhs[i];
        init_solution_rhs_u[i] = this->rhs[i];
      }
  
      for (int i = this->constraints; i < this->constraints+this->variables; i++) {
        init_solution_rhs_l[i] = 0;
        init_solution_rhs_u[i] = INFINITY;
      }
  
      init_solution_objs[this->variables] = 1;
  
      ClpSimplex init_solution_model;
  
      init_solution_model.setOptimizationDirection(-1);
      init_solution_model.loadProblem(init_solution_mat,NULL,
        NULL,init_solution_objs,init_solution_rhs_l,init_solution_rhs_u);
  
      init_solution_model.primal();
      double* init_solution = init_solution_model.primalColumnSolution();
  
      if (init_solution_model.isProvenPrimalInfeasible()) {
  
        double* this_infeasibles = init_solution_model.infeasibilityRay(0);
    
        double* zero_rhs = new double[this->variables+this->constraints];
        memset(zero_rhs, 0, sizeof(double) * (this->constraints + this->variables) );
    
        /*init_solution_mat.times(init_solution,zero_rhs);
  
        cout << "\nInitial solution find reports your problem is infeasible. Now reporting infeasibilities variable.\n";
  
        for (int j = 0; j < 1+variables; j++) {
          if (this_infeasibles[j] > 1.0e-300 || this_infeasibles[j] < -1.0e-300) {
            cout << j << ": " << this_infeasibles[j] << "\n";
            cout << "Corresponding RHS: " << init_solution_rhs_l[j] << "\n";
          }
        }
    
        cout << "\nRHS of zero solution.\n";
  
        for (int j = 0; j < constraints+variables; j++) {
          if (1) {
            cout << j << ": " << zero_rhs[j] << "\n";
            cout << "Corresponding RHS: " << init_solution_rhs_u[j] << "\n";
            cout << "Status: " << init_solution_model.getRowStatus(j) << "\n";
          }
        }
        */
  
        return 1;
      }
  
      if (init_solution_model.isProvenDualInfeasible()) {
  
        cout << "\nInitial solution find reports your problem is unbounded. Let's find out how.\n";
    
        double* feasible_solution_rhs_l = new double[this->constraints];
        memset(feasible_solution_rhs_l, 0, sizeof(double) * (this->constraints) );
        double* feasible_solution_rhs_u = new double[this->constraints];
        memset(feasible_solution_rhs_u, 0, sizeof(double) * (this->constraints) );
        double* feasible_solution_objs = new double[this->variables];
    
        for (int i = 0; i < this->variables; i++) {
          feasible_solution_objs[i] = 1;
        }
  
        CoinPackedMatrix feasible_solution_mat;
        feasible_solution_mat.setDimensions(this->constraints, this->variables);
  
        for (int i = 0; i < this->constraints; i++) {
          for (int j = 0; j < this->constraint_cols[i]; j++) {
           feasible_solution_mat.modifyCoefficient(i,this->problem_matrix[i][j].first,this->problem_matrix[i][j].second);
          }
        }
  
        for (int i = 0; i < constraints; i++) {
          feasible_solution_rhs_l[i] = this->rhs[i];
          feasible_solution_rhs_u[i] = this->rhs[i];
        }
  
        ClpSimplex feasible_solution_model;
  
        feasible_solution_model.setOptimizationDirection(-1);
        feasible_solution_model.loadProblem(feasible_solution_mat,NULL,
          NULL,feasible_solution_objs,feasible_solution_rhs_l,feasible_solution_rhs_u);
    
        feasible_solution_model.primal();
    
        if (!feasible_solution_model.isProvenDualInfeasible()) {
          cout << "\nIt looks bounded or infeasible now. This is a glitch.\n";
          return 1;
        }
    
        double* this_infeasibles = feasible_solution_model.unboundedRay();
        cout << "\nNow reporting unbounded ray by variable.\n";
  
        for (int j = 0; j < constraints+variables; j++) {
          if (this_infeasibles[j] > 1.0e-300 || this_infeasibles[j] < -1.0e-300) {
            cout << j << ": " << this_infeasibles[j] << "\n";
            //cout << "Corresponding RHS: " << feasible_solution_rhs_l[j] << "\n";
          }
        }
    
        return 1;
  
      }
  
      // for (int i = 0; i < this->variables+this->constraints; i++) cout << "  " << rhs[i];

      // Create the initial matrix pair, namely the KKT matrix and a diagonal matrix called inv.

      // Allocating space enough for every possible future pointer.
      double** inv = new double*[this->constraints+2*this->variables];
      double** kkt = new double*[this->constraints+this->variables];
  
      for (int i = 0; i < this->constraints+this->variables; i++) {
        kkt[i] = new double[this->constraints+this->variables];
        inv[i] = new double[this->constraints+this->variables];
        memset(inv[i], 0, sizeof(double) * (this->constraints+this->variables) );
        memset(kkt[i], 0, sizeof(double) * (this->constraints+this->variables) );
        inv[i][i] = 1.0;
      }
  
      // Fill initial KKT matrix.
  
      for (int i = 0; i < variables; i++) {
        if (quadratic_coefficients[i] != 0) {
          kkt[i][i] = -2*quadratic_coefficients[i];
        }
      }

      for (int i = 0; i < this->constraints; i++) {
        for (int j = 0; j < this->constraint_cols[i]; j++) {
          kkt[this->problem_matrix[i][j].first][i + this->variables] = -this->problem_matrix[i][j].second;
          kkt[i + this->variables][this->problem_matrix[i][j].first] = this->problem_matrix[i][j].second;
        }
      }
  
      /* Now we find the "inverse" of the KKT matrix. This allows us to find a solution to the
         initial KKT system *and* to find the solution to the next (one constraint larger) KKT system
         with only linear additional row operations. The whole algorithm should require only O(n^3) 
         scalar operations as a result. 
     
         I put "inverse" in quotes because kkt may be singular, in which a case it will
         only be row reduced echlon at the end of elimination rather than the identity. This still
         allows us to get a solution with matrix-vector multiplication, as long as a solution exists.
      */
  
      compressed_rrmat* kkt_rrmat = new compressed_rrmat;
  
      get_rr_inv_pair(kkt, inv, kkt_rrmat, this->constraints+this->variables);
  
      free(kkt);
      
      // Create the right hand side of the KKT system.
  
      double* constant_kkt_rhs = new double[this->constraints+this->variables];
      memset(constant_kkt_rhs, 0, sizeof(double) * (this->constraints + this->variables) );
  
      for (int i = 0; i < this->variables; i++) {
        if (quadratic_coefficients[i] != 0) {
          constant_kkt_rhs[i] = linear_coefficients[i];
        } else {
          constant_kkt_rhs[i] = 0;
        }
      }
  
      for (int i = this->variables; i < this->variables+this->constraints; i++) {
         constant_kkt_rhs[i] = rhs[i-this->variables];
      }
  
      // Perform the active set algorithm.
  
      for (int i = 0; i < variables; i++) {
        this_solution[i] = init_solution[i];
      }
  
      active_set_alg(kkt_rrmat, this->constraints+this->variables, inv, this_solution,
        constant_kkt_rhs, this->variables, this->constraints);
  
      free(inv);
  
      return 0;
    };
};

#endif