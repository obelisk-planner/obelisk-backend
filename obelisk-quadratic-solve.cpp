#include <map>

#include <cfloat>

#include <tgmath.h>

#include <assert.h>

#include <x86_64-linux-gnu/cblas.h>

#include "obelisk-compressed-rrmat.hpp"

#include "obelisk-quadratic-solve.hpp"

using namespace std;

int swap_pointers(double** x, double** y) {
  double* z = *x;
  *x = *y;
  *y = z;
  return 0;
};

int print_nz(int size, double* array) {
  cout << "\n";
  for (int j = 0; j < size; j++) if (array[j] != 0) cout << "\t" << j << "," << array[j];
  return 0;
}

int column_collapse(compressed_rrmat* lhs, map <int, double*>::iterator lhs_col, double** inv, int inv_cols) {

  // This uses row operations to zero out ("collapse") all the non-zero elements of a column except for a pivot value.

  // We assume the diagonal element is already scaled to one, so there's no need to handle a possible zero value,
  // and a valid iterator is given as input. 
  
  int lhs_col_number = lhs_col->first;
  int matrix_dim = (--(lhs->cols.end()))->first;
  double fac;
  
  if (lhs_col_number != matrix_dim) {
  
    // Zero out the new row entry before proceeding to the loop.
  
    fac = -lhs->get_elem(matrix_dim-1, lhs_col);
  
    if (fac != 0) {
      lhs->rr_daxpy(fac, lhs_col, lhs_col_number, matrix_dim-1);
      cblas_daxpy(inv_cols, fac, inv[lhs_col_number], 1, inv[matrix_dim-1], 1);
    }
    
    for (int i = lhs_col_number-1; i >= 0; i--) {
      fac = -lhs->get_elem(i, lhs_col);
    
      if (fac != 0) {
        lhs->rr_daxpy(fac, lhs_col, lhs_col_number-1, i);
        cblas_daxpy(inv_cols, fac, inv[lhs_col_number-1], 1, inv[i], 1);
      }
    }
    
  } else {
  
    for (int i = lhs_col_number-2; i >= 0; i--) {
      fac = -lhs->get_elem(i, lhs_col);
    
      if (fac != 0) {
        lhs->rr_daxpy(fac, lhs_col, lhs_col_number-1, i);
        cblas_daxpy(inv_cols, fac, inv[lhs_col_number-1], 1, inv[i], 1);
      }
    }
    
  }
  
  return 0;
};

int active_set_alg(compressed_rrmat* lhs, int num_rows, double** inv, double* solution, double* kkt_rhs, int variables, int constraints) {

  double* kkt_solution = new double[variables];

  double* search_dir = new double[variables];
  pair < int, double > first_barrier;

  while (1) {
    
    // Get the values of the subproblem solution corresponding to variables of the original quadratic
    // optimization problem.
    
    for (int i = 0; i < variables; i++) {
      kkt_solution[i] = cblas_ddot(variables+constraints,inv[i],1,kkt_rhs,1);
    }
    
    // Look for the first constraint we hit on the way to the subproblem solution.
    // Because all inactive constraints are positivity contraints this takes a simpler form.
    
    cblas_dcopy(variables,kkt_solution,1,search_dir,1);
    cblas_daxpy(variables,-1,solution,1,search_dir,1);
    
    first_barrier = make_pair(-1, -1.0);
    
    for (int i = 0; i < variables; i++) {
      if (search_dir[i] < 0) {
        double barrier_dist = solution[i]/search_dir[i];
        if (barrier_dist > first_barrier.second) first_barrier = make_pair(i, barrier_dist);
      }
    }
    
    // Advance the solution as far as we've determined is possible;
    
    cblas_daxpy(variables,-first_barrier.second,search_dir,1,solution,1);
    
    if (first_barrier.second == -1.0) {
      return 0; // The subproblem solution is feasible. Return.
    } else {

      // Add the newly active constraint to the problem and restart.
      update_inv(lhs, &num_rows, inv, variables+constraints, first_barrier.first); 
    }
    
  }
};

int update_inv(compressed_rrmat* lhs, int* num_rows, double** inv, int inv_cols, int new_active) {
  
  /* We want to add the an "equals zero" constraint to the tech matrix on linear problem variable new_active.
     This means adding a row and a column onto the KKT matrix. We don't want to have to start the solve
     from scratch, though. The solution is to "simulate" the row operations that have already happened
     so we only must do linearly more row operations. Since the additional row is new and wouldn't have been
     touched we just add the relevant standard basis vector. The operations on the column can be simulated
     by multiplying that vector with the inverse. */
     
  (*num_rows)++;
     
  inv[*num_rows-1] = new double[inv_cols];
  
  // Add the new column to the compressed matrix.
  // Since it's a standard basis vecor, the product is just the relevant column of the inverse matrix.

  double* values = new double[*num_rows+1];
  
  for (int i = 0; i < *num_rows; i++) {
    values[i] = inv[i][new_active];
  }
  values[*num_rows] = 0;
  
  lhs->insert_col(*num_rows, values);
  
  // Add (possibly implicitly) the one element of the new row and perform the first step.
  
  map <int, double*>::iterator current_col = lhs->get_col_exact(new_active);
  
  if (current_col == lhs->cols.end() ) {
  
    current_col = lhs->get_next_col(new_active);
    
    lhs->rr_daxpy(-1, current_col, new_active, *num_rows-1);  // The one in the new row is implicitly erased by non-placement.
    cblas_daxpy(inv_cols, -1, inv[new_active], 1, inv[*num_rows-1], 1);
    
  } else { // The new active constraint corresponds to a singular column. Swap the 1 into it's diagonal value.
  
    lhs->set_elem(*num_rows, current_col, 1);
    
    lhs->rr_dswap(current_col, new_active, *num_rows-1);
    swap_pointers(&inv[new_active], &inv[*num_rows-1]);
  }
  
  map <int, double*>::iterator trailing_col = lhs->cols.end();
  
  for (current_col; current_col != lhs->cols.end(); ++current_col) {
  
    if (trailing_col != lhs->cols.end()) lhs->delete_col(trailing_col);
  
    double leading_value = lhs->get_elem(*num_rows-1, current_col);
    
    if (leading_value != 0) {  //  If this evaluates as false, the row remains (or becomes, in new column) a zero row.
    
      lhs->rr_dscal(1/leading_value, current_col, *num_rows-1);
      cblas_dscal(inv_cols, 1/leading_value, inv[current_col->first-1], 1);
      
      if (current_col->first != *num_rows) {
        // Remember that the row wouldn't be stored if the diagonal was nonzero, so this leaves a zero in the new row.
        lhs->rr_dswap(current_col, current_col->first, *num_rows-1);
        swap_pointers(&inv[new_active], &inv[*num_rows-1]);
      }
      
      column_collapse(lhs, current_col, inv, inv_cols);
      
      trailing_col = current_col;
    } else trailing_col = lhs->cols.end();
  }
  
  if (trailing_col != lhs->cols.end()) lhs->delete_col(trailing_col);
  
  // Get rid of any rounding errors in the row of inv corresponding to the new constraint. It will always be zero in reality.
  
  memset(inv[new_active], 0, sizeof(double) * inv_cols);
  
  return 0;
};

int get_rr_inv_pair(double** input_mat, double** inv, compressed_rrmat* lhs, int dim) {

  for (int i = 0; i < dim; i++) {
  
    // Partial pivoting
    
    pair < int, double > largest_abs = make_pair(-1,0);
    
    for (int j = i; j < dim; j++) {
      if (abs(input_mat[j][i]) > abs(largest_abs.second)){
        largest_abs = make_pair(j,input_mat[j][i]);
      }
    }
    
    // cout << "\n" << largest_abs.first << ", " << largest_abs.second;
    
    if (largest_abs.first != -1) {
      swap_pointers(&input_mat[i],&input_mat[largest_abs.first]);
      swap_pointers(&inv[i],&inv[largest_abs.first]);
    } else {
      continue;  // The column corresponds to a singularity in the matrix. Move on.
    }
    
    // Scale and do row operations.
    
    // cout << "Pivot:" << input_mat[i][i] << "\n";
    
    cblas_dscal(dim,1/input_mat[i][i],inv[i],1);
    cblas_dscal(dim-i,1/input_mat[i][i],(double*) input_mat[i]+i,1);
    
    for (int j = 0; j < dim; j++) {
      if (j != i && input_mat[j][i] != 0) {
        cblas_daxpy(dim,-input_mat[j][i],inv[i],1,inv[j],1);
        cblas_daxpy(dim-i,-input_mat[j][i],(double*) input_mat[i]+i,1,(double*) input_mat[j]+i,1);
          //Problem with input_mat[98][183]?
      }
    }
  }
  
  // Compress the lhs matrix.
  
  for (int i = 0; i < dim; i++) {
    if (input_mat[i][i] != 1) {
    
      double* sing_col = new double[i+2];
      
      for (int j = 0; j < i+1; j++) {
        sing_col[j] = input_mat[j][i];
      }
      sing_col[i+1] = 0;
      
      lhs->insert_col(i,sing_col);
    }
  }
  
  return 0;
};