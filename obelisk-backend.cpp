#include <pqxx/pqxx>

#include <map>

#include <cfloat>

#include <tgmath.h>

#include <assert.h>

#include <coin/ClpSimplex.hpp>

#include <x86_64-linux-gnu/cblas.h>

#include "obelisk-backend.hpp"

#include "obelisk-compressed-rrmat.hpp"

using namespace std;
using namespace pqxx;

int num_recipes = 0;
int variables = 0;
int num_resources = 0;
int constraints = 0;
int num_utils = 0;
double * naturals;
double * max_prices;
double * utility_slopes;
double * activities;
pair < int, double > ** tech_matrix;
int * constraint_cols;

int swap_pointers(double** x, double** y) {
  double* z = *x;
  *x = *y;
  *y = z;
  return 0;
}

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
}

int print_nz(int size, double* array) {
  cout << "\n";
  for (int j = 0; j < size; j++) if (array[j] != 0) cout << "\t" << j << "," << array[j];
  return 0;
}

int column_collapse(compressed_rrmat* lhs, map <int, double*>::iterator lhs_col, double** inv, int inv_cols) {

  // We assume the diagonal element is already scaled to one, so there's no need to handle a possible zero value,
  // and a valid iterator is given as input. 
  
  int num_col = lhs_col->first-1;
  
  for (int i = num_col; i >= 0; i--) {
    int fac = -lhs->get_elem(i, lhs_col);
    
    if (fac != 0) {
      lhs->rr_daxpy(fac, lhs_col, num_col, i);
      cblas_daxpy(inv_cols, fac, inv[num_col], 1, inv[i], 1);
    }
  }
  
  lhs->delete_col(lhs_col);
  
  return 0;
}

int active_set_alg(compressed_rrmat* lhs, int num_rows, double** inv, double* solution, double* kkt_rhs) {

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
}

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
    
    lhs->rr_daxpy(-1, current_col, new_active, *num_rows);
    cblas_daxpy(inv_cols, -1, inv[new_active], 1, inv[*num_rows-1], 1);
    
  } else { // The new active constraint corresponds to a singular column. Swap the 1 into it's diagonal value.
  
    lhs->set_elem(*num_rows, current_col, 1);
    
    lhs->rr_dswap(current_col, new_active, *num_rows);
    swap_pointers(&inv[new_active], &inv[*num_rows-1]);
  }
  
  for (current_col; current_col != lhs->cols.end(); ++current_col) {
  
    double leading_value = lhs->get_elem(new_active, current_col);
    
    if (leading_value != 0) {  //  If this evaluates as false, the row remains (or becomes, in new column) a zero row.
    
      lhs->rr_dscal(1/leading_value, current_col, *num_rows);
      cblas_dscal(inv_cols, 1/leading_value, inv[current_col->first-1], 1);
      
      // Remember that the row wouldn't be stored if the diagonal was nonzero, so this leaves a zero in the new row.
      lhs->rr_dswap(current_col, current_col->first, *num_rows); 
      swap_pointers(&inv[new_active], &inv[*num_rows-1]);
      column_collapse(lhs, current_col, inv, inv_cols);
    }
  }
  
  return 0;
}

int main_algorithm(double* this_solution) {

  // Create the transposed technology matrix.
  
  int* variable_rows = new int[variables];
  memset(variable_rows, 0, sizeof(int) * variables);
  int* variable_rows_filled = new int[variables];
  memset(variable_rows_filled, 0, sizeof(int) * variables);
  pair < int, double > ** t_tech_matrix = new pair < int, double > * [variables];
  
  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < constraint_cols[i]; j++) {
      variable_rows[tech_matrix[i][j].first]++;
    }
  }
  
  for (int i = 0; i < variables; i++) {
    t_tech_matrix[i] = new pair < int, double > [variable_rows[i]];
  }
  
  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < constraint_cols[i]; j++) {
    
      int add_col = tech_matrix[i][j].first;
      
      t_tech_matrix[add_col][variable_rows_filled[add_col]] = make_pair(i,tech_matrix[i][j].second);
      variable_rows_filled[add_col]++;
    }
  }
  
  free(variable_rows_filled);
  
  // Generate a linear programming problem to find an interior point for the quadratic
  // problem solve. It will have a single additional variable t and constraints requiring
  // all other problem and slack variables to be greater than t. Maximize t.
  
  /* If uncommented, print the tech matrix for debugging instead.
  
  double** print_matrix = new double*[constraints];
  for (int i = 0; i < constraints; i++) {
    print_matrix[i] = new double[variables];
    for (int j = 0; j < variables; j++) {
      print_matrix[i][j] = 0;
    }
  }
  
  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < constraint_cols[i]; j++) {
      print_matrix[i][tech_matrix[i][j].first] = tech_matrix[i][j].second;
    }
  }
  
  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < variables; j++) {
      cout << "\t" << print_matrix[i][j] << "\t";
    }
    cout << "\n";
  }
  
  return 0;
  
  */
  
  /* If uncommented, print the transposed tech matrix for debugging instead.
  
  double** print_matrix = new double*[variables];
  for (int i = 0; i < variables; i++) {
    print_matrix[i] = new double[constraints];
    for (int j = 0; j < constraints; j++) {
      print_matrix[i][j] = 0;
    }
  }
  
  for (int i = 0; i < variables; i++) {
    for (int j = 0; j < variable_rows[i]; j++) {
      print_matrix[i][t_tech_matrix[i][j].first] = t_tech_matrix[i][j].second;
    }
  }
  
  for (int i = 0; i < variables; i++) {
    for (int j = 0; j < constraints; j++) {
      cout << "\t" << print_matrix[i][j] << "\t";
    }
    cout << "\n";
  }
  
  return 0;
  
  */
  
  double* init_solution_rhs_l = new double[constraints+variables];
  memset(init_solution_rhs_l, 0, sizeof(double) * (constraints + variables) );
  double* init_solution_rhs_u = new double[constraints+variables];
  memset(init_solution_rhs_u, 0, sizeof(double) * (constraints + variables) );
  double* init_solution_objs = new double[variables+1];
  memset(init_solution_objs, 0, sizeof(double) * (variables + 1) );
  
  CoinPackedMatrix init_solution_mat;
  init_solution_mat.setDimensions(constraints+variables, variables+1);
  
  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < constraint_cols[i]; j++) {
     init_solution_mat.modifyCoefficient(i,tech_matrix[i][j].first,tech_matrix[i][j].second);
    }
  }
    
  for (int i = 0; i < variables; i++) {
    init_solution_mat.modifyCoefficient(constraints+i,i,1);
    init_solution_mat.modifyCoefficient(constraints+i,variables,-1);
  }
  
  for (int i = 0; i < constraints; i++) {
    init_solution_rhs_l[i] = -naturals[i];
    init_solution_rhs_u[i] = -naturals[i];
  }
  
  for (int i = constraints; i < constraints+variables; i++) {
    init_solution_rhs_l[i] = 0;
    init_solution_rhs_u[i] = INFINITY;
  }
  
  init_solution_objs[variables] = 1;
  
  ClpSimplex init_solution_model;
  
  init_solution_model.setOptimizationDirection(-1);
  init_solution_model.loadProblem(init_solution_mat,NULL,
    NULL,init_solution_objs,init_solution_rhs_l,init_solution_rhs_u);
  
  init_solution_model.primal();
  double* init_solution = init_solution_model.primalColumnSolution();
  
  if (init_solution_model.isProvenPrimalInfeasible()) {
  
    double* this_infeasibles = init_solution_model.infeasibilityRay(0);
    
    double* zero_rhs = new double[variables+constraints];
    memset(zero_rhs, 0, sizeof(double) * (constraints + variables) );
    
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
    
    double* feasible_solution_rhs_l = new double[constraints];
    memset(feasible_solution_rhs_l, 0, sizeof(double) * (constraints) );
    double* feasible_solution_rhs_u = new double[constraints];
    memset(feasible_solution_rhs_u, 0, sizeof(double) * (constraints) );
    double* feasible_solution_objs = new double[variables];
    
    for (int i = 0; i < variables; i++) {
      feasible_solution_objs[i] = 1;
    }
  
    CoinPackedMatrix feasible_solution_mat;
    feasible_solution_mat.setDimensions(constraints, variables);
  
    for (int i = 0; i < constraints; i++) {
      for (int j = 0; j < constraint_cols[i]; j++) {
       feasible_solution_mat.modifyCoefficient(i,tech_matrix[i][j].first,tech_matrix[i][j].second);
      }
    }
  
    for (int i = 0; i < constraints; i++) {
      feasible_solution_rhs_l[i] = -naturals[i];
      feasible_solution_rhs_u[i] = -naturals[i];
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
    
  // Create the right hand side of the KKT system.
  
  double* constant_kkt_rhs = new double[constraints+variables];
  memset(constant_kkt_rhs, 0, sizeof(double) * (constraints + variables) );
  
  for (int i = 0; i < variables; i++) {
    if (i < num_recipes && utility_slopes[i] != 0) {
      constant_kkt_rhs[i] = max_prices[i];
    } else {
      constant_kkt_rhs[i] = 0;
    }
  }
  
  for (int i = variables; i < variables+constraints; i++) {
    constant_kkt_rhs[i] = -naturals[i-variables];
  }
  
  // for (int i = 0; i < variables+constraints; i++) cout << "  " << constant_kkt_rhs[i];

  // Create the initial matrix pair, namely the KKT matrix and a diagonal matrix called inv.

  double** inv = new double*[constraints+2*variables];  // Allocating space enough for every possible future pointer.
  double** kkt = new double*[constraints+variables];
  
  for (int i = 0; i < constraints+variables; i++) {
    kkt[i] = new double[constraints+variables];
    inv[i] = new double[constraints+variables];
    memset(inv[i], 0, sizeof(double) * (constraints+variables) );
    memset(kkt[i], 0, sizeof(double) * (constraints+variables) );
    inv[i][i] = 1.0;
  }
  
  // Fill initial KKT matrix.
  
  for (int i = 0; i < variables; i++) {
    if (i < num_recipes) {
      if (utility_slopes[i] != 0) {
        kkt[i][i] = -utility_slopes[i];
      }
    }
  }

  for (int i = 0; i < constraints; i++) {
    for (int j = 0; j < constraint_cols[i]; j++) {
      kkt[tech_matrix[i][j].first][i + variables] = tech_matrix[i][j].second;
      kkt[i + variables][tech_matrix[i][j].first] = tech_matrix[i][j].second;
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
  
  get_rr_inv_pair(kkt, inv, kkt_rrmat, constraints+variables);
  
  free(kkt);
  
  // Perform the active set algorithm.
  
  for (int i = 0; i < variables; i++) {
    this_solution[i] = init_solution[i];
  }
  
  active_set_alg(kkt_rrmat, constraints+variables, inv, this_solution, constant_kkt_rhs);
  
  free(inv);
  
  return 0;
}

int main(int argc, char * argv[]) {

  /*  A sketch of the technology matrix, increasing left to right and top to bottom:
    
       p p p p p p p p p r r
       p p p p p p p p p r r
       e e e e e e e e e 0 0
       e e e e e e e e e 0 0
       t t t t t t t t t r r
       t t t t t t t t t r r
    
      p: Productions directly from the database, with negative values from transport-requiring resources moved to 't'.
      e: "Productions" calculated to limit the amount of chemical elements invested in the entire system.
      t: Productions corresponding to resources after being transported.
      r: Productions added to represent transporting a resource.
    
      Including slacks, which we need to add:
    
       p p p p p p p p p r r 0 0
       p p p p p p p p p r r 0 0
       e e e e e e e e e 0 0 1 0
       e e e e e e e e e 0 0 0 1
       t t t t t t t t t r r 0 0
       t t t t t t t t t r r 0 0
       
  */


  double * elements[83];
  map < string, int > recipe_uuids;
  map < string, int > resource_uuids;
  map < int, int > transports;
  int counter = 0;

  try {
    connection C("dbname = " + (string) argv[1] + " user = " + getenv("USER") +
      " password = " + (string) argv[2] + " hostaddr = 127.0.0.1 port = 5432");
    if (C.is_open()) {
      cout << "Opened database successfully: " << C.dbname() << endl;
    } else {
      cout << "Can't open database" << endl;
      return 1;
    }

    nontransaction N(C);

    result R_rc_count(N.exec("SELECT count(*) FROM obelisk.resources;"));
    result R_rp_count(N.exec("SELECT count(*) FROM obelisk.recipes;"));
    result R_pr_count(N.exec("SELECT count(*) FROM obelisk.production;"));
    result R_rs_rows(N.exec("SELECT * FROM obelisk.resources ORDER BY id ASC;"));
    result R_rp_rows(N.exec("SELECT * FROM obelisk.recipes ORDER BY id ASC;"));
    result R_prod_rows(N.exec("SELECT * FROM obelisk.production ORDER BY resource_id ASC;"));
    result R_sprod_rows(N.exec("SELECT * FROM obelisk.starting_production ORDER BY recipe_id ASC;"));

    for (result::const_iterator c = R_rc_count.begin(); c != R_rc_count.end(); ++c) {
      num_resources = c[0].as < int > ();
    }

    constraints = num_resources;

    for (result::const_iterator c = R_rp_count.begin(); c != R_rp_count.end(); ++c) {
      num_recipes += c[0].as < int > ();
    }
    variables = num_recipes;

    constraint_cols = new int[num_resources * 2 + 83]; // Biggest possible length to avoid yet another counting step.
    int * t_filled = new int[num_resources * 2 + 83];
    memset(constraint_cols, 0, sizeof(int) * (num_resources * 2 + 83) );
    memset(t_filled, 0, sizeof(int) * (num_resources * 2 + 83) );
    max_prices = new double[num_recipes];
    utility_slopes = new double[num_recipes];
    memset(max_prices, 0, sizeof(double) * num_recipes );
    memset(utility_slopes, 0, sizeof(double) * num_recipes );
    activities = new double[num_recipes+num_resources];
    naturals = new double[num_resources * 2 + 83];
    memset(naturals, 0, sizeof(double)  * (num_resources * 2 + 83) );

    for (int i = 0; i < 83; i++) {
      elements[i] = new double[num_resources];
    }

    counter = 0;
    for (result::const_iterator c = R_rs_rows.begin(); c != R_rs_rows.end(); ++c) {
      resource_uuids.insert(pair < string, int > (c[0].as < string > (), counter));
      cout << "id = " << counter << endl;
      cout << "resource_name = " << c[1].as < string > () << endl;
      naturals[counter] = c[2].as < double > ();

      // Parse the elements array and save it in the corresponding elements[][] column.

      array_parser elem_array = c[6].as_array();
      for (int i = 0; i < 83; i++) {
        pair < array_parser::juncture, string > array_pair = elem_array.get_next();
        while (array_pair.first != array_parser::juncture::string_value) {
          array_pair = elem_array.get_next();
        }
        elements[i][counter] = stod(array_pair.second);
      }
      counter++;
    }
    counter = 0;

    double * delays = new double[num_recipes];

    counter = 0;
    for (result::const_iterator c = R_rp_rows.begin(); c != R_rp_rows.end(); ++c) {
      cout << "id = " << counter << endl;
      recipe_uuids.insert(pair < string, int > (c[0].as < string > (), counter));
      cout << "recipe_name = " << c[1].as < string > () << endl;
      if (c[2].as < double > () != 0) {
        num_utils++;
        utility_slopes[counter] = c[3].as < double > () / (c[2].as < double > () - c[5].as < double > ());
        max_prices[counter] =  c[5].as < double > ();
        cout << "Price at supply 0: " << max_prices[counter] << endl;
        cout << "Demand curve slope: " << utility_slopes[counter] << endl;
      }
      delays[counter] = c[4].as < double > ();
      counter++;
    }
    counter = 0;

    // Count transport stuff (a recipe, resource and three productions per transport-requiring resource).
    // This also adds the tranport redirects so sparsity structure can be sussed out.

    counter = 0;
    for (result::const_iterator c = R_rs_rows.begin(); c != R_rs_rows.end(); ++c) {
      if (!c[4].is_null()) {
        constraint_cols[resource_uuids[c[0].as < string > ()]]++;
        constraint_cols[resource_uuids[c[4].as < string > ()]]++;
        constraint_cols[num_resources + 83 + counter]++;

        transports.insert(make_pair(resource_uuids[c[0].as < string > ()], num_resources + 83 + counter));

        variables++;
        constraints++;
        counter++;
      }
    }
    counter = 0;

    constraints += 83;
    variables += 83;

    cout << constraints << " constraints." << endl; // Now we know the full size of the constraints matrix/jacobian.
    cout << variables << " variables." << endl;

    // Go through productions and count the nonzero entries in each constraint.

    for (result::const_iterator c = R_prod_rows.begin(); c != R_prod_rows.end(); ++c) {
      if (transports.find(resource_uuids[c[2].as < string > ()]) == transports.end() || c[3].as < double > () > 0) {
        constraint_cols[resource_uuids[c[2].as < string > ()]]++;
      } else {
        constraint_cols[transports.find(resource_uuids[c[2].as < string > ()]) -> second]++;
      }
    }

    /* Generate "productions" for each recipe corresponding to an element 
       constraint using the starting prods. Use the original resource instead of the transported resource
       because that's where the data is (tokens are intangible and shouldn't have elements). 
       Also, put into a temporary dense matrix first for memory allocation reasons.*/

    double * elem_matrix_dense[83];
    memset(elem_matrix_dense, 0, sizeof(int) * 83);

    for (int elem = 0; elem < 83; elem++) {

      counter = 0;
      double recipe_element = 0;
      elem_matrix_dense[elem] = new double[variables];
      memset(elem_matrix_dense[elem], 0, sizeof(double) * variables);

      for (result::const_iterator c = R_sprod_rows.begin(); c != R_sprod_rows.end(); ++c) {
        if (counter != recipe_uuids[c[1].as < string > ()]) {
          if (recipe_element != 0) {
            elem_matrix_dense[elem][counter] = recipe_element;
            constraint_cols[num_resources + elem]++;
          }
          counter = recipe_uuids[c[1].as < string > ()];
          recipe_element = 0;
        }
        recipe_element += c[3].as < double > () * delays[counter] * elements[elem][resource_uuids[c[2].as < string > ()]];
      }
      
      // Add the slack variable.
      
      constraint_cols[num_resources + elem]++;
      elem_matrix_dense[elem][variables - 83 + elem] = -1;
      
      counter = 0;
    }

    // Allocate technology matrix.

    tech_matrix = new pair < int, double > * [constraints];
    for (int i = 0; i < constraints; i++) {
      tech_matrix[i] = new pair < int, double > [constraint_cols[i]];
    }

    // Add transport stuff (a recipe, resource and three productions per transport-requiring resource).

    counter = 0;
    for (result::const_iterator c = R_rs_rows.begin(); c != R_rs_rows.end(); ++c) {
      if (!c[4].is_null()) {
        tech_matrix[resource_uuids[c[0].as < string > ()]][t_filled[resource_uuids[c[0].as < string > ()]]] =
          make_pair(num_recipes + counter, -1);
        t_filled[resource_uuids[c[0].as < string > ()]]++;
        tech_matrix[resource_uuids[c[4].as < string > ()]][t_filled[resource_uuids[c[4].as < string > ()]]] =
          make_pair(num_recipes + counter, -c[5].as < double > ());
        t_filled[resource_uuids[c[4].as < string > ()]]++;

        tech_matrix[num_resources + 83 + counter][t_filled[num_resources + 83 + counter]] =
          make_pair(num_recipes + counter, 1);
        t_filled[num_resources + 83 + counter]++;

        counter++;
      }
    }
    counter = 0;

    // Add explicit productions, with substitutions made for transport resources.

    for (result::const_iterator c = R_prod_rows.begin(); c != R_prod_rows.end(); ++c) {
      if (transports.find(resource_uuids[c[2].as < string > ()]) == transports.end() || c[3].as < double > () > 0) {
        // There is no transport or this production is positive.
        int row_index = resource_uuids[c[2].as < string > ()];
        int col_index = t_filled[resource_uuids[c[2].as < string > ()]];
        tech_matrix[row_index][col_index] =
          make_pair(recipe_uuids[c[1].as < string > ()], c[3].as < double > ());
        t_filled[resource_uuids[c[2].as < string > ()]]++;
      } else {
        // Otherwise, we put the entries in a post-transport row.
        int row_index = transports[resource_uuids[c[2].as < string > ()]];
        int col_index = t_filled[transports[resource_uuids[c[2].as < string > ()]]];
        tech_matrix[row_index][col_index] =
          make_pair(recipe_uuids[c[1].as < string > ()], c[3].as < double > ());
        t_filled[transports[resource_uuids[c[2].as < string > ()]]]++;
      }
    }

    for (int i=num_resources; i<num_resources+83; i++) {
      naturals[i] = 2e12; // Placeholder.
    }

    // Transfer the nonzero entries in the elements constraints to the main technology matrix.

    for (int elem = 0; elem < 83; elem++) {
      tech_matrix[num_resources + elem] = new pair < int, double > [num_recipes + elem];
      counter = 0;
      for (int j = 0; j < variables; j++) {
        if (elem_matrix_dense[elem][j] != 0) {
          tech_matrix[num_resources + elem][counter] =
            make_pair(j, elem_matrix_dense[elem][j]);
          counter++;
        }
      }
      counter = 0;
    }
   
   double* solution = new double[variables];
   if (main_algorithm(solution)) return 1;
   
   for (result::const_iterator c = R_rp_rows.begin(); c != R_rp_rows.end(); ++c) {
      N.exec("UPDATE obelisk.recipes SET final_activity = " + to_string(solution[recipe_uuids[c[0].as <string>()]]) + 
        " WHERE id = '" + c[0].as <string>() + "';"
      );
   }
      
  } catch (const std::exception & e) {
    cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}
    
    
    