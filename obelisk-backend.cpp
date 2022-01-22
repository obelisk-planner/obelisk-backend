#include <pqxx/pqxx>

#include <map>

#include <cfloat>

#include <tgmath.h>

#include <assert.h>

#include <x86_64-linux-gnu/cblas.h>

#include "obelisk-compressed-rrmat.hpp"

#include "obelisk-quadratic-solve.hpp"

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
    max_prices = new double[num_recipes*2+83];
    utility_slopes = new double[num_recipes*2+83];
    memset(max_prices, 0, sizeof(double) * (num_resources * 2 + 83) );
    memset(utility_slopes, 0, sizeof(double) * (num_resources * 2 + 83) );
    activities = new double[num_recipes+num_resources];
    naturals = new double[num_resources * 2 + 83];
    memset(naturals, 0, sizeof(double)  * (num_resources * 2 + 83) );

    for (int i = 0; i < 83; i++) {
      elements[i] = new double[num_resources];
    }

    counter = 0;
    for (result::const_iterator c = R_rs_rows.begin(); c != R_rs_rows.end(); ++c) {
      resource_uuids.insert(pair < string, int > (c[0].as < string > (), counter));
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
      recipe_uuids.insert(pair < string, int > (c[0].as < string > (), counter));
      if (c[2].as < double > () != 0) {
        num_utils++;
        utility_slopes[counter] = c[3].as < double > () / (c[2].as < double > () - c[5].as < double > ());
        max_prices[counter] =  c[5].as < double > ();
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
   
   // Make the quadratic_programming_problem object.
   
   quadratic_programming_problem plan_economy;
   
   double* problem_rhs = new double[constraints];
   for (int i = 0; i < constraints; i++) {
      problem_rhs[i] = -naturals[i];
   }
   
   double* quadratics = new double[variables];
   for (int i = 0; i < variables; i++) {
      quadratics[i] = 0.5*utility_slopes[i];
   }
   
   plan_economy.problem_matrix = tech_matrix;
   plan_economy.constraint_cols = constraint_cols;
   plan_economy.rhs = problem_rhs;
   plan_economy.quadratic_coefficients = quadratics;
   plan_economy.linear_coefficients = max_prices;
   plan_economy.variables = variables;
   plan_economy.constraints = constraints;
   
   double* solution = new double[variables];
   if (plan_economy.solve(solution)) return 1;
   
   cout << "Plan found succesfully. Entering into the database." << endl;
   
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
