#ifndef _OBELISK_COMPRESSED_RRMAT_
#define _OBELISK_COMPRESSED_RRMAT_
#include <map>

using namespace std;

class compressed_rrmat {

  /*
    This is compressed storage for a square "row reduced echelon" matrix, plus a new column and row.
    Those scare quotes are there because leading ones are envisioned along the main diagonal for convenience.
    
    It's a map of arrays each representing the nonzero entries of one column, with column number (left to right ascending) as
    the key. With the exception of the new column, array length should be key+2, consisting of the nonzeroes from top to bottom
    followed by a zero on the main diagonal, and then a value for the new row. The new column is key+1 in size because the
    new row is on the diagonal, and the final value of the new column is unused unless it becomes an ordinary column later on.
    (The final value of the new column corresponds to some hypothetical future new row)
  */
    
  public:
    map <int, double*> cols;
    
    map <int, double*>::iterator get_col_exact(int col_number) {
      return this->cols.find(col_number);
    }
    
    map <int, double*>::iterator get_next_col(int col_number) {
      map <int, double*>::iterator c;
      for (c = this->cols.begin(); c->first < col_number && c != this->cols.end(); ++c) {}
      return c;
    }
    
    int insert_col(int col_number, double* col_values) {
      this->cols.insert( make_pair(col_number, col_values));
      return 0;
    }
    
    int delete_col(map <int, double*>::iterator col_pointer) {
      
      if (col_pointer == cols.end()) {
        return 1;
      } else {
        free(col_pointer->second);
        this->cols.erase(col_pointer);
      }
      
      return 0;
    }
    
    int get_elem (int row_num, map <int, double*>::iterator col_pointer) {
      
      if (col_pointer == cols.end()) {
        exit(1);
      } else {
        if (col_pointer == this->cols.end()) { // Is this the new column? If so, it's a straight map from row to element.
          return col_pointer->second[row_num];
        } else {
          if (row_num == this->cols.end()->first) {  // Is this the new row? If so, use the devoted final element of the array.
            return col_pointer->second[col_pointer->first];
          } else {
            return col_pointer->second[row_num];
          }
        }
      }
    }
    
    int set_elem (int row_num, map <int, double*>::iterator col_pointer, double new_val) {
      
      if (col_pointer == cols.end()) {
        return 1;
      } else {
        if (col_pointer == this->cols.end()) { // Is this the new column? If so, it's a straight map from row to element.
           col_pointer->second[row_num] = new_val;
        } else {
          if (row_num == this->cols.end()->first) {  // Is this the new row? If so, use the devoted final element of the array.
            col_pointer->second[col_pointer->first] = new_val;
          } else {
            col_pointer->second[row_num] = new_val;
          }
        }
      }
      return 0;
    }
    
    int rr_dscal(double factor, map <int, double*>::iterator start_col, int row_num) {
      
      if (start_col == this->cols.end()) {return 1;}
      
      for (map <int, double*>::iterator c = start_col; c != this->cols.end(); ++c) {
        set_elem(row_num, c, factor*get_elem(row_num,c) );
      }
      return 0;
    }
    
    int rr_daxpy (double factor, map <int, double*>::iterator start_col, int x_row_num, int y_row_num) {
    
      if (start_col == this->cols.end()) {return 1;}
      
      for (map <int, double*>::iterator c = start_col; c != this->cols.end(); ++c) {
        set_elem(y_row_num, c, factor*get_elem(x_row_num, c) + get_elem(y_row_num, c) );
      }
      return 0;
    }
    
    int rr_dswap (map <int, double*>::iterator start_col, int x_row_num, int y_row_num) {
      
      if (start_col == this->cols.end()) {return 1;}
      
      for (map <int, double*>::iterator c = start_col; c != this->cols.end(); ++c) {
        
        double copy_value = get_elem(y_row_num, c);
      
        set_elem(y_row_num, c, get_elem(x_row_num, c) );
        set_elem(x_row_num, c, copy_value);
      }
      
      return 0;
    }
};

#endif