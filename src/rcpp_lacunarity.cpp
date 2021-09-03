#include <Rcpp.h>
using namespace Rcpp;

// Helper function to compute the number of iterations of the main
// for-loop in the rcpp_lacunarity function.
int max_iterations (const IntegerVector w_vec,
                    const int mode,
                    const int mat_width,
                    const int mat_height) {
  int N_r, r, w;
  int iters = w_vec.size();
  for (int i = 0; i < w_vec.size(); i++) {
    w = w_vec[i];
    
    if (mode == 1) {
      r = w;
      N_r = (mat_width - r + 1) * (mat_height - r + 1);
    } else {
      r = (w-1)/2;
      N_r = (mat_width-(2*r))*(mat_height-(2*r));
    }
    iters += N_r;
  }
  
  return(iters);
}


// Calculate Lacunarity based on Plotnik ()
double lacunarity (NumericVector box_masses,
                   const int fun,
                   const int N_r) {
  double lac;
  if (box_masses.size() > 1) {
    if (fun == 1) {
      // 1. Max number of box values
      int max_value = max(box_masses);
      
      // 2. Frequency distribution n(S,r)
      IntegerVector n_S_r(max_value+1, 0);
      for (int j = 0; j < box_masses.size(); j++) {
        n_S_r[box_masses[j]] += 1;
      }
      
      // 3. Probability distribution Q(S,r)
      NumericVector Q_S_r(max_value+1, 0.0);
      for (int k = 0; k < Q_S_r.size(); k++) {
        Q_S_r[k] = n_S_r[k] / double(box_masses.size());
      }
      
      // 4. First and second moments of Q(S,r): S*Q(S,r) and S^2*Q(S,r)
      NumericVector first_moment(max_value+1, 0.0);
      NumericVector second_moment(max_value+1,0.0);
      for (int S = 0; S < Q_S_r.size(); S++) {
        first_moment[S] = S*Q_S_r[S];
        second_moment[S]= S*S*Q_S_r[S];
      }
      double Z_1 = sum(first_moment);
      double Z_2 = sum(second_moment);
      
      // 5. Lacunarity
      lac = Z_2/(Z_1*Z_1);
    } else {
      lac = 1 + (( sd(box_masses)  *  sd(box_masses) ) /
        ( mean(box_masses)*mean(box_masses) ));
      
    }
  } else {
    lac = NA_REAL;
  }
  
  return(lac);
}


#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::export]]
NumericVector rcpp_lacunarity (const NumericMatrix mat,
                               const IntegerVector w_vec,
                               const int fun,
                               const int mode,
                               const int ncores=1,
                               const bool display_progress=false) {
  // Public parameters
  const int mat_width = mat.ncol();
  const int mat_height = mat.nrow();
  int N_r, r, w;
  
  
  // Output
  NumericVector output(w_vec.size());
  
  
  // Progress bar
  int iters = max_iterations(w_vec, mode, mat_width, mat_height);
  Progress p(iters, display_progress);
  
  
  // Begin main loop
  for (int j = 0; j < w_vec.size(); j++) {
    p.increment();
    w = w_vec[j];
    
    
    if (mode == 1) {
      r = w;
      N_r = (mat_width - r + 1) * (mat_height - r + 1);
    } else {
      r = (w-1)/2;
      N_r = (mat_width-(2*r))*(mat_height-(2*r));
    }
    
    
    // Gliding box algorithm
    // Init shared vector for parallel loop
    NumericVector box_masses(N_r);
    
#if defined(_OPENMP)
    omp_set_num_threads(ncores);
#pragma omp parallel for shared(box_masses)
#endif
    for (int i = 0; i < N_r; i++) {
      p.increment();
      
      // private parameters
      double cell_value = 0.0;
      double cell_min = R_PosInf;
      double cell_max = R_NegInf;
      double cell_sum = 0.0;
      int n = 0;
      int x, y;
      
      
      // Get x/y from i
      if (mode == 1) {
        x = floor(i / (mat_width - r + 1));
        y = i % (mat_width - r + 1);
      } else {
        x = floor(i / (mat_width - (2*r))) + r;
        y = (i % (mat_width - (2*r))) + r;
      }
      
      
      // Get box mass (window of r*r):
      // If mode == 1, a square box will be used, else a circular disk.
      // If fun == 1, box-sum will be calculated, else box-range
      
      if (mode == 1) {                  // Square box
        for (int bx = 0; bx < r; bx++) {
          for (int by=0; by < r; by++) {
            cell_value = mat(x+bx,y+by);
            if (!NumericVector::is_na(cell_value)) {
              // Update min
              if (cell_value < cell_min)
                cell_min = cell_value;
              
              // Update max
              if (cell_value > cell_max)
                cell_max = cell_value;
              
              // Update sum
              cell_sum += cell_value;
              
              // Update n
              n += 1;
            }
          }
        }
      } else {                        // Circular disk
        for (int px = x-r; px <= x+r; px++)  {
          for (int py = y-r; py <= y+r; py++){
            cell_value = mat(px,py);
            if (!NumericVector::is_na(cell_value)) {
              if ( (px-x)*(px-x) + (py-y)*(py-y) <= r*r ) {
                // Update min
                if (cell_value < cell_min)
                  cell_min = cell_value;
                
                // Update max
                if (cell_value > cell_max)
                  cell_max = cell_value;
                
                // Update sum
                cell_sum += cell_value;
                
                // Update n
                n += 1;
              }
            }
          }
        }
      }
      
      // Box statistic
      if (n > 0) {
        if (fun == 1) {
          box_masses[i] = cell_sum;
        } else {
          box_masses[i] = (cell_max - cell_min);
        }
      } else {
        box_masses[i] = NA_REAL;
      }
    }
    
    // Remove NA and compute Lacunarity based on
    NumericVector bm_narm = wrap(na_omit(box_masses));
    
    // Compute lacunarity
    double lac = lacunarity(bm_narm, fun, N_r);
    
    output[j] = lac;
  }
  
  return(output);
}