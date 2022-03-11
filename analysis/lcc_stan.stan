// The input data is an adjacency matrix A
data {
  int<lower=1> K;
  int<lower=1> n;
  matrix[n,n] A;
}

// functions{
//   int get_nr(vector g, int r, int K){
//     int res;
//     
//     res = 0;
// 
//     for(k in 1:n){
//       // res += operator==(g[k], r);
//       res += (g[k] == r);
//     }
// 
//     return res;
//   }
// 
//   real l_g(vector g, int K, real g_const){
//     int res;
//     
//     res = 0;
// 
//     for(r in 1:K){
//       res += lgamma(get_nr(g = g, r = r, K = K));
//     }
// 
//     res += g_const;
// 
//     return res;
//   }
// }

transformed data {
  int<lower = 1> K_omega;
  int<lower = 1> n_A;
  real<lower = 1> g_fact_const;
  
  g_fact_const = lchoose(n-1, K-1)*lgamma(n+1);
  K_omega = K*(K-1)/2;
  n_A = n-1;
}

// The parameters accepted by the model.
parameters {
  // vector[K] g;
  int g[K];
  // vector[K_omega] omega;
  matrix[K,K] omega;
}

// The model to be estimated.
model {
  int gi;
  int gj;
  
  // for(r in 1:K_omega){
  //   omega[r] ~ uniform(0,1);
  // }
  
  for(i in 1:n){
    for(j in 1:n){
      omega[i,j] ~ uniform(0,1);
    }
  }
  
  for(i in 1:n){
    for(j in (i+1):n){
      gi = g[i];
      gj = g[j];
      target += bernoulli_lpmf(A[i+1, j+1] | omega[gi, gj]);
    }
  }
  target += l_g(g = g, K = K, g_const = g_fact_const);
}

