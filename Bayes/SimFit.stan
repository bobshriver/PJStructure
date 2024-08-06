
data {
  int tmax;
  int B[tmax];

}


parameters {
  real <lower=0, upper=1> s[tmax];
 real <lower=0, upper=100> b[tmax];
 real <lower=1> Nstart;
}

transformed parameters{
  real N[tmax+1];
  real svec[tmax];
  real svectemp[tmax];
  real Bpred[tmax];
  N[1]=Nstart;
  svectemp=reverse(exp(cumulative_sum(reverse(log(s)))));
  svec[1:tmax]=svectemp[1:tmax];
  for(t in 2:tmax+1){
  N[t]=b[t-1]*N[t-1]+s[t-1]*N[t-1]; ///Note that for simplicity in stan coding indexing for b and s here differs from equations in paper b[t-1]=b_t in paper
  Bpred[t-1]=b[t-1]*N[t-1]*svec[t-1];

  }


}

model {
  Nstart ~normal(50,10);
   s ~ normal(0.95,.1)T[0,1];
  B ~ poisson(Bpred);
  
}

