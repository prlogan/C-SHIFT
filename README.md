# C-SHIFT
Algorithm for recovery of noiseless covariance matrix.
Usage: CShift_fn(Cov_obs, use_trace = T)

The input consists of an observed covariance matrix (Cov_obs) calculated from data with additive noise as described in <insert citation to our paper>.  The default (use_trace = T) is to use the Trace as the objective function.  However, the Frobenius Norm can be used by changing to use_trace = F.  Using the trace is much faster and thus recommended.
  
The output is a covariance matrix that is an estimate of covariance of the underlying data without the additive noise.
