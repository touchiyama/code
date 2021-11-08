data {
    int Na;
    int Np;
    int Nt;
    vector[Nt] T;
    int<lower=1, upper=Np> PersonID[Na];
    int<lower=1, upper=Nt> TimeID[Na];
    vector[Na] Y;
    int new_T;
    real new_Time[new_T];
}
parameters {
    real m_a;
    real m_b;
    vector[Np] log_a;
    vector[Np] log_b;
    real<lower=0> s_a;
    real<lower=0> s_b;
    real <lower=0> s_Y;
}
transformed parameters {
    vector[Np] a;
    vector[Np] b;
    matrix[Np, Nt] mu;
    a = exp(log_a);
    b = exp(log_b);
    for (j in 1:Nt)
        for (i in 1:Np)
            mu[i, j] = a[i]*(1-exp(-b[i]*T[j]));
}
model {
    log_a ~ normal(m_a, s_a);
    log_b ~ normal(m_b, s_b);
    for (i in 1:Na){
        Y[i] ~ normal(mu[PersonID[i], TimeID[i]], s_Y);
    }
}
generated quantities {
    real new_Y[Np, new_T];
    for (i in 1:Np)
        for (j in 1:new_T)
            new_Y[i, j] = normal_rng(a[i]*(1-exp(-b[i]*new_Time[j])), s_Y);
}
