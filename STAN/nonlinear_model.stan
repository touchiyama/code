data {
    int Np;
    int Nt;
    real T[Nt];
    real Y[Np, Nt];
    int new_T;
    real new_Time[new_T];
}
parameters {
    real m_a;
    real m_b;
    real m_c;
    real log_a[Np];
    real log_b[Np];
    real log_c[Np];
    real<lower=0> s_a;
    real<lower=0> s_b;
    real<lower=0> s_c;
    real <lower=0> s_Y;
}
transformed parameters {
    real a[Np];
    real b[Np];
    real c[Np];
    for (i in 1:Np) {
        a[i] = exp(log_a[i]);
        b[i] = exp(log_b[i]);
        c[i] = exp(log_c[i]);
    }

}
model {
    for (i in 1:Np) {
        log_a[i] ~ normal(m_a, s_a);
        log_b[i] ~ normal(m_b, s_b);
        log_c[i] ~ normal(m_c, s_c);
    }

    for (i in 1:Np)
        for (j in 1:Nt)
            Y[i, j] ~ normal(a[i]-b[i]*exp(-c[i]*T[j]), s_Y);
}
generated quantities {
    real new_Y[Np, new_T];
    for (i in 1:Np)
        for (j in 1:new_T)
            new_Y[i, j] = normal_rng(a[i]-b[i]*exp(-c[i]*new_Time[j]), s_Y);
}
