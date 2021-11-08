data {
    int N;
    vector[N] X;
    vector[N] T;
    int new_X;
    real X_new[new_X];
}
parameters {
    real m_a;
    real m_b;
    real m_c;
    real <lower=120, upper=190> a;
    real <lower=0, upper=120> b;
    real <lower=0, upper=1> c;
    real<lower=0> s_a;
    real<lower=0> s_b;
    real<lower=0> s_c;
    real <lower=0> s_T;
}
transformed parameters {
    vector[N] mu;
    mu = a-b*exp(-c*X);
}
model {
    a ~ normal(m_a, s_a);
    b ~ normal(m_b, s_b);
    c ~ normal(m_c, s_c);
    to_vector(T) ~ normal(to_vector(mu), s_T);
}
generated quantities {
    real T_new[new_X];
    for (i in 1:new_X){
        T_new[i] = normal_rng(a-b*exp(-c*X_new[i]), s_T);
    }
}
