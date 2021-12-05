data {
    int N;
    vector[N] Y;
    int new_X;
}
parameters {
    vector[N] a;
    real<lower=0> s_a;
    real <lower=0> s_Y;
}
model {
    a[3:N] ~ normal(2*a[2:(N-1)]-a[1:(N-2)], s_a);
    Y ~ normal(a, s_Y);
}
generated quantities {
    real new_Y[N+new_X];
    real new_a[N+new_X];
    for (i in 1:N){
        new_Y[i] = normal_rng(a[i], s_Y);
        new_a[i] = a[i];
    }
    for (i in N+1:N+new_X){
        new_a[i] = normal_rng(2*new_a[i-1]-new_a[i-2], s_a);
        new_Y[i] = normal_rng(new_a[i], s_Y);
    }
}
