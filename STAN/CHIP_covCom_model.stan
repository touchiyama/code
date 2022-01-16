data {
    int Ni;
    int Nt;
    int Y[Nt, Ni];
}
parameters {
    matrix[Nt, Ni] r;
    real<lower=0> s_r;
}
model {
    for (i in 1:Nt)
        for (j in 3:Ni)
            target += normal_lpdf(
                r[i, j] | 2*r[i, j-1]-r[i, j-2],
                s_r
            );

    for (i in 1:Nt)
        for (j in 1:Ni)
            Y[i, j] ~ poisson_log(r[i, j]);
}
generated quantities {
    matrix[Nt, Ni] Y_mean;
    for (i in 1:Nt)
        for (j in 1:Ni)
            Y_mean[i, j] = exp(r[i, j]);
}
