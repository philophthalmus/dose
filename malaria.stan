functions {
    real[] malaria_ode(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
        /*########## (DECLARE) ##########*/
        /*########## (declare) x_r and x_i ##########*/
        int R_hat;
        int max_iRBC;
        int n_difeq;

        real p;
        real mu;
        real muz;
        real alpha_s;

        /*########## (declare) estimated  parameters ##########*/
        real rho;
        real beta;
        real psi_N1;
        real psi_N2;
        real phi_N1;
        real phi_N2;

        /*########## (declare) variables ##########*/
        real N1;
        real N2;

        real R;

        real I_1;
        real I_2;
        real I_3;
        real I_4;
        real I_5;
        real I_6;
        real I_7;
        real I_8;
        real I_9;
        real I_10;
        real I_11;
        real I_12;

        real M;

        /*########## (declare) dynamics ##########*/
        real dydt[16];

        /*########## (declare) derived  parameters and variables ##########*/
        real R_dev;

        real erythropoiesis;
        real RBCkilling;
        real pRBCkilling;
        real pRBCkilling_total;

        real I_total;

        /*########## (DEFINE) ##########*/
        /*########## (define) fixed  parameters ##########*/
        R_hat = x_i[1];
        max_iRBC  = x_i[2];
	    n_difeq  = x_i[5];

        p = x_r[1];
        mu = x_r[3];
        muz = x_r[4];

        alpha_s  = x_r[7];

        /*########## (define) estimated  parameters ##########*/
        rho = params[1];
        beta = params[2];
        psi_N1 = params[3];
        psi_N2 = params[4];
        phi_N1 = params[5];
        phi_N2 = params[6];

        /*########## (define) variables ##########*/
        N1 = y[1];
        N2 = y[2];

        R = y[3];

        I_1 = y[4];
        I_2 = y[5];
        I_3 = y[6];
        I_4 = y[7];
        I_5 = y[8];
        I_6 = y[9];
        I_7 = y[10];
        I_8 = y[11];
        I_9 = y[12];
        I_10 = y[13];
        I_11 = y[14];
        I_12 = y[15];

        M = y[16];

        /*########## (define) dynamics ##########*/
        I_total = I_1+I_2+I_3+I_4+I_5+I_6+I_7+I_8+I_9+I_10+I_11+I_12;
        if(I_total >= max_iRBC) I_total = max_iRBC;
        if(I_total <= 0) I_total = 0;

        dydt[1] =   psi_N1*(I_total/max_iRBC) * (1-N1) - 1/phi_N1 * N1;
        dydt[2] =   psi_N2*(I_total/max_iRBC) * (1-N2) - 1/phi_N2 * N2;
        
        R_dev = R_hat - (R+I_total);
        if(R_dev<0) {R_dev =0;}
        
        erythropoiesis = mu * R_hat + rho*R_dev;
        RBCkilling = mu + -log(1-N1);
        pRBCkilling = -log(1-N2);

        dydt[3] = erythropoiesis - p*R*M - RBCkilling*R;

        pRBCkilling_total = RBCkilling + pRBCkilling;

        dydt[4] = - pRBCkilling_total*I_1 - 1/alpha_s * I_1 + p*R*M;
        dydt[5] = - pRBCkilling_total*I_2 + 1/alpha_s * (I_1 - I_2);
        dydt[6] = - pRBCkilling_total*I_3 + 1/alpha_s * (I_2 - I_3);
        dydt[7] = - pRBCkilling_total*I_4 + 1/alpha_s * (I_3 - I_4);
        dydt[8] = - pRBCkilling_total*I_5 + 1/alpha_s * (I_4 - I_5);
        dydt[9] = - pRBCkilling_total*I_6 + 1/alpha_s * (I_5 - I_6);
        dydt[10] = - pRBCkilling_total*I_7 + 1/alpha_s * (I_6 - I_7);
        dydt[11] = - pRBCkilling_total*I_8 + 1/alpha_s * (I_7 - I_8);
        dydt[12] = - pRBCkilling_total*I_9 + 1/alpha_s * (I_8 - I_9);
        dydt[13] = - pRBCkilling_total*I_10 + 1/alpha_s * (I_9 - I_10);
        dydt[14] = - pRBCkilling_total*I_11 + 1/alpha_s * (I_10 - I_11);
        dydt[15] = - pRBCkilling_total*I_12 + 1/alpha_s * (I_11 - I_12);

        dydt[16] = beta * 1/alpha_s * I_12 -p*R*M - muz*M;

        return dydt;
    }

    vector malaria_fit(vector phi, vector theta, real[] x_r, int[] x_i) {
        /*########## (DECLARE) ##########*/
        /*########## (declare) fixed  parameters ##########*/
        real sP;
        real R0;

    	real y_hat[x_i[6], x_i[5]];
        real y0[x_i[5]];
        real I_0_cden[x_i[3] + 1];
        real params[x_i[4]];
        real lp_total;

        real sd_RBC;
        real sd_iRBC;

        /*########## (DEFINE) ##########*/
        /*########## (define) fixed  parameters ##########*/
        sP = x_r[2];
        R0  = x_r[5];

        /*########## (define) estimated  parameters ##########*/
        sd_RBC = phi[1];
        sd_iRBC = phi[2];

        params[1] = theta[1];
        params[2] = theta[2];
        params[3] = theta[3];
        params[4] = theta[4];
        params[5] = theta[5];
        params[6] = theta[6];

        /*########## (define) initial conditions ##########*/
        y0[1] = 0;
        y0[2] = 0;

        y0[3] = R0;

        for (iter in 1:size(I_0_cden)){
            I_0_cden[iter] = x_r[6] * beta_cdf((1.0/12) * (iter-1), sP, sP);
        }

        y0[4] = I_0_cden[2] - I_0_cden[1];
        y0[5] = I_0_cden[3] - I_0_cden[2];
        y0[6] = I_0_cden[4] - I_0_cden[3];
        y0[7] = I_0_cden[5] - I_0_cden[4];
        y0[8] = I_0_cden[6] - I_0_cden[5];
        y0[9] = I_0_cden[7] - I_0_cden[6];
        y0[10] = I_0_cden[8] - I_0_cden[7];
        y0[11] = I_0_cden[9] - I_0_cden[8];
        y0[12] = I_0_cden[10] - I_0_cden[9];
        y0[13] = I_0_cden[11] - I_0_cden[10];
        y0[14] = I_0_cden[12] - I_0_cden[11];
        y0[15] = I_0_cden[13] - I_0_cden[12];

        y0[16] = 0;
        
        y_hat[,] = integrate_ode_bdf(malaria_ode, y0, x_r[8], x_r[(8+0*x_i[6]+1):(8+1*x_i[6])], params, x_r, x_i,1e-5,1e-4,1000);

        lp_total = 0;
        for (d in 1:x_i[6]){
            if (x_i[(6+0*x_i[6]+1):(6+1*x_i[6])][d]!=0){
                lp_total += normal_lpdf(x_i[(6+0*x_i[6]+1):(6+1*x_i[6])][d] | sum(y_hat[d,3:15]),sd_RBC);
            }

            if (x_i[(6+1*x_i[6]+1):(6+2*x_i[6])][d]!=0){
                lp_total += normal_lpdf(log10(x_i[(6+1*x_i[6]+1):(6+2*x_i[6])][d]+1) | log10(sum(y_hat[d,4:15])),sd_iRBC);
            }
        }
	
        return [lp_total]';

    }
}

data {
    int<lower = 1> n_odeParams;
    int<lower = 1> n_RndEffs;
    int<lower = 1> n_difeq;

    int<lower = 1> n_obs;

    int<lower = 1> n_mice;

    int<lower = 1> mouse_id[n_mice];

    int<lower = 1> n_doses;
    int dose_slope_id[n_doses];
    int<lower = 1> dose_id[n_mice];

    real R0[n_mice];
    real I0[n_mice];

    int y_RBC[n_obs];
    int y_iRBC[n_obs];

    real t0;
    int<lower = 1> n_ts;
    real ts[n_ts];

}

transformed data {
    real x_r[n_mice,8+n_ts];
    real alpha;
    real p;
    real sP;
    real mu;
    real muz;
    real alpha_s;

    int x_i[n_mice,6+(2*n_ts)];
    int components_I;
    int R_hat;
    int max_iRBC;

    /*########## fixed parameters ##########*/
    p = 8e-06;
    sP = 10.0;
    mu = 0.025;
    muz = 48.0;

    alpha = 1.0;
    components_I = 12;
    alpha_s = alpha/components_I;

    R_hat = 8891286;
    max_iRBC = 2177040;

    for(j in 1:n_mice) {
        x_i[j,1] = R_hat;
        x_i[j,2] = max_iRBC;
        x_i[j,3] = components_I;
        x_i[j,4] = n_odeParams;
        x_i[j,5] = n_difeq;
        x_i[j,6] = n_ts;
        x_i[j,(6+0*n_ts+1):(6+1*n_ts)] = y_RBC[((n_ts*(j-1))+1):((n_ts*(j-1))+n_ts)];
        x_i[j,(6+1*n_ts+1):(6+2*n_ts)] = y_iRBC[((n_ts*(j-1))+1):((n_ts*(j-1))+n_ts)];

        x_r[j,1] = p;
        x_r[j,2] = sP;
        x_r[j,3] = mu;
        x_r[j,4] = muz;
        x_r[j,5] = R0[j];
        x_r[j,6] = I0[j];
        x_r[j,7] = alpha_s;
        x_r[j,8] = t0;
        x_r[j,(8+0*n_ts+1):(8+1*n_ts)] = ts;
    }
}

parameters {
    real rho;
    real beta;
    real psi_N1;
    real psi_N2;
    real phi_N1;
    real phi_N2;

    real dose_slope[n_RndEffs];

    vector<lower = 0>[n_RndEffs*2] sd_u;

    cholesky_factor_corr[n_RndEffs*2] L_u;
    matrix[n_RndEffs*2,n_mice] z_u;

    real<lower = 0> sd_RBC;
    real<lower = 0> sd_iRBC;
}

transformed parameters{
    vector[2] phi;
    matrix[n_RndEffs*2,n_mice] u;
    vector[n_RndEffs] theta[n_mice];

    phi[1] = sd_RBC;
    phi[2] = sd_iRBC;

    u = diag_pre_multiply(sd_u, L_u) * z_u;

    for(j in 1:n_mice) {
        theta[j][1] = exp(rho+     (dose_slope[1]+u[6+1,j])*dose_slope_id[dose_id[j]] + u[1,j])*0.25;
        theta[j][2] = exp(beta+    (dose_slope[2]+u[6+2,j])*dose_slope_id[dose_id[j]] + u[2,j])*7;
        theta[j][3] = exp(psi_N1+  (dose_slope[3]+u[6+3,j])*dose_slope_id[dose_id[j]] + u[3,j]);
        theta[j][4] = exp(psi_N2+  (dose_slope[4]+u[6+4,j])*dose_slope_id[dose_id[j]] + u[4,j]);
        theta[j][5] = exp(phi_N1+  (dose_slope[5]+u[6+5,j])*dose_slope_id[dose_id[j]] + u[5,j]);
        theta[j][6] = exp(phi_N2+  (dose_slope[6]+u[6+6,j])*dose_slope_id[dose_id[j]] + u[6,j]);
    }
}

model {
    rho ~ normal(0,0.25);
    beta ~ normal(0,0.25);
    psi_N1 ~ normal(log(1)+5,sqrt(5));
    psi_N2 ~ normal(log(1)+5,sqrt(5));
    phi_N1 ~ normal(log(1)+5,sqrt(5));
    phi_N2 ~ normal(log(1)+5,sqrt(5));

    dose_slope ~ normal(0,2.5);
    sd_u ~ normal(0,1);

    to_vector(z_u) ~ normal(0,1);
    L_u ~ lkj_corr_cholesky(5);

    sd_RBC ~ normal(5*10^5,(5*10^5)/10);
    sd_iRBC ~ normal(0.13,0.13/10);

    target += map_rect(malaria_fit, phi, theta, x_r, x_i);
}

generated quantities {

}


