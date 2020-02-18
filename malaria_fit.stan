functions {
    real[] malaria_ode(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
        /*########## (DECLARE) ##########*/
        /*########## (declare) x_r and x_i ##########*/
        int max_iRBC;
        int n_odeParams;
        real mu;
        real R0;

        /*########## (declare) estimated  parameters ##########*/
        real phi_N1;
        real phi_N3;
        real psi_N1;
        real psi_N3;
        real n_N1;
        real n_N3;
        
        /*########## (declare) variables ##########*/
        real N1;
        real N3;
        real Ret1;
        real Ret2;
        real Ret3;
        real Nor;
        real I;

        /*########## (declare) dynamics ##########*/
        real dydt[7];

        /*########## (declare) derived  parameters and variables ##########*/
        real RBCkilling;
        real pRBCkilling;
        real I_total;
        real R_dev;

        /*########## (DEFINE) ##########*/
        /*########## (define) fixed  parameters ##########*/
        max_iRBC  = x_i[1];
        n_odeParams  = x_i[2];

        mu = x_r[1];
        R0= x_r[5];

        /*########## (define) estimated  parameters ##########*/
        phi_N1 =  params[1];
        phi_N3 =  params[2];
        psi_N1 =  params[3];
        psi_N3 =  params[4];
        n_N1 =  params[5];
        n_N3 =  params[6];
        
        /*########## (define) variables ##########*/
        N1 =   y[1];
        N3 =   y[2];
        Ret1 = y[3];
        Ret2 = y[4];
        Ret3 = y[5];
        Nor =  y[6];
        I =    y[7];

        /*########## (define) dynamics ##########*/
        // Receptor activation
        R_dev = (R0*exp(-(t-0.5)*mu)) - (Ret1+Ret2+Ret3+Nor);
        if(R_dev<0) R_dev =0;
        
        dydt[1] = psi_N1*(R_dev/(R0*exp(-(t-0.5)*mu)))^n_N1   - 1/phi_N1 * N1; // dN1
        dydt[2] = psi_N3*(I/max_iRBC)^n_N3                    - 1/phi_N3 * N3; // dN3
        
        RBCkilling =  N1;
        pRBCkilling = N3;

        dydt[3] = -mu*Ret1  - RBCkilling*Ret1             ;   // dRet1_dt
        dydt[4] = -mu*Ret2  - RBCkilling*Ret2             ;   // dRet2_dt
        dydt[5] = -mu*Ret3  - RBCkilling*Ret3             ;   // dRet3_dt
        dydt[6] = -mu*Nor   - RBCkilling*Nor              ;   // dNor_dt
        dydt[7] = -mu*I     -(RBCkilling + pRBCkilling)*I ;   // dI_dt

        return dydt;
    }

    vector malaria_fit(vector phi, vector theta, real[] x_r, int[] x_i) {
        /*########## (DECLARE) ##########*/
        /*########## (declare) fixed  parameters ##########*/
        real mu;
        real muz;
        real R0;
        real I0;
        real t0;
        real ds[x_i[4]];

        /*########## (declare) fitted  parameters ##########*/
        real params[x_i[2]]; //x_i[2] == n_odeParams; // Model parameters.
        real rho;
        real g_1;
        real g_2;
        
        real sd_RBC;
        real sd_iRBC;

        real lp_total;

        /*########## (define) derived  parameters ##########*/
        real R_dev;
        real RBC_midn;
        real lambda_ret;
        real lambda_norm;

        /*########## (declare) data storage  ##########*/
        // One for evaluation 12pm starting
        real y_eval[x_i[4]+1,x_i[3]+4]; //x_i[4] == n_ds, x_i[3] == n_difeq plus 1 is for M

        // One for evaluation 12am starting (12h hours delayed)
        real y_midn[x_i[4]+1,x_i[3]+4]; //x_i[4] == n_ds, x_i[3] == n_difeq plus 1 is for M

        // One for storing ODE output: 1st time point for y_eval and 2nd time point for y_midn
        real y_hat[x_i[5],x_i[3]]; //x_i[3]] == n_difeq

        /*########## (DEFINE) ##########*/
        /*########## (define) fixed  parameters ##########*/
        mu = x_r[1];
        muz = x_r[2];
        //p_1 = x_r[3];
        //p_2 = x_r[4];
        R0  = x_r[5];
        I0 = x_r[6];
        t0 = x_r[7];
        
        /*########## (define) estimated  parameters ##########*/
        sd_RBC = phi[1];
        sd_iRBC = phi[2];

        rho =       theta[1];
        g_1 =       theta[2];
        g_2 =       theta[3];
        params[1] = theta[4];
        params[2] = theta[5];
        params[3] = theta[6];
        params[4] = theta[7];
        params[5] = theta[8];
        params[6] = theta[9];
        
        /*########## (define) initial conditions ##########*/
        y_eval[1,1] = 0;                  //N1
        y_eval[1,2] = 0;                  //N3
        y_eval[1,3] = R0*mu*exp(-0.5*mu); //Ret1
        y_eval[1,4] = R0*mu*exp(-1.5*mu); //Ret2
        y_eval[1,5] = R0*mu*exp(-2.5*mu); //Ret3
        y_eval[1,6] = R0 - (R0*mu*exp(-0.5*mu)+R0*mu*exp(-1.5*mu)+R0*mu*exp(-2.5*mu));     //Nor
        y_eval[1,7] = I0;      //I
        y_eval[1,8] = 0;     //M
        y_eval[1,9] = R0*mu;    //Pre1
        y_eval[1,10] = R0*mu;    //Pre2
        y_eval[1,11] = R0*mu;    //Pre3

        //Wind the clock forward by 0.5 day
        y_hat[,] = integrate_ode_bdf(malaria_ode, y_eval[1,1:x_i[3]], x_r[7], x_r[(7+0*x_i[5]+1):(7+1*x_i[5])], params, x_r, x_i,1e-5,1e-4,1000);

        y_midn[1,1:x_i[3]] = y_hat[1,];
        
        // Precusors remain the same
        y_midn[1,9] = y_eval[1,9];
        y_midn[1,10] = y_eval[1,10];
        y_midn[1,11] = y_eval[1,11];
        
        /*########## (MODEL) ##########*/
        for (d in 1:x_i[4]){
            /*########## Model phase 1: growth  ##########*/
            // for now y_midn[d,8] is holding y_midn[after growth,7]
            y_midn[d,8] = sum(y_midn[d,3:5]) * (1-exp(-g_1 * y_midn[d,7])) + y_midn[d,6] * (1-exp(-g_2 * y_midn[d,7]));
        
            /*########## Model phase 2: iRBC bursting (instantaneous)  ##########*/
            y_midn[d,3] = y_midn[d,3] * exp(-g_1 * y_midn[d,7]);
            y_midn[d,4] = y_midn[d,4] * exp(-g_1 * y_midn[d,7]);
            y_midn[d,5] = y_midn[d,5] * exp(-g_1 * y_midn[d,7]);
            y_midn[d,6] = y_midn[d,6] * exp(-g_2 * y_midn[d,7]);
            y_midn[d,7] = y_midn[d,8];

            /*########## Model phase 3: Reticulocyte development (instantaneous) ##########*/
            // The last stage reticulocytes become normocytes
            y_midn[d,6] = y_midn[d,6] + y_midn[d,5];

            // Reticulocytes get one day older
            y_midn[d,5] = y_midn[d,4];
            y_midn[d,4] = y_midn[d,3];
            
            // The last stage precusors become reticulocytes
            y_midn[d,3] = y_midn[d,11];

            // Precusors get one day older
            y_midn[d,11] = y_midn[d,10];
            y_midn[d,10] = y_midn[d,9];

            /*########## Model phase 4: Precusors input (instantaneous) ##########*/
            RBC_midn = sum(y_midn[d,3:7]);
            R_dev = (R0*exp(-0.5*mu)) - RBC_midn;
            if(R_dev<0) {R_dev =0;}
            y_midn[d,9] = mu * R0 + rho*R_dev;
            
            /*########## Model phase 5: RBC and iRBC clearance (continous time) ##########*/
            y_hat[,] = integrate_ode_bdf(malaria_ode, y_midn[d,1:x_i[3]], x_r[7], x_r[(7+0*x_i[5]+1):(7+1*x_i[5])], params, x_r, x_i,1e-5,1e-4,1000);

            y_eval[d+1,1:x_i[3]] = y_hat[1,]; //1st evalution is used for likelihood calculation
            y_midn[d+1,1:x_i[3]] = y_hat[2,]; //2nd evaluation is used for the model loop

            // Precusors remain the same
            y_midn[d+1,9] = y_midn[d,9];
            y_midn[d+1,10] = y_midn[d,10];
            y_midn[d+1,11] = y_midn[d,11];
        }

        /*########## (LIKELIHOOD) ##########*/
        lp_total = 0;
        for (d in 1:x_i[4]){
            // x_i[4] = n_ds
            if (x_i[(5+1*x_i[4]+1):(5+2*x_i[4])][d]!=0){
                // Note: 0 is used as a magic number replacing NA
                lp_total += normal_lpdf(x_i[(5+1*x_i[4]+1):(5+2*x_i[4])][d] | sum(y_eval[d+1,3:7]),sd_RBC);
            }

            if (x_i[(5+2*x_i[4]+1):(5+3*x_i[4])][d]!=0){
                // Note: 0 is used as a magic number replacing NA
                lp_total += normal_lpdf(log10(x_i[(5+2*x_i[4]+1):(5+3*x_i[4])][d]+1) | log10(y_eval[d+1,7]+1),sd_iRBC);
            }
        }
	
        return [lp_total]';

    }
}

data {
    int<lower = 1> n_odeParams; // Number of internal ODE parameters
    int<lower = 1> n_RndEffs; // Number of mouse random effects (correlated)
    int<lower = 1> n_difeq; // Number of differential equations in the system

    int<lower = 1> n_obs; // Number of data points (total)
    int<lower = 1> n_mice; // Number of mice (total)

    int<lower = 1> n_strains; // Number of dose treatments fitted
    int<lower = 1> strain_id[n_mice];

    int<lower = 1> n_treatments; // Number of dose treatments fitted
    int<lower = 1> treatment_id[n_mice];

    real R0[n_mice];
    real I0[n_mice];
    
    int day[n_obs];
    int y_RBC[n_obs];
    int y_iRBC[n_obs];

    real t0; // Initial time point (zero)
    int<lower = 1> n_ts; // Number of time steps per day
    real ts[n_ts]; // Time points (in proportion) that were sampled per day
    int<lower = 1> n_ds; // Number of time steps

}

transformed data {
    // real data
    real x_r[n_mice,7+n_ts];
    real p_1dd;
    real p_2dd;
    real mu;
    real muz;

    // int data
    int x_i[n_mice,5+(3*n_ds)];
    int components_I;
    int R_hat;
    int max_iRBC;

    /*########## fixed parameters ##########*/
    mu = 0.025; // mortality rate of RBC (0.025/day cited in Miller et al. 2010 for mice)
    muz = 48.0; // intrinsic mortality rate of free merozoites
    p_1dd = 1e-06;
    p_2dd = 10e-06;
    max_iRBC = 1350000; //From Pollitt et al.

    for(j in 1:n_mice) {
        x_i[j, 1] = max_iRBC;
        x_i[j, 2] = n_odeParams;
        x_i[j, 3] = n_difeq;
        x_i[j, 4] = n_ds;
        x_i[j, 5] = n_ts;
        x_i[j,(5+0*n_ds+1):(5+1*n_ds)] = day[((n_ds*(j-1))+1):((n_ds*(j-1))+n_ds)];
        x_i[j,(5+1*n_ds+1):(5+2*n_ds)] = y_RBC[((n_ds*(j-1))+1):((n_ds*(j-1))+n_ds)];
        x_i[j,(5+2*n_ds+1):(5+3*n_ds)] = y_iRBC[((n_ds*(j-1))+1):((n_ds*(j-1))+n_ds)];

        x_r[j,1] = mu;
        x_r[j,2] = muz;
        x_r[j,3] = p_1dd;
        x_r[j,4] = p_2dd;
        x_r[j,5] = R0[j];
        x_r[j,6] = I0[j];
        x_r[j,7] = t0;
        x_r[j,(7+0*n_ts+1):(7+1*n_ts)] = ts;
    }
}

parameters {
    real rho;
    real g_1;
    real g_2;
    real phi_N1;
    real phi_N3;
    real psi_N1;
    real psi_N3;
    real n_N1;
    real n_N3;

    vector<lower = 0>[n_RndEffs] sd_s;
    vector<lower = 0>[n_RndEffs] sd_u;

    cholesky_factor_corr[n_RndEffs] L_s;
    cholesky_factor_corr[n_RndEffs] L_u;

    matrix[n_RndEffs,n_strains] z_s;
    matrix[n_RndEffs,n_mice] z_u;

    real<lower = 0> sd_RBC;
    real<lower = 0> sd_iRBC;
}

transformed parameters{
    vector[2] phi; // phi (fixed parameters between subjects: measurement errors for RBC and iRBC; needed to be specified as a vector for map_rect)
    matrix[n_RndEffs,n_strains] s; //strain level random effects
    matrix[n_RndEffs,n_mice] u; //subject level random effects
    vector[n_RndEffs] theta[n_mice]; // boundary indexing at theta[n_mice][n_RndEffs]

    phi[1] = sd_RBC;
    phi[2] = sd_iRBC;

    s = diag_pre_multiply(sd_s, L_s) * z_s;
    u = diag_pre_multiply(sd_u, L_u) * z_u;

    for(j in 1:n_mice) {
        theta[j][1] = exp(rho+           s[1,strain_id[j]] +u[1,j])*0.4;
        theta[j][2] = 1/(10^(exp(g_1+    s[2,strain_id[j]] +u[2,j])*6));
        theta[j][3] = 1/(10^(exp(g_2+    s[3,strain_id[j]] +u[3,j])*6));
        theta[j][4] = exp(phi_N1+        s[4,strain_id[j]] +u[4,j]);
        theta[j][5] = exp(phi_N3+        s[5,strain_id[j]] +u[5,j]);
        theta[j][6] = exp(psi_N1+        s[6,strain_id[j]] +u[6,j]);
        theta[j][7] = exp(psi_N3+        s[7,strain_id[j]] +u[7,j]);
        theta[j][8] = exp(n_N1+          s[8,strain_id[j]] +u[8,j]);
        theta[j][9] = exp(n_N3+          s[9,strain_id[j]] +u[9,j]);
    }
}

model {
    rho ~ normal(0,0.1); // Miller et al.
    g_1 ~ normal(0,0.25); //
    g_2 ~ normal(0,0.25); //
    phi_N1 ~ normal(log(1)+2.5,sqrt(2.5));
    phi_N3 ~ normal(log(1)+2.5,sqrt(2.5));
    psi_N1 ~ normal(log(1)+2.5,sqrt(2.5));
    psi_N3 ~ normal(log(1)+2.5,sqrt(2.5));
    n_N1 ~ normal(0,0.25);
    n_N3 ~ normal(0,0.25);

    sd_s ~ normal(0,2.5);
    sd_u ~ normal(0,2.5);

    to_vector(z_s) ~ normal(0,1);
    to_vector(z_u) ~ normal(0,1);

    L_s ~ lkj_corr_cholesky(5);
    L_u ~ lkj_corr_cholesky(5);

    sd_RBC ~ normal(5*10^5,(5*10^5)/10); // linear-scale sd of 5*10^5 reported by Miller et al. 2010
    sd_iRBC ~ normal(0.13,0.13/10); // log10 sd of 0.2 reported by Mideo et al. 2008

    target += map_rect(malaria_fit, phi, theta, x_r, x_i);
}

generated quantities {

}


