/* file Fishing_Flooding_Model.c */
#include <R.h>
#include <math.h>
static double parms[37];

#define N0_F parms[0]
#define J0_F parms[1]
#define A0_F parms[2]
#define N0_R parms[3]
#define J0_R parms[4]
#define A0_R parms[5]
#define b_M parms[6]
#define comp_d parms[7]
#define comp_b parms[8]
#define comp_m parms[9]
#define cann parms[10]
#define c_N parms[11]
#define c_J parms[12]
#define d_A parms[13]
#define m_J parms[14]
#define d_J parms[15]
#define m_N parms[16]
#define d_N parms[17]
#define naup_catch parms[18]
#define f parms[19]
#define f_J parms[20]
#define f_N parms[21]
#define h parms[22]
#define h_J parms[23]
#define h_N parms[24]
#define i_P parms[25]
#define k_R parms[26]
#define k_F parms[27]
#define latent_stages parms[28]
#define latent_rate parms[29]
#define lambda parms[30]
#define d_W parms[31]
#define d_F parms[32]
#define conversationEff parms[33]
#define ImmigrationRate parms[34]
#define FishingRate parms[35]
#define ImmigrationPeriod parms[36] 

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=37;
odeparms(&N, parms);
}


double sum_array(const double *arr, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += arr[i];
    return s;
} 

/* ^ this is used to compute double sumEs = sum_array(Es, latent_stages_local); */ 

void compute_derivatives(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
    int latent_stages_local = (int)latent_stages;  
    double VOL = 1.0;

    double N = y[0];
    double J = y[1];
    double A = y[2];
    const double *Es = &y[3];
    double I = y[3 + latent_stages_local];
    double Preds = y[4 + latent_stages_local];
    double L3F = y[5 + latent_stages_local];

    // Compute sumEs
    double sumEs = sum_array(Es, latent_stages_local);

    // Compute crowding term
    double crowd = c_N*N + c_J*J + A + sumEs;

    // Density-dependent rates
    double d_A_c = d_A * exp(comp_d / VOL * crowd);
    double d_J_c = d_J * exp(comp_d / VOL * crowd);
    double d_N_c = d_N * exp(comp_d / VOL * crowd);

    double m_N_c = m_N * exp(-comp_m / VOL * crowd);
    double m_J_c = m_J * exp(-comp_m / VOL * crowd);

    // Immigration and fishing
		// || means "or"
		// && means "and"
		// tmod is a floating point double (larger range of numbers)
		// fmod is the floating point remainder from division 

    double tmod = fmod(*t, 365.0); 
	
/*calculates the floating-point remainder of t divided by 365.0 */ 
    
    double immigration = (tmod > 100 && tmod < 100 + ImmigrationPeriod) ? ImmigrationRate : 0; 
	
/*current day of the year (tmod) is greater than 100 and less than 100 plus the duration of the immigration period, and assigns to immigration variable if the condition is true.*/ 
    
    double fishing = (tmod <= 100 || tmod >= 100 + ImmigrationPeriod) ? FishingRate : 0;
	
/*Is the current day of the year less than or equal to 100? OR Is the current day of the year greater than or equal to the day the immigration period ends (100 + ImmigrationPeriod)? FishingRate: This value is assigned to the fishing variable if the condition is true. 
0: This value is assigned to the fishing variable if the condition is false.*/

    // Predation terms
    double Pred_A = f*Preds/(1 + f*h*(A + sumEs + I + f_J*h_J*J + f_N*h_N*N)/VOL + i_P*fmax(Preds-1,0)/VOL);
    double Pred_J = f*f_J*Preds/(1 + f*h*(A + sumEs + I + f_J*h_J*J + f_N*h_N*N)/VOL + i_P*fmax(Preds-1,0)/VOL);
    double Pred_N = f*f_N*Preds/(1 + f*h*(A + sumEs + I + f_J*h_J*J + f_N*h_N*N)/VOL + i_P*fmax(Preds-1,0)/VOL);

    // Latent progression
		// double * tells you where the memory address for the double (latent progression) is stored 
		// malloc is "memory allocation." Requests a block of memory from the system 
		// latent_stages_local * sizeof(double) calculates the total number of memory bytes needed for the array
		// i++ moves onto the next exposed step with each iteration 

    double *latent_progression = (double *)malloc(latent_stages_local * sizeof(double));

    for (int i = 0; i < latent_stages_local; i++)
        latent_progression[i] = latent_rate * Es[i];

/* ^ Iterates from index 0 up to latent_stages_local - 1, calculating a new value for each index in the latent_progression array based on a constant rate and the corresponding value from the Es array.*/ 

    // Compute derivatives
    ydot[0] = b_M*(A + sumEs)/2*exp(-comp_b/VOL * crowd) - (m_N_c + d_N_c)*N - cann*(A + I + sumEs)*N - Pred_N*N;
    ydot[1] = m_N_c*N - (m_J_c + d_J_c)*J - Pred_J*J;
    ydot[2] = m_J_c*J - d_A_c*A - lambda*A - Pred_A*A;

    for (int i = 0; i < latent_stages_local; i++) {
        double gain = (i == 0) ? lambda*A : latent_progression[i-1];     
/* ?: means if else. In the first loop iteration (i = 0),gain is computed as lambda multiplied by A. For later iterations, gain comes from the previous element in the array latent_progression */ 

        ydot[3+i] = -latent_progression[i] - d_A_c*Es[i] + gain - Pred_A*Es[i];
    }

    ydot[3 + latent_stages_local] = latent_progression[latent_stages_local-1] - d_A_c*I - Pred_A*I;
    ydot[4 + latent_stages_local] = Preds*(conversationEff*(Pred_N*N + Pred_J*J + Pred_A*A + Pred_A*sumEs)) - d_F*Preds + immigration - fishing*Preds;
    ydot[5 + latent_stages_local] = Pred_A*I - d_W*L3F;

    free(latent_progression);
} 

/* END file Fishing_Flooding_Model.c */ 