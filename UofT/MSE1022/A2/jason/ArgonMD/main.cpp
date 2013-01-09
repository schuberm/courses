// MD.cpp
// written by Eric Landry
// bare-bones version
// last update May 16, 2006
 
#include <iostream> 
#include <fstream>
#include <math.h>
#include <time.h>
using namespace std;
#include <iomanip>


// declaration of functions/subroutines
double force(double r2); 					// calculates LJ force
double pe(double r2);						// calculates LJ potential energy
double ke(double p[][3], double m[]);		// calculates KE of atom
double partPE(int part, double x[][3]);		// calculates PE of one atom
void chime(int time, double x[][3], char ident_letter[]);		// output for Chime visualization
void atomeye(char name[], int a, int t, double L[3], char ident_letter[], double m[], double x[][3], double p[][3]) ; // output for AtomEye
void gr(int t, double x[][3], double L[3]); // radial distribution function

// declaration of input and output streams
ifstream input("fcc256_40.txt");			// input data
ofstream stats_out;							// energy, temperature, pressure and momentum output (vs. time)
ofstream output;			 				// extra output stream
ofstream chime_out;							// Chime output stream
ofstream atomeye_out;						// Atomeye output stream
ofstream gr_out;							// radial distribution function output

// global variables
const int N = 256;							// number of atoms
double L[3];								// simulation cell size
double pi = atan(1) * 4;					// pi
double f_cutoff, pe_cutoff;					// force and energy at cutoff 
double a = 2.5;								// cutoff
double a2 = a*a;							// cutoff^2		
bool cutoff = true;					     	// if true, use cutoff scheme with discontinuity in force
											// if false, use continuous scheme
		
// LJ values for Argon - used to nondimensionalize values
double epsilon_Ar = 1.67E-21; 				// (J)
double sigma_Ar = 3.4E-10;					// (m)
double mass_Ar = 6.63E-26;					// (kg)
double kB = 1.3806E-23;						// (J/K)

int main() {
	
	int i, j, k, ii;				// counters
	char filename[255];				// dummy text string used to open/close files

	// variable declarations
	double x[N][3];					// position at t and t+delt
	double p[N][3];					// momentum at t and t+delt
	double x_o[N][3];				// initial positions
	double p_half[N][3];			// momentum at t+delt/2
	double m[N];					// mass
	char ident_letter[N];			// particle element letter (C, N, O, etc.)
	double PE, KE, TE;				// potential, kinetic, and total energy
	double r2;						// squared distance between atom pair
	double rij[3];					// component distance between atom pair
	double F[N][3];					// force on atom
	
	double t_step = 0.0002;			// time step (real time)
	int t_total = 500;		 		// simulation time (time steps)
	int t = 0;						// current time (time steps)
	int stats = 10;					// how often energy and momentum statistics are outputted
	int xyz = 1000;					// how often data is outputed for Chime
	int cfg = t_total+1;			// how often data is outputed for Atomeye
	int radial = t_total+1;			// how often the radial distribution function is calculated
		
	// Ensemble control:
	bool thermostat = false;		
	bool barostat = false;			
	bool quench = false;				
	
	
	// Thermostat parameters
	double T;						// instantaneous temperature (ND)
	double Tset_K;					// temperature set (K)
	double Tset;					// temperature set (ND)
	double eta_t = 0.;				// thermostat parameter
	double tau_t = 0.05;			// thermostat time constant
	Tset_K = 40.;
	Tset = Tset_K * kB * (1./epsilon_Ar);
	
	// Barostat parameters
	double P;						// instantaneous pressure (ND)
	double Pset = 0.;				// desired simulation pressure here (ND)
	double P_viral;					// viral pressure term
	double V;						// system volume
	double eps_p = 0.;				// barostat parameter
	double tau_p = 1.0;				// barostat time constant
	
	// Quench parameters
	double eta = 0.;				// quench parameter, remove kinetic energy from the system

	// simulation cell size
	for (k=0; k<3; k++) L[k] = 6.31724;	// desired system size (ND)
	V = L[0]*L[1]*L[2];
			
	// energy and force at cutoff
	pe_cutoff = pe(a2);
	f_cutoff = force(a2);	// includes an extra a in the denominator
				
	// initial positions (read data from input file)
	for (i=0; i<N; i++){
		input >> m[i] >> x_o[i][0] >> x_o[i][1] >> x_o[i][2]>> p[i][0] >> p[i][1] >> p[i][2];
		ident_letter[i] = 'C';
	}
	
	for (ii=0; ii < 1; ii++) { // run multiple simulations
	
		for (i = 0; i<N; i++) {
			for (k = 0; k<3; k++) x[i][k] = x_o[i][k];
		}
	
 	// close and open files
		if (ii != 0) {
			output.close();
			stats_out.close();
			chime_out.close();
		}

		sprintf(filename, "output_%d.txt", ii);
		output.open(filename);
		sprintf(filename, "stats_%d.txt", ii);
		stats_out.open(filename);
		sprintf(filename, "positions_%d.xyz", ii);
		chime_out.open(filename);
		

		// give the atoms a small initial momentum
		double p_sum[3], p_sub[3];
		double p_random;
		srand((unsigned)time(NULL)); 	// initial seed for random function
/*		for (i=0; i<N; i++){
			for (k=0; k<3; k++){
				p_random = rand();
				if (p_random < (RAND_MAX/2)) p_random = p_random*(-1);
				p[i][k] = p_random;
			}
		} */
	
		// scale initial momenta to zero
	
		for (k=0; k<3; k++) p_sum[k] = 0.;			// sum momentum
		for (i=0; i<N; i++){
			for (k=0; k<3; k++) p_sum[k] = p_sum[k] + p[i][k];
		}
		
		for (k=0; k<3; k++) p_sub[k] = p_sum[k]/N;	// find amount to subtract from each atom
	
		for (i=0; i<N; i++){						// subtract momentum
			for (k=0; k<3; k++) p[i][k] = p[i][k] - p_sub[k];
		} 
		
		// scale momentum to get close to desired temperature
		KE = ke(p, m);
		double KE_set, alpha;
		KE_set = 1.5*(N-1.)*Tset;
		alpha = (KE_set)/KE;
		alpha = sqrt(alpha);
		for (i=0; i<N; i++){
			for (k=0; k<3; k++) p[i][k] = p[i][k]*alpha;
		}
		
		for (k=0; k<3; k++) p_sum[k] = 0.;			// sum again and output to screen
		for (i=0; i<N; i++){
			for (k=0; k<3; k++) p_sum[k] = p_sum[k] + p[i][k];
		}
		cout << p_sum[0] << "\t" << p_sum[1] << "\t" << p_sum[2] << endl; 
	

		
		// calculate initial T, PE, KE, and TE
		KE = ke(p, m);
		T = (2.*KE)/(3.*(N-1.));
		PE = 0.;
		for (i=0; i<N; i++){
			for (j=i+1; j<N; j++){
				for (k=0; k<3; k++) {
					rij[k] = x[i][k] - x[j][k];
					if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
					if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
				}
				r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
				if (r2 < a2) {
					if (cutoff == true) PE = PE + pe(r2) - pe_cutoff;
					else {
						PE = PE + pe(r2) - pe_cutoff + (f_cutoff*a*(sqrt(r2) - a));	
					}
				}	
			}	
		}
		TE = KE + PE;
		
		// find initial forces/pressure
		double f;
		for (i=0; i<N; i++){					// clear old values
			for (k=0; k<3; k++) F[i][k] = 0.;
		}
		for (i=0; i<N; i++){
			for (j=i+1; j<N; j++){
				for (k=0; k<3; k++) {
					rij[k] = x[i][k] - x[j][k];
					if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
					if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
				}
				r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
				if (r2 < a2) {
					if (cutoff == true) f = force(r2);
					else f = force(r2) - f_cutoff*a*(1./(sqrt(r2)));
					for (k=0; k<3; k++){
						F[i][k] = F[i][k] + rij[k]*f;
						F[j][k] = F[j][k] - rij[k]*f;
					}
					P_viral = P_viral + r2*f;
				}
			}
		}
		P = (P_viral/(3.*V)) + ((N*T)/V);
		
		// find initial barostat and thermostat parameters
		if (thermostat == true) eta_t = eta_t + t_step*(1./(tau_t*tau_t))*((T/Tset)-1.);
		else eta_t = 0.;
		
		if (barostat == true) eps_p = eps_p + t_step*(1./(tau_p*tau_p))*(P-Pset);
		else eps_p = 0.;
				
		// output "stats" file header and initial data
		stats_out << "Time\tPE/N\tKE/N\tTE/N\tT(K)\tP(MPa)\tpx\tpy\tpz\tL" << endl;
		stats_out << "0\t"<<PE/N<<"\t"<<KE/N<<"\t"<<TE/N<<"\t"<<T*(epsilon_Ar/kB)<<"\t"<<(P*(epsilon_Ar/(sigma_Ar*sigma_Ar*sigma_Ar)))/1000000<<"\t"<<p_sum[0]<<"\t"<<p_sum[1]<<"\t"<<p_sum[2]<<"\t"<<L[0]<<endl;		
		
		/////// simulation start ////////
		for (t=1; t<=t_total; t++) {
			
			// output to screen (to monitor simulation)

			if (t%(stats*10) == 0) cout << t << "\t" << TE << "\t" << T*(epsilon_Ar/kB) << "\t" <<  L[0] << "\t" << sqrt(p_sum[0]*p_sum[0] + p_sum[1]*p_sum[1] + p_sum[2]*p_sum[2]) << endl;
			
			// find momenta at half time step
			for (i=0; i<N; i++){
				for (k=0; k<3; k++) {
					if (quench == true) p_half[i][k] = p[i][k] + (F[i][k] - eta*p[i][k])*0.5*t_step;
					else p_half[i][k] = p[i][k] + (F[i][k] - 1.*(eta_t+eps_p)*p[i][k])*t_step*0.5;
				}
			}
			
			// find new cell size (if barostat is on)
			if (barostat == true) {
				for (k=0; k<3; k++) L[k] = L[k]*pow(1.+(3.*eps_p*t_step),(1./3.));
				V = L[0]*L[1]*L[2];
			}
			
			// find positions at t+delt
			for (i=0; i<N; i++){
				for (k=0; k<3; k++) {
					if (barostat == true) x[i][k] = (1.+ eps_p*t_step + t_step*t_step*(0.5/(tau_p*tau_p))*(P-Pset) + (0.5*t_step*t_step*eps_p*eps_p))*x[i][k] + p_half[i][k]*t_step/m[i] + 0.5*t_step*t_step*eps_p*p[i][k]/m[i];
					else x[i][k] = x[i][k] + p_half[i][k]*t_step*(1/m[i]);
					// periodic boundary conditions
					if (x[i][k] > L[k]) x[i][k] = x[i][k] - L[k];
					if (x[i][k] < 0.) x[i][k] = x[i][k] + L[k];
				}
			}

			// calculate new barostat and thermostat parameters
			if (thermostat == true) eta_t = eta_t + (t_step/(tau_t*tau_t))*((T/Tset)-1.);
			else eta_t = 0.;	
				
			if (barostat == true) eps_p = eps_p + (t_step/(tau_p*tau_p))*(P-Pset);
			else eps_p = 0.;
			

			// find forces/pressure at t+delt
			for (i=0; i<N; i++){					// clear old values
				for (k=0; k<3; k++) F[i][k] = 0.;
			}
			P_viral = 0.;		
			for (i=0; i<N; i++){
				for (j=i+1; j<N; j++){
					for (k=0; k<3; k++) {
						rij[k] = x[i][k] - x[j][k];
						if (rij[k] > (0.5*L[k])) rij[k] = rij[k] - L[k];    // nearest image
						if (rij[k] < (-0.5*L[k])) rij[k] = rij[k] + L[k];
					}
					r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
					if (r2 < a2) {
						if (cutoff == true) f = force(r2);
						else f = force(r2) - f_cutoff*a*(1./(sqrt(r2)));
						for (k=0; k<3; k++){
							F[i][k] = F[i][k] + rij[k]*f;
							F[j][k] = F[j][k] - rij[k]*f;
						}
						P_viral = P_viral + r2*f;
					}
				}
			}
					
			// find momenta at t+delt
			for (i=0; i<N; i++){
				for (k=0; k<3; k++) {
					if (quench == true) p[i][k] = (p_half[i][k] + F[i][k]*0.5*t_step)/(1. + 0.5*eta*t_step);
					else p[i][k] = (p_half[i][k] + F[i][k]*0.5*t_step)/(1. + (eta_t+eps_p)*0.5*t_step);
				}
			}		  
			
			
			// find new temperature
			T = (2.*ke(p,m))/(3.*(N-1.));
			// find new pressure
			P = (P_viral/(3.*V)) + ((N*T)/V);
			
			// scale momentum to desired T (constant kinetic energy ensemble)
			// specify the time range where this is to be applied
	/*		if (t >= 0 && t < 10000) {
				KE = ke(p,m);
				KE_set = 1.5*(N-1)*Tset;	
				alpha = (KE_set)/KE;
				alpha = sqrt(alpha);
				for (i=0; i<N; i++){
					for (k=0; k<3; k++) p[i][k] = p[i][k]*alpha;
				}
			} */
			
			double q = 0.;
			
			for (i = 0; i<N/2; i++) {
				for (j = N/2; j<N; j++) {
					if (i != j) {
						for (k=0; k<3; k++) {
							rij[k] = x[i][k] - x[j][k];
							if (rij[k] > (0.5*L[k])) rij[k] = rij[k] - L[k];    // nearest image
							if (rij[k] < (-0.5*L[k])) rij[k] = rij[k] + L[k];
						}
						r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
						if (r2 < a2) {
							f = force(r2);
							for (k = 0; k<3; k++) q += 0.5*f*rij[k]*((p[j][k]/m[1]) + (p[i][k]/m[0]));
						}
					}
				}
			}

						
			
			

			// output 
			if (t%stats == 0) {			
				// calculate PE, KE and TE
				KE = ke(p,m);
				T = (2.*ke(p,m))/(3.*(N-1.));
				PE = 0.;
				for (i=0; i<N; i++){
					for (j=i+1; j<N; j++){
						for (k=0; k<3; k++) {
							rij[k] = x[i][k] - x[j][k];
							if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
							if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
						}
						r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
						if (r2 < a2) {
							if (cutoff == true) PE = PE + pe(r2) - pe_cutoff;
							else {
								PE = PE + pe(r2) - pe_cutoff + (f_cutoff*a*(sqrt(r2) - a));	
							}
						}
					}	
				}
				TE = KE + PE;
				
				
													
				// calculate total momenta	
				for (k=0; k<3; k++) p_sum[k] = 0.;	
				for (i=0; i<N; i++){
					for (k=0; k<3; k++) p_sum[k] = p_sum[k] + p[i][k];
				}

				//	output data
				
				stats_out<<t<<"\t"<<PE/N<<"\t"<<KE/N<<"\t"<<TE/N<<"\t"<<T*(epsilon_Ar/kB)<<"\t"<<(P*(epsilon_Ar/(sigma_Ar*sigma_Ar*sigma_Ar)))/1000000<<"\t"<<p_sum[0]<<"\t"<<p_sum[1]<<"\t"<<p_sum[2]<<"\t"<<L[0]<<endl;		
			}
			
			// output
			if (t%xyz == 0) chime(t, x, ident_letter);			
			if (t%cfg == 0 && t>0) atomeye(filename, cfg, t, L, ident_letter, m, x, p);
			if (t%radial == 0) gr(t, x, L);
			
		} /////// simulation end ////////
				
	} // multiple simulation loop
	
} // end of main loop
/////////////////////////////////////////////////////////////////////////////////////
double force(double r2){
	// force calculation w/ extra r term in denominator
	double f, r4, r8, inverse_r14;
	
	r4 = r2*r2;
	r8 = r4*r4;
	inverse_r14 = 1./(r8*r4*r2);

	f = 48.*inverse_r14 - 24.*r4*r2*inverse_r14;

	return f;
}
/////////////////////////////////////////////////////////////////////////////////////
double pe(double r2){
	// potential energy calculation
	double inverse_r6, phi;
	
	inverse_r6 = 1./(r2*r2*r2);
	
	phi = 4.*(inverse_r6*inverse_r6 - inverse_r6);

	return phi;
}
////////////////////////////////////////////////////////////////////////////////////
void gr(int t, double x[][3], double L[3]) {
	// calculates the radial distribution function
	int i, j, k;
	double gr[500];
	double rij[3];
	double r2, r;
	char grfilename[255];
	double delr;
	int bin;
	delr = L[0]/500;
	double rho;
	
	for (i=0; i<500; i++) gr[i] = 0.;
	
	rho = N/(L[0]*L[1]*L[2]);
	
	// calculate g(r) for every atom to get a good average
	for (i=0; i<N; i++){
		for (j=0; j<N; j++) {
			if (i!=j){
				for (k=0; k<3; k++) {
					rij[k] = x[i][k] - x[j][k];
					if (rij[k] > (L[k]/2)) rij[k] = rij[k] - L[k];    // nearest image
					if (rij[k] < (-L[k]/2)) rij[k] = rij[k] + L[k];
				}
				r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
				r = sqrt(r2);
				bin = int(floor(r/delr));
				gr[bin] = gr[bin]+1.;
			}
		}
	}

	sprintf(grfilename, "gr%d.txt", t);
	gr_out.open(grfilename);
	
	for (i=0; i<500; i++) {
		gr[i] = gr[i]/N;
		gr[i] = gr[i]*(1/(rho*4*pi*(delr*delr*delr*(i*i + i + 0.25))));
		gr_out << gr[i] << endl;
	}

	gr_out.close ();
}	
/////////////////////////////////////////////////////////////////////////////////////
double partPE(int part, double x[][3]) {
	// calculate the potential energy of one atom
	double part_PHI, r2;
	double rij[3];
	int i, k;
	part_PHI = 0.;
	
	for (i=0; i<N; i++){
		if (i != part) {
			for (k=0; k<3; k++) {
				rij[k] = x[part][k] - x[i][k];
				if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
				if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
			}
			r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
			if (r2 < a2) {
				if (cutoff == true) part_PHI = part_PHI + 0.5*(pe(r2) - pe_cutoff);
				else part_PHI = part_PHI + 0.5*(pe(r2) - pe_cutoff +(f_cutoff*a*(sqrt(r2)-a)));		
			}
		}
	}	
				
	return part_PHI;	
}
/////////////////////////////////////////////////////////////////////////////////////
double ke(double p[][3], double m[]){
	// kinetic energy calculation
	double kinetic = 0.;
	int i;
	
	for (i=0; i<N; i++) {
		kinetic = kinetic + 0.5*(1./m[i])*(p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2]);
	}
	
	return kinetic;
}
/////////////////////////////////////////////////////////////////////////////////////
void chime(int time, double x[][3], char ident_letter[]) {
	// create Chime file

	int i;
	double scale = 2.5; 
	
	chime_out << N << endl<< time << endl;
	
	for(i=0;i<N;i++) {
		chime_out<<ident_letter[i]<<"\t"<<scale*x[i][0]<<"\t"<<scale*x[i][1]<<"\t"<<scale*x[i][2]<<endl;
	}
}
/////////////////////////////////////////////////////////////////////////////////////
void atomeye(char name[], int a, int t, double L[], char ident_letter[], double m[], double x[][3], double p[][3]) { 
	// writes to .cfg file for visualization with AtomEye

	int z;
	int i, j, k;
	z = t/a;
	sprintf (name, "ae_%d.cfg",z);
	atomeye_out.open (name);
	
	// Array to determine which variables are outputted to cfg file (1 = include and 0 = don't include)
	// AuxVariables = {velx, vely, velz, PE, KE, TE}
	int AuxVariables[] = {0, 0, 0, 1, 0, 0};

	int Aux = 0;
	for (i=0;i<6;i++){
		Aux = Aux + AuxVariables[i];
	}
	atomeye_out << "Number of particles = " << N << endl;
	atomeye_out << "A = " << sigma_Ar*1E10 << " Angstrom (basic length-scale)" << endl;
	atomeye_out << "H0(1,1) = " << L[0] << " A" << endl;
	atomeye_out << "H0(1,2) = 0 A" << endl;
	atomeye_out << "H0(1,3) = 0 A" << endl;
	atomeye_out << "H0(2,1) = 0 A" << endl;
	atomeye_out << "H0(2,2) = " << L[1] << " A" << endl;
	atomeye_out << "H0(2,3) = 0 A" << endl;
	atomeye_out << "H0(3,1) = 0 A" << endl;
	atomeye_out << "H0(3,2) = 0 A" << endl;
	atomeye_out << "H0(3,3) = " << L[2] << " A" << endl;
	atomeye_out << ".NO_VELOCITY."<<endl;
	atomeye_out << "entry_count = " << Aux+3 << endl;
	k = 0;
	for (i=0; i<6; i++){
		if (AuxVariables[i]==1) {
			atomeye_out << "auxiliary[" << k << "] = ";
			if (i == 0) atomeye_out << "velx []" << endl;
			if (i == 1) atomeye_out << "vely []" << endl;
			if (i == 2) atomeye_out << "velz []" << endl;
			if (i == 3) atomeye_out << "PE []" << endl;
			if (i == 4) atomeye_out << "KE []" << endl;
			if (i == 5) atomeye_out << "TE []" << endl;
			k=k+1;
		}
	}
	double Kenergy[N], Penergy[N], Tenergy[N];
	for (i=0; i<N; i++) {
		// each atom's KE
		if (AuxVariables[4]==1 || AuxVariables[5]==1){
			Kenergy[i]=0.5*(1/m[i])*((p[i][0]*p[i][0])+(p[i][1]*p[i][1])+(p[i][2]*p[i][2]));
		}		
		// each atom's PE
		if (AuxVariables[3]==1 || AuxVariables[5]==1){
			Penergy[i]= partPE(i, x);
		}			
		// each atom's TE
		if (AuxVariables[5]==1){
			Tenergy[i]=Kenergy[i]+Penergy[i];
		}
	}
	double mass_temp[5];
	for (i=0;i<5;i++){
		mass_temp[i]=0.;
	}
	k=0;
	for (j=0;j<N;j++){
		if (m[j]!=mass_temp[0] && m[j]!=mass_temp[1] && m[j]!=mass_temp[2]&& m[j]!=mass_temp[3]&& m[j]!=mass_temp[4]){
			atomeye_out << m[j] << endl;
			atomeye_out << ident_letter[j] << endl;
			atomeye_out << x[j][0]/L[0]<<" "<<x[j][1]/L[1]<<" "<<x[j][2]/L[2];
			if (AuxVariables[0]==1) atomeye_out << " " << p[j][0]/m[j];
			if (AuxVariables[1]==1) atomeye_out << " " << p[j][1]/m[j];
			if (AuxVariables[2]==1) atomeye_out << " " << p[j][2]/m[j];
			if (AuxVariables[3]==1) atomeye_out << " " << Penergy[j];
			if (AuxVariables[4]==1) atomeye_out << " " << Kenergy[j];
			if (AuxVariables[5]==1) atomeye_out << " " << Tenergy[j];
			atomeye_out << endl;
			mass_temp[k]=m[j];
			k=k+1;
			for (i=j+1; i<N;i++){
				if (m[i]==m[j]){
					atomeye_out << x[i][0]/L[0]<<" "<<x[i][1]/L[1]<<" "<<x[i][2]/L[2];
					if (AuxVariables[0]==1) atomeye_out << " " << p[i][0]/m[i];
					if (AuxVariables[1]==1) atomeye_out << " " << p[i][1]/m[i];
					if (AuxVariables[2]==1) atomeye_out << " " << p[i][2]/m[i];
					if (AuxVariables[3]==1) atomeye_out << " " << Penergy[i];
					if (AuxVariables[4]==1) atomeye_out << " " << Kenergy[i];
					if (AuxVariables[5]==1) atomeye_out << " " << Tenergy[i];
					atomeye_out << endl;
				}
			}	
		}	
	}
	atomeye_out.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////