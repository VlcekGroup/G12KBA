/*
Main file for time stepping of diagonal part of non-equilibrium Green's function
Includes initialisation of Green's function as well as initialisation of system parameters eg system size, Hamiltonian couplings etc
*/

//User defined headers
#include "write_to_file.h"
#include "HF.h" //Holds various versions of  Hartree-Fock Hamiltonian
#include "psi.h"//Holds two particle occupation term
#include "print_matrix.h"//Prints matrices (used for testing)
#include "G1_G2_t_step.h"//Holds update steps for the one and two particle Green's function
#include "commutators.h"//Commutators for one and two particle Green's function with respective HF hamiltonian
#include "collision_integral.h"//Collision integral, describes many body effects on one particle Green's function

//#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <complex>
#include <fstream>
#include <chrono>
using namespace std;
const complex<double> im = {0,1};




class G1_G2_time_step
{
	public:
		double t1,t2,t3,P,U,tmax,dt_fixed;
		int Ns,Ns2,Nb,Nb2;

 		vector<double> epsilon;
		G1_G2_time_step(double interaction_strength = 1, double site_hopping = 1,double band_hopping = 0,double band_site_hopping = 0, double pulse_strength = 1, int sites = 2, int bands = 2,double band_gap = 1)
		{
			dt_fixed = .02;
			tmax = 20;//propagation time
		        P = pulse_strength;//strength of pulse that transistions between bands
			U = interaction_strength;//Strength of interaction term in Hubbard hamiltonian
			t1 = site_hopping;//strength of intersite hopping
			t2 = band_hopping;//Strength of interband hopping
			t3 = band_site_hopping;//Strength of hopping from one band/site to different site and band
			Ns = sites;//Number of sites in the chain
			Ns2 = int(pow(Ns,2));
			Nb = bands;//Number of orbitals per site
			Nb2 = int(pow(Nb,2));
			for(int i = 0; i<Nb;++i)
			{epsilon.push_back(0);}
			epsilon[0] = -band_gap/2-2;//Onsite energy for band 0
			epsilon[1] = band_gap/2-U+2;//Onsite energy for band 1
		}


		int init_G(vector<complex<double>> &G1)
		{
			//Function to initialise the one particle greens function. Set to initialise to half filling, alternating between fully filling a site and half filling a site until there are Ns*Nb/2 particles
			int particle_number = 0;
			for(int i = 0; i < int(Ns); ++i)
	    		{
				if(i%2==0)
				{
					for(int a = 0; a < int(Nb); ++a)
					{
						G1[(i+a*Ns)*Nb*Ns + i + a*Ns] = im;
						particle_number+=1;
						if(particle_number == int(Ns*Nb/2))
							return 0;
					}
				}
				else
				{
					for(int a = 0; a < int(Nb/2); ++a)
                                        {
                                                G1[(i+a*Ns)*Nb*Ns + i + a*Ns] = im;
	                                        particle_number+=1;
                                                if(particle_number == int(Ns*Nb/2))
                                                        return 0;
					}
				}
			}
			return 0;
		}

		void RK_steps(vector<double> &n1, vector<complex<double>> &G1, vector<complex<double>> &G2,vector<vector<complex<double>>> &G_diag, double &avg_wall_t_per_loop)
		{
			double w;
			double t = 0;
			double dt;
			vector<complex<double>> h1G1_comm(Ns*Ns*Nb*Nb);
			vector<complex<double>> I(Ns*Ns*Nb*Nb);
			vector<complex<double>> G1_tmp(Ns*Ns*Nb*Nb);
			vector<complex<double>> h2G2_comm(Ns2*Ns2*Nb2*Nb2);
			vector<complex<double>> psi(Ns2*Ns2*Nb2*Nb2);
    			vector<complex<double>> G2_tmp(Ns2*Ns2*Nb2*Nb2);

			while(t<tmax)
			{
				G_diag.push_back(G1);
				double num_dens = 0;

				for(int i = 0; i < int(Ns/2); ++i)
				{
					for(int a = 0; a < Nb; ++a)
					{
						num_dens+= G1[(i+a*Ns)*Nb*Ns + i + a*Ns].imag();
					}
				}
				n1.push_back(num_dens);
				vector<complex<double>> RK_result_G1(Ns*Ns*Nb*Nb);
				vector<complex<double>> RK_result_G2(Ns2*Ns2*Nb2*Nb2);
				vector<complex<double>> K1(Ns*Ns*Nb*Nb);
				vector<complex<double>> K2(Ns2*Ns2*Nb2*Nb2);
//				float temp_t = t;
//				auto wall_t_start = chrono::steady_clock::now();
				for(int k = 0; k<4; ++k)
				{
					vector<complex<double>> h1(Ns*Ns*Nb*Nb);
					if(k == 1 or k == 2)
				    	{
				        	dt = dt_fixed/2;
				        	w = 2;
				    	}

				   	else
				    	{
				        	dt = dt_fixed;
				        	w = 1;
				    	}

				    	for(int i =0; i<Ns*Ns*Nb*Nb;++i)
				    	{
				        	G1_tmp[i] = complex<double>(G1[i]) + K1[i]*dt;
				    	}

				    	for(int i =0; i<Ns2*Ns2*Nb2*Nb2;++i)
				    	{
				        	G2_tmp[i] = complex<double>(G2[i]) + K2[i]*dt;
				   	}

					h_HF_no_pulse(h1,G1_tmp,Ns,Nb,U,t1,t2,epsilon);

					thread th1(hG_commutator,ref(G1_tmp),ref(h1),ref(h1G1_comm));
                                        thread th2(collision_int,ref(I),ref(G2_tmp),Ns,Nb,U);
                                        thread th3(h2G2_commutator,ref(G2_tmp),ref(h1),ref(h2G2_comm),Ns,Nb);
                        		thread th4(Psi,ref(psi),ref(G1_tmp),Ns,Nb);
                        		th1.join();
                        		th2.join();
					th3.join();
                                        th4.join();

					thread th5(G1_step,ref(h1G1_comm),ref(I),ref(K1),ref(RK_result_G1),w,Ns,Nb);
					thread th6(G2_step,ref(h2G2_comm),ref(psi),ref(K2),ref(RK_result_G2),w,Ns,Nb,U);
					th5.join();
				    	th6.join();
		//			temp_t = t + dt;
				}
//				auto end = chrono::steady_clock::now();
//				avg_wall_t_per_loop += chrono::duration<double>(end -wall_t_start).count();

				for(int i =0; i<Ns*Ns*Nb*Nb;++i)
				{
					G1[i] += double(dt)*RK_result_G1[i]/(complex<double>(6));
				}
				for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i)
				{
					G2[i] += double(dt)*RK_result_G2[i]/(complex<double>(6));
				}
				t+=dt;
				cout<<double(t)/(tmax)*100<<"%\r";
				cout.flush();
			}

		}
};

int main()
{
	auto start = chrono::steady_clock::now();
	G1_G2_time_step SB(0.1,1,0,0,1,2,2,1);

	vector<complex<double>> G1(SB.Ns*SB.Ns*SB.Nb*SB.Nb);
	vector<complex<double>> G2(SB.Ns2*SB.Ns2*SB.Nb2*SB.Nb2);
	vector<double> n1;
	vector<vector<complex<double>>> G_diag;
	double avg_wall_t_per_loop = 0;

	SB.init_G(G1);
	SB.RK_steps(n1,G1,G2,G_diag,avg_wall_t_per_loop);
	write_to_rfile("n1.txt", n1);



	auto end = chrono::steady_clock::now();
	cout<<"Time stepping complete, elapsed time:"<< chrono::duration<double>(end - start).count()<< " sec \n";

	return 0;
}

