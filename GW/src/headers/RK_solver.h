//
//  RK_solver.h
//  
//
//  Created by Cian Reeves on 9/9/22.
//
#include "HF.h" //Holds various versions of  Hartree-Fock Hamiltonian
#include "psi.h"//Holds two particle occupation term
#include "pi.h" //Function for the polarization term in GW
#include "G1_G2_t_step.h"//Holds update steps for the one and two particle Green's function
#include "commutators.h"//Commutators for one and two particle Green's function with respective HF hamiltonian
#include "collision_integral.h"//Collision integral, describes many body effects on one particle Green's function
#include "write_to_file.h"
#include "print_matrix.h"
#include "fermi_func.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <complex>
#include <fstream>
#include <chrono>
#include <string>
using namespace std;
const complex<double> im = {0,1};

class hubbard_RK_solver
{
    public:
        
	bool HF,EHM;
        double gamma,t1,t2,t3,P,U,tmax,dt_fixed,quench_strength;
        int Nq,Ns,Ns2,Nb,Nb2;
        string q_type;
	vector<double> epsilon;
        vector<double> V;
        hubbard_RK_solver(double interaction_strength = 1, double site_hopping = 1,double band_hopping = 0,double band_site_hopping = 0, double pulse_strength = 1, int sites = 2, int bands = 2,double band_gap = 1,double step_size = .02, double prop_time = 20,double w = .2,double alpha = 0,bool model_type = false, bool calc_type = false,int quench_extent=0,string quench_type="full")
        {
	    q_type=quench_type;
            quench_strength = w;
            Nq=quench_extent;
	    dt_fixed = step_size;
            tmax = prop_time;//propagation time
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
            epsilon[0] = -band_gap/2;//Onsite energy for band 0
            epsilon[1] = band_gap/2;//Onsite energy for band 1
	    EHM = model_type;
	    if(EHM == true){
		for(int i = 0; i < Ns-1; ++i){
			V.push_back(U*.5/(i+1));
//			V.push_back(U*exp(-alpha*(i+1)));
		    }
	   }
	   gamma = alpha;
	   HF =calc_type;
        }


        int RK4(vector<complex<double>> &G1u, vector<complex<double>> &G2uuuu,vector<complex<double>> &G2udud,vector<vector<complex<double>>> &G1u_diag)
        {
            double w;
	    double t = 0;
            double dt;
	double t_temp;
            vector<complex<double>> h1uG1u_comm(Ns*Ns*Nb*Nb);
            vector<complex<double>> Iu(Ns*Ns*Nb*Nb);
            vector<complex<double>> G1u_tmp(Ns*Ns*Nb*Nb);

            vector<complex<double>> h2uuG2uuuu_comm(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> h2udG2udud_comm(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> psi_uuuu(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> psi_udud(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> pi_uuuu(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> pi_udud(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> G2uuuu_tmp(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> G2udud_tmp(Ns2*Ns2*Nb2*Nb2);
            while(t<tmax)
            {
                double num_dens = 0;
	        vector<double> progress;
		G1u_diag.push_back(G1u);
                for(int i = 0; i < 1; ++i){
                    for(int a = 0; a < 1; ++a){
                        num_dens+= 2*G1u[(i+a*Ns)*Nb*Ns + i + a*Ns].imag();
                    }
                }
		if(abs(num_dens)>3 || isnan(num_dens)){
		cout<<"ERROR: Number Density "<<num_dens<<" Exceeds initial Filling\n";
		return 1;}

                vector<complex<double>> RK_result_G1u(Ns*Ns*Nb*Nb);
                vector<complex<double>> K1u(Ns*Ns*Nb*Nb);

                vector<complex<double>> RK_result_G2uuuu(Ns2*Ns2*Nb2*Nb2);
                vector<complex<double>> RK_result_G2udud(Ns2*Ns2*Nb2*Nb2);

                vector<complex<double>> K2uuuu(Ns2*Ns2*Nb2*Nb2);
                vector<complex<double>> K2udud(Ns2*Ns2*Nb2*Nb2);
                for(int k = 0; k<4; ++k)
                {
                    vector<complex<double>> h1u(Ns*Ns*Nb*Nb);
                    if(k == 1 or k == 2){
                        dt = dt_fixed/2;
                        w = 2;
                    }
                    else{
                        dt = dt_fixed;
                        w = 1;
                    }
		    if(k!=0){
			t_temp =t+dt;
		    }

                    for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
                        G1u_tmp[i] = G1u[i] + K1u[i]*dt;
                    }
              //#pragma omp parallel for
                    for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
                        G2uuuu_tmp[i] = G2uuuu[i] + K2uuuu[i]*dt;
                        G2udud_tmp[i] = G2udud[i] + K2udud[i]*dt;
                   }
	           if(t_temp<50){
		            h_HFx_quench(h1u,G1u_tmp,G1u_tmp,Ns,Nb,U,t1,t2,0,epsilon,V,0);
			}
		   else {
			    if(q_type=="full"){
		            	h_HFx_quench(h1u,G1u_tmp,G1u_tmp,Ns,Nb,U,t1,t2,fermi(t_temp,50,.2)*quench_strength,epsilon,V,Nq);
			    }
			    if(q_type=="pulse"){
		        	h_HFx_quench(h1u,G1u_tmp,G1u_tmp,Ns,Nb,U,t1,t2,(fermi(t_temp,50,.2)*(1-fermi(t_temp,55,.2)))*quench_strength,epsilon,V,Nq);
		            }
                       }
#if NOTHREADS
                    hG_commutator(G1u_tmp,h1u,h1uG1u_comm);
                    if(U!=0 && HF == false){
                        collision_int(Iu,G2udud_tmp,G2uuuu_tmp,Ns,Nb,U,V);

                        h2xyG2xyxy_commutator(G2uuuu_tmp,h1u,h1u,h2uuG2uuuu_comm,Ns,Nb);
                        h2xyG2xyxy_commutator(G2udud_tmp,h1u,h1u,h2udG2udud_comm,Ns,Nb);

                        Psi_xxxx(psi_uuuu,G1u_tmp,Ns,Nb,U,V);
                        Psi_xyxy(psi_udud,G1u_tmp,Ns,Nb,U,V);

                        Pi_xyxy(pi_uuuu,G1u_tmp,G2udud_tmp,G2uuuu_tmp,U,Ns,Nb,V);
                        Pi_xyxy(pi_udud,G1u_tmp,G2uuuu_tmp,G2udud_tmp,U,Ns,Nb,V);

                        G2_step(h2uuG2uuuu_comm,psi_uuuu,pi_uuuu,K2uuuu,RK_result_G2uuuu,w,Ns,Nb);
                        G2_step(h2udG2udud_comm,psi_udud,pi_udud,K2udud,RK_result_G2udud,w,Ns,Nb);
                    }
                    G1_step(h1uG1u_comm,Iu,K1u,RK_result_G1u,w,Ns,Nb);
#else
                    thread th_G1u_comm(hG_commutator,ref(G1u_tmp),ref(h1u),ref(h1uG1u_comm));

                    thread th_collision_int_u(collision_int,ref(Iu),ref(G2udud_tmp),ref(G2uuuu_tmp),Ns,Nb,U,ref(V));

                    thread th_h2uuG2uuuu_comm(h2xyG2xyxy_commutator,ref(G2uuuu_tmp),ref(h1u),ref(h1u),ref(h2uuG2uuuu_comm),Ns,Nb);
                    thread th_h2udG2udud_comm(h2xyG2xyxy_commutator,ref(G2udud_tmp),ref(h1u),ref(h1u),ref(h2udG2udud_comm),Ns,Nb);

                    thread th_psi_uuuu(Psi_xxxx,ref(psi_uuuu),ref(G1u_tmp),U,Ns,Nb,ref(V));
                    thread th_psi_uuuu(Psi_xyxy,ref(psi_udud),ref(G1u_tmp),ref(G1u_tmp),U,Ns,Nb,ref(V));

                    thread th_pi_uuuu(Pi_xyxy,ref(pi_uuuu),ref(G1u_tmp),ref(G2udud_tmp),ref(G2uuuu_tmp),U,Ns,Nb,ref(V));
                    thread th_pi_udud(Pi_xyxy,ref(pi_udud),ref(G1u_tmp),ref(G2uuuu_tmp),ref(G2udud_tmp),U,Ns,Nb,ref(V));

                    th_G1u_comm.join();
                    th_collision_int_u.join();
                    th_h2uuG2uuuu_comm.join();
                    th_h2udG2udud_comm.join();
                    th_psi_uuuu.join()
                    th_psi_udud.join()
                    th_pi_uuuu.join()
                    th_pi_udud.join()

                    thread th_G1u_step(G1_step,ref(h1uG1u_comm),ref(Iu),ref(K1u),ref(RK_result_G1u),w,Ns,Nb);

                    thread th_G2uuuu_step(G2_step,ref(h2uuG2uuuu_comm),ref(psi_uuuu),ref(pi),ref(K2_uuuu),ref(RK_result_G2uuuu),w,Ns,Nb);
                    thread th_G2udud_step(G2_step,ref(h2udG2udud_comm),ref(psi_udud),ref(pi),ref(K2_udud),ref(RK_result_G2udud),w,Ns,Nb);

                    th_G1u_step.join();
                    th_G2uuuu_step.join();
                    th_G2udud_step.join();
#endif
		}

                for(int i =0; i<Ns*Ns*Nb*Nb;++i)
                {
                    G1u[i] += RK_result_G1u[i]*dt/6.0;

                }
                for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i)
                {
                    G2uuuu[i] += RK_result_G2uuuu[i]*dt/6.0;
                    G2udud[i] += RK_result_G2udud[i]*dt/6.0;

                }

                t+=dt;
                cout<<double(t)/(tmax)*100<<"%\r";
                cout.flush();
		progress.push_back(double(t)/(tmax)*100);
		write_to_rfile("progress.txt",progress);

            }
	return 0;
        }

        int Adiabatic_switching(vector<complex<double>> &G1u, vector<complex<double>> &G2uuuu,vector<complex<double>> &G2udud,vector<vector<complex<double>>> &G1u_diag)
        {
            double w;
       	    double t_temp;
            double U_temp;
	    double t = 0;
            double dt;
            vector<complex<double>> h1uG1u_comm(Ns*Ns*Nb*Nb);
            vector<complex<double>> Iu(Ns*Ns*Nb*Nb);
            vector<complex<double>> G1u_tmp(Ns*Ns*Nb*Nb);

            vector<complex<double>> h2uuG2uuuu_comm(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> h2udG2udud_comm(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> psi_uuuu(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> psi_udud(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> pi_uuuu(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> pi_udud(Ns2*Ns2*Nb2*Nb2);

            vector<complex<double>> G2uuuu_tmp(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> G2udud_tmp(Ns2*Ns2*Nb2*Nb2);
            int count = 0;
	    while(count < 20)
            {
		vector<double> progress;
		if(abs(U_temp-U)<=0){
		count+=1;}

                vector<complex<double>> RK_result_G1u(Ns*Ns*Nb*Nb);
                vector<complex<double>> K1u(Ns*Ns*Nb*Nb);

                vector<complex<double>> RK_result_G2uuuu(Ns2*Ns2*Nb2*Nb2);
                vector<complex<double>> RK_result_G2udud(Ns2*Ns2*Nb2*Nb2);

                vector<complex<double>> K2uuuu(Ns2*Ns2*Nb2*Nb2);
                vector<complex<double>> K2udud(Ns2*Ns2*Nb2*Nb2);
		t_temp = t;
                for(int k = 0; k<4; ++k)
                {
                    vector<complex<double>> h1u(Ns*Ns*Nb*Nb);
                    if(k == 1 or k == 2){
                        dt = .1/2;
                        w = 2;
                    }
                    else{
                        dt = .1;
                        w = 1;
                    }
		    if(k!=0){
			t_temp =t+dt;
		    }

                    for(int i = 0; i < Ns*Ns*Nb*Nb; ++i){
                        G1u_tmp[i] = complex<double>(G1u[i]) + K1u[i]*dt;
                    }

                    for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
                        G2uuuu_tmp[i] = complex<double>(G2uuuu[i]) + K2uuuu[i]*dt;
                        G2udud_tmp[i] = complex<double>(G2udud[i]) + K2udud[i]*dt;
                    }
		    U_temp = fermi(t_temp,25,3)*U;
		    if(EHM == true){
                    for(int i = 0; i < Ns-1; ++i){
                        V[i]=(U_temp*.5/(i+1));
			//V[i]=(U_temp*exp(-gamma*(i+1)));
                    }
           }
	            h_HFx_quench(h1u,G1u_tmp,G1u_tmp,Ns,Nb,U_temp,t1,t2,0,epsilon,V,0);
#if NOTHREADS
                    hG_commutator(G1u_tmp,h1u,h1uG1u_comm);
                    if(U!=0 && HF == false){
                        collision_int(Iu,G2udud_tmp,G2uuuu_tmp,Ns,Nb,U_temp,V);

                        h2xyG2xyxy_commutator(G2uuuu_tmp,h1u,h1u,h2uuG2uuuu_comm,Ns,Nb);
                        h2xyG2xyxy_commutator(G2udud_tmp,h1u,h1u,h2udG2udud_comm,Ns,Nb);

                        Psi_xxxx(psi_uuuu,G1u_tmp,Ns,Nb,U_temp,V);
                        Psi_xyxy(psi_udud,G1u_tmp,Ns,Nb,U_temp,V);

                        Pi_xyxy(pi_uuuu,G1u_tmp,G2udud_tmp,G2uuuu_tmp,U_temp,Ns,Nb,V);
                        Pi_xyxy(pi_udud,G1u_tmp,G2uuuu_tmp,G2udud_tmp,U_temp,Ns,Nb,V);

                        G2_step(h2uuG2uuuu_comm,psi_uuuu,pi_uuuu,K2uuuu,RK_result_G2uuuu,w,Ns,Nb);
                        G2_step(h2udG2udud_comm,psi_udud,pi_udud,K2udud,RK_result_G2udud,w,Ns,Nb);
                    }
                    G1_step(h1uG1u_comm,Iu,K1u,RK_result_G1u,w,Ns,Nb);
#else
                    thread th_G1u_comm(hG_commutator,ref(G1u_tmp),ref(h1u),ref(h1uG1u_comm));

                    thread th_collision_int_u(collision_int,ref(Iu),ref(G2udud_tmp),ref(G2uuuu_tmp),Ns,Nb,U_temp,ref(V));

                    thread th_h2uuG2uuuu_comm(h2xyG2xyxy_commutator,ref(G2uuuu_tmp),ref(h1u),ref(h1u),ref(h2uuG2uuuu_comm),Ns,Nb);
                    thread th_h2udG2udud_comm(h2xyG2xyxy_commutator,ref(G2udud_tmp),ref(h1u),ref(h1u),ref(h2udG2udud_comm),Ns,Nb);

                    thread th_psi_uuuu(Psi_xxxx,ref(psi_uuuu),ref(G1u_tmp),U_temp,Ns,Nb,ref(V));
                    thread th_psi_uuuu(Psi_xyxy,ref(psi_udud),ref(G1u_tmp),ref(G1u_tmp),U_temp,Ns,Nb,ref(V));

                    thread th_pi_uuuu(Pi_xyxy,ref(pi_uuuu),ref(G1u_tmp),ref(G2udud_tmp),ref(G2uuuu_tmp),U_temp,Ns,Nb,ref(V));
                    thread th_pi_udud(Pi_xyxy,ref(pi_udud),ref(G1u_tmp),ref(G2uuuu_tmp),ref(G2udud_tmp),U_temp,Ns,Nb,ref(V));

                    th_G1u_comm.join();
                    th_collision_int_u.join();
                    th_h2uuG2uuuu_comm.join();
                    th_h2udG2udud_comm.join();
                    th_psi_uuuu.join()
                    th_psi_udud.join()
                    th_pi_uuuu.join()
                    th_pi_udud.join()

                    thread th_G1u_step(G1_step,ref(h1uG1u_comm),ref(Iu),ref(K1u),ref(RK_result_G1u),w,Ns,Nb);

                    thread th_G2uuuu_step(G2_step,ref(h2uuG2uuuu_comm),ref(psi_uuuu),ref(pi),ref(K2_uuuu),ref(RK_result_G2uuuu),w,Ns,Nb);
                    thread th_G2udud_step(G2_step,ref(h2udG2udud_comm),ref(psi_udud),ref(pi),ref(K2_udud),ref(RK_result_G2udud),w,Ns,Nb);

                    th_G1u_step.join();
                    th_G2uuuu_step.join();
                    th_G2udud_step.join();
#endif
		}

                for(int i =0; i<Ns*Ns*Nb*Nb;++i)
                {
                    G1u[i] += RK_result_G1u[i]*dt/6.0;

                }
                for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i)
                {
                    G2uuuu[i] += RK_result_G2uuuu[i]*dt/6.0;
                    G2udud[i] += RK_result_G2udud[i]*dt/6.0;

                }
		t+=dt;
	        cout<<"U = "<<U_temp<<"\r";
                cout.flush();
		progress.push_back(double(U_temp)/(U)*100);
		write_to_rfile("progress.txt",progress);

            }
	return 0;
        }



        int RK6(vector<double> &n1, vector<complex<double>> &G1u, vector<complex<double>> &G2uuuu,vector<complex<double>> &G2udud,vector<vector<complex<double>>> &Gu_diag, double &avg_wall_t_per_loop)

    {
        double w;
        double t = 0;
        double dt = dt_fixed;
        vector<complex<double>> h1uG1u_comm(Ns*Ns*Nb*Nb);
        vector<complex<double>> Iu(Ns*Ns*Nb*Nb);
        vector<complex<double>> G1u_tmp(Ns*Ns*Nb*Nb);

        vector<complex<double>> h2uuG2uuuu_comm(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> h2udG2udud_comm(Ns2*Ns2*Nb2*Nb2);

        vector<complex<double>> psi_uuuu(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> psi_udud(Ns2*Ns2*Nb2*Nb2);

        vector<complex<double>> pi_uuuu(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> pi_udud(Ns2*Ns2*Nb2*Nb2);

        vector<complex<double>> G2uuuu_tmp(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> G2udud_tmp(Ns2*Ns2*Nb2*Nb2);


        while(t<tmax)
        {
//                G_diag.push_back(G1);
            double num_dens = 0;

            for(int i = 0; i < 1; ++i)
            {
                for(int a = 0; a <1; ++a)
                {
                    num_dens+= 2*G1u[(i+a*Ns)*Nb*Ns + i + a*Ns].imag();
                }
            }
		if(abs(num_dens)>10){
		cout<<"ERROR: Number Density "<<num_dens<<" Exceeds initial Filling\n";
		return 1;}

            n1.push_back(num_dens);

            vector<complex<double>> RK_result_G1u(Ns*Ns*Nb*Nb);
            vector<vector<complex<double>>> K1u(6, vector<complex<double>> (Ns*Ns*Nb*Nb));

            vector<complex<double>> RK_result_G2uuuu(Ns2*Ns2*Nb2*Nb2);
            vector<complex<double>> RK_result_G2udud(Ns2*Ns2*Nb2*Nb2);

            vector<vector<complex<double>>> K2uuuu(6, vector<complex<double>> (Ns2*Ns2*Nb2*Nb2));
            vector<vector<complex<double>>> K2udud(6, vector<complex<double>> (Ns2*Ns2*Nb2*Nb2));

//                auto wall_t_start = chrono::steady_clock::now();
            for(int k = 0; k<6; ++k)
            {
                vector<complex<double>> h1u(Ns*Ns*Nb*Nb);
                vector<int> indices;
                vector<double> weights;
                if(k == 0){
                    for(int i = 0; i<Ns*Ns*Nb*Nb; ++i){
                        G1u_tmp[i] = complex<double>(G1u[i]);
                    }
              //#pragma omp parallel for
                    for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i){
                        G2uuuu_tmp[i] = complex<double>(G2uuuu[i]);
                        G2udud_tmp[i] = complex<double>(G2udud[i]);
                    }
                }
                else if(k == 1){
                    indices.push_back(0);
                    weights.push_back(dt_fixed/4);
                    }
                else if(k==2){
                    indices.push_back(0);
                    indices.push_back(1);
                    weights.push_back(dt_fixed/8);
                    weights.push_back(dt_fixed/8);
                }
                else if(k==3){
                    indices.push_back(1);
                    indices.push_back(2);
                    weights.push_back(-dt_fixed/2);
                    weights.push_back(dt_fixed);
                }
                
                else if(k==4){
                    indices.push_back(0);
                    indices.push_back(3);
                    weights.push_back(3*dt_fixed/16);
                    weights.push_back(9*dt_fixed/16);
                }
                
                else if(k==5){
                    indices.push_back(0);
                    indices.push_back(1);
                    indices.push_back(2);
                    indices.push_back(3);
                    indices.push_back(4);
                    weights.push_back(-3*dt_fixed/7);
                    weights.push_back(2*dt_fixed/7);
                    weights.push_back(12*dt_fixed/7);
                    weights.push_back(-12*dt_fixed/7);
                    weights.push_back(8*dt_fixed/7);

                }

                for(int i = 0; i<Ns*Ns*Nb*Nb; ++i){
                    for(int j = 0; j<indices.size(); ++j){
                        G1u_tmp[i] = complex<double>(G1u[i]) + K1u[indices[j]][i]*weights[j];
                    }
                }
          //#pragma omp parallel for
                for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i){
                    for(int j = 0; j<indices.size(); ++j){
                        G2uuuu_tmp[i] = complex<double>(G2uuuu[i]) + K2uuuu[indices[j]][i]*weights[j];
                        G2udud_tmp[i] = complex<double>(G2udud[i]) + K2udud[indices[j]][i]*weights[j];
                    }
                }


	            h_HFx_quench(h1u,G1u_tmp,G1u_tmp,Ns,Nb,U,t1,t2,quench_strength,epsilon,V,0);
#if NOTHREADS
                    hG_commutator(G1u_tmp,h1u,h1uG1u_comm);
                    if(U!=0 && HF == false){
                        collision_int(Iu,G2udud_tmp,G2uuuu_tmp,Ns,Nb,U,V);

                        h2xyG2xyxy_commutator(G2uuuu_tmp,h1u,h1u,h2uuG2uuuu_comm,Ns,Nb);
                        h2xyG2xyxy_commutator(G2udud_tmp,h1u,h1u,h2udG2udud_comm,Ns,Nb);

                        Psi_xxxx(psi_uuuu,G1u_tmp,Ns,Nb,U,V);
                        Psi_xyxy(psi_udud,G1u_tmp,Ns,Nb,U,V);

                        Pi_xyxy(pi_uuuu,G1u_tmp,G2udud_tmp,G2uuuu_tmp,U,Ns,Nb,V);
                        Pi_xyxy(pi_udud,G1u_tmp,G2uuuu_tmp,G2udud_tmp,U,Ns,Nb,V);

                        G2_step(h2uuG2uuuu_comm,psi_uuuu,pi_uuuu,K2uuuu[k],RK_result_G2uuuu,w,Ns,Nb);
                        G2_step(h2udG2udud_comm,psi_udud,pi_udud,K2udud[k],RK_result_G2udud,w,Ns,Nb);
                    }
                    G1_step(h1uG1u_comm,Iu,K1u[k],RK_result_G1u,w,Ns,Nb);
#else
                    thread th_G1u_comm(hG_commutator,ref(G1u_tmp),ref(h1u),ref(h1uG1u_comm));

                    thread th_collision_int_u(collision_int,ref(Iu),ref(G2udud_tmp),ref(G2uuuu_tmp),Ns,Nb,U,ref(V));

                    thread th_h2uuG2uuuu_comm(h2xyG2xyxy_commutator,ref(G2uuuu_tmp),ref(h1u),ref(h1u),ref(h2uuG2uuuu_comm),Ns,Nb);
                    thread th_h2udG2udud_comm(h2xyG2xyxy_commutator,ref(G2udud_tmp),ref(h1u),ref(h1u),ref(h2udG2udud_comm),Ns,Nb);

                    thread th_psi_uuuu(Psi_xxxx,ref(psi_uuuu),ref(G1u_tmp),U,Ns,Nb,ref(V));
                    thread th_psi_uuuu(Psi_xyxy,ref(psi_udud),ref(G1u_tmp),ref(G1u_tmp),U,Ns,Nb,ref(V));

                    thread th_pi_uuuu(Pi_xyxy,ref(pi_uuuu),ref(G1u_tmp),ref(G2udud_tmp),ref(G2uuuu_tmp),U,Ns,Nb,ref(V));
                    thread th_pi_udud(Pi_xyxy,ref(pi_udud),ref(G1u_tmp),ref(G2uuuu_tmp),ref(G2udud_tmp),U,Ns,Nb,ref(V));

                    th_G1u_comm.join();
                    th_collision_int_u.join();
                    th_h2uuG2uuuu_comm.join();
                    th_h2udG2udud_comm.join();
                    th_psi_uuuu.join()
                    th_psi_udud.join()
                    th_pi_uuuu.join()
                    th_pi_udud.join()

                    thread th_G1u_step(G1_step,ref(h1uG1u_comm),ref(Iu),ref(K1u),ref(RK_result_G1u),w,Ns,Nb);

                    thread th_G2uuuu_step(G2_step,ref(h2uuG2uuuu_comm),ref(psi_uuuu),ref(pi),ref(K2_uuuu),ref(RK_result_G2uuuu),w,Ns,Nb);
                    thread th_G2udud_step(G2_step,ref(h2udG2udud_comm),ref(psi_udud),ref(pi),ref(K2_udud),ref(RK_result_G2udud),w,Ns,Nb);

                    th_G1u_step.join();
                    th_G2uuuu_step.join();
                    th_G2udud_step.join();
#endif
            }
//                auto end = chrono::steady_clock::now();
//                avg_wall_t_per_loop += chrono::duration<double>(end -wall_t_start).count();

            for(int i =0; i<Ns*Ns*Nb*Nb;++i){
                G1u[i] += double(dt)*(7.0*K1u[0][i] + 32.0*K1u[2][i] + 12.0*K1u[3][i] + 32.0*K1u[4][i] + 7.0*K1u[5][i])/(complex<double>(90));

            }
	    
            for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i){
                G2uuuu[i] += double(dt)*(7.0*K2uuuu[0][i] + 32.0*K2uuuu[2][i] + 12.0*K2uuuu[3][i] + 32.0*K2uuuu[4][i] + 7.0*K2uuuu[5][i])/(complex<double>(90));
                G2udud[i] += double(dt)*(7.0*K2udud[0][i] + 32.0*K2udud[2][i] + 12.0*K2udud[3][i] + 32.0*K2udud[4][i] + 7.0*K2udud[5][i])/(complex<double>(90));

            }
            t+=dt;
            cout<<double(t)/(tmax)*100<<"%\r";
            cout.flush();
        }

    }
};
