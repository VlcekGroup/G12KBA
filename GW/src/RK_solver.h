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

class hubbard_RK_solver
{
    public:
        
        double t1,t2,t3,P,U,tmax,dt_fixed;
        int Ns,Ns2,Nb,Nb2;
        vector<double> epsilon;
        hubbard_RK_solver(double interaction_strength = 1, double site_hopping = 1,double band_hopping = 0,double band_site_hopping = 0, double pulse_strength = 1, int sites = 2, int bands = 2,double band_gap = 1)
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
            epsilon[0] = -band_gap/2;//Onsite energy for band 0
            epsilon[1] = band_gap/2;//Onsite energy for band 1
        }


        void init_G(vector<complex<double>> &G1)
        {
            int particle_number = 0;
            for(int i = 0; i < int(Ns/2); ++i){
                for(int a = 0; a < Nb; ++a){
                    G1[(i+a*Ns)*Nb*Ns + i + a*Ns] = im;
                }
            }
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
            vector<complex<double>> pi(Ns2*Ns2*Nb2*Nb2);
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
                vector<complex<double>> RK_result_G1(Ns*Ns*Nb*Nb);
                vector<complex<double>> RK_result_G2(Ns2*Ns2*Nb2*Nb2);
                vector<complex<double>> K1(Ns*Ns*Nb*Nb);
                vector<complex<double>> K2(Ns2*Ns2*Nb2*Nb2);

                n1.push_back(num_dens);
//                auto wall_t_start = chrono::steady_clock::now();
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
              //#pragma omp parallel for
                        for(int i =0; i<Ns2*Ns2*Nb2*Nb2;++i)
                        {
                            G2_tmp[i] = complex<double>(G2[i]) + K2[i]*dt;
                        }
                    h_HF_no_pulse(h1,G1_tmp,Ns,Nb,U,t1,t2,epsilon);
#if NOTHREADS
                    hG_commutator(G1_tmp,h1,h1G1_comm);
                    collision_int(I,G2_tmp,Ns,Nb,U);
                    h2G2_commutator(G2_tmp,h1,h2G2_comm,Ns,Nb);
                    Psi(psi,G1_tmp,Ns,Nb);
                    Pi(pi,G1_tmp,G2_tmp,U,Ns,Nb);
                    G1_step(h1G1_comm,I,K1,RK_result_G1,w,Ns,Nb);
                    G2_step(h2G2_comm,psi,pi,K2,RK_result_G2,w,Ns,Nb,U);
#else
                    thread th1(hG_commutator,ref(G1_tmp),ref(h1),ref(h1G1_comm));
                    thread th2(collision_int,ref(I),ref(G2_tmp),Ns,Nb,U);
                    thread th3(h2G2_commutator,ref(G2_tmp),ref(h1),ref(h2G2_comm),Ns,Nb);
                    thread th4(Psi,ref(psi),ref(G1_tmp),Ns,Nb);
                    thread th5(Pi,ref(pi),ref(G1_tmp),ref(G2_tmp),U,Ns,Nb);
                    th1.join();
                    th2.join();
                    th3.join();
                    th4.join();
                    th5.join();
                    
                    thread th6(G1_step,ref(h1G1_comm),ref(I),ref(K1),ref(RK_result_G1),w,Ns,Nb);
                    thread th7(G2_step,ref(h2G2_comm),ref(psi),ref(pi),ref(K2),ref(RK_result_G2),w,Ns,Nb,U);
                    th6.join();
                    th7.join();
#endif
                }
//                auto end = chrono::steady_clock::now();
//                avg_wall_t_per_loop += chrono::duration<double>(end -wall_t_start).count();

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
