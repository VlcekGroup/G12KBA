//
//
//   The routine(s) in this file are a part of the
//                     G12KBA
//   suite, developed 2022, and copyrighted
//   to the authors: Cian Reeves and Vojtech Vlcek
//   at the University of California, Santa Barbara
//   and Khaled Ibrahim
//   at Lawrence Berkeley National Lab, Berkeley.
//
//
//  If you use or modify any part of this routine
//  the header should be kept and unmodified.
//
//
//

#ifndef rk_solver_h
#define rk_solver_h
#include "HF.h" //Holds various versions of  Hartree-Fock Hamiltonian
#include "psi.h"//Holds two particle occupation term
#include "pi.h" //Function for the polarization term in GW
#include "G1_G2_t_step.h"//Holds update steps for the one and two particle Green's function
#include "commutators.h"//Commutators for one and two particle Green's function with respective HF hamiltonian
#include "collision_integral.h"//Collision integral, describes many body effects on one particle Green's function
#include "write_to_file.h"
#include "print_matrix.h"
#include "fermi_func.h"
#include "read_input.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <thread>
#include <complex>
#include <chrono>
#include <string>

using namespace std;

extern complex<double> im;
extern bool HF;
extern bool EHM;
extern double decay_rate;
extern double t1;
extern double t2;
extern double t3;
extern double P;
extern double U;
extern int time_steps;
extern double dt_fixed;
extern double quench_strength;
extern int Nq;
extern int Ns;
extern int Ns2;
extern int Nb;
extern int Nb2;
extern string q_type;
extern vector<double> epsilon;
extern double quench_rate;
extern double quench_on_midpoint;
extern double quench_off_midpoint;
int RK4(vector<complex<double>> &G1x, vector<complex<double>> &G2xxxx,vector<complex<double>> &G2xyxy,vector<vector<complex<double>>> &G1x_diag)
{
    double tmax=double(time_steps)*dt_fixed;
    cout<<"Beginning Time Stepping Procedure\n";
    auto start = chrono::steady_clock::now();

    if(EHM == true){
        vector<double> V;
        for(int i = 0; i < Ns-1; ++i){
//            V.push_back(U*.5/(i+1));
            V.push_back(U*exp(-decay_rate*(i+1)));
        }
    }
    double w;
    double t = 0;
    double dt;
    double t_temp;
    vector<complex<double>> h1xG1x_comm(Ns*Ns*Nb*Nb);
    vector<complex<double>> Ix(Ns*Ns*Nb*Nb);
    vector<complex<double>> G1x_tmp(Ns*Ns*Nb*Nb);

    vector<complex<double>> h2xxG2xxxx_comm(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> h2xyG2xyxy_comm(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> psi_xxxx(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> psi_xyxy(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> pi_xxxx(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> pi_xyxy(Ns2*Ns2*Nb2*Nb2);

    vector<complex<double>> G2xxxx_tmp(Ns2*Ns2*Nb2*Nb2);
    vector<complex<double>> G2xyxy_tmp(Ns2*Ns2*Nb2*Nb2);
    while(t<tmax)
    {
        double num_dens = 0;
        vector<double> progress;
        G1x_diag.push_back(G1x);
        for(int i = 0; i < 1; ++i){
            for(int a = 0; a < 1; ++a){
                num_dens+= 2*G1x[(i+a*Ns)*Nb*Ns + i + a*Ns].imag();
            }
        }
        if(abs(num_dens)>3 || isnan(num_dens)){
            auto end_fail = chrono::steady_clock::now();
            cout<<"ERROR: Unphysical Number Density on site 1:"<<num_dens<<" elapsed time:"<< chrono::duration<double>(end_fail - start).count()<< " sec \n";
            return 1;
        }

        vector<complex<double>> RK_result_G1x(Ns*Ns*Nb*Nb);
        vector<complex<double>> K1x(Ns*Ns*Nb*Nb);

        vector<complex<double>> RK_result_G2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> RK_result_G2xyxy(Ns2*Ns2*Nb2*Nb2);

        vector<complex<double>> K2xxxx(Ns2*Ns2*Nb2*Nb2);
        vector<complex<double>> K2xyxy(Ns2*Ns2*Nb2*Nb2);
        for(int k = 0; k<4; ++k){
            vector<complex<double>> h1x(Ns*Ns*Nb*Nb);
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
                G1x_tmp[i] = G1x[i] + K1x[i]*dt;
            }
            for(int i = 0; i < Ns2*Ns2*Nb2*Nb2; ++i){
                G2xxxx_tmp[i] = G2xxxx[i] + K2xxxx[i]*dt;
                G2xyxy_tmp[i] = G2xyxy[i] + K2xyxy[i]*dt;
            }
            if(t_temp<50){
	                h_HFx_quench(h1x,G1x_tmp,G1x_tmp,Ns,Nb,U,t1,t2,0*fermi(t_temp,quench_on_midpoint,quench_rate)*quench_strength,epsilon,V,0);
            }
            else{
                if(q_type=="full"){
                    h_HFx_quench(h1x,G1x_tmp,G1x_tmp,Ns,Nb,U,t1,t2,fermi(t_temp,quench_on_midpoint,quench_rate)*quench_strength,epsilon,V,Nq);
                }
                if(q_type=="pulse"){
                    h_HFx_quench(h1x,G1x_tmp,G1x_tmp,Ns,Nb,U,t1,t2,(fermi(t_temp,quench_on_midpoint,quench_rate)*(1-fermi(t_temp,quench_off_midpoint,quench_rate)))*quench_strength,epsilon,V,Nq);
                }
            }
#if NOTHREADS
            hG_commutator(G1x_tmp,h1x,h1xG1x_comm);
            if(U!=0 && HF == false){
                collision_int(Ix,G2xyxy_tmp,G2xxxx_tmp,Ns,Nb,U,V);

                h2xyG2xyxy_commutator(G2xxxx_tmp,h1x,h1x,h2xxG2xxxx_comm,Ns,Nb);
                h2xyG2xyxy_commutator(G2xyxy_tmp,h1x,h1x,h2xyG2xyxy_comm,Ns,Nb);

                Psi_xxxx(psi_xxxx,G1x_tmp,Ns,Nb,U,V);
                Psi_xyxy(psi_xyxy,G1x_tmp,Ns,Nb,U,V);

                Pi_xyxy(pi_xxxx,G1x_tmp,G2xyxy_tmp,G2xxxx_tmp,U,Ns,Nb,V);
                Pi_xyxy(pi_xyxy,G1x_tmp,G2xxxx_tmp,G2xyxy_tmp,U,Ns,Nb,V);

                G2_step(h2xxG2xxxx_comm,psi_xxxx,pi_xxxx,K2xxxx,RK_result_G2xxxx,w,Ns,Nb);
                G2_step(h2xyG2xyxy_comm,psi_xyxy,pi_xyxy,K2xyxy,RK_result_G2xyxy,w,Ns,Nb);
            }
            G1_step(h1xG1x_comm,Ix,K1x,RK_result_G1x,w,Ns,Nb);
#else
            thread th_G1x_comm(hG_commutator,ref(G1x_tmp),ref(h1x),ref(h1xG1x_comm));

            thread th_collision_int(collision_int,ref(Ix),ref(G2xyxy_tmp),ref(G2xxxx_tmp),Ns,Nb,U,ref(V));

            thread th_h2xxG2xxxx_comm(h2xyG2xyxy_commutator,ref(G2xxxx_tmp),ref(h1x),ref(h1x),ref(h2xxG2xxxx_comm),Ns,Nb);
            thread th_h2xyG2xyxy_comm(h2xyG2xyxy_commutator,ref(G2xyxy_tmp),ref(h1x),ref(h1x),ref(h2xyG2xyxy_comm),Ns,Nb);

            thread th_psi_xxxx(Psi_xxxx,ref(psi_xxxx),ref(G1x_tmp),Ns,Nb,U,ref(V));
            thread th_psi_xyxy(Psi_xyxy,ref(psi_xyxy),ref(G1x_tmp),Ns,Nb,U,ref(V));

            thread th_pi_xxxx(Pi_xyxy,ref(pi_xxxx),ref(G1x_tmp),ref(G2xyxy_tmp),ref(G2xxxx_tmp),U,Ns,Nb,ref(V));
            thread th_pi_xyxy(Pi_xyxy,ref(pi_xyxy),ref(G1x_tmp),ref(G2xxxx_tmp),ref(G2xyxy_tmp),U,Ns,Nb,ref(V));

            th_G1x_comm.join();
            th_collision_int.join();
            th_h2xxG2xxxx_comm.join();
            th_h2xyG2xyxy_comm.join();
            th_psi_xxxx.join();
            th_psi_xyxy.join();
            th_pi_xxxx.join();
            th_pi_xyxy.join();

            thread th_G1x_step(G1_step,ref(h1xG1x_comm),ref(Ix),ref(K1x),ref(RK_result_G1x),w,Ns,Nb);

            thread th_G2xxxx_step(G2_step,ref(h2xxG2xxxx_comm),ref(psi_xxxx),ref(pi_xxxx),ref(K2xxxx),ref(RK_result_G2xxxx),w,Ns,Nb);
            thread th_G2xyxy_step(G2_step,ref(h2xyG2xyxy_comm),ref(psi_xyxy),ref(pi_xyxy),ref(K2xyxy),ref(RK_result_G2xyxy),w,Ns,Nb);

            th_G1x_step.join();
            th_G2xxxx_step.join();
            th_G2xyxy_step.join();
#endif
        }

        for(int i =0; i<Ns*Ns*Nb*Nb;++i)
        {
            G1x[i] += RK_result_G1x[i]*dt/6.0;

        }
        for(int i = 0; i<Ns2*Ns2*Nb2*Nb2;++i)
        {
            G2xxxx[i] += RK_result_G2xxxx[i]*dt/6.0;
            G2xyxy[i] += RK_result_G2xyxy[i]*dt/6.0;

        }

        t+=dt;
        cout<<double(t)/(tmax)*100<<"%\r";
        cout.flush();
        progress.push_back(double(t)/(tmax)*100);
        write_to_rfile("progress.txt",progress);

    }
    auto end_suc = chrono::steady_clock::now();
    cout<<"\n Execution Succesful, elapsed time:"<< chrono::duration<double>(end_suc - start).count()<< " sec \n";
    return 0;
}

#endif

 

