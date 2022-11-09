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

#ifndef read_input_h
#define read_input_h
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
using namespace std;


bool HF=false;
bool EHM=false;
double decay_rate=.7;
double t1=1;
double t2=0;
double t3=0;
double U=.5;
int time_steps=1000;
double dt_fixed=.02;
double quench_strength=1;
int Ns=2;
int Nq=0;
int Ns2=int(pow(Ns,2));
int Nb=1;
int Nb2=int(pow(Nb,2));
string q_type="full";
vector<double> epsilon(Nb);
vector<double> V;
double AS_rate=3;
double AS_midpoint=25;
double quench_rate=.2;
double quench_on_midpoint=50;
double quench_off_midpoint=55;



string read_var_val(string filename,string var)
{
    
    double j = 0;
    std::fstream myfile1;
    myfile1.open(filename);
    
    int line_counter = 0;
    string line;
    while (getline(myfile1,line)){
        line_counter++;
        if (line.find(var) != string::npos){
            break;
        }
    }
    myfile1.close();
    
    std::fstream myfile2;
    myfile2.open(filename);
    string s;
    for(int i = 0; i < line_counter + 1; ++i){
        getline(myfile2,s);
    }
    if(s==""){
        cout<<"Missing value for " <<var<< " initialising to 0\n";
        return "0";
    }
    return s;
    
    
}

void assign_vals()
{

    t1=stod(read_var_val("input","t1"));
    t2=stod(read_var_val("input","t2"));
    t3=stod(read_var_val("input","t3"));
    U=stod(read_var_val("input","U"));
    time_steps=stoi(read_var_val("input","time_steps"));
    dt_fixed=stod(read_var_val("input","dt"));
    quench_strength=stod(read_var_val("input","quench_strength"));
    Ns=stoi(read_var_val("input","Ns"));
    Ns2=int(pow(Ns,2));
    Nb=stoi(read_var_val("input","Nb"));
    Nb2=int(pow(Nb,2));

    Nq=stoi(read_var_val("input","Nq"));

    AS_rate=stod(read_var_val("input","AS_rate"));
    AS_midpoint=stod(read_var_val("input","AS_midpoint"));
    quench_on_midpoint=stod(read_var_val("input","quench_on_midpoint"));
    quench_off_midpoint=stod(read_var_val("input","quench_off_midpoint"));
    quench_rate=stod(read_var_val("input","quench_rate"));


    if(read_var_val("input","HF") == "true"){
        HF=true;
    }
    else if(read_var_val("input","HF") == "false"){
        HF=false;
    }
    else{
        cout<<"Invalid specification of HF calculation parameter, performing full HF-GKBA calculation.\n";
        HF=false;
    }

    if(read_var_val("input","EHM") == "true"){
        EHM=true;
    }
    else if(read_var_val("input","EHM") == "false"){
        EHM=false;
    }
    else{
        cout<<"Invalid specification for EHM parameter. performing full calculation for onsite model.\n";
        EHM=false;
    }


    if(read_var_val("input","q_type")=="full"){
	q_type="full";
    }
    else if(read_var_val("input","q_type")=="pulse"){
        q_type="pulse";
    }
    else if(read_var_val("input","q_type")=="none"){
	quench_strength=0;
	Nq=0;
    }
    else{
	cout<<"Invalid specification of quench type, no quench initiated\n";
    }

}


#endif /* read_input_h */
