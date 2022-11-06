#ifndef write_to_file_h
#define write_to_file_h

#include <fstream>
#include <vector>
#include <complex>
using namespace std;

void write_to_rfile(string filename, vector<double> data)
{
        std::ofstream myfile;
        myfile.open(filename);

        for(int i = 0; i<(data.size()); i++)
        {
                myfile << data[i] << '\n';
        }
        myfile.close();
}

void write_to_cfile(string filename,vector<complex<double>> data)
{
        std::ofstream myfile1;
        std::ofstream myfile2;
        myfile1.open("imag_" + filename);
        myfile2.open("real_" + filename);
        for(int i = 0; i<(data.size()); i++){
		myfile1 <<data[i].imag()<<" ";
                myfile2 <<data[i].real()<<" ";
        }
        myfile1<<"\n";
        myfile2<<"\n";
        
        

        myfile1.close();
        myfile2.close();
}
void write_to_cfile2D(string filename,vector<vector<complex<double>>> data)
{
        std::ofstream myfile1;
        std::ofstream myfile2;
        myfile1.open("imag_" + filename);
        myfile2.open("real_" + filename);

	for(int j = 0; j < data[0].size(); ++j){
                for(int i = 0; i<(data.size()); i++){
	        myfile1 <<data[i][j].imag()<<" ";
	        myfile2 <<data[i][j].real()<<" ";

        }
	myfile1<<"\n";
	myfile2<<"\n";
	}
	

        myfile1.close();
        myfile2.close();
}


void read_in_cfile(string filename,vector<complex<double>> &data)
{
        std::fstream myfile1;
        std::fstream myfile2;
        myfile1.open("imag_" + filename);
        myfile2.open("real_" + filename);
	double temp_i;
	double temp_r;
	for(int i = 0; i < data.size(); ++i){
	        myfile1 >> temp_i;
	        myfile2 >> temp_r;
		data[i] = {temp_r,temp_i};
        }


        myfile1.close();
        myfile2.close();
}



#endif
