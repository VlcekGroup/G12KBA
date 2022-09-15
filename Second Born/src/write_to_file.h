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
void write_to_cfile(string filename, vector<complex<double>> data)
{
        std::ofstream myfile;
        myfile.open(filename);

        for(int i = 0; i<(data.size()); i++)
        {
                myfile << data[i].real()<<"," <<data[i].imag()<< '\n';
        }
        myfile.close();
}

