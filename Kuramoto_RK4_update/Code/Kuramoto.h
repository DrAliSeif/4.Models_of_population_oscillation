/************************************************************************************************/
/*** Topic: Kuramoto model with Runge-Kutta 4th Order Method for 500 nodes in fully connection  ***/
/***   network with cupling                                                                   ***/
/*** Version Release 17.12 rev 11256                                                Ali-Seif  ***/
/*** Date: 5/8/2021                                                                           ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include<iostream>                                                              //for cout
#include<vector>                                                                //for vector
#include <fstream>                                                              //infile
#include <random>                                                               //random_device & ...
#include <chrono>                                                               //chrono
using namespace std;                                                            //using standard 
typedef vector<int> vector_1d_int;                                              //define type of vector integer one dimensional
using vec2i = vector<vector<int>>;                                              //using of vector integer two dimensional
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                class Kuramoto                                  @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class Kuramoto {                                                                //create and define class 
private:                                                                        //private values
    const double Pi=3.14159265359;                                              //pi number
    const int N;                                                                //number of node
    const double dt;                                                            //step lenth
    double tfinal;                                                              //time final
    string address_matrix;                                                      //name matrix file
    double* phi = new double[N];                                                //create 1d pointer
    double* W = new double[N];                                                  //create 1d pointer
    vec2i connection;                                                           //connection in [row][column]
public:                                                                         //public values
                                                                                //
    Kuramoto                                                                    //
    (int N_, double dt_, double t_final_, double coupling_, string address_matrix_) ://input data in class
        N(N_), dt(dt_), tfinal(t_final_), coupling(coupling_), address_matrix(address_matrix_){}//change name input to privet
    double syncer = 0;                                                          //sync variable to print and type in the int main
    double coupling;                                                            //stringh cupling
    void matrix_rider();                                                        //read matrix connection
    void matrix_fully();                                                        //create matrix fully connection
    void theta_print();                                                         //crate rand between 0 to 1 uniform_real_distribution for initial phi and print this
    void theta_rider();                                                         //read rand between 0 to 1 uniform_real_distribution for initial phi
    void Omega_print();                                                         //crate rand mean 1 and standard deviation 0.01 normal_distribution and print this 
    void Omega_rider();                                                         //read rand mean 1 and standard deviation 0.01 normal_distribution
    void Run();                                                                 //run program
    void scale_2_pi(double*);                                                   //put variable phi in 0 to 2pi
    double sync(double*);                                                       //calculate syncrony
    void runge_kutta4_integrator(double*);                                      //calculate runge kutta4
    void dydt(double*, const double*);                                          //calculate dphi/dt
};                                                                              //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                          Read matrix connection                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*void Kuramoto::matrix_rider() {                                               //read matrix connection     
    std::ifstream ifile(address_matrix + ".txt");                               //read address of file .txt
    if (ifile.is_open()) {                                                      //if file was available 
        vector_1d_int matrix_1d;                                                //create vector integer one dimensional(matrix_1d)
        int num;                                                                //create integer number(num)
        while (ifile >> num) {                                                  //Set the read number of the file to the defined integer(num)
            matrix_1d.push_back(num);                                           //set defined integer(num) into the One after the last cell vector integer one dimensional(matrix_1d)
            if (matrix_1d.size() == N) {                                        //if size of vector is full
                connection.push_back(matrix_1d);                                //set vector integer one dimensional(matrix_1d) into the One after the last cell vector integer tow dimensional(matrix_2d)
                matrix_1d.clear();                                              //clean all cels of vector integer one dimensional(matrix_1d) and delete it
            }                                                                   //
        }                                                                       //
    }                                                                           //
    else {                                                                      //else if file wasn't available
        std::cout << "There was an error opening the input file!\n";            //print "..."
        exit(1);                                                                //means program(process) terminate normally unsuccessfully..
    }                                                                           //
} */                                                                            //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                          Read matrix connection                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::matrix_fully() {                                                       //read matrix connection     
    float** conn;
    conn = new float* [N];
    for (int i = 0; i < N; i++)
    {
        conn[i] = new float[N];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                conn[i][j] = 1;
            }
            else
            {
                conn[i][j] = 0;
            }

        }
    }
    vector_1d_int matrix_1d;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_1d.push_back(conn[i][j]);
        }
        connection.push_back(matrix_1d);
        matrix_1d.clear();
    }
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     theta                                      @@@@
//@@@              crate rand between 0 to 1 uniform_real_distribution               @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::theta_print() {                                                  //calculate initial theta
    ofstream theta("theta.txt", ios::out | ios::trunc);                         //create file .txt for save sync volt
    random_device rd;                                                           //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());                                                          //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> uniform_distribution_teta(0.0, 1.0);      //Uniform real distribution from 0 to 1
    for (int i = 0; i < N; i++) {                                               //for loop all node
        phi[i] = uniform_distribution_teta(gen) * 2 * Pi;
        theta << phi[i] << endl;
    }
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                  theta rider                                   @@@@
//@@@              crate rand between 0 to 1 uniform_real_distribution               @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::theta_rider() {                                                  //calculate initial theta
    ifstream infile;                                                            //create space for read file
    infile.open("theta.txt");                                                   //read address and file .txt
    if (infile.fail())                                                          //if could not find file
    {                                                                           //
        cout << "Could not open file numbers." << "\n";                         //type this paragraph
    }                                                                           //
    else {                                                                      //if find file
        double data;                                                            //define double for transport .txt file
        infile >> data;                                                         //transport first data from file to variable 
        int i = 0;                                                              //create conter for cont 
        while (!infile.eof()) {                                                 //while ine the file was data do this  
            phi[i] = data;                                                      //save variable in to the pointer1D randborg
            i++;                                                                //i=i+1
            infile >> data;                                                     //transport next data from file to variable 
        }                                                                       //
    }                                                                           //
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     Omega                                      @@@@
//@@@         rand mean 1 and standard deviation 0.01 normal_distribution            @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::Omega_print(){                                                   //calculate w
    ofstream Omega("Omega.txt", ios::out | ios::trunc);                         //create file .txt for save sync volt
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();     //using seed time
    default_random_engine generator(seed);                                      //random generator
    normal_distribution<double> dist(0.0, 1.0);                                 //creat random that This distribution produces random numbers around the distribution mean 1 with a specific standard deviation 0.01.
    for (int i = 0; i < N; i++) {                                               //for loop all node
        W[i] = dist(generator);                                                 //create w random for each node
        Omega << W[i] << endl;
    }
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     Omega                                      @@@@
//@@@         rand mean 1 and standard deviation 0.01 normal_distribution            @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::Omega_rider() {                                                  //calculate w
    ifstream infile;                                                            //create space for read file
    infile.open("Omega.txt");                                                   //read address and file .txt
    if (infile.fail())                                                          //if could not find file
    {                                                                           //
        cout << "Could not open file numbers." << "\n";                         //type this paragraph
    }                                                                           //
    else {                                                                      //if find file
        double data;                                                            //define double for transport .txt file
        infile >> data;                                                         //transport first data from file to variable 
        int i = 0;                                                              //create conter for cont 
        while (!infile.eof()) {                                                 //while ine the file was data do this  
            W[i] = data;                                                        //save variable in to the pointer1D randborg
            i++;                                                                //i=i+1
            infile >> data;                                                     //transport next data from file to variable 
        }                                                                       //
    }                                                                           //
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                     RUN                                        @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::Run() {                                                          //call Run function
    //ofstream flactuated("syncrony_flactuated.txt", ios::out | ios::trunc);    //create file .txt for save sync volt
    int conter = 0;
    double sum_sync = 0;
    for (double t = 0; t < tfinal; t += dt) {                                   //loop over time step
        runge_kutta4_integrator(phi);                                           //calculate runge_kutta for all node
        scale_2_pi(phi);                                                        //change values to be in range 0 to 2*Pi
        if (t > 4) {
            conter++;
            sum_sync += sync(phi);
        }
        //flactuated << t<<'\t'<<sync(phi) << endl;
    }   
    syncer = sum_sync / conter;
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                         scaled in 2_pi if Cross 2pi                            @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::scale_2_pi(double* phi){                                         //put all phi node in 2pi scale
    for (int i = 0; i < N; i++)                                                 //for loop all node
        phi[i] = phi[i] - 2 * Pi * static_cast<int>(phi[i] / (2 * Pi));         //if each node was upper 2pi put on in the 0 to 2pi
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     sync                                       @@@@
//@@@                              Order parameter synchrony                         @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double Kuramoto::sync(double* phi){                                             //calculate syncrony
    float rc = 0.0, rs = 0.0;                                                   //create variable for //r_data << t * dt << "\t" << order_param(phi) << endl;
    for (int j = 0; j < N; j++) {                                               //for loop to change each node
        rc += cos(phi[j]);                                                      //calculate sum all cos phi for j node to calculate rc
        rs += sin(phi[j]);                                                      //calculate sum all sin phi for j node to calculate rs
    }                                                                           //
    return sqrt(pow(rc, 2) + pow(rs, 2)) / (1.0 * N);                           //Order parameter sync
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                            runge_kutta4_integrator                             @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::runge_kutta4_integrator(double* y) {                             //with have IC(V h n s) we calculate runge kutta for all parameter
    double* k1 = new double[N];                                                 //created 1d pointer for k1
    double* k2 = new double[N];                                                 //created 1d pointer
    double* k3 = new double[N];                                                 //created 1d pointer
    double* k4 = new double[N];                                                 //created 1d pointer
    double* f = new double[N];                                                  //created 1d pointer for change to runge
    for (int i = 0; i < N; i++)                                                 //for loop from 0 to 4*200 V h n s
        f[i] = y[i];                                                            //each phi[i] save to f[i]
    dydt(k1, f);                                                                //calculate for all V h n s one step k1
    for (int i = 0; i < N; i++)                                                 //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + 0.5 * dt * k1[i];                                         //calculate gradient k1
    dydt(k2, f);                                                                //calculate for all V h n s one step k2
    for (int i = 0; i < N; i++)                                                 //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + 0.5 * dt * k2[i];                                         //calculate gradient k2
    dydt(k3, f);                                                                //calculate for all V h n s one step k3
    for (int i = 0; i < N; i++)                                                 //for loop from 0 to 4*200 V h n s
        f[i] = y[i] + dt * k3[i];                                               //calculate gradient k3
    dydt(k4, f);                                                                //calculate for all V h n s one step k4
    for (int i = 0; i < N; i++)                                                 //for loop from 0 to 4*200 V h n s
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;             //calculate gradient total
    delete f;                                                                   //remove pointer
    delete k1;                                                                  //remove pointer
    delete k2;                                                                  //remove pointer
    delete k3;                                                                  //remove pointer
    delete k4;                                                                  //remove pointer
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     dydt                                       @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::dydt(double* ff, const double* x0) {                             //calculate one step for parameter phi
    for (int i = 0; i < N; i++) {                                               //for loop to change each node
        double a = 0.0;                                                         //define initional variable
        for (int j = 0; j < N; j++) {                                           //for loop to change each node
           // a += (  (connection[i][j] * sin((x0[j] - x0[i]))) + ((1-connection[i][j]) * (1-cos((x0[j] - x0[i])))/2)  );//
            a += (connection[i][j] * sin((x0[j] - x0[i]))) ;                    //simple kuramato 
        }                                                                       //
        ff[i] = W[i] + (coupling / N) * a;                                      //calculate phi for i node
    }                                                                           //
}                                                                               //