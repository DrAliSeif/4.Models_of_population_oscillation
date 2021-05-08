/************************************************************************************************/
/*** Topic: Kuramoto model with Runge-Kutta 4th Order Method for 9 nodes in fully connection  ***/
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
    double coupling;                                                            //stringh cupling
    string address_matrix;                                                      //name matrix file
    double* phi = new double[N];                                                //create 1d pointer
    double* W = new double[N];                                                  //create 1d pointer
    vec2i connection;                                                           //connection in [row][column]
public:                                                                         //public values
                                                                                //
    Kuramoto                                                                    //
    (int N_, double dt_, double t_final_, double coupling_, string address_matrix_) ://input data in class
        N(N_), dt(dt_), tfinal(t_final_), coupling(coupling_), address_matrix(address_matrix_){}//change name input to privet
    void initial();                                                             //run initial variable for example read matrix connection and initial rand for w and phi
    void matrix();                                                              //read matrix connection
    void theta();                                                               //crate rand between 0 to 1 uniform_real_distribution for initial phi
    void Omega();                                                               //rand mean 1 and standard deviation 0.01 normal_distribution 
    void Run();                                                                 //run program
    void scale_2_pi(double*);                                                   //put variable phi in 0 to 2pi
    double sync(double*);                                                       //calculate syncrony
    void runge_kutta4_integrator(double*);                                      //calculate runge kutta4
    void dydt(double*, const double*);                                          //calculate dphi/dt
};                                                                              //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                    initial                                     @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::initial() {                                                      //call initial in wangbuzsaki class
    matrix();                                                                   //read connection matrix
    theta();                                                                    //create random between 0 to 2pi normal
    Omega();                                                                    //create random between 0 to 1 ufiform
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                          Read matrix connection                                @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::matrix() {                                                       //read matrix connection     
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
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     theta                                      @@@@
//@@@              crate rand between 0 to 1 uniform_real_distribution               @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::theta() {                                                        //calculate initial theta
    random_device rd;                                                           //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());                                                          //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<float> uniform_distribution_teta(0.0, 1.0);       //Uniform real distribution from 0 to 1
    for (int i = 0; i < N; i++)                                                 //for loop all node
        phi[i] = uniform_distribution_teta(gen) * 2 * Pi;                       //create phi random for each node
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                     Omega                                      @@@@
//@@@         rand mean 1 and standard deviation 0.01 normal_distribution            @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::Omega(){                                                         //calculate w
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();     //using seed time
    default_random_engine generator(seed);                                      //random generator
    normal_distribution<double> dist(0, 0.2);                                //creat random that This distribution produces random numbers around the distribution mean 1 with a specific standard deviation 0.01.
    for (int i = 0; i < N; i++)                                                 //for loop all node
        W[i] = dist(generator);                                                 //create w random for each node
}                                                                               //
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                     RUN                                        @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Kuramoto::Run() {                                                          //call Run function
    ofstream syncrony("syncrony.txt", ios::out | ios::trunc);                   //create file .txt for save sync volt
    for (double t = 0; t < tfinal; t += dt) {                                   //loop over time step
        runge_kutta4_integrator(phi);                                           //calculate runge_kutta for all node
        scale_2_pi(phi);                                                        //change values to be in range 0 to 2*Pi
        syncrony << t << "\t" << sync(phi) << endl;                             //print time and sync in .txt
    }                                                                           //
    syncrony.close();                                                           //close file .txt
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