# 4.Models of population oscillation
Models of population oscillation

__________________________________________________________________
--------------------------------------------------------------
# Kuramoto_RK4
## Kuramoto model with Runge-Kutta 4th Order Method for 9 nodes in fully connection network with cupling
<p align="center">
 <img src="https://github.com/aliseif321/4.Models_of_population_oscillation/blob/main/Kuramoto_RK4/Pic/Untitled.png?raw=true" >
 </p>
 
_______________________________________________________
```ruby
/************************************************************************************************/
/*** Topic: Kuramoto model with Runge-Kutta 4th Order Method for 9 nodes in fully connection  ***/
/***   network with cupling                                                                   ***/
/*** Version Release 17.12 rev 11256                                                Ali-Seif  ***/
/*** Date: 5/8/2021                                                                           ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>                                                             //import library for cout and endl and etc.
#include"Kuramoto.h"                                                            //import library Kuramoto
using namespace std;                                                            //using standard 
int main() {                                                                    //
    clock_t t = clock();                                                        //time start record
    const string address_matrix = "fully50";                                      //name matrix file
    const int size_matrix = 50;                                                //number of node size of matrix clom and row
    const double dt = 0.05;                                                     //step lenth
    const double t_final = 500;                                                 //time final
    const double coupling = 0.8;                                                //stringh cupling
    Kuramoto kR(size_matrix, dt, t_final, coupling, address_matrix);            //call class and input initial variable
    kR.initial();                                                               //run initial cundition
    kR.Run();                                                                   //run kuramato model
    cout << "\nTime taken by program is :\t" << ((double)clock() - t) / CLOCKS_PER_SEC << " sec\nfinish\n";
    return 0;                                                                   //run program was correct
}


``` 
__________________________________________________________________
--------------------------------------------------------------
