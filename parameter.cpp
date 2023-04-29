#include "parameter.h"

double via_C[9] = {0.000026296 , 0.000026296 , 0.0000252959 , 0.0000255218 , 0.0000288807 , 0.0000264854 , 0.0000263049 , 0.0000219274 , 0.0000307698};
double via_R[9] = {0.00271234 , 0.00271234 , 0.00271234 , 0.00271234 , 0.00055652 , 0.00055652 , 0.00055652 ,  0.0000993014 , 0.0000993014};
double Lambda = 100;
double Miu = 1;
double Beta = 0.01;
double driver_R = 1;
double load_C = 0.001;
TABLE RCtable[9];
int mainmode = 0;
int weightmode =0;
int parallelmode =0;
wirePropertyTable wpt = wpTableLoader().load("data\\wireTable");
int tmode = 1;
int nodeweightingthreshold = 1;
double nodeweightingparamter = 1.0;
double historyparamter = 1.0;
double parallelLevelthreshold = 0.01;
double parallelpercentthreshold = 1.0;
int fcfsthreshold = 1;

void readparameter( std::string filename){
    std::ifstream infile(filename.c_str());

    if(!infile){
        std::cout << "can not open parameter file!\n";
        exit(1);
    }

    std::string stmp;
    int ntmp;
    double dtmp;
    while(!infile.eof()){
        infile >> stmp;
        if(stmp == "Lambda"){
            infile >> dtmp;
            Lambda = dtmp;
        }
        else if(stmp == "Beta"){
            infile >> dtmp;
            Beta = dtmp;
        }
        else if(stmp == "Miu"){
            infile >> dtmp;
            Miu = dtmp;
        }
        else if(stmp == "driver_R"){
            infile >> dtmp;
            driver_R = dtmp;
        }
        else if(stmp == "load_C"){
            infile >> dtmp;
            load_C = dtmp;
        }
        else if(stmp == "mainmode"){
            infile >> ntmp;
            mainmode = ntmp;
        }
        else if(stmp == "weightmode"){
            infile >> ntmp;
            weightmode = ntmp;
        }
        else if(stmp == "parallelmode"){
            infile >> ntmp;
            parallelmode = ntmp;
        }
        else if(stmp== "nodeweightingthreshold"){
            infile >> ntmp;
            nodeweightingthreshold = ntmp;
        }
        else if(stmp == "nodeweightingparamter"){
            infile >> dtmp;
            nodeweightingparamter = dtmp;
        }
        else if(stmp == "historyparamter"){
            infile >> dtmp;
            historyparamter = dtmp;
        }
        else if(stmp == "parallelLevelthreshold"){
            infile >> dtmp;
            parallelLevelthreshold = dtmp;
        }
        else if(stmp == "parallelpercentthreshold"){
            infile >> dtmp;
            parallelpercentthreshold = dtmp;
        }
        else if(stmp == "fcfsthreshold"){
            infile >> ntmp;
            fcfsthreshold = ntmp;
        }
    }
}
