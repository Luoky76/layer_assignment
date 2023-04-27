#ifndef PARAMETER_H
#define PARAMETER_H

#include "wirePropertyTable.h"
struct TABLE {
    std::vector<std::vector<double> > UC;
    std::vector<std::vector<double> > UR;
};
extern double via_C[9];
extern double via_R[9];
extern double Lambda;
extern double Miu;
extern double Beta;
extern double driver_R;
extern double load_C;
extern wirePropertyTable wpt;
extern int mainmode;
extern int weightmode;
extern int parallelmode;
extern int tmode;
extern TABLE RCtable[9];
extern int nodeweightingthreshold;
extern double nodeweightingparamter;
extern double historyparamter;
extern double parallelLevelthreshold;
extern double parallelpercentthreshold;
extern int fcfsthreshold;
//void tablecreating();
void readparameter(std::string ); 
#endif 
