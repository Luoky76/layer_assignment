#include "wirePropertyTable.h"
#include <math.h>

inline int rounding2(double a){ 
    return (int) (a+0.5);
}


wirePropertyTable::wirePropertyTable()
  : _d1Max(10), _d2Max(10), _d3Max(10), _layerMax(10), _fatMax(2), _ctable(10*10*10*10*2*2,-1), _rtable(10*10*10*10*2*2,-1)
{
  _wireSpacing.push_back(-1); //poly
  _wireSpacing.push_back(32); // M1
  _wireSpacing.push_back(32); // M2
  _wireSpacing.push_back(32); // M3
  _wireSpacing.push_back(32); // M4
  _wireSpacing.push_back(64); // M5
  _wireSpacing.push_back(64); // M6
  _wireSpacing.push_back(64); // M7
  _wireSpacing.push_back(128); // M8
  _wireSpacing.push_back(128); // M9

  _wireWidth.push_back(-1); //poly
  _wireWidth.push_back(32); // M1
  _wireWidth.push_back(32); // M2
  _wireWidth.push_back(32); // M3
  _wireWidth.push_back(32); // M4
  _wireWidth.push_back(64); // M5
  _wireWidth.push_back(64); // M6
  _wireWidth.push_back(64); // M7
  _wireWidth.push_back(128); // M8
  _wireWidth.push_back(128); // M9

  _numMasks.push_back(-1); //poly
  _numMasks.push_back(2); // M1
  _numMasks.push_back(2); // M2
  _numMasks.push_back(2); // M3
  _numMasks.push_back(2); // M4
  _numMasks.push_back(1); // M5
  _numMasks.push_back(1); // M6
  _numMasks.push_back(1); // M7
  _numMasks.push_back(1); // M8
  _numMasks.push_back(1); // M9
}


double wirePropertyTable::getUnitLengthCapacitance(int layerIndex, int wireDensity, int useFatWire){
  return getUnitLengthCapacitance(layerIndex,wireDensity,wireDensity,wireDensity,useFatWire);
}
 
 double wirePropertyTable::getUnitLengthCapacitance(int layerIndex, int wireDensity, int wireDensityUnder, int wireDensityOver, int useFatWire){
    if(!(1 <= wireDensity && wireDensity <= 10))
        std::cout << wireDensity << std::endl;
  assert(2 <= layerIndex  && layerIndex <= 9);
  assert(1 <= wireDensity && wireDensity <= 10);
  assert(useFatWire == 0 || useFatWire == 1);
  int ind = index(layerIndex,wireDensity,wireDensityUnder,wireDensityOver,useFatWire);
  double unitCap = getCapEntry(ind);
  return unitCap;
}

double wirePropertyTable::getUC(int layer, double wireD,double wireDU,double wireDO,int fat){
    /*if(  rounding2(wireD) %10==0){
        int d = rounding2(wireD);
        int ch = d/10;
        if((ch != 10 && fat!=1)&&ch!=0)
            return getUnitLengthCapacitance(layer,d,d,d,fat);
     }*/
    int a = floor(wireD);
    int b = ceil(wireD);
    if(a==b){
        return getUnitLengthCapacitance(layer,a,a,a,fat);
    }

    if(wireD <1){
        a++;
        b++;
        if(rounding2(wireD)==0)
            b++;
    }
    /*else if(wireD >= 9&& fat==1){
        a--;
        b--;
        if(rounding2(wireD)==10)
            a--;
    }*/


    double x = getUnitLengthCapacitance(layer,a,a,a,fat);
    double y = getUnitLengthCapacitance(layer,b,b,b,fat);

    return x + double( double( wireD - a)/(b-a) * (y-x));
}

int wirePropertyTable::numMasks(int layer){return _numMasks[layer];}
int wirePropertyTable::wireWidth(int layer){return _wireWidth[layer];}
int wirePropertyTable::wireSpacing(int layer){return _wireSpacing[layer];}


double wirePropertyTable::getUnitLengthResistance(int layerIndex, int wireDensity, int useFatWire){
  return getUnitLengthResistance(layerIndex,wireDensity,wireDensity,wireDensity,useFatWire);
}
 
 double wirePropertyTable::getUnitLengthResistance(int layerIndex, int wireDensity, int wireDensityUnder, int wireDensityOver, int useFatWire){
  assert(2 <= layerIndex  && layerIndex <= 9);
  assert(1 <= wireDensity && wireDensity <= 10);
  assert(useFatWire == 0 || useFatWire == 1);
  int ind = index(layerIndex,wireDensity,wireDensityUnder,wireDensityOver,useFatWire);
  double unitRes = getResEntry(ind);
  return unitRes;
}

double wirePropertyTable::getUR(int layer, double wireD,double wireDU,double wireDO,int fat){
    /*if(  rounding2(wireD) %10==0){
        int d = rounding2(wireD);
        int ch = d/10;
        if((ch != 10 && fat!=1)&&ch!=0)
            return getUnitLengthResistance(layer,d,d,d,fat);
    }*/

    int a = floor(wireD);
    int b = ceil(wireD);
    if(a==b){
      //  if(!(a==10&&fat==1))
            return getUnitLengthResistance(layer,a,a,a,fat);
    }

    if(wireD <1){
        a++;
        b++;
        if(rounding2(wireD)==0)
            b++;
    }
    /*else if(wireD >= 9&& fat==1){
        a--;
        b--;
        if(rounding2(wireD)==10)
            a--;
    }*/
    double x = getUnitLengthResistance(layer,a,a,a,fat);
    double y = getUnitLengthResistance(layer,b,b,b,fat);
    
    return x + double( double( wireD - a)/(b-a) * (y-x));
}



int wirePropertyTable::index(int lay, int d1, int d2, int d3, int fat){
  return fat + _fatMax*lay + d1*_fatMax*_layerMax + d2*_fatMax*_layerMax*_d1Max + d3*_fatMax*_layerMax*_d1Max*_d2Max; 
}

double wirePropertyTable::getCapEntry(int idx){
  if(!(idx <_ctable.size()) || !(_ctable[idx] != -1)){
    std::cout << "Error: Requested cap entry out of bounds or an non initialized entry" << std::endl;
  }

  assert(idx <_ctable.size());
  assert(_ctable[idx] != -1);
  return _ctable[idx];
}

double wirePropertyTable::getResEntry(int idx){
  if(!(idx <_rtable.size()) || !(_rtable[idx] != -1)){
    std::cout << "Error: Requested res entry out of bounds or an non initialized entry" << std::endl;
  }
  assert(idx <_rtable.size());
  assert(_rtable[idx] != -1);
  return _rtable[idx];
}


wpTableLoader::wpTableLoader(){}


wirePropertyTable wpTableLoader::load(std::string path){
  //wirePropertyTable(){}
  std::ifstream in;
  in.open(path.data(),std::ifstream::in);
  assert(in.is_open()); 
  load(in);
  in.close();
  return _wpt;
}


void wpTableLoader::load(std::ifstream & in){
  std::string word;
  in >> word;
  std::cout << word << std::endl;
  assert(word.compare("wpTable") == 0);
  int numEntries = -1;
  in >> word; in >> numEntries;
  assert(word.compare("numEntries") == 0);
  assert(numEntries > 0);
  int layer; int fat; int d1; int d2; int d3; double unitCap; double unitRes;
  
   std::string dummy1; std::string dummy2; std::string dummy3; std::string dummy4;
   std::string dummy5; std::string dummy6; std::string dummy7; std::string dummy8;
//      lay     fat       d1        d2        d3
  in >> dummy1 >> dummy2 >> dummy3 >> dummy4 >> dummy5 >> dummy6 >> dummy7;
 
  // std::cout << dummy1 << dummy2 << dummy3 << dummy4 << dummy5 << dummy6 << dummy7 << std::endl; 
  
  for(int i = 0; i < numEntries;i++){
    in >> layer >> fat >> d1 >> d2 >> d3 >> unitCap >> unitRes;
    //std::cout << "Layer " << layer << " fat " << fat << " d1 " << d1 << " d2 " << d2 << " d3 " << d3 << " unitCap " << unitCap << " unitRes " << unitRes << std::endl; 
    assert(0 <= fat);
    assert(fat <= 1);
    assert(2 <= layer);
    assert(layer <= 9);
    assert(1 <= d1/10);
    assert(d1/10 <= 10);
    assert(1 <= d2/10);
    assert(d2/10 <= 10);
    assert(1 <= d3/10);
    assert(d3/10 <= 10);
    
    int idx = wpt().index(layer,d2/10,d1/10,d3/10,fat);
    //std::cout << "idx " << idx << " " << wpt()._ctable.size() << " " <<  wpt()._rtable.size() << std::endl;
    
    wpt()._ctable[idx] = unitCap;
    wpt()._rtable[idx] = unitRes;
  }
  
}
