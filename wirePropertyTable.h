#ifndef WIRE_PROPERTY_TABLE_H
#define WIRE_PROPERTY_TABLE_H
#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <assert.h>
#include <iostream>

/* Created by Rickard Ewetz at Purdue University in 2015*/
/* This file contains a wirePropertyTable class and loader class that is used to load the wirePropertyTable from file*/ 


class wirePropertyTable{
  friend class wpTableLoader;
 public:
  // layerIndex 2 => wire on M2
  // wireDensity 1 => wire density is 10%
  // fat=0 => single normal wire
  // fat=1 on M1 to M4 => is parallel wire
  // fat=1 on M5 to M9 => is a wide wire.
  
  //CAUTION: There is no fat wire for wireDensity=100 

  // the capacitance units are fF/nm and ohm/nm
  double getUnitLengthCapacitance(int layerIndex, int wireDensity, int useFatWire);
  double getUnitLengthResistance(int layerIndex, int wireDensity, int useFatWire);
 
 
  
  int numMasks(int layer);
  int wireWidth(int layer);
  int wireSpacing(int layer);




  // There two functions are under development
  double getUnitLengthCapacitance(int layerIndex, int wireDensity, int wireDensityUnder, int wireDensityOver, int useFatWire);
  double getUC(int,double,double,double,int);
  double getUnitLengthResistance(int layerIndex, int wireDensity, int wireDensityUnder, int wireDensityOver, int useFatWire);
  double getUR(int,double,double,double,int);


  
 private:
  wirePropertyTable();
  int index(int layer, int wireDensity, int wireDensityUnder, int wireDensityOver, int useFat);
  
  double getCapEntry(int);
  double getResEntry(int);
  
  int _d1Max;
  int _d2Max;
  int _d3Max;
  int _layerMax;
  int _fatMax;

 private:
  std::vector<double> _ctable;
  std::vector<double> _rtable;
  std::vector<int> _numMasks;
  std::vector<int> _wireSpacing;
  std::vector<int> _wireWidth;

};


class wpTableLoader{
 public:
  wpTableLoader();
  wirePropertyTable load(std::string tablePath);
  

 private:
  void load(std::ifstream & in);
  wirePropertyTable & wpt(){return _wpt;}

 private:
  wirePropertyTable _wpt; 
};


#endif
