#include <iostream>
#include "string.h"
#include <vector>
#include <map>
#include <math.h>
#include <queue>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include <utility>
#include <algorithm>
#include "parameter.h"
//#include <sys/time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include "taskflow.hpp"
using namespace std;

class SORTEL{//sort element
public:
    double eldelay;
    int index;
};

class RECORD_DP{//the dynamic programming tree of a net in layer assignment
public:
    double cost[8], C[8], delay[8];
    int parallelwirenumber;// record the wire type used by the assignment (0 is default, 1 is parallel or NDR wire)
    int toplayerindex; //in layer assignment, record the top layer index which can use parallel or NDR wires.
    vector<vector<int> > child_combine;//all sub-trees of child nodes
    bool done;//true=visited  false=unvisit
};

class TERMINAL{//sink or source 
public:
    int grid_index_x;
    int grid_index_y;
    int layer;
    double elmoredelay;
    int RTree_index;//in which tree node (TREE_NODE) of a 2D routing net
};

class TREE_NODE {//routed net tree, (input = 2D,  output = 3D)
public:
    int x,y;
    int layer;
    vector<int> pinlayer;//If the node do not have pin(source or sink), the vector[0]= -1,  otherwise, vector[0~n] = sink pins . In root node , the vector[0] = source pin, other vector[i] = sink pins.
    int parent_index; 
    vector<int> child_index;
    double child_C;// downstream capacitance 
    int wirenumber;//= track number (1= default, 2=parallel wire, 3=NDR wire)
    double weight;//delay cost = weight*delay*Lambda
    double level;//= node level / max node level
    bool underlevelthreshould;//determine use parallel or NDR or not (true= can use, false=cannot use)

    TREE_NODE()
    {
        wirenumber = 0;
    }
};

class SEGMENT_POINT{
public:
    int layer;
    int x;
    int y;
    SEGMENT_POINT ()
    {
        x = y = layer = 0;
    }
};

class SEGMENT{//layer assignment result of a net
public:
    SEGMENT_POINT p1,p2;
    bool via;//true = via
    int VORH;// 0=via,1=V,2=H
    bool done;//using in create routing tree
    int wirenumber;//= track number (1= default, 2=parallel wire, 3=NDR wire)
    bool needripup;//true = need rip up(reduce capacity of 2d grid edge)
};

class EDGE_C{
public:
    int layernumber;
    vector<int> C_V;// capacity of 3D grid edge (Vertical)
    vector<int> C_H;//           ....           (Horiziontal)
    vector<double> historic;// historic cost 
    int total_c_v,total_c_h;// capacity of 2D grid edge
    int aveg_c_v,aveg_c_h;// average capacity of 2D grid edge
    map<int,map<int,int> > netid;//which net routed at this grid edge, 2D map, first index is layer, second index is netlist index (In order to quickly modify, "map(red-black tree struct)" is used.) 
};

class NET{
public:
    NET(){
        havesegment = false;//if the net has routing tree, value = true.
        needreassign = false;//using in ripup and reassign
        maxlevel=0;// the max level of nodes of 2D routed tree. the level of root node is 1.
        fatwirelevelthreshold=0.0;//the node level threshould
        parallel=-1;
    };
    void setname(string _name){
        name = _name;
    };
    string getname(){
        return name;
    };
    double eldelay;//elmore delay of the net
    bool needreassign;//true= need to ripup and reassignment
    double fatwirelevelthreshold;//the threshold of parallel an NDR, if < fatwirelevelthreshold, can use parallel and NDR. if >= fatwirelevelthreshold, cannot.
    double maxlevel;// max level of the 2D routed tree
    vector<TERMINAL> sinks_marge;// sinks vector (input is the marge file), vector[0] is source pin
    vector<SEGMENT> segments;// the input 3D segment of the net (need to modify to 2D grid edge model)
    vector<SEGMENT> twoDsegments;// the modified 2D segment
    vector<TREE_NODE> RTree;// the routed tree of the net.
    vector<SEGMENT> resultsegments;// the layer assignment result
    int parallel;
    bool havesegment;
    void createroutingtree(vector<SEGMENT>& ,int,vector<vector<EDGE_C> >& ,vector<int>& , vector<int>& ,vector<int>&,vector<int>& );
    bool findcrosspoint(SEGMENT , SEGMENT, SEGMENT_POINT &);
    void splitsegment(SEGMENT_POINT, vector<SEGMENT> & , int );
    void splitsegmentto2pinsegment(vector<SEGMENT>&,int);
    void modifyroutingtree();
    void ct(vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);//create timing information (optimize slack)
    double elmoredelay(int,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);//compute a segment elmore delay (using in timing model)
    double computechild_C(int,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);//comput a sgement downstream capacitance (using in timing mode)
    double elmoredelay_assign(TERMINAL&,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>& );//using in layer assignment
    double computechild_C_assign(int,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);//using in layer assignmnet
    double viaDelay(int,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);
    void ThreeDprojectto2D(int,int,int,int); 
    void updateRTreeweight();//node weighting (recursive function) 
    void passupdate(int,double);// using in updateRTreeweight
    void nodelevel(int);//compute the level of each node (percentage)
    void computelevelthreshold();// determine each node whether can use parallel and NDR wires or not
    void computeinitialRTreeweight(int);
    int viacount();
    inline double computedensity(TREE_NODE,TREE_NODE,int,vector<vector<EDGE_C> >&,vector<int>&,vector<int>&,vector<int>&,vector<int>&);
private:
    string name;
};


class SHAPE{
public:
    SHAPE(){
        x=0;y=0;w=0;h=0;
    };
    SHAPE(int _x,int _y,int _w,int _h){
        x=_x;y=_y;w=_w;h=_h;
    };
    int x,y;
    int w,h;
};

class NODE{
public:
    NODE() {
        name = "";
        type = 0;
        width = 0;
        height = 0;
        nonrectangular = false;
        pinlayer = 1;
    };

    NODE(string _name, int _type, int w, int h) {
        name = _name;
        type = _type;
        width = w;
        height = h;
        nonrectangular = true;
        pinlayer =1;
    };

    void setname(string _name) {
        name = _name;
    };
    string getname() {
        return name;
    };
    void settype(int _type) {
        type = 0;
    };
    int gettype() {
        return type;
    };
    void setpos_x(double x) {
        pos_x = x;
    };
    int getpos_x() {
        return pos_x;
    };
    void setpos_y(double y) {
        pos_y = y;
    };
    int getpos_y() {
        return pos_y;
    };
    void setpinlayer(int _layer){
        pinlayer = _layer;
    };
    int getpinlayer(){
        return pinlayer;
    };
    void setblocklayer(int _layer){
        blocklayer.push_back(_layer);
    };
    vector<int> getblocklayer(){
        return blocklayer;
    };
    void setnonrect(){
        nonrectangular=true;
    };
    bool getnonrect(){
        return nonrectangular;
    };
    void setshape(vector<SHAPE> _shape){
        shapes = _shape;
    };
    vector<SHAPE> getshape(){
        return shapes;
    };
    int getw() {
        return width;
    };
    int geth() {
        return height;
    };
private:
    string name;
    int type; //0: normal node  1: terminal node  2: terminal_N node
    bool nonrectangular; // blockage is nonrectangular
    vector<SHAPE> shapes;
    int width;
    int height;
    int pos_x;
    int pos_y;
    int pinlayer;
    vector<int> blocklayer;
};

class GLOBALCELL{
public:
    int pos_x,pos_y;
};


class CIRCUIT {
public:
    CIRCUIT(std::vector<string> filenames) {
        ablenet=0;
        parser_route(filenames[1]);
        cout << "done route\n";
        creatroutegraph();
        parser_net(filenames[1]);
        cout << "done net\n";
        parser_nctugr(filenames[2]);
        cout << "done nctugr\n";
        readmarge(filenames[3]);
        cout << "done marge\n";
        //int mode = atoi(argv[5]);//mode 0: read nctugr global routing result and create timing data(require time)
                                 //mode 1: read nctugr global routing result and compress to 2d routing
        
        //compress 3D to 2D and create routed tree
        for(int i=0;i<nets.size();i++){
            if(!nets[i].havesegment)
               continue;
            if(mainmode==0){
                nets[i].splitsegmentto2pinsegment(nets[i].segments,0);
                nets[i].createroutingtree(nets[i].segments,0,edge_c,C_V,C_H,minwidth,minspacing);//for 3D routing tree
            }
            else if(mainmode==1){
                nets[i].ThreeDprojectto2D(grid_x_number,grid_y_number,gridsize_x,gridsize_y);
                nets[i].splitsegmentto2pinsegment(nets[i].twoDsegments,1);
                nets[i].createroutingtree(nets[i].twoDsegments,1,edge_c,C_V,C_H,minwidth,minspacing);
            }

        }
        cout << "done routing tree creating\n";
        
        createlookuptable();//the coupling effect look up table

		for (int i = 0; i < layernumber; i++){
			maxcap[i] = 1.0*(C_V[i] + C_H[i]) / (minwidth[i] + minspacing[i]);
		}
        if(mainmode==0){
            for(int i=0;i<nets.size();i++){
                if(!nets[i].havesegment)
                    continue;
                computegecllcapacityfornctugr(nets[i].RTree);
            }
            for(int i=0;i<nets.size();i++){
                if(!nets[i].havesegment)
                    continue;
                nets[i].ct(edge_c,C_V,C_H,minwidth,minspacing);
            }
            cout << "start timing data creating\n";
            createtimingdata(filenames[1]);
            cout << "done timing data creating\n";
        }
        else if(mainmode==1){

            computeaveragedensity();
			
			cout << "start initial dp layer assignment\n";
            //double tmp = Miu;
            //Miu=0;//variable of congestion cost
			Lambda = 10; Beta = 1; Miu = 1; parallelmode = 0;//注意：Miu=1
            initialLA();
			//Miu = tmp; 
						
			cout << "start RRA\n";
			Lambda = 10; Beta = 1; Miu = 1; parallelmode = 0;
            RRA_main();//NDLA
			computeparallelnumber();
            emptyedge_c();
            computecritical();
            computeviacount();
            
			
			
			cout << "start refinemnet\n";
			Lambda = 10; Beta = 1; Miu = 1; parallelmode = 1;
            greedyrefinement();
			
			cout << "start RRA2\n";
			Lambda = 10; Beta = 3.5; Miu = 1; parallelmode = 1;
			RRA_main();//NDLA
			computeparallelnumber();
			emptyedge_c();
			computecritical();
			computeviacount();



			Lambda = 10; Beta = 1.5; Miu = 1; parallelmode = 1;
			postop();
            computeparallelnumber();
            emptyedge_c();
            computecritical();
            computeviacount();
            outputroutingresult(filenames[4]);
            cout << "output routing result: "<<filenames[4]<<endl;
        }

    };
    void initialLA();
    void parser_node(string filename);
    void parser_pl(string filename);
    void parser_route(string filename);
    void parser_shape(string filename);
    void parser_net(string filename);
    void parser_nctugr(string filename);
    void readmarge(string filename);
    void createlookuptable();
    void creatroutegraph();
    void computegecllcapacityfornctugr(vector<TREE_NODE>&);
    void computeaveragedensity();
    void createtimingdata(string filename);
    void readtimingdata(string filename);
    void outputroutingresult(string filename);
    void computinggcellcapacity();
    void dynamic_program_main(NET& ,int,int,int,double);
    void dynamic_program(vector<TREE_NODE>& ,int , vector<RECORD_DP>&, vector<SEGMENT>& ,int,int,int,double);
    //void findallcombine(vector<int>& ,int ,int , TREE_NODE& rt, vector<RECORD_DP>& , int , vector<int>& );
    void findallcombinefor( int ,vector<RECORD_DP>& , vector<TREE_NODE>& ,int,int,int);
    void createsolution(vector<int>& ,int ,vector<RECORD_DP>& ,vector<SEGMENT>& , vector<TREE_NODE>& ,int);//topdown, after root assigned
    vector<SORTEL> checkillegalnet(int);
    void ripup(NET&,int);
    void RRA_main();
    void computedelay();
    void greedyrefinement();
    //void speedup(vector<RECORD_DP>&,vector<TREE_NODE>&,int);
    double computedensity(TREE_NODE,TREE_NODE,int,int);
    void singlenetdelay(int);
    void computeoverflow();
    void recoverpathtoedge(NET&,int);
    void computeviacount();
    void computeparallelnumber();
    void emptyedge_c();
    void computecritical();
    void postop();
    map<string,int> netsid2index;
    vector<NET> nets;

    //新增函数
    //建立依赖
    void buildDependency(vector<TREE_NODE>& _rt, vector<RECORD_DP>& record_dp, 
    vector<SEGMENT>& resultsegments,int netindex,int mode,int greedy, double q, 
    tf::Taskflow &taskflow);
private:
    map<string,int> nodesid2index;
    vector<NODE> nodes;
    int grid_x_number,grid_y_number;  
    int layernumber;
    vector<int> C_V;
    vector<int> C_H;
    vector<int> minwidth;
    vector<int> minspacing;
    int gridsize_x;
    int gridsize_y;    
    int grid_oring_x;
    int grid_oring_y;
    vector<vector<GLOBALCELL> > gcell;
    vector<vector<EDGE_C> > edge_c;
	vector<vector<EDGE_C> > edge_cap;
	double maxcap[9];
	int ablenet;
    double ableparadelay;
};
