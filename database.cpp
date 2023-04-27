#include "database.h"
#include "parameter.h"
int rounding(double);

void CIRCUIT::parser_node(string _filename){

    string _filename_nodes = _filename + ".nodes";
    ifstream infile(_filename_nodes.c_str());
    if(!infile){
        cout << "cant open file: nodes!\n";
        exit(1);
    }
    else{
    int nodenumber=0;
    int terminalnumber=0;
    while(!infile.eof()){
        string token;
        infile >> token;
        
        if(token == "NumNodes") {
            infile >>token >> nodenumber;
        }
        else if(token == "NumTerminals") {
            infile >> token >> terminalnumber;
            nodenumber = nodenumber - terminalnumber;
            for(int i=0;i<nodenumber;i++){
                string stmp;
                int xtmp,ytmp;
                infile >> stmp >> xtmp >> ytmp;
                NODE nodetmp(stmp,0,xtmp,ytmp);
                nodes.push_back(nodetmp);
                nodesid2index[stmp] = nodes.size()-1;
            }
            
            for(int i=0;i<terminalnumber;i++){
                string stmp, stmp2;
                int typetmp;
                int xtmp,ytmp;
                infile >> stmp >> xtmp >> ytmp >> stmp2;
                if(stmp2 == "terminal")
                    typetmp=1;
                else
                    typetmp=2;
                NODE nodetmp(stmp,typetmp,xtmp,ytmp);
                nodes.push_back(nodetmp);
                nodesid2index[stmp] = nodes.size()-1;
            }
        }
    }
    infile.close();
    }
}

void CIRCUIT::parser_pl(string _filename) {

    string _filename_pl = _filename + ".pl";
    ifstream infile(_filename_pl.c_str());
    if(!infile){
        cout << "cant open file: pl!\n";
        exit(1);
    }
    else{
    while(!infile.eof()){
        string token;
        infile >> token >> token >> token;
        for(int i=0;i<nodes.size();i++){
            string stmp;
            infile >> stmp;
            string name = stmp;
            int nodeindex = nodesid2index[stmp];
            int x,y;
            infile >> x >> y >>stmp >> stmp;
            nodes[nodeindex].setpos_x(x);
            nodes[nodeindex].setpos_y(y);
            int type = nodes[nodeindex].gettype();
            if(type!=0)
                infile >> stmp;
        }
        break;
    }
    infile.close();
    }
    
}

void CIRCUIT::parser_route(string _filename) {
    string _filename_route = _filename + ".route";
    ifstream infile(_filename_route.c_str());
    if(!infile){
        cout << "cant open file : route!\n";
        exit(1);
    }
    else{
        string token;
        infile >> token;
    while(!infile.eof()){
        if(token=="Grid"){
            infile >>token >>grid_x_number  >> grid_y_number >> layernumber;
        }
        else if(token=="VerticalCapacity"){
            infile >> token;
            for(int i=0;i<layernumber;i++){
                int tmp;
                infile >> tmp;
                C_V.push_back(tmp);
            }
        }
        else if(token=="HorizontalCapacity"){
            infile >> token;
            for(int i=0;i<layernumber;i++){
                 int tmp;
                 infile >> tmp;
                 C_H.push_back(tmp);
            }
        }
        else if(token=="MinWireWidth"){
            infile >> token;
            for(int i=0;i<layernumber;i++){
                int tmp;
                infile >> tmp;
                minwidth.push_back(tmp);
            }
        }
        else if(token=="MinWireSpacing"){
            infile >> token;
            for(int i=0;i<layernumber;i++){
                int tmp;
                infile >> tmp;
                minspacing.push_back(tmp);
            }
        }
        else if(token=="GridOrigin"){
            infile >>token >> grid_oring_x >> grid_oring_y;
        }
        else if(token=="TileSize"){
            infile >>token >>gridsize_x >> gridsize_y;
        }
        else if(token=="NumNiTerminals"){
            int ntmp;
            infile >> token >> ntmp;
            for(int i=0;i<ntmp;i++){
                string stmp;
                int ntmp2;
                infile >> stmp >> ntmp2;
            }
        }
        else if(token=="NumBlockageNodes"){
            int ntmp;
            infile >> token >> ntmp;
            for(int i=0;i<ntmp;i++){
                string stmp;
                int ntmp2;
                infile >> stmp >> ntmp2;
                for(int j=0;j<ntmp2;j++){
                    int ntmp3;
                    infile >> ntmp3;
                }
            }
        }
        infile >>token;
    }
    infile.close();
    }
    
}




void CIRCUIT::parser_shape(string _filename) {
    string _filename_pl = _filename + ".pl";
    ifstream infile(_filename_pl.c_str());
    if(!infile){
        cout << "cant open file: shape!\n";
        exit(1);
    }
    else{
        string token;
        infile >> token;
        while(!infile.eof()){
            if(token=="NumNonRectangularNodes"){
                int nodenumber;
                infile >> token >> nodenumber;
                for(int i=0;i<nodenumber;i++){
                    string stmp;
                    int shapenumber;
                    infile >> stmp >>token >> shapenumber;
                    int index = nodesid2index[stmp];
                    vector<SHAPE> vectortmp;
                    for(int j=0;j<shapenumber;j++){
                        int n1,n2,n3,n4;
                        infile >> token >> n1 >> n2 >> n3 >> n4;
                        SHAPE shtmp(n1,n2,n3,n4);
                        vectortmp.push_back(shtmp);
                    }
                    nodes[index].setnonrect();
                    nodes[index].setshape(vectortmp);
                }
            }
            infile >> token;
        }
    }
}

void CIRCUIT::creatroutegraph(){
    for(int i=0;i<grid_x_number;i++){
        vector<GLOBALCELL> vectorGtmp;
        vector<EDGE_C> vectorEtmp;
        for(int j=0;j<grid_y_number;j++){
            GLOBALCELL Gtmp;
            EDGE_C Etmp;
            Gtmp.pos_x = i*gridsize_x + grid_oring_x;
            Gtmp.pos_y = j*gridsize_y + grid_oring_y;
            vectorGtmp.push_back(Gtmp);
            Etmp.layernumber = layernumber;
            for(int k=0;k<layernumber;k++){
                int V = C_V[k];
                int H = C_H[k];
                Etmp.C_V.push_back(V);
                Etmp.C_H.push_back(H);
                Etmp.historic.push_back(1.0);
            }
            vectorEtmp.push_back(Etmp);
        }
        gcell.push_back(vectorGtmp);
        edge_c.push_back(vectorEtmp);
		edge_cap.push_back(vectorEtmp);
    }
};

void CIRCUIT::parser_net(string _filename){
    string _filename_pl = _filename + ".nets";
    ifstream infile(_filename_pl.c_str());
    if(!infile){
        cout << "cant open file: net!\n";
        exit(1);
    }
    else{
        string token;
        infile >> token;
        while(!infile.eof()){
            if(token == "NumNets"){
                int netnumber;
                infile >>token >> netnumber >> token>>token>>token;
                for(int i=0;i<netnumber;i++){
                    NET nettmp;
                    int degree;
                    string netname;
                    infile >> token >> token >> degree >> netname;
                    nettmp.setname(netname);
                    for(int j=0;j<degree;j++){
                        string stmp;
                        infile >> stmp;
                        infile >> stmp;
                        infile >> stmp;
                        double xtmp,ytmp;
                        infile >> xtmp >> ytmp;
                    }
                    nets.push_back(nettmp);
                    netsid2index[netname] = i;
                }
            }
            infile >> token;
        }
        infile.close();
    }

}


void CIRCUIT::parser_nctugr(string _filename){


    ifstream infile( _filename.c_str());
    if(!infile){
        cout << "cant open file: net!\n";
        exit(1);
    }
    else{
        string netname;
        int ntmp;
        char ctmp[100];
        infile >> netname >>ntmp>>ctmp;
        int index = netsid2index[netname];
        nets[index].havesegment = true;
        while(!infile.eof()){
            if(strcmp(ctmp,"!")){
                char *pch = NULL;
                pch = strtok(ctmp,"(-,)");
                int arraytmp[6];
                arraytmp[0] = atoi(pch);
                for(int i=0;i<5;i++){
                    pch = strtok(NULL,"(-,)");
                    arraytmp[i+1] = atoi(pch);
                }
                SEGMENT segtmp;
                segtmp.p1.x = arraytmp[0];
                segtmp.p1.y = arraytmp[1];
                segtmp.p1.layer = arraytmp[2];
                segtmp.p2.x = arraytmp[3];
                segtmp.p2.y = arraytmp[4];
                segtmp.p2.layer = arraytmp[5];
                if(segtmp.p2.layer != segtmp.p1.layer)
                    segtmp.via = true;
                else
                    segtmp.via = false;
                if(segtmp.via)
                    segtmp.VORH=0;
                else if(segtmp.p1.x == segtmp.p2.x)
                    segtmp.VORH=1;
                else
                    segtmp.VORH=2;
                segtmp.done = false;//for create routing tree, it mean that the segment is finded.
                nets[index].segments.push_back(segtmp);
            }
            else{
                infile >> netname;
                infile >> ntmp;
                index = netsid2index[netname];
                nets[index].havesegment = true;
            }
            infile >> ctmp;
        }
    }
    infile.close();
};

void CIRCUIT::readmarge(string _filename) {
    ifstream infile( _filename.c_str());
    if(!infile){
        cout << "cant open file: net!\n";
        exit(1);
    }
    else{
        string stmp;
        infile >> stmp;
        while(!infile.eof()){
            if(stmp == "num"){
                int netnum;
                infile >> stmp >> netnum;
                for(int i=0;i<netnum;i++){
                    string netname;
                    int netpinnum;
                    int ntmp2;
                    infile >> netname >> ntmp2 >> netpinnum >> ntmp2;
                    int netindex = netsid2index[netname];
                    for(int j=0;j<netpinnum;j++){
                        TERMINAL TERtmp;
                        int tmpx,tmpy;
                        infile >> tmpx >> tmpy >> TERtmp.layer;
                        TERtmp.layer--;
                        TERtmp.grid_index_x = tmpx/gridsize_x;
                        TERtmp.grid_index_x *= gridsize_x;
                        TERtmp.grid_index_y = tmpy/gridsize_y;
                        TERtmp.grid_index_y *= gridsize_y;
                        TERtmp.elmoredelay = 0.0;
                        bool repeat = false;
                        for(int k=0;k<nets[netindex].sinks_marge.size();k++){
                            if(nets[netindex].sinks_marge[k].grid_index_x == TERtmp.grid_index_x && nets[netindex].sinks_marge[k].grid_index_y == TERtmp.grid_index_y && nets[netindex].sinks_marge[k].layer == TERtmp.layer ){
                                repeat = true;
                                break;
                            }
                        }
                        if(!repeat)
                            nets[netindex].sinks_marge.push_back(TERtmp);
                    }

                }

                //read routing resource
                int tmpnum;
                infile >> tmpnum;
                for(int i=0;i<tmpnum;i++){
                    int x1,x2,y1,y2,layer,capacity;
                    infile >> x1 >> y1 >> layer >> x2 >> y2 >> layer >> capacity;
                    if(x1 == x2){
                        edge_c[x1][y1].C_V[layer-1] = capacity;
						edge_cap[x1][y1].C_V[layer - 1] = capacity;
                    }
                    else if(y1==y2){
                        edge_c[x1][y1].C_H[layer-1] = capacity;
						edge_cap[x1][y1].C_H[layer - 1] = capacity;

                    }
                }
            }
            else{
                infile >> stmp;
            }
        }
    }
    infile.close();
    
    for(int i=0;i<grid_x_number;i++){
        for(int j=0;j<grid_y_number;j++){
            edge_c[i][j].total_c_v = 0;
            edge_c[i][j].total_c_h = 0;
            for(int k=0;k<layernumber;k++){
                edge_c[i][j].C_V[k] = edge_c[i][j].C_V[k] / (minwidth[k]+minspacing[k]);
				edge_cap[i][j].C_V[k] = edge_cap[i][j].C_V[k] / (minwidth[k] + minspacing[k]);
                edge_c[i][j].C_H[k] = edge_c[i][j].C_H[k] / (minwidth[k]+minspacing[k]);
				edge_cap[i][j].C_H[k] = edge_cap[i][j].C_H[k] / (minwidth[k] + minspacing[k]);
                edge_c[i][j].total_c_v += edge_c[i][j].C_V[k];
                edge_c[i][j].total_c_h += edge_c[i][j].C_H[k];
            }
        }
    }

}

void NET::splitsegmentto2pinsegment(vector<SEGMENT> &_segments, int mode){


    int m,n;
    m=1000;n=1000;
    int **data;
    //if segments size > 1000, the crosspint is too much, so that the runtime will be to long. Use map to record crosspoint which was used
    if(_segments.size()>1000){
        data = new int*[m];
        for(int i=0;i<m;i++){
            data[i] = new int[n];
        }
    

        for(int i=0;i<m;i++)
            for(int j=0;j<n;j++){
                data[i][j] =0;
            }
    }

    //find all cross point
    vector<SEGMENT_POINT> crosspoint;
    for(int i=0;i<_segments.size()-1;i++){
        for(int j=i+1;j<_segments.size();j++){
            SEGMENT_POINT segtmp;
            bool havecross;
            havecross = false;
            if(_segments[i].VORH == _segments[j].VORH){//same direction
                if(_segments[i].VORH ==1 && _segments[i].p1.x == _segments[j].p1.x && _segments[i].p1.layer == _segments[j].p1.layer){//Vertical 
                    if(_segments[i].p1.y == _segments[j].p2.y){
                        segtmp.x = _segments[i].p1.x;
                        segtmp.layer = _segments[i].p1.layer;
                        segtmp.y = _segments[i].p1.y;
                        havecross = true;
                    }
                    else if(_segments[i].p2.y == _segments[j].p1.y){
                        segtmp.x = _segments[i].p2.x;
                        segtmp.layer = _segments[i].p2.layer;
                        segtmp.y = _segments[i].p2.y;
                        havecross = true;
                    }
                } 
                else if(_segments[i].VORH ==2 && _segments[i].p1.y == _segments[j].p1.y && _segments[i].p1.layer == _segments[j].p1.layer){//Horizontal
                    if(_segments[i].p1.x == _segments[j].p2.x){
                        segtmp.x = _segments[i].p1.x;
                        segtmp.layer = _segments[i].p1.layer;
                        segtmp.y = _segments[i].p1.y;
                        havecross = true;
                    }
                    else if(_segments[i].p2.x == _segments[j].p1.x){
                        segtmp.x = _segments[i].p2.x;
                        segtmp.layer = _segments[i].p2.layer;
                        segtmp.y = _segments[i].p2.y;
                        havecross = true;
                    }
                }
            }
            else if(_segments[i].VORH != 0 && _segments[j].VORH != 0 && _segments[i].p1.layer == _segments[j].p1.layer){// direction = V or H
                int cross_x,cross_y;
                if( _segments[i].VORH == 1){// i = V, j = H
                    cross_x = _segments[i].p1.x;
                    cross_y = _segments[j].p1.y;
                    havecross = true;
                }
                else if( _segments[i].VORH ==2){
                    cross_x = _segments[i].p1.y;
                    cross_y = _segments[j].p1.x;
                    havecross = true;
                }

                if(_segments[i].p1.x <= cross_x && cross_x <= _segments[i].p2.x && _segments[i].p1.y <= cross_y && cross_y <= _segments[i].p2.y){
                    if(_segments[j].p1.x <= cross_x && cross_x <= _segments[j].p2.x && _segments[j].p1.y <= cross_y && cross_y <= _segments[j].p2.y){
                        segtmp.x = cross_x;
                        segtmp.y = cross_y;
                        segtmp.layer = _segments[i].p1.layer;
                        havecross = true;
                    }
                }
            }
            else if(_segments[i].VORH == 0 && _segments[j].VORH != 0){// i = via, j = V or H
                int cross_x,cross_y;
                cross_x = _segments[i].p1.x;
                cross_y = _segments[i].p1.y;
                if(_segments[j].p1.x <= cross_x && cross_x <= _segments[j].p2.x && _segments[j].p1.y <= cross_y && cross_y <= _segments[j].p2.y){
                    if(_segments[i].p1.layer <= _segments[j].p1.layer && _segments[j].p1.layer <= _segments[i].p2.layer){
                        segtmp.x = cross_x;
                        segtmp.y = cross_y;
                        segtmp.layer = _segments[j].p1.layer;
                        havecross = true;
                    }
                }
            }
            else if(_segments[i].VORH != 0 && _segments[j].VORH == 0){// i = V or H, j = via
                int cross_x,cross_y;
                cross_x = _segments[j].p1.x;
                cross_y = _segments[j].p1.y;
                if(_segments[i].p1.x <= cross_x && cross_x <= _segments[i].p2.x && _segments[i].p1.y <= cross_y && cross_y <= _segments[i].p2.y){
                    if(_segments[j].p1.layer <= _segments[i].p1.layer && _segments[i].p1.layer <= _segments[j].p2.layer){
                        segtmp.x = cross_x;
                        segtmp.y = cross_y;
                        segtmp.layer = _segments[i].p1.layer;
                        havecross = true;
                    }
                }
            }
            if(havecross){
                if(_segments.size()>1000){
                    int level = segtmp.layer;
                    int tmp = data[segtmp.x/32][segtmp.y/40]%rounding(pow(10,level));
                    if(tmp < rounding(pow(10,level-1))){
                        crosspoint.push_back(segtmp);
                        data[segtmp.x/32][segtmp.y/40] += rounding(pow(10,level-1));

                    }
                }
                else{
                    crosspoint.push_back(segtmp);
                }

            }
            
        }
    }

    if(_segments.size()>1000){
        for(int i = 0; i < m; i++){
            delete [] data[i];
        }
        delete [] data;
    }

    for(int i=0;i<sinks_marge.size();i++){
        SEGMENT_POINT Stmp;
        Stmp.x = sinks_marge[i].grid_index_x;
        Stmp.y = sinks_marge[i].grid_index_y;
        if(mode ==0)
            Stmp.layer = sinks_marge[i].layer+1;
        else{
            Stmp.layer = 1;
        }
        crosspoint.push_back(Stmp);
    }
    for(int i=0;i<crosspoint.size();i++){
        for(int j=0;j<_segments.size();j++){
            if(_segments[j].VORH==0){
                if(crosspoint[i].x == _segments[j].p1.x && crosspoint[i].y == _segments[j].p1.y){
                    if(_segments[j].p1.layer < crosspoint[i].layer && crosspoint[i].layer  < _segments[j].p2.layer){
                        splitsegment(crosspoint[i], _segments, j);
                    }
                }
            }
            else if(_segments[j].VORH == 1){
                if(crosspoint[i].x == _segments[j].p1.x && crosspoint[i].layer == _segments[j].p1.layer){
                    if(_segments[j].p1.y < crosspoint[i].y && crosspoint[i].y  < _segments[j].p2.y){
                        splitsegment(crosspoint[i], _segments, j);
                    }
                }
            }
            else if(_segments[j].VORH ==2){
                if(crosspoint[i].y == _segments[j].p1.y && crosspoint[i].layer == _segments[j].p1.layer){
                    if(_segments[j].p1.x < crosspoint[i].x && crosspoint[i].x  < _segments[j].p2.x){
                        splitsegment(crosspoint[i], _segments, j);
                    }
                }
            }
            else{
                cout << _segments[j].VORH << endl;
                cout << "bug in XXX !\n0";
                exit(1);
            }
        }
    }
}

void NET::createroutingtree(vector<SEGMENT>& _segments,int mode,vector<vector<EDGE_C> >& edge,vector<int>& cv, vector<int>& ch,vector<int>& mw,vector<int>& ms ){
    int segnum = _segments.size(); 
    TREE_NODE currentpoint;
    int current_index=0;
    currentpoint.x = sinks_marge[0].grid_index_x;
    currentpoint.y = sinks_marge[0].grid_index_y;
    if(mode==0)
        currentpoint.layer = sinks_marge[0].layer+1;
    else
        currentpoint.layer = 1;

    currentpoint.parent_index=-1;
    currentpoint.child_C =-1;
    //currentpoint.nodeel = -1;
    currentpoint.weight =0;
    currentpoint.pinlayer.push_back(-1);
    currentpoint.level=0;
    RTree.push_back(currentpoint);
    
    while(segnum!=0){
        for(int i=0;i<_segments.size();i++){
            if(_segments[i].done)
                continue;

            if(_segments[i].p1.x == currentpoint.x && _segments[i].p1.y == currentpoint.y && _segments[i].p1.layer == currentpoint.layer){
                TREE_NODE Ttmp;
                Ttmp.x = _segments[i].p2.x;
                Ttmp.y = _segments[i].p2.y;
                Ttmp.layer = _segments[i].p2.layer;
                Ttmp.parent_index = current_index;
                Ttmp.child_C=-1;
                //Ttmp.nodeel = -1;
                Ttmp.pinlayer.push_back(-1);// = -1;
                RTree[current_index].child_index.push_back(RTree.size());
                Ttmp.weight = 0;
                Ttmp.level=0;
                RTree.push_back(Ttmp);
                _segments[i].done = true;
                segnum--;
            }
            else if(_segments[i].p2.x == currentpoint.x && _segments[i].p2.y == currentpoint.y && _segments[i].p2.layer == currentpoint.layer){
                TREE_NODE Ttmp;
                Ttmp.x = _segments[i].p1.x;
                Ttmp.y = _segments[i].p1.y;
                Ttmp.layer = _segments[i].p1.layer;
                Ttmp.parent_index = current_index;
                Ttmp.child_C =-1;
                //Ttmp.nodeel = -1;
                Ttmp.pinlayer.push_back(-1);// = -1;
                RTree[current_index].child_index.push_back(RTree.size());
                Ttmp.weight = 0;
                Ttmp.level=0;
                RTree.push_back(Ttmp);
                _segments[i].done = true;
                segnum--;
            }
        }
        current_index++;
        currentpoint.x = RTree[current_index].x;
        currentpoint.y = RTree[current_index].y;
        currentpoint.layer = RTree[current_index].layer;
    }

    modifyroutingtree();
    nodelevel(0); 
    
    //find the index which is the sink or the source pin
    for(int i=0;i<RTree.size();i++){
        for(int j=0;j<sinks_marge.size();j++){
            if(mode==0){
                if(sinks_marge[j].grid_index_x == RTree[i].x && sinks_marge[j].grid_index_y == RTree[i].y && sinks_marge[j].layer+1 == RTree[i].layer){
                    sinks_marge[j].RTree_index = i;
                    RTree[i].pinlayer[0] = 1;
                }
            }
            else{
                if(sinks_marge[j].grid_index_x == RTree[i].x && sinks_marge[j].grid_index_y == RTree[i].y){
                    sinks_marge[j].RTree_index = i;
                    if(RTree[i].pinlayer[0]==-1)
                        RTree[i].pinlayer[0] = sinks_marge[j].layer;
                    else{
                        RTree[i].pinlayer.push_back(sinks_marge[j].layer);
                    }
                }
            }
        }
    }

    computeinitialRTreeweight(0);

}

void NET::ct(vector<vector<EDGE_C> >& edge,vector<int>& cv, vector<int>& ch,vector<int>& mw,vector<int>& ms){
    for(int i=1;i<sinks_marge.size();i++){
        double delay;
        delay = elmoredelay(sinks_marge[i].RTree_index,edge,cv,ch,mw,ms);
        sinks_marge[i].elmoredelay = delay;
    }

}

void NET::splitsegment(SEGMENT_POINT crosspoint, vector<SEGMENT> & _segments, int index){
    SEGMENT Stmp;
    Stmp.p1.x = crosspoint.x;
    Stmp.p1.y = crosspoint.y;
    Stmp.p1.layer = crosspoint.layer;
    Stmp.p2.x = _segments[index].p2.x;
    Stmp.p2.y = _segments[index].p2.y;
    Stmp.p2.layer = _segments[index].p2.layer;
    Stmp.via = _segments[index].via;
    Stmp.VORH = _segments[index].VORH;
    Stmp.done = false;

    _segments[index].p2.x = crosspoint.x;
    _segments[index].p2.y = crosspoint.y;
    _segments[index].p2.layer = crosspoint.layer;

    _segments.push_back(Stmp);

}

void NET::ThreeDprojectto2D(int width,int high,int w_size,int h_size){
    vector<SEGMENT> twoDstmp_V;
    vector<SEGMENT> twoDstmp_H;
    for(int i=0;i<segments.size();i++){
        SEGMENT Stmp;
        if(segments[i].VORH==0)
            continue;
        
        Stmp.VORH = segments[i].VORH;
        Stmp.p1 = segments[i].p1;
        Stmp.p2 = segments[i].p2;

        if(Stmp.VORH ==1)
            twoDstmp_V.push_back(Stmp);
        else
            twoDstmp_H.push_back(Stmp);
    }
   
    vector<vector<int> > indextmp;//grouping segment with x or y position
    //marge two overlap segment to one segment with vertical segment
    for(int i=0;i<twoDstmp_V.size();i++){//use x position to group this segment
        int targ_x = twoDstmp_V[i].p1.x;
        int idtmp =-1;
        for(int j=0;j<indextmp.size();j++){
            if(targ_x == indextmp[j][0]){
                idtmp = j;
                indextmp[j].push_back(i);                
                break;
            }
        }
        if(idtmp==-1){
            vector<int> inttmp;
            inttmp.push_back(targ_x);
            inttmp.push_back(i);
            indextmp.push_back(inttmp);
        }
    }
    
    for(int i=0;i<indextmp.size();i++){
        vector<int> havesegment;
        havesegment.resize(high+1,0);
        for(int j=1;j<indextmp[i].size();j++){
            int y1 = twoDstmp_V[indextmp[i][j]].p1.y/h_size;
            int y2 = twoDstmp_V[indextmp[i][j]].p2.y/h_size;
            if(y1==y2)
                continue;
            for(int k=y1;k<y2;k++){ // debug change <= -> <
                havesegment[k]=1;
            }
        }

        int mode =0;
        int y1,y2;
        for(int j=0;j<high;j++){
            if(mode==0){
                if(havesegment[j]==1){
                    mode =1;
                    y1 = j;
                }
            }
            else{
                if(havesegment[j+1]==0||havesegment[j]==0){
                    mode=0;
                    if(havesegment[j]==0)
                        y2 = j;
                    else
                        y2 = j+1;
                    //y2 = j; // debug change  j -> j+1
                    SEGMENT sstmp;
                    sstmp.p1.x = indextmp[i][0];
                    sstmp.p1.y = y1*h_size;
                    sstmp.p1.layer = 1;
                    sstmp.p2.x = sstmp.p1.x;
                    sstmp.p2.y = y2*h_size;
                    sstmp.p2.layer = 1;
                    sstmp.done = false;
                    sstmp.via = false;
                    sstmp.VORH = 1;
                    twoDsegments.push_back(sstmp);
         //           tttmp.push_back(sstmp);
                }
            }
        }
    }

    indextmp.clear();
    
    //for horiziontal
    for(int i=0;i<twoDstmp_H.size();i++){//use x position to group this segment
        int targ_y = twoDstmp_H[i].p1.y;
        int idtmp =-1;
        for(int j=0;j<indextmp.size();j++){
            if(targ_y == indextmp[j][0]){
                idtmp = j;
                indextmp[j].push_back(i);
                break;
            }
        }
        if(idtmp==-1){
            vector<int> inttmp;
            inttmp.push_back(targ_y);
            inttmp.push_back(i);
            indextmp.push_back(inttmp);
        }
    }   

    for(int i=0;i<indextmp.size();i++){
        vector<int> havesegment;
        havesegment.resize(width+1,0);
        for(int j=1;j<indextmp[i].size();j++){
            int x1 = twoDstmp_H[indextmp[i][j]].p1.x/w_size;
            int x2 = twoDstmp_H[indextmp[i][j]].p2.x/w_size;
            if(x1==x2)
                continue;
            for(int k=x1;k<x2;k++){ // debug change <= -> <
                havesegment[k]=1;
            }
        }
        int mode =0;
        int x1,x2;
        //vector<SEGMENT> tttmp;
        for(int j=0;j<width;j++){
            if(mode==0){
                if(havesegment[j]==1){
                    mode =1;
                    x1 = j;
                }
            }
            else{
                if(havesegment[j+1]==0||havesegment[j]==0){
                    mode=0;
                    if(havesegment[j]==0)
                        x2 = j;
                    else
                        x2 = j+1;
                    //x2 = j; //debug change j >= j+1
                    SEGMENT sstmp;
                    sstmp.p1.x = x1*w_size;
                    sstmp.p1.y = indextmp[i][0];
                    sstmp.p1.layer =1;
                    sstmp.p2.x = x2*w_size;
                    sstmp.p2.y = sstmp.p1.y;
                    sstmp.p2.layer =1;
                    sstmp.done = false;
                    sstmp.via = false;
                    sstmp.VORH = 2;
                    twoDsegments.push_back(sstmp);
                }
            }
        }
        havesegment.clear();
    }
    indextmp.clear();
    
    
}

double NET::elmoredelay(int index ,vector<vector<EDGE_C> >& edge,vector<int>& cv, vector<int>& ch,vector<int>&mw,vector<int>& ms){
    //here layer number is really layer number(nctugr using)
    double delay=0.0;
    while(index != -1){
        if(RTree[index].child_C == -1){
            RTree[index].child_C = computechild_C(index,edge,cv,ch,mw,ms);
        }

        if(RTree[index].parent_index == -1){
            delay += driver_R*RTree[index].child_C;
        }
        else if(RTree[index].layer!= RTree[RTree[index].parent_index].layer){
            int low = RTree[index].layer;
            if(low > RTree[RTree[index].parent_index].layer)
                low = RTree[RTree[index].parent_index].layer;
            delay += (via_C[low-1]/2+RTree[index].child_C) * (via_R[low-1]);
        }
        else if(RTree[index].x == RTree[RTree[index].parent_index].x){
            double density = computedensity(RTree[index],RTree[RTree[index].parent_index],0,edge,cv,ch,mw,ms);
            delay += ((20)*wpt.getUC(RTree[index].layer,density,density,density,0)+RTree[index].child_C) * (40*wpt.getUR(RTree[index].layer,density,density,density,0));
        }
        else if(RTree[index].y == RTree[RTree[index].parent_index].y){
            double density = computedensity(RTree[index],RTree[RTree[index].parent_index],0,edge,cv,ch,mw,ms);
            delay += ((16)*wpt.getUC(RTree[index].layer,density,density,density,0) +RTree[index].child_C) * (32*wpt.getUR(RTree[index].layer,density,density,density,0));
        }
        else{
            cout << "bug in elmoredelay\n!";
            exit(1);
        }

        index = RTree[index].parent_index;
    }
    return delay; 
}

void NET::computeinitialRTreeweight(int index){
    if(RTree[index].child_index.size()==0){
        if(RTree[index].pinlayer[0] !=-1)
            RTree[index].weight = (RTree[index].pinlayer.size())*double(1.0/(sinks_marge.size()-1));
    }

    for(int i=0;i<RTree[index].child_index.size();i++)
        computeinitialRTreeweight(RTree[index].child_index[i]);

    double w=0;
    for(int i=0;i<RTree[index].child_index.size();i++)
        w += RTree[RTree[index].child_index[i]].weight;
    if(RTree[index].pinlayer[0]!=-1){
        if(index !=0)
            w += (RTree[index].pinlayer.size())*double(1.0/(sinks_marge.size()-1));
        else{
            w += (RTree[index].pinlayer.size()-1)*double(1.0/(sinks_marge.size()-1));
        }
    }
    RTree[index].weight = w;
}


double NET::computechild_C(int index,vector<vector<EDGE_C> >& edge,vector<int>& cv, vector<int>& ch,vector<int>&mw,vector<int>& ms){
    //include self C
    if(RTree[index].child_index.size()==0){
        if(RTree[index].pinlayer[0]==1)
            return load_C;
        else
            return 0;
    }
    else{
        double Ctmp =0.0;
        for(int i=0;i<RTree[index].child_index.size();i++){
            if(RTree[RTree[index].child_index[i]].child_C == -1){
                RTree[RTree[index].child_index[i]].child_C = computechild_C(RTree[index].child_index[i],edge,cv,ch,mw,ms);
            }
            Ctmp += RTree[RTree[index].child_index[i]].child_C;
            
            if(RTree[index].layer != RTree[RTree[index].child_index[i]].layer){
                if(RTree[index].layer > RTree[RTree[index].child_index[i]].layer)
                    Ctmp += via_C[RTree[index].layer-2];
                else
                    Ctmp += via_C[RTree[index].layer-1];
            }
            else if(RTree[index].x == RTree[RTree[index].child_index[i]].x){
                double d = computedensity(RTree[RTree[index].child_index[i]],RTree[index],0,edge,cv,ch,mw,ms);
                Ctmp += 40*wpt.getUC(RTree[index].layer,d,d,d,0);
            }
            else if (RTree[index].y == RTree[RTree[index].child_index[i]].y){
                double d = computedensity(RTree[RTree[index].child_index[i]],RTree[index],0,edge,cv,ch,mw,ms);
                Ctmp += 32*wpt.getUC(RTree[index].layer,d,d,d,0);
            }
            else{
                cout << "bug in computechild_C\n";
                exit(1);
            }
        }

        if(RTree[index].pinlayer[0]==1&&RTree[index].parent_index != -1){
            Ctmp +=load_C;
        }
        return Ctmp;
    }
}

double NET::elmoredelay_assign(TERMINAL& sink,vector<vector<EDGE_C> >& edge_c,vector<int>& cv, vector<int>& ch,vector<int>&mw,vector<int>& ms){
    //here layer number is layer index (layer assignment using)
    //tmp int Cvalue and Rvaule
    double Cvalue;
    double Rvalue;
    double delay=0.0;
    int index = sink.RTree_index;

    //compute the node which is connect sink first, aviod other elmore delay
    double total_C=0.0;
    double ctmp = 0.0;
    int layertmp = RTree[index].layer;
    if(RTree[index].parent_index == -1)
        layertmp = RTree[index].pinlayer[0];


    //high via delay
    if( sink.layer > layertmp){

        
        int top = sink.layer;

        //add load and check that pin layer > top
        for(int i=0;i < RTree[index].pinlayer.size();i++){
            if(RTree[index].pinlayer[i] > sink.layer)
                ctmp += load_C;
            if(top < RTree[index].pinlayer[i])
                top = RTree[index].pinlayer[i];
        }        

        //add child C(include wire segment from child to self)
        for(int i=0;i< RTree[index].child_index.size();i++){
            if( RTree[RTree[index].child_index[i]].layer > sink.layer){
                if(RTree[RTree[index].child_index[i]].child_C == -1)
                    RTree[RTree[index].child_index[i]].child_C = computechild_C_assign(RTree[index].child_index[i],edge_c,cv,ch,mw,ms);
                ctmp += RTree[RTree[index].child_index[i]].child_C;

                if(RTree[index].x == RTree[RTree[index].child_index[i]].x){
                    int x = RTree[index].x/32;
                    int y = (RTree[index].y > RTree[RTree[index].child_index[i]].y ? RTree[RTree[index].child_index[i]].y/40  : RTree[index].y/40);
                    int dindex = edge_c[x][y].C_V[RTree[ RTree[index].child_index[i]].layer];
                    if(dindex<0)
                        dindex=0;
                    int fatornot=0;
                    if(RTree[RTree[index].child_index[i]].wirenumber > 1)
                        fatornot =1;
                    ctmp += 40*RCtable[RTree[ RTree[index].child_index[i]].layer].UC[dindex][fatornot];
                }
                else if (RTree[index].y == RTree[RTree[index].child_index[i]].y){
                    int x = (RTree[index].x > RTree[RTree[index].child_index[i]].x ? RTree[RTree[index].child_index[i]].x/32  : RTree[index].x/32);
                    int y = RTree[index].y/40;
                    int dindex = edge_c[x][y].C_H[RTree[ RTree[index].child_index[i]].layer];
                    if(dindex<0)
                        dindex=0;
                    int fatornot=0;
                    if(RTree[RTree[index].child_index[i]].wirenumber > 1)
                       fatornot=1; 
                    ctmp += 32*RCtable[RTree[ RTree[index].child_index[i]].layer].UC[dindex][fatornot];
                }

                //find top layer
                if(top < RTree[RTree[index].child_index[i]].layer)
                    top = RTree[RTree[index].child_index[i]].layer;

            }
                
        }


        //add via C
        for(int i=top; i>sink.layer;i--)
            ctmp += via_C[i-1];

        //compute via delay from sink layer to self layer
        for(int i=sink.layer;i > layertmp;i--){
            int jtmp=0;
            if(RTree[index].parent_index==-1)
                jtmp=1;
            for(int j=jtmp;j<RTree[index].pinlayer.size();j++){
                if(RTree[index].pinlayer[j]==i){
                    ctmp += load_C;
                }
            }
            for(int j=0;j<RTree[index].child_index.size();j++){
                if(RTree[RTree[index].child_index[j]].layer == i){
                    if(RTree[RTree[index].child_index[j]].child_C == -1){
                        RTree[RTree[index].child_index[j]].child_C = computechild_C_assign(RTree[index].child_index[j],edge_c,cv,ch,mw,ms);
                    }
                    ctmp += RTree[RTree[index].child_index[j]].child_C;

                    if(RTree[index].x == RTree[RTree[index].child_index[j]].x){
                        int x = RTree[index].x/32;
                        int y = (RTree[index].y > RTree[RTree[index].child_index[j]].y ?  RTree[RTree[index].child_index[j]].y/40 : RTree[index].y/40 );
                        int dindex = edge_c[x][y].C_V[RTree[ RTree[index].child_index[j]].layer];
                        if(dindex<0)
                            dindex=0;
                        int fatornot = 0;
                        if(RTree[RTree[index].child_index[j]].wirenumber >1)
                            fatornot=1;
                        ctmp += 40*RCtable[RTree[ RTree[index].child_index[j]].layer].UC[dindex][fatornot];
                    }
                    else if (RTree[index].y == RTree[RTree[index].child_index[j]].y){
                        int x =  (RTree[index].x > RTree[RTree[index].child_index[j]].x ?  RTree[RTree[index].child_index[j]].x/32 : RTree[index].x/32 );
                        int y = RTree[index].y/40;
                        int dindex = edge_c[x][y].C_H[RTree[ RTree[index].child_index[j]].layer ];
                        if(dindex<0)
                            dindex=0;
                        int fatornot=0;
                        if(RTree[RTree[index].child_index[j]].wirenumber>1)
                            fatornot=1;
                        ctmp += 32*RCtable[RTree[ RTree[index].child_index[j] ].layer].UC[dindex][fatornot];
                    }

                }

            }
            delay += (via_C[i-1]/2+ctmp) * via_R[i-1];
            ctmp += via_C[i-1];
        }

    }
    else{
        int low = sink.layer;
        for(int i=0;i < RTree[index].pinlayer.size();i++){
            if(RTree[index].pinlayer[i] < sink.layer)
                ctmp += load_C;

            if(RTree[index].pinlayer[i] < low)
                low = RTree[index].pinlayer[i];
        }

        for(int i=0;i< RTree[index].child_index.size();i++){
            int fatornot =0;
            if(RTree[ RTree[index].child_index[i] ].wirenumber>1)
                fatornot=1;
            if(RTree[RTree[index].child_index[i]].layer < sink.layer){
                if(RTree[RTree[index].child_index[i]].child_C == -1)
                    RTree[RTree[index].child_index[i]].child_C = computechild_C_assign(RTree[index].child_index[i],edge_c,cv,ch,mw,ms);
                ctmp += RTree[RTree[index].child_index[i]].child_C;
                if(RTree[index].x == RTree[RTree[index].child_index[i]].x){
                    int x = RTree[index].x/32;
                    int y = (RTree[index].y > RTree[RTree[index].child_index[i]].y ?  RTree[RTree[index].child_index[i]].y/40 : RTree[index].y/40 );
                    int dindex = edge_c[x][y].C_V[RTree[ RTree[index].child_index[i]].layer];
                    if(dindex<0)
                        dindex=0;
                    ctmp += 40*RCtable[RTree[ RTree[index].child_index[i] ].layer].UC[dindex][fatornot];
                }
                else if (RTree[index].y == RTree[RTree[index].child_index[i]].y){
                    int x = (RTree[index].x > RTree[RTree[index].child_index[i]].x ?  RTree[RTree[index].child_index[i]].x/32 : RTree[index].x/32 );
                    int y = RTree[index].y/40;
                    int dindex = edge_c[x][y].C_H[RTree[ RTree[index].child_index[i]].layer];
                    if(dindex<0)
                        dindex=0;

                    ctmp += 32*RCtable[RTree[ RTree[index].child_index[i] ].layer].UC[dindex][fatornot];
                }

                if(RTree[RTree[index].child_index[i]].layer<low)
                    low = RTree[RTree[index].child_index[i]].layer;
            }

        }

        for(int i=low; i< sink.layer;i++)
            ctmp += via_C[i];
        for(int i=sink.layer; i < layertmp;i++){
            int jtmp=0;
            if(RTree[index].parent_index==-1)
                jtmp=1;
            for(int j=jtmp;j<RTree[index].pinlayer.size();j++){
                if(RTree[index].pinlayer[j]==i){
                    ctmp += load_C;
                }
            }
            for(int j=0;j<RTree[index].child_index.size();j++){
                int fatornot=0;
                if(RTree[ RTree[index].child_index[j] ].wirenumber>1)
                    fatornot=1;
                if(RTree[RTree[index].child_index[j]].layer == i){
                    if(RTree[RTree[index].child_index[j]].child_C == -1){
                        RTree[RTree[index].child_index[j]].child_C = computechild_C_assign(RTree[index].child_index[j],edge_c,cv,ch,mw,ms);
                    }
                    ctmp += RTree[RTree[index].child_index[j]].child_C;
                    if(RTree[index].x == RTree[RTree[index].child_index[j]].x){
                        int x = RTree[index].x/32;
                        int y = (RTree[index].y > RTree[RTree[index].child_index[j]].y? RTree[RTree[index].child_index[j]].y/40 : RTree[index].y/40);
                        int dindex = edge_c[x][y].C_V[RTree[ RTree[index].child_index[j] ].layer];
                        if(dindex<0)
                            dindex=0;
                        ctmp += 40*RCtable[RTree[ RTree[index].child_index[j] ].layer].UC[dindex][fatornot];
                    }
                    else if (RTree[index].y == RTree[RTree[index].child_index[j]].y){
                        int x = (RTree[index].x > RTree[RTree[index].child_index[j]].x? RTree[RTree[index].child_index[j]].x/32 : RTree[index].x/32);
                        int y = RTree[index].y/40;
                        int dindex = edge_c[x][y].C_H[RTree[ RTree[index].child_index[j] ].layer];
                        if(dindex<0)
                            dindex=0;
                        ctmp += 32*RCtable[RTree[ RTree[index].child_index[j] ].layer].UC[dindex][fatornot];
                    }
                }
            }

            delay += (via_C[i]/2+ctmp) * via_R[i];
            ctmp += via_C[i];
        }
    }

    while(index != -1){
        if(RTree[index].child_C == -1){
            RTree[index].child_C = computechild_C_assign(index,edge_c,cv,ch,mw,ms);
        }

        if(RTree[index].parent_index==-1){
            delay += driver_R * RTree[index].child_C;
        }
        else if(RTree[index].x == RTree[RTree[index].parent_index].x){

            int x = RTree[index].x/32;
            int y = (RTree[index].y > RTree[RTree[index].parent_index].y? RTree[RTree[index].parent_index].y/40 : RTree[index].y/40);
            int dindex = edge_c[x][y].C_V[RTree[index].layer];
            if(dindex<0)
                dindex=0;
            int fatornot=0;
            if(RTree[index].wirenumber>1)
                fatornot=1;
            delay += (20*RCtable[RTree[index].layer].UC[dindex][fatornot] + RTree[index].child_C) * (40*RCtable[RTree[index].layer].UR[dindex][fatornot]);
        }
        else if(RTree[index].y == RTree[RTree[index].parent_index].y){
            int x = (RTree[index].x > RTree[RTree[index].parent_index].x ? RTree[RTree[index].parent_index].x/32 : RTree[index].x/32);
            int y = RTree[index].y/40;
            int dindex = edge_c[x][y].C_H[RTree[index].layer];
            if(dindex<0)
                dindex=0;
            int fatornot=0;
            if(RTree[index].wirenumber>1)
                fatornot=1;
            delay += (16*RCtable[RTree[index].layer].UC[dindex][fatornot] + RTree[index].child_C) * (32*RCtable[RTree[index].layer].UR[dindex][fatornot]);
        }
        else{
            cout << "bug in elmoredelay_assign\n";
            exit(1);
        }
      

        double dd=0.0;
        if(RTree[index].parent_index != -1){
            dd = viaDelay(index,edge_c,cv,ch,mw,ms);
            delay += dd;
        }

        index = RTree[index].parent_index;  
    }
    return delay; 
}

double NET::viaDelay(int index,vector<vector<EDGE_C> >& edge_c,vector<int>& cv, vector<int>& ch,vector<int>&mw,vector<int>& ms){
// there layer is layertableindex
    int layer = RTree[index].layer;
    int parent_layer = RTree[RTree[index].parent_index].layer;
    int p_index = RTree[index].parent_index;
    if(RTree[p_index].parent_index==-1)
        parent_layer = RTree[p_index].pinlayer[0];

    if(layer > parent_layer){
        double Dtmp =0.0;
        double Ctmp =0.0;

        int top=layer;
        for(int i=0;i<RTree[p_index].child_index.size();i++){
            if(RTree[RTree[p_index].child_index[i]].layer > layer){
                if(RTree[RTree[p_index].child_index[i]].child_C ==-1)
                    RTree[RTree[p_index].child_index[i]].child_C = computechild_C_assign(RTree[p_index].child_index[i],edge_c,cv,ch,mw,ms);
                Ctmp += RTree[RTree[p_index].child_index[i]].child_C;
                
                int fatornot=0;
                if(RTree[RTree[p_index].child_index[i] ].wirenumber>1)
                    fatornot=1;

                if(RTree[ RTree[p_index].child_index[i] ].x == RTree[p_index].x){
                    int x = RTree[p_index].x/32;
                    int y = (RTree[ RTree[p_index].child_index[i] ].y > RTree[p_index].y ? RTree[p_index].y/40  : RTree[ RTree[p_index].child_index[i] ].y/40  );
                    int dindex = edge_c[x][y].C_V[RTree[RTree[p_index].child_index[i] ].layer];
                    if(dindex<0)
                        dindex = 0;
                    Ctmp += 40* RCtable[RTree[RTree[p_index].child_index[i] ].layer ].UC[dindex][fatornot];
                }
                else if(RTree[ RTree[p_index].child_index[i] ].y == RTree[p_index].y){
                    int x = (RTree[ RTree[p_index].child_index[i] ].x > RTree[p_index].x ? RTree[p_index].x/32 : RTree[ RTree[p_index].child_index[i] ].x/32);
                    int y = RTree[p_index].y/40;
                    int dindex = edge_c[x][y].C_H[RTree[RTree[p_index].child_index[i] ].layer];
                    if(dindex<0)
                        dindex=0;
                    Ctmp += 32*RCtable[RTree[RTree[p_index].child_index[i] ].layer].UC[dindex][fatornot];
                }
                if(RTree[RTree[p_index].child_index[i]].layer>top)
                    top = RTree[RTree[p_index].child_index[i]].layer;
            }
        }

        if(RTree[p_index].pinlayer[0]!=-1)
            for(int i=0;i<RTree[p_index].pinlayer.size();i++){
                if(RTree[p_index].pinlayer[i] > layer)
                    Ctmp += load_C;

                if(RTree[p_index].pinlayer[i] > top)
                    top = RTree[p_index].pinlayer[i];
            }

        for(int i=top;i>layer;i--)
            Ctmp += via_C[i-1];

        for(int i=layer;i>parent_layer;i--){
            for(int j=0;j<RTree[p_index].child_index.size();j++){
                if( RTree[RTree[p_index].child_index[j]].layer==i){
                    if(RTree[RTree[p_index].child_index[j]].child_C==-1){
                        RTree[RTree[p_index].child_index[j]].child_C = computechild_C_assign(RTree[p_index].child_index[j],edge_c,cv,ch,mw,ms);
                    }
                    Ctmp += RTree[RTree[p_index].child_index[j]].child_C;

                    int fatornot=0;
                    if(RTree[RTree[p_index].child_index[j]].wirenumber>1)
                        fatornot=1;


                    if(RTree[ RTree[p_index].child_index[j] ].x == RTree[p_index].x){
                        int x = RTree[p_index].x/32;
                        int y = (RTree[ RTree[p_index].child_index[j] ].y > RTree[p_index].y ? RTree[p_index].y/40  : RTree[ RTree[p_index].child_index[j] ].y/40  );
                        int dindex = edge_c[x][y].C_V[RTree[RTree[p_index].child_index[j] ].layer];
                        if(dindex<0)
                            dindex=0;
                        Ctmp += 40*RCtable[RTree[RTree[p_index].child_index[j]].layer].UC[dindex][fatornot];
                    }
                    else if(RTree[ RTree[p_index].child_index[j] ].y == RTree[p_index].y){
                        int x = (RTree[ RTree[p_index].child_index[j] ].x > RTree[p_index].x ? RTree[p_index].x/32  : RTree[ RTree[p_index].child_index[j] ].x/32  );
                        int y = RTree[p_index].y/40;
                        int dindex = edge_c[x][y].C_H[RTree[RTree[p_index].child_index[j]].layer];
                        if(dindex<0)
                            dindex=0;
                        Ctmp += 32*RCtable[RTree[RTree[p_index].child_index[j]].layer].UC[dindex][fatornot];
                    }
                }
            }

            int jtmp=0;
            if(RTree[p_index].parent_index == -1)
                jtmp = 1;
            for(int j=jtmp;j<RTree[p_index].pinlayer.size();j++){
                if(RTree[p_index].pinlayer[j]==i)
                    Ctmp += load_C;
            }

            Dtmp += (via_C[i-1]/2+Ctmp)*via_R[i-1];

            Ctmp += via_C[i-1];
        }

        return Dtmp;

    }
    else{
        double Dtmp =0.0;
        double Ctmp =0.0;
        int low=layer;
        for(int i=0;i<RTree[p_index].child_index.size();i++){
            if(RTree[RTree[p_index].child_index[i]].layer < layer){
                if(RTree[RTree[p_index].child_index[i]].child_C ==-1)
                    RTree[RTree[p_index].child_index[i]].child_C = computechild_C_assign(RTree[p_index].child_index[i],edge_c,cv,ch,mw,ms);
                Ctmp += RTree[RTree[p_index].child_index[i]].child_C;
                int fatornot=0;
                if(RTree[RTree[p_index].child_index[i] ].wirenumber>1)
                    fatornot=1;
                if(RTree[ RTree[p_index].child_index[i] ].x == RTree[p_index].x){
                    int x = RTree[p_index].x/32;
                    int y = (RTree[ RTree[p_index].child_index[i] ].y > RTree[p_index].y ? RTree[p_index].y/40  : RTree[ RTree[p_index].child_index[i] ].y/40  );
                    int dindex = edge_c[x][y].C_V[RTree[RTree[p_index].child_index[i] ].layer];
                    if(dindex<0)
                        dindex = 0;
                    Ctmp += 40* RCtable[RTree[RTree[p_index].child_index[i] ].layer ].UC[dindex][fatornot];
                }
                else if(RTree[ RTree[p_index].child_index[i] ].y == RTree[p_index].y){
                    int x = (RTree[ RTree[p_index].child_index[i] ].x > RTree[p_index].x ? RTree[p_index].x/32 : RTree[ RTree[p_index].child_index[i] ].x/32);
                    int y = RTree[p_index].y/40;
                    int dindex = edge_c[x][y].C_H[RTree[RTree[p_index].child_index[i] ].layer];
                    if(dindex<0)
                        dindex=0;
                    Ctmp += 32*RCtable[RTree[RTree[p_index].child_index[i] ].layer].UC[dindex][fatornot];
                }

                if(RTree[RTree[p_index].child_index[i]].layer < low)
                    low = RTree[RTree[p_index].child_index[i]].layer;
            }
        }

        if(RTree[p_index].pinlayer[0]!=-1)
            for(int i=0;i<RTree[p_index].pinlayer.size();i++){
                if(RTree[p_index].pinlayer[i] < layer)
                    Ctmp += load_C;

                if(RTree[p_index].pinlayer[i] < low)
                    low = RTree[p_index].pinlayer[i];
            }

        for(int i=low;i <layer;i++)
            Ctmp += via_C[i];

        for(int i=layer;i < parent_layer;i++){
            for(int j=0;j<RTree[p_index].child_index.size();j++){
                if( RTree[RTree[p_index].child_index[j]].layer==i){
                    if(RTree[RTree[p_index].child_index[j]].child_C==-1){
                        RTree[RTree[p_index].child_index[j]].child_C = computechild_C_assign(RTree[p_index].child_index[j],edge_c,cv,ch,mw,ms);
                    }
                    Ctmp += RTree[RTree[p_index].child_index[j]].child_C;

                    int fatornot=0;
                    if(RTree[RTree[p_index].child_index[j]].wirenumber>1)
                        fatornot=1;

                    if(RTree[ RTree[p_index].child_index[j] ].x == RTree[p_index].x){
                        int x = RTree[p_index].x/32;
                        int y = (RTree[ RTree[p_index].child_index[j] ].y > RTree[p_index].y ? RTree[p_index].y/40  : RTree[ RTree[p_index].child_index[j] ].y/40  );
                        int dindex = edge_c[x][y].C_V[RTree[RTree[p_index].child_index[j] ].layer];
                        if(dindex<0)
                            dindex=0;
                        Ctmp += 40*RCtable[RTree[RTree[p_index].child_index[j]].layer].UC[dindex][fatornot];
                    }
                    else if(RTree[ RTree[p_index].child_index[j] ].y == RTree[p_index].y){
                        int x = (RTree[ RTree[p_index].child_index[j] ].x > RTree[p_index].x ? RTree[p_index].x/32  : RTree[ RTree[p_index].child_index[j] ].x/32  );
                        int y = RTree[p_index].y/40;
                        int dindex = edge_c[x][y].C_H[RTree[RTree[p_index].child_index[j]].layer];
                        if(dindex<0)
                            dindex=0;
                        Ctmp += 32*RCtable[RTree[RTree[p_index].child_index[j]].layer].UC[dindex][fatornot];
                    }
                }
            }

            int jtmp=0;
            if(RTree[p_index].parent_index == -1)
                jtmp = 1;
            for(int j=jtmp;j<RTree[p_index].pinlayer.size();j++){
                if(RTree[p_index].pinlayer[j]==i)
                    Ctmp += load_C;
            }

            Dtmp += (via_C[i]/2+Ctmp)*via_R[i];
            Ctmp += via_C[i];
        }
        
        return Dtmp; 
    }

}

double NET::computechild_C_assign(int index,vector<vector<EDGE_C> >& edge_c,vector<int>& cv, vector<int>& ch,vector<int>&mw,vector<int>& ms){
    //include self C, it mean the wire C connect index with index's all child
    //there layer is layertableindex 
    if(RTree[index].child_index.size()==0){
        double total_C=0.0;
        double dtmp=0.0;
        if(RTree[index].pinlayer[0] != -1){
            int top,low;
            top = RTree[index].layer;
            if(RTree[index].parent_index==-1)
                top = RTree[index].pinlayer[0];
            low = top;
                
            for(int i=0;i<RTree[index].pinlayer.size();i++){
                if(top < RTree[index].pinlayer[i])
                    top = RTree[index].pinlayer[i];
                if(low > RTree[index].pinlayer[i])
                    low = RTree[index].pinlayer[i];
            }

            for(int i=low;i<top;i++){
                total_C += via_C[i];
            }
            return (RTree[index].parent_index !=-1 ? total_C + RTree[index].pinlayer.size()*load_C:total_C + (RTree[index].pinlayer.size()-1)*load_C);

        }
        else{
            return 0.0;
        }
    }
    else{
        double Ctmp =0.0;
        int low,top,present;
        low = RTree[index].layer;
        if(RTree[index].parent_index == -1)
            low = RTree[index].pinlayer[0];
        top = low;
        present = low;

        //add total child C
        for(int i=0;i<RTree[index].child_index.size();i++){
            if(RTree[RTree[index].child_index[i]].child_C == -1){
                RTree[RTree[index].child_index[i]].child_C = computechild_C_assign(RTree[index].child_index[i],edge_c,cv,ch,mw,ms);
            }
            Ctmp += RTree[RTree[index].child_index[i]].child_C;

            int fatornot=0;
            if(RTree[RTree[index].child_index[i] ].wirenumber>1)
                fatornot=1;

            if(RTree[index].x == RTree[RTree[index].child_index[i]].x){
                int x = RTree[index].x/32;
                int y = (RTree[index].y > RTree[RTree[index].child_index[i]].y? RTree[RTree[index].child_index[i]].y/40 :RTree[index].y/40);
                int dindex = edge_c[x][y].C_V[RTree[RTree[index].child_index[i] ].layer];
                if(dindex<0)
                    dindex=0;
                Ctmp += 40*RCtable[RTree[RTree[index].child_index[i] ].layer].UC[dindex][fatornot];
            }
            else if (RTree[index].y == RTree[RTree[index].child_index[i]].y){
                int x = (RTree[index].x > RTree[RTree[index].child_index[i]].x? RTree[RTree[index].child_index[i]].x/32 :RTree[index].x/32);
                int y = RTree[index].y/40;
                int dindex = edge_c[x][y].C_H[RTree[RTree[index].child_index[i] ].layer];
                if(dindex<0)
                    dindex=0;
                Ctmp += 32*RCtable[RTree[RTree[index].child_index[i] ].layer].UC[dindex][fatornot];
            }
            else{
                cout << name << "  " << RTree[index].x << "  " << RTree[index].y << "  " << RTree[index].layer << "  " << index << endl;
                cout << name << "  " << RTree[RTree[index].parent_index].x << "  " << RTree[RTree[index].parent_index].y << "  " << RTree[RTree[index].parent_index].layer << "  " << RTree[index].parent_index << endl;
                cout << "bug in computechild_C_assign\n";
                exit(1);
            }

            if(RTree[RTree[index].child_index[i]].layer > top)
                top = RTree[RTree[index].child_index[i]].layer;
            else if(RTree[RTree[index].child_index[i]].layer < low)
                low = RTree[RTree[index].child_index[i]].layer;

        }

        //find top and low layer. If the node has sink pin, add the load C
        if(RTree[index].pinlayer[0]!=-1){
            for(int i=0;i<RTree[index].pinlayer.size();i++){
                if(RTree[index].pinlayer[i] > top)
                    top = RTree[index].pinlayer[i];
                if(RTree[index].pinlayer[i] < low)
                    low = RTree[index].pinlayer[i];

            }
            Ctmp += load_C*RTree[index].pinlayer.size();
            // sub the load C which is driver R 
            if(RTree[index].parent_index ==-1)
                Ctmp -= load_C;
        }

        //add total via C
        for(int i=low;i<top;i++){
            Ctmp += via_C[i];
        }

        return Ctmp;
    }
}

int NET::viacount(){
    int count=0;
    for(int i=0;i<RTree.size();i++){
        int top,low;
        if(i==0){
            top = RTree[i].pinlayer[0];
            low = top;
        }
        else{
            top = RTree[i].layer;
            low = top;
        }
        if(RTree[i].pinlayer[0]!=-1){
            for(int j=0;j<RTree[i].pinlayer.size();j++){
                if(RTree[i].pinlayer[j] > top)
                    top = RTree[i].pinlayer[j];
                if(RTree[i].pinlayer[j] < low)
                    low = RTree[i].pinlayer[j];
            }
        }
        
        for(int j=0;j<RTree[i].child_index.size();j++){
            if(RTree[RTree[i].child_index[j]].layer > top)
                top = RTree[RTree[i].child_index[j]].layer;
            if(RTree[RTree[i].child_index[j]].layer < low)
                low = RTree[RTree[i].child_index[j]].layer;
        }
        
        count += top-low;
    }
    return count;
}

void CIRCUIT::createtimingdata(string filename){
    filename = filename + ".timingdata";
    ofstream outfile(filename.c_str());
    int netnum=0;
    for(int i=0;i<nets.size();i++){
        if(nets[i].havesegment)
            netnum++;
    }
    double totalslack=0;
    double minslack=0;
    srand(time(NULL));
    outfile << "netnumber: " << netnum << endl;
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        
        outfile << nets[i].getname() << "  " << nets[i].sinks_marge.size() << endl;
        outfile << nets[i].sinks_marge[0].grid_index_x << "  " << nets[i].sinks_marge[0].grid_index_y << "  " << nets[i].sinks_marge[0].layer << "  -1\n";
        for(int j=1;j<nets[i].sinks_marge.size();j++){
            int random_range = ((rand()%20)+1)-10;
            double random_percent = (random_range *1.0) / 100;
            double timingdata = nets[i].sinks_marge[j].elmoredelay * (1.0+random_percent);
            if(timingdata - nets[i].sinks_marge[j].elmoredelay < 0)
                totalslack += timingdata - nets[i].sinks_marge[j].elmoredelay;
            if(timingdata - nets[i].sinks_marge[j].elmoredelay < minslack)
                minslack = timingdata - nets[i].sinks_marge[j].elmoredelay ;
            outfile << nets[i].sinks_marge[j].grid_index_x << "  " << nets[i].sinks_marge[j].grid_index_y << "  " << nets[i].sinks_marge[j].layer << "  " << timingdata << endl;

        }
    }
    outfile.close();
    
    cout << "total nslack: " << totalslack << endl << "wns: " << minslack << endl;
}

void CIRCUIT::readtimingdata(string filename){
    filename = filename + ".timingdata";
    ifstream infile(filename.c_str());
    int netnum;
    string stmp;
    infile >> stmp >> netnum;
    for(int i=0;i<netnum;i++){
        string netname;
        int pinnum;
        infile >> netname >> pinnum;
        int index = netsid2index[netname];
        infile >> stmp >> stmp >> stmp >> stmp;
        for(int j=1;j<nets[index].sinks_marge.size();j++){
            infile >> stmp >> stmp >> stmp >> nets[index].sinks_marge[j].elmoredelay;
        }
    }
    infile.close();
}

void CIRCUIT::outputroutingresult(string filename){
    ofstream outfile(filename.c_str());
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        outfile << nets[i].getname() << " 0 " << nets[i].resultsegments.size()<<endl;
        for(int j=0;j<nets[i].resultsegments.size();j++){     
            outfile << "(" << nets[i].resultsegments[j].p1.x <<","<<nets[i].resultsegments[j].p1.y<<","<<nets[i].resultsegments[j].p1.layer+1<<")-("<<nets[i].resultsegments[j].p2.x<<","<<nets[i].resultsegments[j].p2.y<<","<<nets[i].resultsegments[j].p2.layer+1<<") " << nets[i].resultsegments[j].wirenumber << endl;
        }
        outfile << "!\n";
    }
    outfile.close();
}
void NET::modifyroutingtree(){
    int l,h;
    l=32;
    h=40;
    int number = RTree.size();
    
    for(int i=0;i<number;i++){
        if(RTree[i].child_index.size()==0)
            continue;
        for(int j=0;j<RTree[i].child_index.size();j++){
            int parent_index = i;
            int firstnewnodeindex = RTree.size();
            bool addornot = false;
            if(RTree[i].layer != RTree[RTree[i].child_index[j]].layer){
                if(RTree[i].layer < RTree[RTree[i].child_index[j]].layer){
                    for(int k=RTree[i].layer+1; k< RTree[RTree[i].child_index[j]].layer ; k=k+1){
                        TREE_NODE  T_tmp;
                        T_tmp.x = RTree[i].x;
                        T_tmp.y = RTree[i].y;
                        T_tmp.layer = k;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }
                }
                else{
                    for(int k=RTree[i].layer-1; k > RTree[RTree[i].child_index[j]].layer ; k=k-1){
                        TREE_NODE  T_tmp;
                        T_tmp.x = RTree[i].x;
                        T_tmp.y = RTree[i].y;
                        T_tmp.layer = k;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }
                }
            }
            else if(RTree[i].x == RTree[RTree[i].child_index[j]].x){//vertical
                if(RTree[i].y <  RTree[RTree[i].child_index[j]].y){
                    for(int k=RTree[i].y+h; k < RTree[RTree[i].child_index[j]].y;k=k+h){
                        TREE_NODE  T_tmp;
                        T_tmp.x = RTree[i].x;
                        T_tmp.y = k;               
                        T_tmp.layer = RTree[i].layer;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);// = -1;
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }

                }
                else{
                    for(int k=RTree[i].y-h; k > RTree[RTree[i].child_index[j]].y;k=k-h){
                        TREE_NODE  T_tmp;
                        T_tmp.x = RTree[i].x;
                        T_tmp.y = k;
                        T_tmp.layer = RTree[i].layer;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);// = -1;
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }

                }
            }
            else if(RTree[i].y == RTree[RTree[i].child_index[j]].y){
                if(RTree[i].x <  RTree[RTree[i].child_index[j]].x){
                    for(int k=RTree[i].x+l; k < RTree[RTree[i].child_index[j]].x;k=k+l){
                        TREE_NODE  T_tmp;
                        T_tmp.x = k;
                        T_tmp.y = RTree[i].y;
                        T_tmp.layer = RTree[i].layer;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);// = -1;
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }
                }
                else{
                    for(int k=RTree[i].x-l; k > RTree[RTree[i].child_index[j]].x;k=k-l){
                        TREE_NODE  T_tmp;
                        T_tmp.x = k;
                        T_tmp.y = RTree[i].y;
                        T_tmp.layer = RTree[i].layer;
                        T_tmp.parent_index = parent_index;
                        T_tmp.child_index.push_back(RTree.size()+1);
                        parent_index = RTree.size();
                        T_tmp.child_C = -1;
                        T_tmp.pinlayer.push_back(-1);// = -1;
                        T_tmp.weight = 1;
                        T_tmp.level=0;
                        RTree.push_back(T_tmp);
                        addornot = true;
                    }
                }
            }

            if(!addornot)
                continue;
            RTree[RTree[i].child_index[j]].parent_index = parent_index;
            RTree[parent_index].child_index[0] = RTree[i].child_index[j];
            RTree[i].child_index[j]= firstnewnodeindex;
            
        }
    }

}

inline double NET::computedensity(TREE_NODE n, TREE_NODE p,int wirenumber,vector<vector<EDGE_C> >& edge_c,vector<int>& C_V, vector<int>& C_H,vector<int>&minwidth,vector<int>& minspacing){
    double density=0;
    if(n.x == p.x){
        int y = (n.y>p.y ? p.y : n.y);
        if(edge_c[n.x/32][y/40].C_V[n.layer]-wirenumber <= 0)
            density = 1;
        else{
            int totalc = C_V[n.layer] / (minspacing[n.layer] + minwidth[n.layer]);
            double a = totalc-edge_c[n.x/32][y/40].C_V[n.layer]+wirenumber;
            double b = totalc*1.0;
            density = a/b;
        }

    }
    else if(n.y==p.y){
        int x = (n.x>p.x ? p.x  : n.x);
        if(edge_c[x/32][n.y/40].C_H[n.layer]-wirenumber <= 0)
            density = 1;
        else{
            int totalc = C_H[n.layer] / (minspacing[n.layer] + minwidth[n.layer]);
            double a =  (totalc-edge_c[x/32][n.y/40].C_H[n.layer]+wirenumber);
            double b = totalc*1.0;

            density = a/b;

        }
    }
    else{
        cout << "bug in computedensity!\n";
        exit(1);
     
    }
    

    return density*10;
}

double CIRCUIT::computedensity(TREE_NODE n, TREE_NODE p,int wirenumber,int layer){
    double density=0;
    if(n.x == p.x){
        int y = (n.y>p.y ? p.y : n.y);
        if(edge_c[n.x/32][y/40].C_V[layer]-wirenumber <= 0)
            density = 1;
        else{
            int totalc = C_V[layer] / (minspacing[layer] + minwidth[layer]);
            double a = totalc-edge_c[n.x/32][y/40].C_V[layer]+wirenumber;
            double b = totalc;
            density = a/b;
        }
    }
    else if(n.y==p.y){
        int x = (n.x>p.x ? p.x  : n.x);
        if(edge_c[x/32][n.y/40].C_H[layer]-wirenumber <= 0)
            density = 1;
        else{
            int totalc = C_H[layer] / (minspacing[layer] + minwidth[layer]);
            double a =  (totalc-edge_c[x/32][n.y/40].C_H[layer]+wirenumber);
            double b = totalc;
            density = a/b;

        }
    }
    else{
        cout << "bug in computedensity!\n";
        exit(1);
    }
    return density*10;
}


void CIRCUIT::computegecllcapacityfornctugr(vector<TREE_NODE>& rt){
    for(int i=1;i<rt.size();i++){
        if(rt[i].layer != rt[rt[i].parent_index].layer)
            continue;
        else if(rt[i].x == rt[rt[i].parent_index].x){
            int y = (rt[i].y > rt[rt[i].parent_index].y ? rt[rt[i].parent_index].y : rt[i].y);
            edge_c[rt[i].x/32][y/40].C_V[rt[i].layer-1]--;
        }
        else{
            int x = (rt[i].x > rt[rt[i].parent_index].x ? rt[rt[i].parent_index].x : rt[i].x);
            edge_c[x/32][rt[i].y/40].C_H[rt[i].layer-1]--;
        }
    }
}

void CIRCUIT::computeoverflow(){
    int total=0;
    int max=0;
    for(int i=0;i<edge_c.size()-1;i++){
        for(int j=0;j<edge_c[i].size();j++){
            for(int k=1;k<=8;k++){
                int ctmp;
                if(k%2==0)
                    ctmp = edge_c[i][j].C_H[k];
                else
                    ctmp = edge_c[i][j].C_V[k];
                if(ctmp<0){
                    ctmp = ctmp *-1 * (minwidth[k]+minspacing[k]);
                    total += ctmp;
                    if(ctmp > max)
                        max = ctmp;
                }
            }
        }
    }
    cout << "total :" << total << "  max:" << max <<endl;
}

void CIRCUIT::computeaveragedensity(){
    for(int i=0;i<edge_c.size();i++){
        for(int j=0;j<edge_c[i].size();j++){
            edge_c[i][j].aveg_c_v = edge_c[i][j].total_c_v;
            edge_c[i][j].aveg_c_h = edge_c[i][j].total_c_h;
        }
    }

    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        for(int j=1;j<nets[i].RTree.size();j++){
            if(nets[i].RTree[j].x== nets[i].RTree[nets[i].RTree[j].parent_index].x){
                int x = nets[i].RTree[j].x/32;
                int y = (nets[i].RTree[j].y > nets[i].RTree[nets[i].RTree[j].parent_index].y ? nets[i].RTree[nets[i].RTree[j].parent_index].y/40 : nets[i].RTree[j].y/40);
                edge_c[x][y].aveg_c_v--;
            }
            else{
                int x = (nets[i].RTree[j].x > nets[i].RTree[nets[i].RTree[j].parent_index].x ? nets[i].RTree[nets[i].RTree[j].parent_index].x/32 : nets[i].RTree[j].x/32);
                int y = nets[i].RTree[j].y/40;
                edge_c[x][y].aveg_c_h--;
            }
        }
    }

    for(int i=0;i<edge_c.size();i++){
        for(int j=0;j<edge_c[i].size();j++){
            edge_c[i][j].aveg_c_v = edge_c[i][j].aveg_c_v/4;
            edge_c[i][j].aveg_c_h = edge_c[i][j].aveg_c_h/4;
        }
    }
}

void CIRCUIT::createlookuptable(){

    for(int i=1;i<9;i++){
        int totalc = (C_V[i]+C_H[i]) / (minspacing[i] + minwidth[i]);
        for(int j=totalc;j >= 1;j--){
            double a = j*1.0;
            double b = totalc*1.0;
            double density = (a/b)*10;
            double C0 = wpt.getUC(i+1,density,density,density,0);
            double C1 = wpt.getUC(i+1,density,density,density,1);
            double R0 = wpt.getUR(i+1,density,density,density,0);
            double R1 = wpt.getUR(i+1,density,density,density,1);
            
            vector<double> cdtmp;
            vector<double> rdtmp;
            cdtmp.push_back(C0);
            cdtmp.push_back(C1);
            rdtmp.push_back(R0);
            rdtmp.push_back(R1);
            RCtable[i].UC.push_back(cdtmp);
            RCtable[i].UR.push_back(rdtmp);
        }
    }
}

inline int rounding(double val) {  
      return (int) (val + 0.5);
}


