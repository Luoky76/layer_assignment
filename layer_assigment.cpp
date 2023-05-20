#include "database.h"
#include <limits>
// #include <omp.h>

inline int rounding(double val) {  
      return (int) (val + 0.5);
}

double k1 = 1;
int globaliteration;
inline double astar(TREE_NODE rt,double child_c){
    //return RCtable[8].UR[0][0]*(rt.level-1)*(0.5*RCtable[8].UC[0][0]+child_c*(rt.level-1)+RCtable[8].UC[0][0]*0.5*(rt.level));
    return 0;
}

bool struct_sortel_cmp(SORTEL a, SORTEL b){
    return a.eldelay > b.eldelay;
}

bool struct_sortel_cmp2(SORTEL a, SORTEL b){
    return a.eldelay < b.eldelay;
}

void CIRCUIT::initialLA(){
    clock_t start,end;

    //double time;
    //struct timeval s_t, e_t;
    //gettimeofday(&s_t, NULL);
    chrono::steady_clock::time_point STc = chrono::steady_clock::now();
    
    computinggcellcapacity();
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        ablenet++;
        dynamic_program_main(nets[i],i,0,0,0);//initial layer assignment without parallel wires and NDR wires
    } 
    computeoverflow();

    //gettimeofday(&e_t, NULL);
    chrono::steady_clock::time_point ETc = chrono::steady_clock::now();

    //time = ((e_t.tv_sec*1000000)+e_t.tv_usec)-((s_t.tv_sec*1000000)+s_t.tv_usec);
    //time /= 1000000;
    //cout << "done initalLA , runtime: " << time << endl;
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double> >(ETc - STc);
    cout << "done initalLA , runtime: " << time_span.count() << " (s)" << endl;
    
    computedelay();
    computecritical();

}

void CIRCUIT::computinggcellcapacity(){//the capacity of 2D grid edges are reduced by the 2D routed net 
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        for(int j=0;j<nets[i].twoDsegments.size();j++){
            int frond,end;
            if(nets[i].twoDsegments[j].VORH==1){
                frond = nets[i].twoDsegments[j].p1.y/gridsize_y;
                end = nets[i].twoDsegments[j].p2.y/gridsize_y;
                for(int k=frond;k<end;k++){
                    edge_c[nets[i].twoDsegments[j].p1.x/gridsize_x][k].total_c_v--;
                }
            }
            else if(nets[i].twoDsegments[j].VORH==2){
                frond = nets[i].twoDsegments[j].p1.x/gridsize_x;
                end = nets[i].twoDsegments[j].p2.x/gridsize_x;
                for(int k=frond;k<end;k++){
                    edge_c[k][nets[i].twoDsegments[j].p1.y/gridsize_y].total_c_h--;
                }
            }
        }
            
    }
}

vector<SORTEL> CIRCUIT::checkillegalnet(int itreation){
    vector<SORTEL> _netvector;
    int converge = itreation/5;
    int le = 2;
    double wns=0.0;
    int wnsindex=0;
    double base = 0.05;
    double max=0;
    double increasecost= base*pow(le,itreation);
    for(int i=0;i<edge_c.size();i++){
        for(int j=0;j<edge_c[i].size();j++){
            for(int k=1;k<9;k++){
                if(k%2==1){//vertical
                    if(edge_c[i][j].C_V[k] < 0){
                        map<int, int>::iterator it;
                        for(it = edge_c[i][j].netid[k].begin(); it != edge_c[i][j].netid[k].end(); it++){
                            if(!nets[it->second].needreassign){
                                SORTEL Stmp;
                                Stmp.eldelay = nets[it->second].eldelay;
                                Stmp.index = it->second;
                                _netvector.push_back(Stmp);
                                nets[it->second].needreassign = true;
                            }
                        }
                        edge_c[i][j].historic[k] +=increasecost;
                    }
                }
                else{
                    if(edge_c[i][j].C_H[k] < 0){
                        map<int,int>::iterator it;
                        for(it = edge_c[i][j].netid[k].begin(); it != edge_c[i][j].netid[k].end(); it++){
                            if(!nets[it->second].needreassign){
                                SORTEL Stmp;
                                Stmp.eldelay = nets[it->second].eldelay;
                                Stmp.index = it->second;
                                _netvector.push_back(Stmp);
                                nets[it->second].needreassign = true;
                            }
                        }

                        edge_c[i][j].historic[k] +=increasecost;
                    }
                }
            }
        }
    }

    
    if(globaliteration>fcfsthreshold+1 && parallelmode == 1)
        //sort(_netvector.begin(),_netvector.end(),struct_sortel_cmp);
        stable_sort(_netvector.begin(),_netvector.end(),struct_sortel_cmp);
    return _netvector;  
};

void CIRCUIT::ripup(NET& _net,int netindex){
    for(int i=0;i<_net.resultsegments.size();i++){
        if( _net.resultsegments[i].VORH == 1){
            edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].C_V[_net.resultsegments[i].p1.layer]+= _net.resultsegments[i].wirenumber;
            if(_net.resultsegments[i].wirenumber>1&&_net.resultsegments[i].needripup==true)
                edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].total_c_v += (_net.resultsegments[i].wirenumber-1);
        }
        else if( _net.resultsegments[i].VORH == 2 ){
            edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].C_H[_net.resultsegments[i].p1.layer]+= _net.resultsegments[i].wirenumber;

            if(_net.resultsegments[i].wirenumber>1&&_net.resultsegments[i].needripup==true)
                edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].total_c_h += (_net.resultsegments[i].wirenumber-1);
        }

        if( _net.resultsegments[i].VORH!=0){
            map<int,int>::iterator it;
            it = edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].netid[_net.resultsegments[i].p1.layer].find( netindex);
        
            if(it != edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].netid[_net.resultsegments[i].p1.layer].end()){
                edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].netid[_net.resultsegments[i].p1.layer].erase(it);
            }
            else{
                cout << "bug in ripup\n!";
                exit(1);
            }
        }

    }
    _net.resultsegments.clear();
};

void CIRCUIT::computeviacount(){
    int count=0;
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        count += nets[i].viacount();
    }

    cout << "totalviacount : " << count << endl;
}

void CIRCUIT::recoverpathtoedge(NET& _net,int netindex){

    for(int i=0;i<_net.resultsegments.size();i++){
        if( _net.resultsegments[i].VORH == 1){
            edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].C_V[_net.resultsegments[i].p1.layer] -= _net.resultsegments[i].wirenumber;

            if(_net.resultsegments[i].wirenumber>1&&_net.resultsegments[i].needripup==true)
                edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].total_c_v += (_net.resultsegments[i].wirenumber-1);
        }
        else if( _net.resultsegments[i].VORH == 2 ){
            edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].C_H[_net.resultsegments[i].p1.layer] -= _net.resultsegments[i].wirenumber;
            if(_net.resultsegments[i].wirenumber>1&&_net.resultsegments[i].needripup==true)
                edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].total_c_h += (_net.resultsegments[i].wirenumber-1);
        }

        if( _net.resultsegments[i].VORH!=0){
            edge_c[ _net.resultsegments[i].p1.x/32][_net.resultsegments[i].p1.y/40].netid[_net.resultsegments[i].p1.layer][netindex] = netindex;
        }
    }
}

void CIRCUIT::RRA_main(){
    int itnumber=1;

    //double tdtime = 0.0;
    //double ttime =0.0;
    //struct timeval s_t, e_t;
    //gettimeofday(&s_t, NULL);
    chrono::duration<double> RRA_time_span, DP_time_span;
    chrono::steady_clock::time_point STc = chrono::steady_clock::now();

    globaliteration=1;
    while(1){
        for(int i=0;i<nets.size();i++){
            if(!nets[i].havesegment)
                continue;
            nets[i].needreassign = false;
        }
        //struct timeval s_t2, e_t2;
        //double time;
        //gettimeofday(&s_t2, NULL);
        chrono::steady_clock::time_point STc2 = chrono::steady_clock::now();

        vector<SORTEL> illegalnetid = checkillegalnet(itnumber);
        
        //check the wire congestion constraint
        if(illegalnetid.empty())
            break;
        
        //node weight threshold 
        if(itnumber<7&&weightmode==1){
            for(int i=0;i< nets.size();i++){
                if(!nets[i].havesegment)
                    continue;
                nets[i].updateRTreeweight();
            }
        }
        
        cout << "start ripup and reassign:\n ";

        int illnumber = illegalnetid.size();
        for(int j=0;j<illnumber;j++){
            ripup(nets[illegalnetid[j].index], illegalnetid[j].index);
            dynamic_program_main(nets[illegalnetid[j].index],illegalnetid[j].index,parallelmode,0,0);
        }

        computeoverflow();

        //gettimeofday(&e_t2, NULL);
        //time = ((e_t2.tv_sec*1000000)+e_t2.tv_usec)-((s_t2.tv_sec*1000000)+s_t2.tv_usec);
        //time /= 1000000;
        //tdtime += time;//dptime
        chrono::steady_clock::time_point ETc2 = chrono::steady_clock::now();
        chrono::duration<double> time_span_single_cycle = chrono::duration_cast<chrono::duration<double> >(ETc2 - STc2);

        cout << "done ripup and reaaaign iteration: " << itnumber;
        cout << "  runtime: " << time_span_single_cycle.count() << " (s)" << endl << endl;
        
        computedelay();
        itnumber ++;
        globaliteration++;
        
    }

    cout << "\ndone RRALA\n";
    computedelay();
    
    //gettimeofday(&e_t, NULL);
    chrono::steady_clock::time_point ETc = chrono::steady_clock::now();
    //ttime = ((e_t.tv_sec*1000000)+e_t.tv_usec)-((s_t.tv_sec*1000000)+s_t.tv_usec);
    //ttime = ttime/1000000;//double(et-st)/CLOCKS_PER_SEC;
    RRA_time_span = chrono::duration_cast<chrono::duration<double> >(ETc - STc);

    //cout << "totalDPtime : " << tdtime << endl;

    //cout << "totalRRAtime : " << ttime << endl;
    cout << "totalRRAtime : " << RRA_time_span.count() << " (s)" << endl;

    cout << "totalitreation : " << itnumber << endl << endl;

};

void CIRCUIT::postop(){
    cout << "start postop:\n";
    vector<SORTEL> nlist;

    //struct timeval s_t2, e_t2;
    //double time;
    //gettimeofday(&s_t2, NULL);
    chrono::steady_clock::time_point STc = chrono::steady_clock::now();

    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        SORTEL Stmp;
        Stmp.eldelay = nets[i].eldelay;
        Stmp.index = i;
        nlist.push_back(Stmp);
    }
    
    //sort(nlist.begin(),nlist.end(),struct_sortel_cmp);
    stable_sort(nlist.begin(),nlist.end(),struct_sortel_cmp);
    

    int limit = nlist.size()/20;
    for(int i=0;i<nlist.size();i++ ){
        NET oldnet = nets[nlist[i].index];
        double od = nets[nlist[i].index].eldelay;
        ripup(nets[nlist[i].index] , nlist[i].index);
        if(i<limit)
            dynamic_program_main(nets[nlist[i].index],nlist[i].index,parallelmode,1,0);
        else
            dynamic_program_main(nets[nlist[i].index],nlist[i].index,0,1,0);
        singlenetdelay(nlist[i].index);
        if(nets[nlist[i].index].eldelay > od){
            ripup(nets[nlist[i].index], nlist[i].index);
            nets[nlist[i].index] = oldnet;
            recoverpathtoedge(nets[nlist[i].index] , nlist[i].index);
        }
    }

    computedelay();
    computeoverflow();

    //gettimeofday(&e_t2, NULL);
    //time = ((e_t2.tv_sec*1000000)+e_t2.tv_usec)-((s_t2.tv_sec*1000000)+s_t2.tv_usec);
    //time /= 1000000;
    //cout << "postop time: " << time << "\ndone postop\n";
    chrono::steady_clock::time_point ETc = chrono::steady_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double> >(ETc - STc);
    cout << "postop time: " << time_span.count() << " (s)\ndone postop\n";

}

void CIRCUIT::greedyrefinement(){

    //struct timeval s_t, e_t;
    //gettimeofday(&s_t, NULL);
    //double ttime;
    chrono::steady_clock::time_point STc = chrono::steady_clock::now();

    int maxindex=0;
    double maxd=0;
    //find wns
    for(int i=0;i< nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        nets[i].needreassign = false;
        if(nets[i].eldelay > maxd){
            maxd = nets[i].eldelay;
            maxindex=i;
        }
    }
    bool flag=true;

	int ccount = 0;
	double mark = 0.0;
	int tagtotal = 1;
	int tagmax = 1;

    int itreation=1;
    //itreation greedy
    while(flag){
        if(weightmode==1){
            for(int i=0;i< nets.size();i++){
                if(!nets[i].havesegment)
                    continue;
                nets[i].updateRTreeweight();
            }
        }
        cout <<"itreation: " << itreation  <<"  maxdelay: " << maxd << endl;
        //Miu=0;
        //Lambda = 1000;

        NET oldnet;
        oldnet = nets[maxindex];//backup net
        ripup(nets[maxindex],maxindex);
        dynamic_program_main(nets[maxindex],maxindex,parallelmode,0,0);


        //computeslack
        singlenetdelay(maxindex);
        //slack is worse
        if(nets[maxindex].eldelay >= maxd || ((maxd-nets[maxindex].eldelay)/maxd)<0.01){
            cout << "break " << nets[maxindex].eldelay << "  " << maxd << endl;
            ripup(nets[maxindex],maxindex);
            nets[maxindex] = oldnet;
            recoverpathtoedge(nets[maxindex],maxindex);
            break;
        }

        //find overflow net
        vector<SORTEL> illnet = checkillegalnet(0);
        cout <<"overflow net number : " <<illnet.size() << endl;
        for(int i=illnet.size()-1;i >0;i--){
            if(illnet[i].index == maxindex)
                continue;
            double oldmaxd = nets[illnet[i].index].eldelay;
            NET oldnet2=nets[illnet[i].index];//backup
            ripup(nets[illnet[i].index],illnet[i].index);
            dynamic_program_main(nets[illnet[i].index],illnet[i].index,parallelmode,1,0);
            singlenetdelay(illnet[i].index);

            if(nets[illnet[i].index].eldelay > maxd){
                ripup(nets[illnet[i].index],illnet[i].index);
                nets[illnet[i].index] = oldnet2;//recovernet
                recoverpathtoedge(nets[illnet[i].index],illnet[i].index);

                oldnet2 = nets[maxindex];
                ripup(nets[maxindex],maxindex);
                dynamic_program_main(nets[maxindex],maxindex,parallelmode,1,0);
                singlenetdelay(maxindex);
                if(nets[maxindex].eldelay > maxd){
                    flag=false;
                    ripup(nets[maxindex],maxindex);
                    nets[maxindex] = oldnet2;
                    recoverpathtoedge(nets[maxindex],maxindex);
                    break;
                }
            }

        }
        
        computedelay();
		int total = 0;
		int max = 0;
		for (int i = 0; i<edge_c.size() - 1; i++){
			for (int j = 0; j<edge_c[i].size(); j++){
				for (int k = 1; k <= 8; k++){
					int ctmp;
					if (k % 2 == 0)
						ctmp = edge_c[i][j].C_H[k];
					else
						ctmp = edge_c[i][j].C_V[k];
					if (ctmp<0){
						ctmp = ctmp *-1 * (minwidth[k] + minspacing[k]);
						total += ctmp;
						if (ctmp > max)
							max = ctmp;
					}
				}
			}
		}
		cout << "total :" << total << "  max:" << max << endl;
		tagtotal = total;
		tagmax = max;
        maxd=0;
        for(int i=0;i< nets.size();i++){
            if(!nets[i].havesegment)
                continue;
            nets[i].needreassign = false;
            if(nets[i].eldelay > maxd){
                maxd = nets[i].eldelay;
                maxindex=i;
            }
        }
        itreation++;
		
		if (mark == maxd) ccount++;
		else { mark = maxd; ccount = 0; }
		if (ccount == 5) { break; }
		if (itreation - 1 > 11 && tagtotal == 0 && tagmax == 0){ break; }
    }
    cout << "\ndone refinement\n";
    computedelay();

    //gettimeofday(&e_t, NULL);
    //ttime = ((e_t.tv_sec*1000000)+e_t.tv_usec)-((s_t.tv_sec*1000000)+s_t.tv_usec);
    //ttime /= 1000000;
    chrono::steady_clock::time_point ETc = chrono::steady_clock::now();
    chrono::duration<double> total_time_span = chrono::duration_cast<chrono::duration<double> >(ETc - STc);

    cout << "totalitreation : " << itreation-1<<endl;
    //cout << "totaltime : " << ttime <<endl<< endl;
    cout << "totaltime : " << total_time_span.count() << " (s)" << endl << endl;

}

void CIRCUIT::singlenetdelay(int index){
    nets[index].eldelay=0;

    for(int j=0;j<nets[index].RTree.size();j++){
        nets[index].RTree[j].child_C=-1;
    }
    for(int i=1;i<nets[index].sinks_marge.size();i++){
        double etmp = nets[index].elmoredelay_assign(nets[index].sinks_marge[i],edge_c,C_V,C_H,minwidth,minspacing);
        nets[index].eldelay += etmp * (1.0/(nets[index].sinks_marge.size()-1));
        nets[index].sinks_marge[i].elmoredelay = etmp;
    }

}

void CIRCUIT::computedelay(){

    double totaled = 0.0;
    double maxd = 0.0;
    vector<SORTEL> netvector;
    for(int i=0;i< nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        nets[i].eldelay = 0.0;

        for(int j=0;j<nets[i].RTree.size();j++){
            nets[i].RTree[j].child_C=-1;
        }
        for(int j=1;j<nets[i].sinks_marge.size();j++){
            double etmp = nets[i].elmoredelay_assign(nets[i].sinks_marge[j],edge_c,C_V,C_H,minwidth,minspacing);
            nets[i].eldelay += etmp*(1.0/(nets[i].sinks_marge.size()-1));
            nets[i].sinks_marge[j].elmoredelay = etmp;
        }
        
        if(nets[i].eldelay>maxd)
            maxd = nets[i].eldelay;
        
        totaled += nets[i].eldelay;

        if(globaliteration==0){
            SORTEL Stmp;
            Stmp.eldelay = nets[i].eldelay;
            Stmp.index = i;
            netvector.push_back(Stmp);       
        }
    }
     
    if(globaliteration==0){
        //sort(netvector.begin(),netvector.end(),struct_sortel_cmp);
        stable_sort(netvector.begin(),netvector.end(),struct_sortel_cmp);
        int total = netvector.size();
        int unit = total/10000;
        double reduce=0.02;
        int times=0;
        for(int i=0;i<netvector.size();i++){
            nets[netvector[i].index].fatwirelevelthreshold = \
            (parallelLevelthreshold - double(reduce*times)>0 ? parallelLevelthreshold - double(reduce*times) : 0);
            nets[netvector[i].index].computelevelthreshold();
            if(i%unit==unit-1)
                times++;
            if(i==0)
                cout << nets[netvector[i].index].fatwirelevelthreshold << " " << nets[netvector[i].index].maxlevel << endl;
        }
    }
    cout << "totalelmoredelay : " << totaled <<"  " << maxd << endl;
}

void NET::updateRTreeweight(){
    for(int i=1;i<sinks_marge.size();i++){
        passupdate(sinks_marge[i].RTree_index,sinks_marge[i].elmoredelay/(sinks_marge.size()-1));
    }
}

void NET::nodelevel(int index){
    if(index==0)
        RTree[index].level=1;
    else{
        RTree[index].level =  RTree[ RTree[index].parent_index ].level+1;
    }
    if(RTree[index].level>maxlevel)
        maxlevel=RTree[index].level;
    for(int i=0;i< RTree[index].child_index.size();i++)
        nodelevel(RTree[index].child_index[i]);
}

void NET::computelevelthreshold() {
    for(int i=0;i<RTree.size();i++){
        if( double(RTree[i].level/double(maxlevel*1.0)) > fatwirelevelthreshold)
            RTree[i].underlevelthreshould = false;
        else
            RTree[i].underlevelthreshould = true;
    } 
}

void NET::passupdate(int index,double slack){

    RTree[index].weight += slack/100;// / RTree[index].level;
    if(RTree[index].parent_index != -1)
        passupdate(RTree[index].parent_index,slack);

}

static int calcIndicators(
    vector<TREE_NODE>& _rt, int _rt_index, // base
    int _rt_node_height, map<int, int> &_rt_node_height_dis, // distribution of height of _rt's node 
    map<int, int> &_rt_node_re_height_dis // distribution of re height of _rt's node
) {
    _rt_node_height_dis[_rt_node_height]++;
    int _rt_node_re_height = -1;
    for(int i = 0; i < _rt[_rt_index].child_index.size(); i++) {
       int _rt_son_node_re_height =  calcIndicators(_rt, _rt[_rt_index].child_index[i], _rt_node_height+1, _rt_node_height_dis, _rt_node_re_height_dis);
       _rt_node_re_height = max(_rt_node_re_height, _rt_son_node_re_height); 
    }
    _rt_node_re_height++;
    _rt_node_re_height_dis[_rt_node_re_height]++;
    return _rt_node_re_height;
}

void CIRCUIT::dynamic_program_main(NET& _net,int netindex,int mode,int greedy,double q){
    // static int dynamic_program_main_access_number = 0;
    // dynamic_program_main_access_number++;
    // static ofstream out_rt_max_height_dis("analyze/rt_max_height_dis.txt");
    // static ofstream out_rt_node_height_dis("analyze/rt_node_height_dis.txt");
    // static ofstream out_rt_max_re_height_dis("analyze/rt_max_re_height_dis.txt");
    // static ofstream out_rt_node_re_height_dis("analyze/rt_node_re_height_dis.txt");
    // // 计算每一次与树相关的指标
    // map<int, int> _rt_node_height_dis;
    // map<int, int> _rt_node_re_height_dis;
    // calcIndicators(_net.RTree, 0, 0, _rt_node_height_dis, _rt_node_re_height_dis); 
    // out_rt_max_height_dis << dynamic_program_main_access_number << "," << _rt_node_height_dis.rbegin()->first << endl;
    // for(auto &it: _rt_node_height_dis) {
    //     out_rt_node_height_dis << dynamic_program_main_access_number << "," << it.first << "," << it.second << endl;
    // }
    // out_rt_max_re_height_dis << dynamic_program_main_access_number << "," << _rt_node_height_dis.rbegin()->first << endl;
    // for(auto &it: _rt_node_re_height_dis) {
    //     out_rt_node_re_height_dis << dynamic_program_main_access_number << "," << it.first << "," << it.second << endl;
    // }
    // dynamic_program(_net.RTree,0, record_dp,_net.resultsegments,netindex,mode,greedy,_net.fatwirelevelthreshold);
    vector<RECORD_DP>  record_dp(_net.RTree.size());
    vector<TREE_NODE>& _rt = _net.RTree;
    vector<SEGMENT>& resultsegments = _net.resultsegments;
    int _rt_i = 0;
    for(auto &record: record_dp) {
        record.done = false;
        if(_rt[_rt_i].child_index.size()) {
            record.child_combine.resize(8);
            for(auto &child_combine: record.child_combine) {
                child_combine.resize( _rt[_rt_i].child_index.size()  );
            }
        }
        ++_rt_i;
    }
    // 原来的初始化
    // for(int i=0;i<record_dp.size();i++){
    //     record_dp[i].done = false;
    //     if (_rt[i].child_index.size()) {
    //         record_dp[i].child_combine.resize(8);
    //         for(int j=0;j < 8; j++)
    //             record_dp[i].child_combine[j].resize( _net.RTree[i].child_index.size()  );
    //     }
    // }
    // BFS
    // if (record_dp.size() <= 50) {
    //     vector<int> index_list(_rt.size());
    //     int l=0, r = 0;
    //     while(l <= r) {
    //         for(auto &it : _rt[index_list[l]].child_index) {
    //             index_list[++r] = it;
    //         }
    //         l++;
    //     }
    //     for(int i = r; i>=0; i--) {
    //         dynamic_program(_rt, index_list[i], record_dp, resultsegments, netindex, mode, greedy, q);
    //     }
    if (_net.parallel == -1) {
        int leaf_number = 0;
        for(auto &it: _rt) leaf_number += int(it.child_index.size() == 0);
        // _net.parallel = min(12, (leaf_number+19)/20);
        _net.parallel = 1; //min(12, (leaf_number+9)/10);
    }
    // } else {
        tf::Taskflow taskflow;
        buildDependency(_rt, record_dp, _net.resultsegments, netindex, mode, greedy, _net.fatwirelevelthreshold, taskflow);
        //自适应线程数
        tf::Executor executor(_net.parallel);
        executor.run(taskflow).wait();
        //taskflow.dump(std::cout);
        //exit(0);
    // }
}

void CIRCUIT::buildDependency(vector<TREE_NODE>& _rt, vector<RECORD_DP>& record_dp, 
    vector<SEGMENT>& resultsegments,int netindex,int mode,int greedy, double q, 
    tf::Taskflow &taskflow)
{
    //利用队列构造任务依赖

    //树节点队列
    vector <int> nodes(_rt.size());
    //任务节点队列
    vector <tf::Task> tasks(_rt.size());
    //队列头尾指针
    int l = 0, r = 0;
    //将根节点0加入队列
    nodes[0] = 0;
    //将根节点任务加入队列
    tasks[0] = taskflow.emplace([=, &_rt, &record_dp, &resultsegments]()
    {
        dynamic_program(_rt, 0, record_dp, resultsegments, netindex, mode, greedy, q);
    });
    while (l<=r)
    {
        for (auto &ch:_rt[nodes[l]].child_index)
        {
            ++r;
            nodes[r] = ch;
            tasks[r] = taskflow.emplace([=, &_rt, &record_dp, &resultsegments]()
            {
                dynamic_program(_rt, ch, record_dp, resultsegments, netindex, mode, greedy, 0);
            });
            //父节点在子节点之后执行
            tasks[l].succeed(tasks[r]);
        }
        ++l;
    }
}

void CIRCUIT::dynamic_program(vector<TREE_NODE>& _rt,int index, vector<RECORD_DP>& record_dp ,vector<SEGMENT>& resultsegments,int netindex,int mode,int greedy, double q){
    clock_t start,end;
    double time;
    if(_rt[index].child_index.size()==0){
        int p_index = _rt[index].parent_index;

        if( _rt[p_index].x == _rt[index].x ){//vertical layer 1 3 5 7
            int y,x;
            int parallelnumber;
            x = _rt[index].x/32;
            if( _rt[p_index].y < _rt[index].y){
                y = _rt[p_index].y/40;
            }
            else{
                y = _rt[index].y/40;
            }
            
            if(edge_c[x][y].total_c_v >0){
                parallelnumber = 1;
                record_dp[index].parallelwirenumber = 1;
            }
            else{
                parallelnumber = 0;
                record_dp[index].parallelwirenumber = 0;
            }
           

            //unable using parallel wire 
            if(mode==0 || ( !_rt[index].underlevelthreshould && !greedy)){
                parallelnumber = 0;
                record_dp[index].parallelwirenumber = 0;
            }
      
            int toplayerindex = (parallelnumber+1)*4-1;
            if(edge_c[x][y].total_c_v ==1 && mode==1){
                toplayerindex = 5;
            }
            record_dp[index].toplayerindex = toplayerindex;
            for(int i=0;i<=toplayerindex;i++){
                int presentlayer = (2*i+1)%8;
                int tableindex = (2*i+1)-1;
                double via_c_tmp = 0.0;
                double via_delay_tmp = 0.0;

                if(greedy){
                    int wirenumber = (i/4)+1;
                    if(presentlayer>3 && wirenumber >1)
                        wirenumber++;
                    if(edge_c[x][y].C_V[presentlayer] - wirenumber < 0){
                        record_dp[index].cost[i] = -1; //std::numeric_limits<double>::max();
                        record_dp[index].C[i] = -2;//std::numeric_limits<double>::max();
                        continue;
                    }
                        
                }
                //find top and low pin to compute via delay
                int top,low;
                top = _rt[index].pinlayer[0];
                low = top;
                if(_rt[index].pinlayer[0]!=-1){

                    for(int j=0;j<_rt[index].pinlayer.size();j++){
                        if(_rt[index].pinlayer[j]> top ){
                            top = _rt[index].pinlayer[j];
                        }
                        else if(_rt[index].pinlayer[j]<low){
                            low = _rt[index].pinlayer[j];
                        }
                    }

                    //lowviadalay
                    if( low < presentlayer){
                        for(int j= low ;j<presentlayer ;j++){
                            for(int k=0;k<_rt[index].pinlayer.size();k++){
                                if(_rt[index].pinlayer[k] == j)
                                    via_c_tmp += load_C;
                            }
                            via_delay_tmp += (via_c_tmp + via_C[j]/2)*via_R[j];
                            via_c_tmp += via_C[j];
                        }
                    }
                    //highviadelay
                    via_c_tmp=0;
                    if(presentlayer < top){
                        for(int j= top ; j > presentlayer ;j--){
                            for(int k=0;k<_rt[index].pinlayer.size();k++){
                                if(_rt[index].pinlayer[k] == j)
                                    via_c_tmp += load_C;
                            }
                            if (j-1>=0)
                            {
                                via_delay_tmp += (via_c_tmp + via_C[j-1]/2)*via_R[j-1];
                                via_c_tmp += via_C[j-1];
                            }
                        }
                    }
                    
                    //total capacity
                    via_c_tmp=0;
                    if(presentlayer > top)
                        top = presentlayer;
                    if(presentlayer < low)
                        low = presentlayer;
                    for(int j=low;j<top;j++){
                        if (j-1>=0)
                        {
                            via_c_tmp += via_C[j-1];
                        }
                    }
                    if(_rt[index].pinlayer[0]!=-1)
                        for(int j=0;j<_rt[index].pinlayer.size();j++)
                            via_c_tmp += load_C;
                }
                else{
                    via_c_tmp = 0.0;
                    low=0;
                    top=0;
                }
                //compute total dalay and capacity
                double dtmp;
                int fatornot = i/4;
                int wirenumber = fatornot+1;
                //use NDR
                if(presentlayer >3 && fatornot==1)
                    wirenumber++;
                int dindex = edge_c[x][y].C_V[presentlayer]-wirenumber;
                if(dindex <0)
                    dindex=0;
                else if(dindex > edge_c[x][y].aveg_c_v)
                    dindex = edge_c[x][y].aveg_c_v;
                dtmp = (20*RCtable[presentlayer].UC[dindex][fatornot] + via_c_tmp) * 40 * RCtable[presentlayer].UR[dindex][fatornot];
                double ctmp;
                ctmp = 40*RCtable[presentlayer].UC[dindex][fatornot]; 
                

                record_dp[index].delay[i] = dtmp+via_delay_tmp;
                record_dp[index].C[i] = ctmp+via_c_tmp;
                

                //via number of cost function
                int vianumber = abs( top-low);

                //congestion of cost function
				double congestion_cost = 1.0*(edge_cap[x][y].C_V[presentlayer] - edge_c[x][y].C_V[presentlayer]) / edge_cap[x][y].C_V[presentlayer] + 1.0*wirenumber / edge_cap[x][y].C_V[presentlayer] + 1.0*(maxcap[presentlayer] - edge_cap[x][y].C_V[presentlayer]) / maxcap[presentlayer];
                if(edge_c[x][y].C_V[presentlayer]-wirenumber < 0)
                    congestion_cost += k1*(0-(edge_c[x][y].C_V[presentlayer]-wirenumber))*(minwidth[presentlayer]+minspacing[presentlayer])* (edge_c[x][y].historic[presentlayer]);

                record_dp[index].cost[i] =  _rt[index].weight*Lambda*record_dp[index].delay[i] + Beta*vianumber + Miu*congestion_cost;
            }
        }
        else{//horizontal layer 2 4 6 8
            int x,y;
            y = _rt[index].y/40;
            if( _rt[index].x < _rt[p_index].x)
                x = _rt[index].x/32;
            else
                x = _rt[p_index].x/32;

            int parallelnumber;
            if(edge_c[x][y].total_c_h > 0){
                parallelnumber = 1;
                record_dp[index].parallelwirenumber = 1;
            }
            else{
                parallelnumber = 0;
                record_dp[index].parallelwirenumber = 0;
            }
            
            if(mode==0 || ( !_rt[index].underlevelthreshould && !greedy)){
                parallelnumber = 0;
                record_dp[index].parallelwirenumber=0;
            }

            int toplayerindex = (parallelnumber+1)*4-1;
            if(edge_c[x][y].total_c_h==1&&mode==1)
                toplayerindex = 4;
            
            record_dp[index].toplayerindex = toplayerindex;

            for(int i=0;i<=toplayerindex;i++){
                //compute via number
                int presentlayer = ((i+1)*2)%8;
                int tableindex = ((i+1)*2)-1;
                if( ((i+1)*2)%8==0)
                    presentlayer =8;

                double via_c_tmp=0.0;
                double via_delay_tmp=0.0;

                if(greedy){
                    int wirenumber = (i/4)+1;
                    if(presentlayer >3 && wirenumber>1)
                        wirenumber++;
                    if(edge_c[x][y].C_H[presentlayer]- wirenumber < 0){
                        record_dp[index].cost[i] = -1;
                        record_dp[index].C[i] = -2;
                        continue;
                    }
                }
                int top,low;
                top = _rt[index].pinlayer[0];
                low = top;

                if(_rt[index].pinlayer[0]!=-1){

                    for(int j=0;j<_rt[index].pinlayer.size();j++){
                        if( _rt[index].pinlayer[j] > top)
                            top = _rt[index].pinlayer[j];
                        if( _rt[index].pinlayer[j] < low)
                            low = _rt[index].pinlayer[j];
                    }
                    
                    //lowviadelay
                    if( low < presentlayer){
                        for(int j= low ;j<presentlayer ;j++){
                            for(int k=0;k< _rt[index].pinlayer.size();k++){
                                if(j== _rt[index].pinlayer[k])
                                    via_c_tmp += load_C;
                            }
                            via_delay_tmp += (via_c_tmp + via_C[j]/2)*via_R[j];
                            via_c_tmp += via_C[j];
                        }
                    }
                    
                    //highviadelay
                    via_c_tmp=0.0;
                    if(top>presentlayer){
                        for(int j= top ; j > presentlayer ;j--){
                            for(int k=0;k< _rt[index].pinlayer.size();k++){
                                if(j== _rt[index].pinlayer[k])
                                    via_c_tmp += load_C;
                            }
                            if (j-1>=0)
                            {
                                via_delay_tmp += (via_c_tmp + via_C[j-1]/2)*via_R[j-1];
                                via_c_tmp += via_C[j-1];
                            }
                        }
                    }

                    //total capacity
                    via_c_tmp=0;
                    if(presentlayer > top)
                         top = presentlayer;
                    if(presentlayer < low)
                         low = presentlayer;
                    for(int j=low;j<top;j++){
                        if (j-1>=0)
                        {
                            via_c_tmp += via_C[j-1];
                        }
                    }
                    if(_rt[index].pinlayer[0]!=-1)
                        for(int j=0;j<_rt[index].pinlayer.size();j++)
                            via_c_tmp += load_C;
                }
                else{
                    via_c_tmp = 0.0;
                    top =0;
                    low =0;
                }
                //need C_value for different i
                double dtmp;
                int fatornot = i/4;
                int wirenumber = fatornot+1;
                //use NDR
                if(fatornot==1 && presentlayer>3)
                    wirenumber++;
                int dindex = edge_c[x][y].C_H[presentlayer]-wirenumber;
                if(dindex <0)
                    dindex =0;
                else if(dindex > edge_c[x][y].aveg_c_h)
                   dindex = edge_c[x][y].aveg_c_h; 
                dtmp = ( 16*RCtable[presentlayer].UC[dindex][fatornot] + via_c_tmp) * 32*RCtable[presentlayer].UR[dindex][fatornot];
                double ctmp;
                ctmp = 32*RCtable[presentlayer].UC[dindex][fatornot];

                record_dp[index].delay[i] = dtmp+ via_delay_tmp;//delay only
                record_dp[index].C[i] = dtmp+via_c_tmp;
                

                //via number of cost function
                int vianumber = abs( top-low);
                //congestion of cost function
				double congestion_cost = 1.0*(edge_cap[x][y].C_H[presentlayer] - edge_c[x][y].C_H[presentlayer]) / edge_cap[x][y].C_H[presentlayer] + 1.0*wirenumber / edge_cap[x][y].C_H[presentlayer] + 1.0*(maxcap[presentlayer] - edge_cap[x][y].C_H[presentlayer]) / maxcap[presentlayer];
                if(edge_c[x][y].C_H[presentlayer]-wirenumber < 0)
                    congestion_cost += k1*(0-(edge_c[x][y].C_H[presentlayer]-wirenumber))*(minwidth[presentlayer]+minspacing[presentlayer])* (edge_c[x][y].historic[presentlayer]);
                record_dp[index].cost[i] =  _rt[index].weight*Lambda*record_dp[index].delay[i] + Beta*vianumber + Miu*congestion_cost;
                if(congestion_cost < 0)
                   cout << "fuck"; 
            }

        }
        record_dp[index].done = true;
    }
    else if( _rt[index].parent_index == -1){//routing tree root

        /* 由于已经建立了依赖图，故动规函数不再需要递归
        *
        for(int i=0;i< _rt[index].child_index.size();i++){
            if( !record_dp[_rt[index].child_index[i]].done ){
                dynamic_program( _rt ,_rt[index].child_index[i], record_dp, resultsegments,netindex,mode,greedy,0);

            }
        }
        */

        findallcombinefor(_rt[index].child_index.size(),record_dp, _rt, index,_rt[index].pinlayer[0],0);
        record_dp[index].done = true;
        vector<int> dp_layer;
        dp_layer.reserve( _rt.size());
        dp_layer[0] = _rt[index].pinlayer[0];
        if(record_dp[0].cost[_rt[index].pinlayer[0]]<0)
            cout << netindex << "  " << record_dp[0].cost[_rt[index].pinlayer[0]] << endl;
        createsolution(dp_layer,0, record_dp, resultsegments, _rt,netindex);
    }
    else{

        int p_index = _rt[index].parent_index;
        /* 由于已经建立了依赖图，故动规函数不再需要递归
        *
        for(int i=0;i< _rt[index].child_index.size();i++){
            if( !record_dp[_rt[index].child_index[i]].done ){
                dynamic_program( _rt ,_rt[index].child_index[i], record_dp, resultsegments,netindex,mode,greedy,0);
            }
        
        }
        */
        int x,y,parallelnumber,VORH;
        int toplayerindex = 0;
        if( _rt[index].x == _rt[p_index].x){//vertical
            x = _rt[index].x/32;
            if( _rt[index].y< _rt[p_index].y)
                y = _rt[index].y/40;
            else
                y = _rt[p_index].y/40;

            if(edge_c[x][y].total_c_v> 0){
                parallelnumber = 1;
                record_dp[index].parallelwirenumber = 1;
            }
            else{
                parallelnumber = 0;
                record_dp[index].parallelwirenumber = 0;
            }

            if(mode==0 || ( !_rt[index].underlevelthreshould && !greedy) ){
                record_dp[index].parallelwirenumber = 0;
                parallelnumber=0;
            }
            VORH=1;
            toplayerindex = (parallelnumber+1)*4-1;
            if(edge_c[x][y].total_c_v ==1&&mode==1)
                toplayerindex = 5;
        }
        else{
            y = _rt[index].y/40;
            if(_rt[index].x < _rt[p_index].x)
                x = _rt[index].x/32;
            else
                x = _rt[p_index].x/32;

            if(edge_c[x][y].total_c_h>0){
                parallelnumber = 1;
                record_dp[index].parallelwirenumber = 1;
            }
            else{
                parallelnumber = 0;
                record_dp[index].parallelwirenumber = 0;
            }
            
            if(mode==0 || ( !_rt[index].underlevelthreshould && !greedy) ){
                record_dp[index].parallelwirenumber = 0;
                parallelnumber=0;
            }
            
            VORH=2;
            toplayerindex = (parallelnumber+1)*4-1;
            if(edge_c[x][y].total_c_h==1&&mode==1)
                toplayerindex = 4;
        }
        
        record_dp[index].toplayerindex = toplayerindex;

        for(int j=0;j<=toplayerindex;j++){//for different layer and parallel wire with corrent node
            if(greedy){
                if(VORH==1){
                    int presentlayer = (2*j+1)%8;
                    int wirenumber = ((j/4)+1);
                    if (presentlayer > 3&& wirenumber>1)
                        wirenumber++;
                    if(edge_c[x][y].C_V[presentlayer]-wirenumber < 0){
                        record_dp[index].cost[j] = -1;
                        record_dp[index].C[j] = -2; 
                        continue;
                    }
                }
                else if(VORH==2){
                    int presentlayer = ((j+1)*2)%8;
                    int wirenumber = ((j/4)+1);
                    if( ((j+1)*2)%8==0)
                        presentlayer =8;
                    if (presentlayer >3&&wirenumber>1)
                       wirenumber++; 
                    if(edge_c[x][y].C_H[presentlayer]- wirenumber < 0){
                        record_dp[index].cost[j] = -1;//std::numeric_limits<double>::max();
                        record_dp[index].C[j] = -2;//std::numeric_limits<double>::max();
                        continue;
                    }
                }
            }
            findallcombinefor( _rt[index].child_index.size(),record_dp, _rt, index,j,1);
        }

        record_dp[index].done = true;
    }

}

void CIRCUIT::findallcombinefor( int n,vector<RECORD_DP>& record_dp, vector<TREE_NODE>& rt,int treeindex,int presentlayerindex,int mode){//mode 0 = routing tree root / mode 1 = others
    // n = child size;
    // presentlayerinde is the record_dp[treeindex] solution index(ex cost[x] C[x] combine[x])
    int presentlayer,tableindex;
    vector<int> child_parallel;
    for(int i=0;i<rt[treeindex].child_index.size();i++){
        child_parallel.push_back(record_dp[rt[treeindex].child_index[i]].toplayerindex+1);
    }
    if( rt[treeindex].parent_index==-1)
        presentlayer = rt[treeindex].pinlayer[0];
    else{
        if(rt[treeindex].x == rt[rt[treeindex].parent_index].x){
            presentlayer = (presentlayerindex*2+1)%8;
            tableindex = (presentlayerindex*2+1)-1;
        }
        else{
            presentlayer = ((presentlayerindex+1)*2)%8;
            if(presentlayer==0)
                presentlayer = 8;
            tableindex = ((presentlayerindex+1)*2)-1;
        }
    }

    vector<int> min_combine;
    min_combine.reserve(n);
    int base = 8;
    int *combine = new int [n];
    for(int i=0;i<n;i++)
        combine[i] = 0;

	for (int i = 0; i<n; i++)
		min_combine[i] = 0;

    int l = 0;
    double min_cost = numeric_limits<double>::max();
    double min_delay=0.0;
    double min_C=0.0;
    double ast_tmp=0.0;
    while(l<n){
        for(int i=0;i<child_parallel[0];i++){
            combine[0] = i;
            //get one combine
            
            if(mode==0){
                double highviadelay=0.0;
                double lowviadelay=0.0;
                double delay_tmp=0.0;
                double cost_tmp=0.0;
                int top,low;

                //find top and low pin 
                top = rt[treeindex].pinlayer[0];
                low = rt[treeindex].pinlayer[0];
                for(int j=1;j<rt[treeindex].pinlayer.size();j++){
                    if(top<rt[treeindex].pinlayer[j])
                        top = rt[treeindex].pinlayer[j];
                    if(low>rt[treeindex].pinlayer[j])
                        low = rt[treeindex].pinlayer[j];
                }
                bool ill = false;
                for(int j=0;j<n;j++){
                    if(record_dp[rt[treeindex].child_index[j]].cost[combine[j]]==-1){
                        ill=true;
                        break;
                    }
                    cost_tmp += record_dp[ rt[treeindex].child_index[j] ].cost[combine[j]];
                    int layertmp;
                    if( rt[treeindex].x == rt[rt[treeindex].child_index[j]].x){
                        layertmp = (combine[j]*2+1)%8;
                    }
                    else{
                        layertmp = ((combine[j]+1)*2)%8;
                        if(layertmp==0)
                            layertmp=8;
                    }
                    if(layertmp>top)
                        top=layertmp;
                    if(layertmp<low)
                        low=layertmp;
                }
                
                if(ill){
                    continue;
                }


                //highviadelay
                double ctmp=0.0;
                for(int j=top;j>rt[treeindex].pinlayer[0];j--){
                    for(int k=0;k<rt[treeindex].child_index.size();k++){
                        int tmplayer;
                        if(rt[treeindex].x == rt[rt[treeindex].child_index[k]].x){
                            tmplayer = (combine[k]*2+1)%8;
                        }
                        else{
                            tmplayer = ((combine[k]+1)*2)%8;
                            if(tmplayer==0)
                                tmplayer=8;
                        }

                        if(j==tmplayer){
                            ctmp += record_dp[ rt[treeindex].child_index[k]  ].C[combine[k]];
                        }
                    }
                    for(int k=1;k<rt[treeindex].pinlayer.size();k++){
                        if(rt[treeindex].pinlayer[k]==j)
                            ctmp += load_C;
                    }
                    if (j-1>=0)
                    {
                        highviadelay += (via_C[j-1]/2+ctmp)*via_R[j-1];
                        ctmp += via_C[j-1];
                    }
                }

                //lowviadelay
                ctmp=0.0;
                for(int j=low;j<rt[treeindex].pinlayer[0];j++){
                    for(int k=0;k<rt[treeindex].child_index.size();k++){
                        int tmplayer;
                        if(rt[treeindex].x == rt[rt[treeindex].child_index[k]].x){
                            tmplayer = (combine[k]*2+1)%8;
                        }
                        else{
                            tmplayer = ((combine[k]+1)*2)%8;
                            if(tmplayer==0)
                                tmplayer=8;
                        }

                        if(j==tmplayer){
                            ctmp += record_dp[ rt[treeindex].child_index[k]  ].C[combine[k]];
                        }
                    }
                    for(int k=1;k<rt[treeindex].pinlayer.size();k++){
                        if(rt[treeindex].pinlayer[k]==j)
                            ctmp += load_C;
                    }

                    lowviadelay += (via_C[j]/2+ctmp)*via_R[j];
                    ctmp += via_C[j];

                }

                //driverdelay
                double driverdelay=0.0;
                ctmp=0.0;
                if(top < rt[treeindex].pinlayer[0])
                    top = rt[treeindex].pinlayer[0];
                if(low > rt[treeindex].pinlayer[0])
                    low = rt[treeindex].pinlayer[0];

                for(int j=low;j<top;j++){
                    ctmp += via_C[j];
                }
                for(int j=0;j<rt[treeindex].child_index.size();j++){
                    ctmp += record_dp[ rt[treeindex].child_index[j]].C[combine[j]];
                }
                for(int j=1;j<rt[treeindex].pinlayer.size();j++)
                    ctmp += load_C;
                driverdelay = driver_R*ctmp;

                //total delay
                delay_tmp = driverdelay+highviadelay+lowviadelay;
                
                int vianumber = top-low;
                cost_tmp += Lambda*delay_tmp + Beta*vianumber; 

                if(cost_tmp<min_cost){
                    min_cost = cost_tmp;
                    record_dp[treeindex].cost[presentlayerindex] = min_cost;
                    for(int j=0;j<n;j++)
                        min_combine[j] = combine[j];
                }

                continue;     
            }
            
            //find max and min layer of child's layer and compute the highviadelay, lowviadelay,wiredelay
            double cost_tmp=0.0;
            double icost_tmp=0.0;
            vector<int> layervector;
            vector<int> indexvector;
            double highviadelay = 0.0;
            double lowviadelay = 0.0;
            double wiredelay = 0.0;
            double C_tmp=0.0;
            double delay_tmp=0.0;
            for(int j=0;j<rt[treeindex].child_index.size();j++){
                int layertmp;
                if( rt[treeindex].x == rt[rt[treeindex].child_index[j]].x  ){//vertical 1 3 5 7
                    layertmp = (combine[j]*2+1)%8;
                }
                else{//horizontal
                    layertmp = ((combine[j]+1)*2)%8;
                    if( layertmp ==0)
                        layertmp=8;
                }
                layervector.push_back(layertmp);
                indexvector.push_back(j);
            }

            //sort layer
            for(int j=0;j<n-1;j++){
                int max=layervector[j];
                int maxindex=j;
                for(int k=j+1;k<n;k++){
                    if(layervector[k] > max){
                        max = layervector[k];
                        maxindex = k;
                    }
                }
                
                if(maxindex!=n-1){
                    int tmp= layervector[maxindex];
                    layervector[maxindex] =  layervector[n-1];
                    layervector[n-1] = tmp;

                    tmp = indexvector[maxindex];
                    indexvector[maxindex] = indexvector[n-1];
                    indexvector[n-1] = tmp;
                }
            }
            
            int top=layervector[n-1];
            int low=layervector[0];
            int retmp = top;
            
            if(rt[treeindex].pinlayer[0]!=-1){
                for(int j=0;j<rt[treeindex].pinlayer.size();j++){
                    if(rt[treeindex].pinlayer[j] > top)
                        top = rt[treeindex].pinlayer[j];
                    else if(rt[treeindex].pinlayer[j] < low)
                        low = rt[treeindex].pinlayer[j];
                }
            }
            int retmp2=top;
            if(presentlayer > top)
               top = presentlayer;
            if(presentlayer < low)
               low = presentlayer; 
            double via_c_tmp=0.0;
            double via_delay_tmp=0.0;
            int index_tmp=n-1;//sorted layer vector index
            //highviadelay (top to current layer)
            for(int j=top;j > presentlayer;j--){
                while(layervector[index_tmp]==j&& index_tmp >= 0){
                    int child_index = rt[treeindex].child_index[indexvector[index_tmp]];
                    via_c_tmp += record_dp[ child_index ].C[ combine[ indexvector[index_tmp]] ];
                    index_tmp--;
                    
                    if(index_tmp < 0){
                        break;
                    }
                }
                for(int k=0;k<rt[treeindex].pinlayer.size();k++){
                    if(rt[treeindex].pinlayer[k] == j)
                        via_c_tmp += load_C;
                }
                if (j-1>=0)
                {
                    via_delay_tmp += ((via_C[j-1]/2)+via_c_tmp)*via_R[j-1];
                    via_c_tmp += via_C[j-1];
                }
            }
            highviadelay = via_delay_tmp;

            //lowviadelay (low to current layer)
            index_tmp=0;
            via_c_tmp=0.0;
            via_delay_tmp = 0.0;
            for(int j=low;j<presentlayer;j++){
                while(layervector[index_tmp]==j&& index_tmp < n){
                    int child_index = rt[treeindex].child_index[indexvector[index_tmp]];
                    via_c_tmp += record_dp[ child_index ].C[ combine[ indexvector[index_tmp]] ];
                    index_tmp++;
                    if(index_tmp > n-1)
                        break;

                }
                for(int k=0;k<rt[treeindex].pinlayer.size();k++){
                    if(rt[treeindex].pinlayer[k]==j)
                        via_c_tmp += load_C;
                }
                via_delay_tmp += (via_C[j]/2 + via_c_tmp)*via_R[j];
                via_c_tmp += via_C[j];
                
            }
            lowviadelay = via_delay_tmp;


            //cost_tmp is record the child cost of child delay with this combine.  C_tmp = total child C with this combine and wire C
            for(int j=low;j<top;j++)
                C_tmp += via_C[j];

            if(rt[treeindex].pinlayer[0] !=-1)
                C_tmp += load_C*rt[treeindex].pinlayer.size();


            bool ill=false;

            for(int j=0;j<rt[treeindex].child_index.size();j++){
                int child_index = rt[treeindex].child_index[j];
                if(record_dp[ child_index  ].cost[ combine[j] ]==-1){
                    ill=true;
                    break;
                }
                C_tmp += record_dp[ child_index  ].C[ combine[j] ];
                cost_tmp += record_dp[ child_index  ].cost[ combine[j] ];
            }


            if(ill)
                continue;

            int xindex,yindex;
            int fatornot = presentlayerindex/4;
            int wirenumber = fatornot+1;
            //use NDR
            if(fatornot==1 && presentlayer >3)
                wirenumber++;
            if(rt[treeindex].x == rt[rt[treeindex].parent_index].x){//vertical
                xindex = rt[treeindex].x/32;
                if(rt[treeindex].y > rt[ rt[treeindex].parent_index ].y){
                    yindex = rt[ rt[treeindex].parent_index ].y/40;
                }
                else
                    yindex = rt[treeindex].y/40;

                int dindex = edge_c[xindex][yindex].C_V[presentlayer]-wirenumber;
                if(dindex<0)
                    dindex =0;
                else if(dindex > edge_c[xindex][yindex].aveg_c_v)
                    dindex = edge_c[xindex][yindex].aveg_c_v;
                wiredelay = ( 20*RCtable[presentlayer].UC[dindex][fatornot]  + C_tmp) * 40*RCtable[presentlayer].UR[dindex][fatornot];
                C_tmp += 40*RCtable[presentlayer].UC[dindex][fatornot];
            }
            else if(rt[treeindex].y == rt[rt[treeindex].parent_index].y){//horziontal
                yindex = rt[treeindex].y/40;
                if(rt[treeindex].x > rt[ rt[treeindex].parent_index ].x){
                    xindex = rt[ rt[treeindex].parent_index ].x/32;
                }
                else{
                    xindex = rt[treeindex].x/32;
                }
                int dindex = edge_c[xindex][yindex].C_H[presentlayer]-wirenumber;
                if(dindex<0)
                    dindex=0;
                else if(dindex > edge_c[xindex][yindex].aveg_c_h)
                    dindex = edge_c[xindex][yindex].aveg_c_h;
                wiredelay = ( 16*RCtable[presentlayer].UC[dindex][fatornot] + C_tmp) * 32* RCtable[presentlayer].UR[dindex][fatornot];
                C_tmp += 32*RCtable[presentlayer].UC[dindex][fatornot];

            }
            else{
                cout << "bug!!: in layer_assigment.cpp (findallcombinefor: compute wiredelay, wiresegment is via)\n";
                exit(1);
            }


            //delay_tmp is record the delay(high/low viadelay and wire delay) with this layer and child delay with this combine
            delay_tmp = highviadelay + lowviadelay + wiredelay;


            //via number of cost function
            int vianumber = abs( top - low);

            //congestion of cost function
            double congestion_cost=0.0;
            if( rt[treeindex].x == rt[ rt[treeindex].parent_index ].x){
				congestion_cost = 1.0*(edge_cap[xindex][yindex].C_V[presentlayer] - edge_c[xindex][yindex].C_V[presentlayer]) / edge_cap[xindex][yindex].C_V[presentlayer] + 1.0*wirenumber / edge_cap[xindex][yindex].C_V[presentlayer] + 1.0*(maxcap[presentlayer] - edge_cap[xindex][yindex].C_V[presentlayer]) / maxcap[presentlayer];
                if(edge_c[xindex][yindex].C_V[presentlayer]-wirenumber < 0)
                    congestion_cost += k1*(0-(edge_c[xindex][yindex].C_V[presentlayer]-wirenumber))*(minwidth[presentlayer]+minspacing[presentlayer])*(edge_c[xindex][yindex].historic[presentlayer]);
            }
            else{
                //congestion_cost = (1.0+k1/(1.0+exp(k2*(edge_c[xindex][yindex].C_H[presentlayer] - wirenumber) ) ) ) * (edge_c[xindex][yindex].historic[presentlayer]);
				congestion_cost = 1.0*(edge_cap[xindex][yindex].C_H[presentlayer] - edge_c[xindex][yindex].C_H[presentlayer]) / edge_cap[xindex][yindex].C_H[presentlayer] + 1.0*wirenumber / edge_cap[xindex][yindex].C_H[presentlayer] + 1.0*(maxcap[presentlayer] - edge_cap[xindex][yindex].C_H[presentlayer]) / maxcap[presentlayer];
                if(edge_c[xindex][yindex].C_H[presentlayer]-wirenumber < 0)
                    //congestion_cost = k1*(0-(edge_c[xindex][yindex].C_H[presentlayer]-wirenumber))*(edge_c[xindex][yindex].historic[presentlayer]);
                    congestion_cost += k1*(0-(edge_c[xindex][yindex].C_H[presentlayer]-wirenumber))*(minwidth[presentlayer]+minspacing[presentlayer])*(edge_c[xindex][yindex].historic[presentlayer]);
            }
            icost_tmp = rt[treeindex].weight*Lambda*delay_tmp + Beta*vianumber + Miu*congestion_cost;

                if(congestion_cost < 0)
                   cout << "fuck"; 
            cost_tmp += icost_tmp;


            if(cost_tmp < min_cost){ 
                min_cost = cost_tmp;
                min_C = C_tmp;
                min_delay = delay_tmp;
                for(int j=0;j<n;j++)
                    min_combine[j] = combine[j];
            }
        }
        while( (++l<n) && (combine[l] == child_parallel[l]-1)  );
        if(l<n){
            ++combine[l];
            while(l)
                combine[--l]=0;
        }
    }
    record_dp[treeindex].done = true;
    record_dp[treeindex].C[presentlayerindex] = min_C;
    record_dp[treeindex].cost[presentlayerindex] = min_cost;
    for(int i=0;i<n;i++){
        record_dp[treeindex].child_combine[presentlayerindex][i] = min_combine[i];
    }
    delete [] combine;

}

void CIRCUIT::createsolution(vector<int>& dp_layer,int index,vector<RECORD_DP>& record_dp,vector<SEGMENT>& resultsegments, vector<TREE_NODE>& _rt, int netindex){
    int toplayer,lowlayer;
    _rt[index].child_C=-1;//for reset child_C which is used by compute elmore
    if( _rt[index].parent_index==-1){
        toplayer = dp_layer[index];
        lowlayer = dp_layer[index];
    }
    else if(_rt[index].x == _rt[_rt[index].parent_index].x){
        toplayer = (dp_layer[index]*2+1)%8;
        lowlayer = toplayer;
    }
    else{
        toplayer = ( ((dp_layer[index]+1)*2)%8==0 ? 8 :((dp_layer[index]+1)*2)%8);
        lowlayer = toplayer;
    }
    //create best child combine layer(wire segment) for node i and compute the top/low layer to create via segement
    for(int j=0;j< _rt[index].child_index.size();j++){
        int wirenumber;
        SEGMENT S_tmp;
        int child_layer = record_dp[index].child_combine[dp_layer[index]][j];
        dp_layer[ _rt[index].child_index[j] ] = child_layer;

        //the child_layer is layer "index"
        wirenumber = (child_layer/4)+1;
        if(_rt[index].x==_rt[_rt[index].child_index[j]].x){
            child_layer = (child_layer*2+1)%8;
        }
        else{
            if(((child_layer+1)*2)%8==0)
                child_layer=8;
            else
                child_layer = ((child_layer+1)*2)%8 ;
        }
        //after here, the child_layer is real layer 

        //if use NDR, wirenumber++
        if(wirenumber > 1 && child_layer >3)
            wirenumber++;

        //record the edge layer (u,v) u = index, v = child of index , it mean that RTree[i].layer is edge layer connect i and i's parent
        _rt[ _rt[index].child_index[j]].layer = child_layer;
        _rt[ _rt[index].child_index[j]].wirenumber = wirenumber;
        if( _rt[index].x == _rt[ _rt[index].child_index[j]].x  ){//vertical
            S_tmp.p1.x = _rt[index].x;
            S_tmp.p2.x = _rt[index].x;
            if(_rt[index].y < _rt[ _rt[index].child_index[j]].y){
                S_tmp.p1.y = _rt[index].y;
                S_tmp.p2.y = _rt[ _rt[index].child_index[j] ].y;
            }
            else{
                S_tmp.p2.y = _rt[index].y;
                S_tmp.p1.y = _rt[ _rt[index].child_index[j] ].y;
            }
            S_tmp.via = false;
            S_tmp.VORH = 1;
        }
        else{//horiziontal
            S_tmp.p1.y = _rt[index].y;
            S_tmp.p2.y = _rt[index].y;
            if(_rt[index].x < _rt[ _rt[index].child_index[j]].x){
                S_tmp.p1.x = _rt[index].x;
                S_tmp.p2.x = _rt[ _rt[index].child_index[j] ].x;
            }
            else{
                S_tmp.p2.x = _rt[index].x;
                S_tmp.p1.x = _rt[ _rt[index].child_index[j] ].x;
            }
            S_tmp.via = false;
            S_tmp.VORH = 2;
        }
        S_tmp.done = false;
        S_tmp.needripup=false;
        S_tmp.wirenumber = wirenumber;
        S_tmp.p1.layer = child_layer;
        S_tmp.p2.layer = child_layer;
        //congestion update 
        if(_rt[index].x == _rt[ _rt[index].child_index[j]].x){
            edge_c[S_tmp.p1.x/32][S_tmp.p1.y/40].C_V[S_tmp.p1.layer] -= wirenumber;
            if(wirenumber > 1&&globaliteration>fcfsthreshold)
                edge_c[S_tmp.p1.x/32][S_tmp.p1.y/40].total_c_v-=(wirenumber-1);
        }
        else{
            edge_c[S_tmp.p1.x/32][S_tmp.p1.y/40].C_H[S_tmp.p1.layer] -= wirenumber;
            if(wirenumber > 1&&globaliteration>fcfsthreshold)
                edge_c[S_tmp.p1.x/32][S_tmp.p1.y/40].total_c_h-=(wirenumber-1);
        }

        //updata grid edge's net 
        if(globaliteration>fcfsthreshold)
            S_tmp.needripup=true;
        edge_c[S_tmp.p1.x/32][S_tmp.p1.y/40].netid[S_tmp.p1.layer][netindex] = netindex;
        resultsegments.push_back(S_tmp);
        if(child_layer > toplayer)
            toplayer = child_layer;
        if(child_layer < lowlayer)
            lowlayer = child_layer;
    }
    if(_rt[index].pinlayer[0]!=-1){
        for(int j=0;j< _rt[index].pinlayer.size();j++){
            if(_rt[index].pinlayer[j] > toplayer)
                toplayer = _rt[index].pinlayer[j];
            if( _rt[index].pinlayer[j] < lowlayer)
                lowlayer = _rt[index].pinlayer[j];
        }
    }

    for(int j=lowlayer;j<toplayer;j++){
        SEGMENT S_tmp;
        S_tmp.p1.x = _rt[index].x;
        S_tmp.p2.x = S_tmp.p1.x;
        S_tmp.p1.y = _rt[index].y;
        S_tmp.p2.y = S_tmp.p1.y;
        S_tmp.p1.layer = j;
        S_tmp.p2.layer = j+1;
        S_tmp.via = true;
        S_tmp.VORH = 0;
        S_tmp.done = false;
        S_tmp.needripup=false;
        S_tmp.wirenumber = 1;
        resultsegments.push_back(S_tmp);
    }
    for(int i=0;i<_rt[index].child_index.size();i++)
        createsolution(dp_layer,_rt[index].child_index[i],record_dp,resultsegments, _rt,netindex);
}

void CIRCUIT::computeparallelnumber(){
    int number=0;
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        for(int j=0;j<nets[i].resultsegments.size();j++){
            if(nets[i].resultsegments[j].wirenumber==2)
                number++;
        }
    }

    cout << "parallel wire number: " << number << endl;
}

void CIRCUIT::emptyedge_c(){
    int t=0;
    for(int i=0;i<edge_c.size();i++){
        for(int j=0;j<edge_c[i].size();j++){
            for(int k=1;k<9;k++){
                if(k%2==1)
                    t += edge_c[i][j].C_V[k];
                else
                    t += edge_c[i][j].C_H[k];
            }
        }
    }

    cout << "empty : " << t << endl;
}

void CIRCUIT::computecritical(){
    vector<SORTEL> vtmp;
    int number=0;
    double total=0;
    double max=0;
    for(int i=0;i<nets.size();i++){
        if(!nets[i].havesegment)
            continue;
        number++;
        total += nets[i].eldelay;
        if(max < nets[i].eldelay)
            max = nets[i].eldelay;
        SORTEL stmp;
        stmp.index = i;
        stmp.eldelay = nets[i].eldelay;
        vtmp.push_back(stmp);
    }

    double aveg = total / number;

    //sort(vtmp.begin(),vtmp.end(),struct_sortel_cmp2);
    stable_sort(vtmp.begin(),vtmp.end(),struct_sortel_cmp2);
    
    number=0;
    double op5=0;
    double tt=0;
    for(int i= vtmp.size()-1;i >= vtmp.size()*0.995;i--){
        tt += vtmp[i].eldelay;
        number++;
    }
    double avegop5 = tt / number;

   
    double onetotal=0;
    number=0;
    for(int i= vtmp.size()-1;i >= vtmp.size()*0.99;i--){
        onetotal += vtmp[i].eldelay;
        number++;
    }
    double avegone = onetotal / number;


    number=0;
    double five=0;
    tt=0;
    for(int i= vtmp.size()-1;i >= vtmp.size()*0.95;i--){
        tt += vtmp[i].eldelay;
        number++;
    }
    double avegfive = tt / number;

    double tentotal=0;
    number=0;
    for(int i= vtmp.size()-1;i >= vtmp.size()*0.9;i--){
        tentotal += vtmp[i].eldelay;
        number++;
    }
    double avegten = tentotal / number;

    number=0;

    cout << "critical 0.5 %: " << avegop5 << endl;
    cout << "critical 1 %: " << avegone << "  num: " << vtmp.size()*0.01 << endl;
    cout << "critical 5 %: " << avegfive << endl;
    cout << "critical 10 %: " << avegten << " num: " << vtmp.size()*0.1 << endl;
    cout << "max :" << vtmp[vtmp.size()-1].eldelay << endl;
}


