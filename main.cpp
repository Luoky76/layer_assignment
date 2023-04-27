/*
 * usage : ./main benchmark(superblueXX) globla_routing_result marge_file outputfile mode
 *
*/

//手动开栈
// 16777216 = 1024*1024*16 -> 16m // 33554432 = 1024*1024*32 -> 32m
// 1024*1024*128 = 134217728 -> 128m // 268435456 , 1073741824 = 1Gb
#pragma comment(linker,"/STACK:16777216")

#include "parameter.h"
#include "database.h"
//#include <sys/time.h>
#include <ctime>
#include <ratio>
#include <chrono>
//void readparameter(string );
int main(int argc, char* argv[]) {
    vector <string> filenames;
    filenames.resize(6);
    filenames[0] = "-----error-----";
    string cir_number; //数据编号
    cout << "input the superblue:";
    cin >> cir_number;
    filenames[1] = "data\\superblue" + cir_number;
    filenames[2] = "data\\HQ_GRsb" + cir_number;
    filenames[3] = "data\\sb" + cir_number + ".marge";
    filenames[4] = "data\\output"; // FIXME can add cir_number.
    filenames[5] = "data\\my_wp_parameter_v2.txt"; // alined all benchmark

    //struct timeval s_t;
    //gettimeofday(&s_t, NULL);
    //开始时间
    chrono::steady_clock::time_point STc = chrono::steady_clock::now();

    readparameter(filenames[5]);
    CIRCUIT cir(filenames);

    //struct timeval e_t;
    //gettimeofday(&e_t, NULL);
    //结束时间
    chrono::steady_clock::time_point ETc = chrono::steady_clock::now();

    //double time;
    //time = ((e_t.tv_sec*1000000)+e_t.tv_usec)-((s_t.tv_sec*1000000)+s_t.tv_usec);
    //间隔时间
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double> >(ETc - STc);
    
    //printf("time: %lf \n" ,time/1000000);
    //打印间隔时间
    cout << "Ctime: " << time_span.count() << " (s)" << endl;

    return 0;
}
