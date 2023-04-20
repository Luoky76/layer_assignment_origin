/*
 * usage : ./main benchmark(superblueXX) globla_routing_result marge_file outputfile mode
 *
*/

#pragma comment(linker,"/STACK:16777216")   //手动开栈

#include "parameter.h"
#include "database.h"
//#include <sys/time.h>
#include <ctime>
#include <ratio>
#include <chrono>
//void readparameter(string );
int main(int argc, char* argv[]) {
    std::vector <string> filenames;
    filenames.resize(6);
    filenames[0] = "-----error-----";
    std::string cir_number; //数据编号
    std::cout << "input the superblue:";
    std::cin >> cir_number;
    filenames[1] = "superblue" + cir_number;
    filenames[2] = "HQ_GRsb" + cir_number;
    filenames[3] = "sb" + cir_number + ".marge";
    filenames[4] = "output"; // FIXME can add cir_number.
    filenames[5] = "my_wp_parameter_v2.txt"; // alined all benchmark

    //struct timeval s_t;
    //gettimeofday(&s_t, NULL);
    //开始时间
    std::chrono::steady_clock::time_point STc = std::chrono::steady_clock::now();

    readparameter(filenames[5]);
    CIRCUIT cir(filenames);

    //struct timeval e_t;
    //gettimeofday(&e_t, NULL);
    //结束时间
    std::chrono::steady_clock::time_point ETc = std::chrono::steady_clock::now();

    //double time;
    //time = ((e_t.tv_sec*1000000)+e_t.tv_usec)-((s_t.tv_sec*1000000)+s_t.tv_usec);
    //间隔时间
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(ETc - STc);
    
    //printf("time: %lf \n" ,time/1000000);
    //打印间隔时间
    std::cout << "Ctime: " << time_span.count() << " (s)" << std::endl;

    return 0;
}
