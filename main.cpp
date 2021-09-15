//
// Created by Xandru Mifsud on 09/09/2021.
// Please cite as: Mifsud X., "\[Lambda]--Core distance partitions", Linear Algebra Appl. (2021), https://doi.org/10.1016/j.laa.2020.12.012
//

#include <fstream>
#include <ctime>

#include "LaplacianLambdaCDP.h"

int main(int argc, char *argv[]){
    ull N; // dimension of the input adjacency matrix
    string adj_fpath; // file location of the CSV representation of the adjacency matrix

    // option checking...
    if(argc <= 2){
        throw std::runtime_error("At least one of the adjacency matrix CSV or the dimension was not specified...exiting...");
    }
    else if(argc > 3){
        throw std::runtime_error("Too many parameters specified...exiting...");
    }
    else{
        N = stoll(argv[1], nullptr);
        adj_fpath = argv[2];
    }

    std::ifstream adj_f(adj_fpath);

    MatrixXd laplacian(N, N);

    // Logging...
    auto timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tLoading adjacency matrix..." << endl;

    // From CSV file, load data as a MatrixXd instance
    if(adj_f.is_open()){
        std::string line;
        int v, u; v = u = 0;

        while(std::getline(adj_f, line)){
            char* ptr = (char*) line.c_str();
            ull len = line.length();
            u = 0;

            char* start = ptr;
            for(ull i = 0; i < len; i++){
                if(ptr[i] == ','){
                    laplacian(v, u++) = (-1)*strtod(start, nullptr);
                    start = ptr + i + 1;
                }
            }

            laplacian(v, u) = (-1)*strtod(start, nullptr);
            v++;
        }

        adj_f.close();
    }

    // Compute the degrees and set the diagonal entries...
    for(ull v = 0; v < N; v++){
        double d = 0;

        for(ull u = 0; u < N; u++){
            d += laplacian(v, u); // aggregate the degrees
        }

        laplacian(v, v) = (-1)*d; // set the diagonal entry to the sum of degrees across the v^th row
    }

    // Find the 'optimal' Lambda CDP for the Laplacian matrix...
    getLaplacianLambdaCDP(&laplacian, laplacian.rows());

    // Logging...
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tCompleted." << endl;

    return 0;
}