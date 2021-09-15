//
// Created by Xandru Mifsud on 09/09/2021.
// Please cite as: Mifsud X., "\[Lambda]--Core distance partitions", Linear Algebra Appl. (2021), https://doi.org/10.1016/j.laa.2020.12.012
//

#ifndef LAPLACIANLAMBDACDP_H
#define LAPLACIANLAMBDACDP_H

// Macros for determining if two values are approximately equal to each other for a given epsilon
#define EPSILON 1.0e-15
#define approx_eig(l1, l2, N) ((l1 < (l2 + (EPSILON*N))) && (l1 > (l2 - (EPSILON*N))))
#define approx_zero(l, N) ((l < (EPSILON*N)) && (l > (-1*EPSILON*N)))

#include <iostream>
#include <vector>
#include <ctime>
#include <set>

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

typedef unsigned long long ull;
typedef vector<set<ull>*>* cdp; // CDPs are given as a pointer to a vectors of pointers to sets, each set being a block

cdp getLaplacianLambdaCDP(MatrixXd* laplacian, ull N);
cdp nestedNeighbourhoods(MatrixXd* laplacian, set<ull>* cores);
double LambdaCDPIndex(cdp lambdaCDP, ull N);

/* Computes the 'optimal' Lambda-CDP (LCDP) for the Laplacian of a given graph. The optimal LCDP is defined as the LCDP
 * with the largest LCDP Index; this is generally the LCDP with the largest number of blocks, of approximately equal size
 * (or the closest possible from all generated LCDPs).
 *
 * Parameters:
 *  MatrixXd* laplacian: A pointer to a MatrixXd instance representing an Laplacian matrix of a graph. It is assumed that
 *                       laplacian is a valid N x N symmetric matrix satisfying Laplace's equation L = D - A.
 *  ull N              : The number of vertices in the graph.
 *
 * Outline:
 *  (1) Compute the eigen-decomposition of L, given by L = USU^T, there the columns of U are the eigenvectors of L and
 *      the diagonal entries of S are the corresponding eigenvalues.
 *  (2) Compute the set of core vertices for each eigenvalue of L, via the eigen-decomposition.
 *  (3) Compute the LCDP associated with each eigenvalue, via the corresponding set of core vertices.
 *  (4) Return the LCDP with the largest LCDP Index, while freeing any memory associated with the un-returned LCDPs.
 */
cdp getLaplacianLambdaCDP(MatrixXd* laplacian, ull N){
    // Logging...
    auto timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tComputing optimal CDP..." << endl;

    // Logging...
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tFinding the eigensystem of the Laplacian matrix..." << endl;

    /* Theorem: For a real-symmetric positive semi-definite matrix, its singular value decomposition (SVD) and eigen-
     *          decomposition correspond to each other.
     *
     * Proposition: The Laplacian of a graph is a real-symmetric positive semi-definite matrix.
     *
     * Consequently, we can use the divide-and-conquer SVD implementation of the Eigen package, which albeit slower then
     * the eigen-decomposition implementation for small dimensionality, the SVD implementation scales better for larger
     * dimensionality while retaining superior accuracy.
     */
    auto svd = laplacian->bdcSvd(ComputeThinU); // compute the SVD for the Laplacian L
    const auto& eigs = svd.singularValues(); // get the eigenvalues of L; note that these are sorted in decreasing order

    // Logging...
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tComputing the core vertex sets..." << endl;

    // Compute the set of core vertices for each eigenvalues of L
    auto core_sets = new vector<set<ull>*>;
    for(ull eig_min = 0; eig_min < N; eig_min++){
        auto cores = new set<ull>;

        // Find the maximum range of indices (eig_min to eig_max) of the eigs vector such that all entries are (approx)
        // the same eigenvalue. We use an approximate, up to some epsilon, to account for numerical errors.
        ull eig_max = eig_min + 1;
        for(; eig_max < N; eig_max++){
            if(!approx_eig(eigs(eig_max), eigs(eig_min), N)){
                break;
            }
        }
        eig_max--;

        /* Hence the columns of U in the SVD of L with indices between eig_min and eig_max for an eigen-basis for the
         * same eigenvalue. From these columns we can then construct the core vertex set of this eigenvalue; a vertex v
         * is in the core vertex set if there is one such column vector with the v^th entry (approx.) not zero. We use
         * an approximate, up to some epsilon, to account for numerical errors.
         */
        for(ull v = 0; v < N; v++){ // for each vertex v
            for(ull i = eig_min; i <= eig_max; i++){ // for each column vector of U for this eigenvalue
                if(!approx_zero((svd.matrixU())(v, i), N)){ // if the v^th entry is non-zero, v is a core vertex
                    cores->insert(v); // hence we insert it into the core vertex set
                    break; // and there is no need to check the v^th entry of the remaining columns
                }
            }
        }

        core_sets->push_back(cores);
        eig_min = eig_max;
    }

    // Logging...
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\t# Core Vertex Sets = " << core_sets->size() << endl;
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Log " << ctime(&timenow) <<"\tComputing and finding the optimal LCDP..." << endl;

    /* For each set of core vertices associated with a distinct eigenvalue, we find the full LCDP and its LCDP Index.
     * Note that we maintain the largest (optimal) LCDP Index and the associated LCDP; if the LCDP associated with a
     * core vertex set has and LCDP Index less than the current maximum, then it is not 'optimal' and hence we discard
     * it (in particular we carry out memory management on the go, especially useful for large numbers of vertices).
     */
    cdp ret_cdp, curr_cdp;
    double max_index = 0;
    for(auto cv: *core_sets){
        curr_cdp = nestedNeighbourhoods(laplacian, cv); // calculate the full LCDP
        double index = LambdaCDPIndex(curr_cdp, N); // and find the corresponding index
        if(index >= max_index){ // if index exceed the current maximum, this LCDP is more optimal
            max_index = index;
            ret_cdp = curr_cdp;
        }
        else{ // otherwise we discard it by freeing any associated memory
            for(auto p: *curr_cdp){ delete p;}
            delete curr_cdp;
        }
    }

    delete core_sets;

    return ret_cdp; // lastly, we return the optimally computed LCDP
}

/* Given a Laplacian matrix of a graph and a subset of the vertex set, this function computes the nested neighbourhoods
 * of this subset, which is a partition of the vertex set when the graph is connected. For a subset S of the vertex set,
 * we define its neighbourhood N(S) as the set of vertices not in S that are adjacent to some vertex in S.
 */
cdp nestedNeighbourhoods(MatrixXd* laplacian, set<ull>* cores){
    // Vector structure to to hold the nested neighbourhood set instances, in order of distance.
    auto nestedN = new vector<set<ull>*>;
    nestedN->push_back(cores);

    // The set 'remaining_vertices' initially has all the vertices except those in the initial subset passed. Termination
    // follows when this set is empty. Termination is guaranteed if the graph is connected.
    set<ull> remaining_vertices;
    for(ull i = 0; i < laplacian->rows(); i++){ // for each vertex in the graph
        if(cores->count(i) == 0){ // insert it into the set 'remaining_vertices' if it is not in the set 'cores'
            remaining_vertices.insert(i);
        }
    }

    // Construct nested neighbourhoods until we have visited all the vertices in the graph...
    while(!remaining_vertices.empty()){
        auto curr_neighbourhood = new set<ull>; // new set instance for the next nested neighbourhood

        set<ull>::iterator v;
        for(v = remaining_vertices.begin(); v != remaining_vertices.end();){ // Iterate through the unvisited vertices
            bool v_deleted = false;
            for(auto& u: *(nestedN->back())){ // for every vertex in the last neighbourhood...
                if((*laplacian)(*v, u) == -1){ // if adjacent to some vertex in the last neighbourhood
                    curr_neighbourhood->insert(*v); // then insert into new neighbourhood
                    v = remaining_vertices.erase(v); // and remove from the list of remaining vertices
                    v_deleted = true;

                    break;
                }
            }

            // If not deleted, i.e. not a neighbour to any vertex in the last neighbourhood, increment iterator.
            if(!v_deleted){ ++v;}
        }

        nestedN->push_back(curr_neighbourhood);
    }

    return nestedN;
}

/* Computes the LCDP Index for a given LCDP for a graph on N vertices. The formula is given by:
 *                                (-1/N)*Sum(|p| * Log10(|p| / N), p in LCDP)
 */
double LambdaCDPIndex(cdp lambdaCDP, ull N){
    double sum = 0;

    for(auto p: *lambdaCDP){ // Sum over p in LCDP
        sum += ((p->size()) * log10((double) (p->size())/N)); // |p| * Log10(|p| / N)
    }

    return (-1.0) * (sum / N);
}

#endif //LAPLACIANLAMBDACDP_H
