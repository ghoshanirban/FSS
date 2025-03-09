#include "Utilities.h"
#include "PointGenerators.h"
#include "FSS/FSS.h"
#include "FSS/StretchFactorCalculator.h"
#include "measurements/diameter.h"

int main() {

    std::vector<Point> P;

    // choose a distribution and generate a random pointset
    // unsigned n = 10E3; // set the number of points to be generated ( n > 1 )
    // std::cout << timeStamp() << "Generating points randomly ... " << std::endl;
    // RandomPointGenerator g;
    // g.uni_square(n, 1000, P);
    // std::cout << timeStamp() << "Point generation is complete." << std::endl;

    /*
     Other available generators:
        g.normal_clustered(10, P);
        g.grid_random(n, P);
        g.annulus(n, 500, 400, P);
        g.galaxy(n,5,P);
        g.convex(n, 100000, P);
        g.spokes(n,100000,P);
    */

    // Points can also be loaded from files
    readPointsFromFile("RealWorldPointsets/monalisa.txt",P);

    std::vector<Edge> E; // every edge is represented using a pair of unsigned integers
    std::cout << timeStamp() << "Generating spanner ... " << std::endl;
    FSSEuclideanGraph G(P);
    bool mayHaveDuplicates = true; // If P can have duplicates, then 'mayHaveDuplicate' must be set to true
    int threadCount = 4;
    G.fastSparseSpanner(E, 1.1, mayHaveDuplicates, threadCount, false);
    std::cout << timeStamp() << "Spanner generation is complete. " << std::endl;
    std::cout << timeStamp() << "Computing stretch-factor ... " << std::endl;
    std::cout << timeStamp() << "Stretch-factor: " <<  G.fastStretchFactor(threadCount) << std::endl;
    std::cout << timeStamp() << "|E|: " << E.size() << std::endl;
    std::cout << timeStamp() << "Average degree: " << (2.00 * (double)E.size()) / (double)P.size() << std::endl;
    std::cout << timeStamp() << "Diameter: " << diameter(P, E, threadCount) << std::endl;
    return 0;
}
