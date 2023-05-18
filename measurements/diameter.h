#ifndef FASTSPARSESPANNER_DIAMETER_H
#define FASTSPARSESPANNER_DIAMETER_H

#include "../Utilities.h"


unsigned bfs(unsigned src,  std::vector<std::unordered_set<unsigned >> &adjGraph){
    std::queue<unsigned > Q;
    std::vector<bool> visited(adjGraph.size(),false);
    std::vector<unsigned > dist(adjGraph.size(),0);
    unsigned maxDistance = 0;
    Q.push(src);
    visited[src] = true;
    while(!Q.empty()){
        unsigned node = Q.front();
        Q.pop();
        for(auto neighbor: adjGraph[node]){
            if(!visited[neighbor]){
                visited[neighbor] = true;
                dist[neighbor] = 1 + dist[node];
                maxDistance = std::max(maxDistance,dist[neighbor]);
                Q.push(neighbor);
            }
        }
    }
    return maxDistance;
}

unsigned diameter(const std::vector<Point> &P, const std::vector<Edge> &E, unsigned threadCount = 1) {
    std::vector<std::unordered_set<unsigned >> adjGraph(P.size(),std::unordered_set<unsigned >());

    for( auto &edge : E ) {
        adjGraph[edge.first].insert(edge.second);
        adjGraph[edge.second].insert(edge.first);
    }

    unsigned diameter = 0;
    std::vector<unsigned > diameters(threadCount,0);

    #pragma omp parallel for default(none) shared(P,adjGraph,diameters) schedule(dynamic) num_threads(threadCount)
    for( unsigned i = 0; i<P.size(); i++ )
        diameters[omp_get_thread_num()] = std::max(diameters[omp_get_thread_num()],bfs(i,adjGraph));

    for( auto nextDiameter: diameters )
        diameter = std::max(diameter,nextDiameter);

    return diameter;
}

#endif //FASTSPARSESPANNER_DIAMETER_H

