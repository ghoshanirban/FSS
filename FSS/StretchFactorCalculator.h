#ifndef FASTSPARSESPANNER_STRETCHFACTORCALCULATOR_H
#define FASTSPARSESPANNER_STRETCHFACTORCALCULATOR_H

#include "omp.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <CGAL/Memory_sizer.h>
#include "../Utilities.h"


struct priorityComparator {
    bool operator()(const std::pair<double, unsigned> &p1, const std::pair<double, unsigned> &p2) const {
        return p1.first > p2.first;
    }
};

typedef boost::heap::fibonacci_heap<std::pair<double, unsigned>, boost::heap::compare<priorityComparator>> FibonacciHeap;

double aStar(const std::vector<Point> &P, const std::vector<std::unordered_set<unsigned>> &G, const unsigned startVertex,
      const unsigned goalVertex) {

    FibonacciHeap openSet;
    std::vector<int> cameFrom(P.size(), -1);
    std::vector<bool> isInOpenSet(P.size(), false);
    std::vector<FibonacciHeap::handle_type> handleOfVertex(P.size());
    std::vector<double> g(P.size(), std::numeric_limits<double>::infinity());

    handleOfVertex[startVertex] = openSet.push(std::make_pair(L2Distance(P[startVertex], P[goalVertex]), startVertex));
    isInOpenSet[startVertex] = true;
    g[startVertex] = 0;

    while (!openSet.empty()) {
        std::pair<double, unsigned> current = openSet.top();
        openSet.pop();
        isInOpenSet[current.second] = false;
        if (current.second == goalVertex) {
            double pathLength = 0;
            unsigned currentVertex = current.second;
            while (currentVertex != startVertex) {
                pathLength += L2Distance(P[currentVertex], P[cameFrom[currentVertex]]);
                currentVertex = cameFrom[currentVertex];
            }
            return pathLength;
        }

        for (unsigned neighbor: G[current.second]) {
            double tentativeGscore = g[current.second] + L2Distance(P[current.second], P[neighbor]);
            if (tentativeGscore < g[neighbor]) {
                cameFrom[neighbor] = current.second;
                g[neighbor] = tentativeGscore;
                double fOfNeighbor = g[neighbor] + L2Distance(P[neighbor], P[goalVertex]);

                if (!isInOpenSet[neighbor]) {
                    handleOfVertex[neighbor] = openSet.push(std::make_pair(fOfNeighbor, neighbor));
                    isInOpenSet[neighbor] = true;
                } else
                    openSet.decrease(handleOfVertex[neighbor], std::make_pair(fOfNeighbor, neighbor));
            }
        }
    }
    return std::numeric_limits<double>::infinity();
}

auto exactStretchFactorByDijkstraParallel(const std::vector<Point> &P, const std::vector<Edge> &E, const double &t,
                                          std::vector<std::pair<Point, Point>> &pointPairsWithoutTPaths, unsigned threadCount = 1) {
    try {
        if (P.empty())
            throw std::runtime_error("Point set is empty! Aborting.");

        if(t < 1)
            throw std::runtime_error("t must be > 1. Aborting!");
    }
    catch(std::runtime_error &e) {
        std::cout << e.what() << std::endl;
        return ;
    }
#if LOG == ON
    std::cout << timeStamp() << "--------------------------------------Parallel Dijkstra has started--------------------------------"<<std::endl;
#endif
    CGAL::Real_timer parallelDijkstraTime;
    parallelDijkstraTime.start();

    using namespace boost;
    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double, property<edge_weight2_t, double>>> Graph;
    typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    Graph g(&E[0], &E[E.size()], P.size());
    graph_traits<Graph>::edge_iterator e, e_end;
    property_map<Graph, edge_weight_t>::type w = get(edge_weight, g);

    std::vector<double> weights(E.size());
    for (unsigned i = 0; i < E.size(); i++)
        weights[i] = L2Distance(P[E[i].first], P[E[i].second]);

    double *wp = &weights[0];

    for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
        w[*e] = *wp++;

    std::vector<double> stretchFactorPerThread(threadCount, 0.0);

    std::vector<std::vector<std::pair<Point, Point>>> threadSpecificBadPairs(threadCount);

#pragma omp parallel for default(none) shared(P, threadSpecificBadPairs, stretchFactorPerThread, t, g) schedule(dynamic) num_threads(threadCount)
    for (unsigned s = 0; s < P.size() - 1; s++) {
        double stretchFactor = 0.0;
        std::vector<double> d(num_vertices(g));
        std::vector<vertex_descriptor> p(num_vertices(g));

        dijkstra_shortest_paths_no_color_map(g, vertex(s, g), boost::predecessor_map(&p[0]).distance_map(&d[0]));

        for (unsigned i = s + 1; i < P.size(); i++) {
            double currentStretchFactor = d[i] / L2Distance(P[s], P[i]);
            if (currentStretchFactor > stretchFactor) {
                stretchFactor = currentStretchFactor;
                if (stretchFactor > t)
                    threadSpecificBadPairs[omp_get_thread_num()].emplace_back(P.at(s), P.at(i));
            }
        }
        stretchFactorPerThread[omp_get_thread_num()] = std::max(stretchFactor,
                                                                stretchFactorPerThread[omp_get_thread_num()]);
    }
    for (auto &v: threadSpecificBadPairs)
        for (auto pair: v)
            pointPairsWithoutTPaths.emplace_back(pair);

    double stretchFactorofG = 0.0;
    for (auto sf: stretchFactorPerThread)
        stretchFactorofG = std::max(stretchFactorofG, sf);
    parallelDijkstraTime.stop();

#if LOG == ON
    std::cout << timeStamp() << "Done. Took " << parallelDijkstraTime.time() << "s = "
              << parallelDijkstraTime.time() / 60.0 << "m = "
              << parallelDijkstraTime.time() / 60.0 / 60.0 << "h " << std::endl;
    std::cout << timeStamp() << "Stretch factor: " << stretchFactorofG << std::endl;
    std::cout << timeStamp() << "--------------------------------------Parallel Dijkstra ended-----------------------------------"<<std::endl<<std::endl;
#endif
}


double greedyHeuristic(const std::vector<Point> &P, const std::vector<std::unordered_set<unsigned>> &G,
                       const unsigned startVertex, const unsigned goalVertex){

    std::vector<bool> isInOpenSet(P.size(), false);

    unsigned currentNode = startVertex;
    isInOpenSet[currentNode] = true;
    double pathLength = 0.0;

    while( currentNode != goalVertex ){
        auto localMinHeuristic = DBL_MAX;
        unsigned bestNeighbor = 0;
        double bestNeighborDist = 0.0;

        for( auto neighbor: G[currentNode] ){
            if( !isInOpenSet[neighbor] ){
                isInOpenSet[neighbor] = true;

                double currentNodeToNeighborDist = L2Distance(P[currentNode],P[neighbor]);
                double neighborHeuristic = pathLength + currentNodeToNeighborDist + L2Distance(P[neighbor],P[goalVertex]);

                if( localMinHeuristic >= neighborHeuristic ){
                    localMinHeuristic = neighborHeuristic;
                    bestNeighbor = neighbor;
                    bestNeighborDist = currentNodeToNeighborDist;
                }
            }
            if( neighbor == goalVertex ) {
                pathLength += bestNeighborDist;
                return pathLength;
            }
        }
        if( localMinHeuristic == DBL_MAX ) break;
        pathLength += bestNeighborDist;
        currentNode = bestNeighbor;
    }
    return std::numeric_limits<double>::infinity();
}

#endif //FASTSPARSESPANNER_STRETCHFACTORCALCULATOR_H
