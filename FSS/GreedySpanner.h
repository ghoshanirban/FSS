
#ifndef FASTSPARSESPANNER_GREEDYSPANNER_H
#define FASTSPARSESPANNER_GREEDYSPANNER_H

#include "../Utilities.h"
#include "boost/graph/dijkstra_shortest_paths_no_color_map.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <CGAL/Quadtree.h>

typedef std::vector<Point> Point_vector;
typedef CGAL::Quadtree<K, Point_vector> Quadtree;
typedef Quadtree::Node QuadTreeLeaf;

struct QuadTreeLeafInfo {
    unsigned leaderPId = UINT_MAX;
    std::vector<unsigned> internalPIDs;
    unsigned leafID = UINT_MAX;
    CGAL::Bbox_2 bboxOfLeaf;
    CGAL::Bbox_2 extendedBboxOfLeaf;
    Point centerOfBbox = Point(DBL_MAX, DBL_MAX);
    std::vector<Point> bboxOfVertices;
};

// The first point in the leaf is considered to be its key. If Leaf is empty then compute from bbox.
inline Point keyOf(const Quadtree &T, const QuadTreeLeaf &leaf) {
    if (leaf.empty()) {
        const CGAL::Bbox_2 &B = T.bbox(leaf);
        return Point{(B.xmax() - B.xmin()) / 2 + B.xmin(),
                     (B.ymax() - B.ymin()) / 2 + B.ymin()};
    }
    return *leaf.begin();
}

std::vector<Edge> constructFG_GreedySpanner(const std::vector<Point> &P,
                                            const QuadTreeLeaf &leaf,
                                            const Quadtree &T,
                                            const std::unordered_map<Point, QuadTreeLeafInfo> &THelper,
                                            const double t,
                                            unsigned &countE) {
    std::unordered_map<unsigned, unsigned> globalToLocalPIDs;
    std::vector<unsigned> localToGlobalPIDs;
    globalToLocalPIDs.reserve(leaf.size());
    localToGlobalPIDs.reserve(leaf.size());

    const QuadTreeLeafInfo &leafInfo = THelper.at(keyOf(T, leaf));
    for (unsigned i = 0; i < leafInfo.internalPIDs.size(); ++i) {
        globalToLocalPIDs.insert({leafInfo.internalPIDs[i], i});
        localToGlobalPIDs.emplace_back(leafInfo.internalPIDs[i]);
    }
    std::vector<Edge> localLeafEdges;
    std::vector<std::vector<double>> distances(leafInfo.internalPIDs.size(),
                                               std::vector<double>(leafInfo.internalPIDs.size(), 0)),
            shortestPath(leafInfo.internalPIDs.size(),
                         std::vector<double>(leafInfo.internalPIDs.size(), DBL_MAX));

    std::vector<Edge> completeGraphE;
    completeGraphE.reserve((leafInfo.internalPIDs.size() * (leafInfo.internalPIDs.size() - 1)) / 2);

    for (unsigned i = 0; i < globalToLocalPIDs.size() - 1; i++)
        for (unsigned j = i + 1; j < globalToLocalPIDs.size(); j++) {
            completeGraphE.emplace_back(i, j);
            distances[i][j] = distances[j][i] =
                    std::sqrt(squared_distance(P[localToGlobalPIDs[i]],
                                               P[localToGlobalPIDs[j]]));
        }
    std::sort(completeGraphE.begin(), completeGraphE.end(),
              [&distances](const Edge &e1, const Edge &e2) {
                  return distances[e1.first][e1.second] < distances[e2.first][e2.second];
              });

    /////////// boost's stuff
    using namespace boost;
    typedef
    adjacency_list<
            vecS, vecS, undirectedS, no_property,
            property<edge_weight_t, double, property<edge_weight2_t, double>>>
            Graph;
    typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    Graph g(&localLeafEdges[0], &localLeafEdges[localLeafEdges.size()], leafInfo.internalPIDs.size());
    std::vector<double> d(num_vertices(g));
    std::vector<vertex_descriptor> p(num_vertices(g));
    property_map<Graph, edge_weight_t>::type w = get(edge_weight, g);
    std::vector<double> weightsForDijkstraGraph(localLeafEdges.size());
    weightsForDijkstraGraph.reserve((size_t)4.5 * leafInfo.internalPIDs.size());
    graph_traits<Graph>::edge_iterator edge, e_end;
    boost::tie(edge, e_end) = edges(g);

    for (unsigned i = 0; i < localLeafEdges.size(); i++)
        weightsForDijkstraGraph[i] =
                distances[localLeafEdges[i].first][localLeafEdges[i].second];

    double *wp = &weightsForDijkstraGraph[0];
    //////////////////////

    std::vector<Edge> E;
    for (const Edge &e: completeGraphE) {
        double tTimesDistance = t * distances[e.first][e.second];
        if (shortestPath[e.first][e.second] <= tTimesDistance)
            continue;
        /////////////////// using boost's dijkstra
        for (; edge != e_end; ++edge)
            w[*edge] = *wp++;

        dijkstra_shortest_paths_no_color_map(
                g, vertex(e.first, g),
                boost::predecessor_map(&p[0]).distance_map(&d[0]));

        for (unsigned i = 0; i < d.size(); i++) {
            shortestPath[e.first][i] = shortestPath[i][e.first] = std::min(
                    shortestPath[e.first][i], d[i]);
        }
        ////////////////////
        if (shortestPath[e.first][e.second] > tTimesDistance) {
            countE++;
            E.emplace_back(localToGlobalPIDs[e.first],localToGlobalPIDs[e.second]);
            localLeafEdges.emplace_back(e.first, e.second);
            boost::add_edge(e.first, e.second, distances[e.first][e.second], g);
            weightsForDijkstraGraph.emplace_back(distances[e.first][e.second]);
            shortestPath[e.first][e.second] = shortestPath[e.second][e.first] = distances[e.first][e.second]; /// it was not present in main code
        }
    }
    return E;
}

#endif // FASTSPARSESPANNER_GREEDYSPANNER_H
