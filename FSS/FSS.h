#ifndef FASTSPARSESPANNER_FSS_H
#define FASTSPARSESPANNER_FSS_H

#include <boost/dynamic_bitset.hpp>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Quadtree.h>
#include <CGAL/Real_timer.h>

#include "GreedySpanner.h"
#include "WSPDSpanner.h"
#include "StretchFactorCalculator.h"
#include "../Utilities.h"


/**
* @param points it is assumed that there are no duplicates in 'points' otherwise the constructed spanner may have INF stretch factor
* @param t needs to be greater than 1.0
* @param E should be empty
*/

class FSSEuclideanGraph {
private:
    typedef std::tuple<unsigned, unsigned, double> IDsDistanceTuple;
    Quadtree *T = nullptr;
    std::unordered_map<Point, unsigned> PointToIndexMap;
    std::vector<std::unordered_set<unsigned>> H;
    unsigned totalNumOfLeaf = UINT_MAX;
    std::unordered_map<Point, QuadTreeLeafInfo> PointToLeafMap;
    std::vector<QuadTreeLeaf> leafIdToLeaf;
    std::unordered_set<index_tPair, index_tPairHash> FSSMergedLeafPairs;
    std::vector<Point> P;
    std::vector<Point> copyOfP;
    double expectedStretchFactor = DBL_MIN;
    bool isFSSCompleted = false;
    bool isFSFCompleted = false;
    double fastStretchFactor_t = DBL_MIN;
    double t = DBL_MAX;
    double sizeOfSmallestLeaf = DBL_MAX;
public:
    explicit FSSEuclideanGraph(const std::vector<Point> &P) {
        this->P = P;
        this->copyOfP = P;
        if (!P.empty())
            this->T = new Quadtree(copyOfP);
        copyOfP.clear();
    }

    inline Point keyOf(const QuadTreeLeaf &Leaf) {
        if (Leaf.empty()) {
            const CGAL::Bbox_2 &B = T->bbox(Leaf);
            return Point{(B.xmax() - B.xmin()) / 2 + B.xmin(),
                         (B.ymax() - B.ymin()) / 2 + B.ymin()};
        }
        return *Leaf.begin();
    }

    static inline CGAL::Bbox_2 enlargeBbox(CGAL::Bbox_2 &bboxOfLeaf, double enlargeFactor) {
        return {bboxOfLeaf.xmin() - enlargeFactor, bboxOfLeaf.ymin() - enlargeFactor,
                bboxOfLeaf.xmax() + enlargeFactor, bboxOfLeaf.ymax() + enlargeFactor};
    }

    bool isDiagonal(const QuadTreeLeaf &firstLeaf, const QuadTreeLeaf &secondLeaf) {
        QuadTreeLeafInfo &firstLeafInfo = PointToLeafMap[keyOf(firstLeaf)];
        QuadTreeLeafInfo &secondLeafInfo = PointToLeafMap[keyOf(secondLeaf)];
        for (Point &firstVertex: firstLeafInfo.bboxOfVertices) {
            if (CGAL::collinear(firstLeafInfo.centerOfBbox, secondLeafInfo.centerOfBbox, firstVertex)) {
                for (Point &secondVertex: secondLeafInfo.bboxOfVertices) {
                    if (CGAL::collinear(firstLeafInfo.centerOfBbox, secondLeafInfo.centerOfBbox, secondVertex)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    void getNextLayerNeighbors(std::vector<std::unordered_set<unsigned >> &neighborsGraph,
                               std::vector<QuadTreeLeaf> &LeafIdToLeaf,
                               boost::dynamic_bitset<> &processedLeafs,
                               std::vector<QuadTreeLeaf> &curLeaves,
                               std::vector<QuadTreeLeaf> &expectedLayerNeighbors,
                               unsigned curLayerNo, unsigned expectedLayerNo) {
        std::vector<QuadTreeLeaf> nextLayerNeighbors;
        for (const auto &curLeaf: curLeaves) {
            if (!curLeaf.is_null()) {
                QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(curLeaf));
                for (auto &nextNeighborId: neighborsGraph[curLeafInfo.leafID]) {
                    const QuadTreeLeaf &nextNeighborLeaf = LeafIdToLeaf[nextNeighborId];
                    QuadTreeLeafInfo &nextNeighborInfo = PointToLeafMap.at(keyOf(nextNeighborLeaf));
                    if (!nextNeighborLeaf.is_null() && !processedLeafs.test(nextNeighborId)) {
                        processedLeafs[nextNeighborId] = true;
                        nextLayerNeighbors.emplace_back(nextNeighborLeaf);
                        if (curLayerNo == expectedLayerNo) {
                            expectedLayerNeighbors.emplace_back(nextNeighborLeaf);
                        }
                    }
                }
            }
        }
        curLeaves.clear();
        curLeaves = nextLayerNeighbors;
    }


    void intersectedLeafs(const CGAL::Bbox_2 &B, const QuadTreeLeaf &Leaf, std::vector<QuadTreeLeaf> &output) {
        if (B != T->bbox(Leaf) && CGAL::do_overlap(B, T->bbox(Leaf))) {
            if (Leaf.is_leaf()) {
                output.emplace_back(Leaf);
                return;
            }
            for (int i = 0; i < QuadTreeLeaf::Degree::value; ++i) {
                intersectedLeafs(B, Leaf[i], output);
            }
        }
    }


    Point computeLeader(const QuadTreeLeaf &leaf) {
        QuadTreeLeafInfo &leafOfQuadTree = PointToLeafMap[keyOf(leaf)];
        /// leader point computation of individual quad tree leaf Leaf
        CGAL::Bbox_2 bboxOfLeaf = CGAL::bbox_2(leaf.begin(), leaf.end());
        Point centerOfBoundingBox((bboxOfLeaf.xmax() - bboxOfLeaf.xmin()) / 2 + bboxOfLeaf.xmin(),
                                  (bboxOfLeaf.ymax() - bboxOfLeaf.ymin()) / 2 + bboxOfLeaf.ymin());
        double minDistance = std::numeric_limits<double>::infinity();
        Point leafLeaderP;
        for (Point const &point: leaf) {
            double distance = (CGAL::squared_distance(point, centerOfBoundingBox));
            if (distance < minDistance) {
                minDistance = distance;
                leafLeaderP = point;
            }
        }
        leafOfQuadTree.leaderPId = PointToIndexMap[leafLeaderP];
        return leafLeaderP;
    }

    void greedyMerge(
            const QuadTreeLeaf &firstLeaf,
            const QuadTreeLeaf &secondLeaf,
            unsigned &countE,
            std::vector<std::unordered_set<unsigned>> &G,
            std::vector<Edge> &edgesByThread) {

        QuadTreeLeafInfo &firstLeafInfo = PointToLeafMap[keyOf(firstLeaf)], &secondLeafInfo = PointToLeafMap[keyOf(
                secondLeaf)];
        std::vector<unsigned> &firstLeafInfoPIds = firstLeafInfo.internalPIDs, &secondLeafInfoPIds = secondLeafInfo.internalPIDs;
        unsigned totalLeavesPsSize = firstLeafInfoPIds.size() + secondLeafInfoPIds.size();

        std::unordered_map<unsigned, unsigned> firstLeafGlobalToLocalPIDs(totalLeavesPsSize);
        std::vector<unsigned> firstLeafLocalToGlobalPIDs(totalLeavesPsSize);

        // First Leaf points mapping
        for (unsigned i = 0; i < firstLeafInfoPIds.size(); ++i) {
            firstLeafGlobalToLocalPIDs.insert({firstLeafInfoPIds[i], i});
            firstLeafLocalToGlobalPIDs[i] = firstLeafInfoPIds[i];
        }
        // Second Leaf points mapping
        std::unordered_map<unsigned, unsigned> secondLeafGlobalToLocalPIDs(totalLeavesPsSize);
        std::vector<unsigned> secondLeafLocalToGlobalPIDs(totalLeavesPsSize);
        for (unsigned i = 0; i < secondLeafInfoPIds.size(); ++i) {
            secondLeafGlobalToLocalPIDs.insert({secondLeafInfoPIds[i], i});
            secondLeafLocalToGlobalPIDs[i] = secondLeafInfoPIds[i];
        }
        // First and second Leaf leader points info
        unsigned &firstLeafLeaderP = firstLeafInfo.leaderPId, &secondLeafLeaderP = secondLeafInfo.leaderPId,
                firstLeafLeaderLocalId = firstLeafGlobalToLocalPIDs[firstLeafLeaderP],
                secondLeafLeaderLocalId = secondLeafGlobalToLocalPIDs[secondLeafLeaderP];
        double euclideanDistanceOfLeaderPs = std::sqrt(
                CGAL::squared_distance(P[firstLeafLeaderP], P[secondLeafLeaderP]));

        std::vector<std::vector<double>> firstLeafPointsIntraDist(firstLeaf.size(),
                                                                  std::vector<double>(firstLeaf.size(), 0.0));
        std::vector<std::vector<double>> secondLeafDistance(secondLeaf.size(),
                                                            std::vector<double>(secondLeaf.size(), 0.0));

        // If direct bride exist then compute all distances based on leader point
        bool isDirectBridgePresent = false;
        if (H[firstLeafLeaderP].find(secondLeafLeaderP) != H[firstLeafLeaderP].end()) {
            isDirectBridgePresent = true;
            for (unsigned i = 0; i < firstLeafInfoPIds.size(); i++) {
                firstLeafPointsIntraDist[i][firstLeafLeaderLocalId] = firstLeafPointsIntraDist[firstLeafLeaderLocalId][i] =
                        (std::sqrt(CGAL::squared_distance(P[firstLeafInfoPIds[i]], P[firstLeafLeaderP])));
            }
            for (unsigned i = 0; i < secondLeafInfoPIds.size(); i++) {
                secondLeafDistance[i][secondLeafLeaderLocalId] =
                secondLeafDistance[secondLeafLeaderLocalId][i] =
                        (std::sqrt(CGAL::squared_distance(P[secondLeafInfoPIds[i]], P[secondLeafLeaderP])));
            }
        }

        std::vector<IDsDistanceTuple> interLeafPairs;
        interLeafPairs.reserve((totalLeavesPsSize * (totalLeavesPsSize - 1)) / 2);

        // If direct bridge exists then prune otherwise compute all pairs for further merging.
        for (unsigned firstLeafInfoPId: firstLeafInfoPIds) {
            unsigned uGlobalId = firstLeafInfoPId, u = firstLeafGlobalToLocalPIDs[uGlobalId];
            Point uPoint = P[uGlobalId];
            for (unsigned secondLeafInfoPId: secondLeafInfoPIds) {
                unsigned vGlobalId = secondLeafInfoPId, v = secondLeafGlobalToLocalPIDs[vGlobalId];
                Point vPoint = P[vGlobalId];
                double euclideanDistanceOfUV = std::sqrt(CGAL::squared_distance(uPoint, vPoint));
                if (isDirectBridgePresent) {
                    double leftSide = t * (firstLeafPointsIntraDist[u][firstLeafLeaderLocalId] +
                                           secondLeafDistance[secondLeafLeaderLocalId][v]
                                           - euclideanDistanceOfUV);
                    if (leftSide > -euclideanDistanceOfLeaderPs)
                        interLeafPairs.emplace_back(u, v, euclideanDistanceOfUV);
                } else
                    interLeafPairs.emplace_back(u, v, euclideanDistanceOfUV);
            }
        }
        std::sort(interLeafPairs.begin(), interLeafPairs.end(),
                  [](const IDsDistanceTuple &e1, const IDsDistanceTuple &e2) {
                      return get<2>(e1) < get<2>(e2);
                  });

        // Main processing inter Leaf pairs for merging of two Leafs
        std::vector<Edge> leavesBridges;
        std::vector<double> bridgesEuclideanDistance;

        for (const IDsDistanceTuple &interLeafPair: interLeafPairs) {
            unsigned u = get<0>(interLeafPair), v = get<1>(
                    interLeafPair), uGlobalPId = firstLeafLocalToGlobalPIDs[u], vGlobalPId = secondLeafLocalToGlobalPIDs[v];
            double euclideanDistanceUV = get<2>(interLeafPair);
            bool isNewBridgeNeeded = true;
            for (int b = (int) leavesBridges.size() - 1; b >= 0; b--) {
                unsigned leavesBridgesFirstPointLocalId = leavesBridges[b].first,
                        leavesBridgesSecondPointLocalId = leavesBridges[b].second;
                if (CGAL::is_zero(firstLeafPointsIntraDist[u][leavesBridgesFirstPointLocalId])) {
                    firstLeafPointsIntraDist[u][leavesBridgesFirstPointLocalId] = firstLeafPointsIntraDist[leavesBridgesFirstPointLocalId][u] =
                            std::sqrt(CGAL::squared_distance(P[uGlobalPId],
                                                             P[firstLeafLocalToGlobalPIDs[leavesBridgesFirstPointLocalId]]));
                }
                if (CGAL::is_zero(secondLeafDistance[leavesBridgesSecondPointLocalId][v])) {
                    secondLeafDistance[leavesBridgesSecondPointLocalId][v] = secondLeafDistance[v][leavesBridgesSecondPointLocalId] =
                            std::sqrt(CGAL::squared_distance(P[vGlobalPId],
                                                             P[secondLeafLocalToGlobalPIDs[leavesBridgesSecondPointLocalId]]));
                }
                double leftSide = t * (firstLeafPointsIntraDist[u][leavesBridgesFirstPointLocalId] +
                                       secondLeafDistance[leavesBridgesSecondPointLocalId][v] -
                                       euclideanDistanceUV);
                if (leftSide <= -bridgesEuclideanDistance[b]) {
                    isNewBridgeNeeded = false;
                    break;
                }
            }
            if (isNewBridgeNeeded) {
                leavesBridges.emplace_back(u, v);
                auto greedyHeuristicDist = DBL_MAX;
                greedyHeuristicDist = greedyHeuristic(P, G, uGlobalPId, vGlobalPId);
                if (greedyHeuristicDist <= t * euclideanDistanceUV) {
                    bridgesEuclideanDistance.emplace_back(greedyHeuristicDist);
                } else {
                    double aStarShortestPathDist = aStar(P, G, uGlobalPId, vGlobalPId);
                    if (aStarShortestPathDist <= t * euclideanDistanceUV) {
                        bridgesEuclideanDistance.emplace_back(aStarShortestPathDist);
                    } else {
                        bridgesEuclideanDistance.emplace_back(euclideanDistanceUV);
                        G[uGlobalPId].insert(vGlobalPId);
                        G[vGlobalPId].insert(uGlobalPId);
                        edgesByThread.emplace_back(uGlobalPId, vGlobalPId);
                        countE++;
                    }
                }
            }
        }
    }

    void greedyMergeLight(
            const QuadTreeLeaf &firstLeaf,
            const QuadTreeLeaf &secondLeaf,
            std::vector<std::unordered_set<unsigned>> &G,
            std::vector<Edge> &edgesByThread) {
        QuadTreeLeafInfo &firstLeafInfo = PointToLeafMap[keyOf(firstLeaf)], &secondLeafInfo = PointToLeafMap[keyOf(
                secondLeaf)];

        const std::vector<unsigned> &firstLeafInfoPIds = firstLeafInfo.internalPIDs,
                &secondLeafInfoPIds = secondLeafInfo.internalPIDs;

        std::unordered_map<unsigned, unsigned> firstLeafGlobalToLocalPId(firstLeafInfoPIds.size());
        std::vector<unsigned> firstLeafLocalToGlobalPId(firstLeafInfoPIds.size());

        std::unordered_map<unsigned, unsigned> secondLeafGlobalToLocal(secondLeafInfoPIds.size());
        std::vector<unsigned> secondLeafLocalToGlobal(secondLeafInfoPIds.size());

        // First Leaf points mapping
        for (unsigned i = 0; i < firstLeafInfoPIds.size(); ++i) {
            firstLeafGlobalToLocalPId.insert({firstLeafInfoPIds[i], i});
            firstLeafLocalToGlobalPId[i] = firstLeafInfoPIds[i];
        }
        // Second Leaf points mapping
        for (unsigned i = 0; i < secondLeafInfoPIds.size(); ++i) {
            secondLeafGlobalToLocal.insert({secondLeafInfoPIds[i], i});
            secondLeafLocalToGlobal[i] = secondLeafInfoPIds[i];
        }
        std::vector<std::vector<double>> firstLeafPointsIntraDist(firstLeafInfoPIds.size(),
                                                                  std::vector<double>(firstLeafInfoPIds.size(), 0.0));
        std::vector<std::vector<double>> secondLeafPointsIntraDist(secondLeafInfoPIds.size(),
                                                                   std::vector<double>(secondLeafInfoPIds.size(), 0.0));

        // First and second Leaf leader points info
        const unsigned &firstLeafLeaderP = firstLeafInfo.leaderPId, &secondLeafLeaderP = secondLeafInfo.leaderPId;
        unsigned firstLeafLeaderLocalPId = firstLeafGlobalToLocalPId[firstLeafLeaderP],
                secondLeafLeaderLocalPId = secondLeafGlobalToLocal[secondLeafLeaderP];

        std::vector<unsigned> firstLeafLocalPIds(firstLeafInfoPIds.size()), secondLeafLocalPIds(
                secondLeafInfoPIds.size());

        std::vector<double> firstLeafPointsToSecondLeaderDist(firstLeafInfoPIds.size(), DBL_MAX);

        for (unsigned i = 0; i < firstLeafInfoPIds.size(); i++) {
            firstLeafLocalPIds[i] = i;
            firstLeafPointsToSecondLeaderDist[i] = std::sqrt(
                    CGAL::squared_distance(P[firstLeafInfoPIds[i]], P[secondLeafLeaderP]));
        }

        std::vector<double> secondLeafPointsToFirstLeaderDist(secondLeafInfoPIds.size(), DBL_MAX);
        for (unsigned i = 0; i < secondLeafInfoPIds.size(); i++) {
            secondLeafLocalPIds[i] = i;
            secondLeafPointsToFirstLeaderDist[i] = std::sqrt(squared_distance(P[secondLeafInfoPIds[i]],
                                                                              P[firstLeafLeaderP]));
        }
        // Sort the first Leaf based on the second Leaf leader point
        std::sort(firstLeafLocalPIds.begin(), firstLeafLocalPIds.end(),
                  [&firstLeafPointsToSecondLeaderDist](const unsigned fId, const unsigned sId) {
                      return firstLeafPointsToSecondLeaderDist[fId] < firstLeafPointsToSecondLeaderDist[sId];
                  });
        // Sort the second Leaf based on the first Leaf leader point
        std::sort(secondLeafLocalPIds.begin(), secondLeafLocalPIds.end(),
                  [&secondLeafPointsToFirstLeaderDist](const unsigned fId, const unsigned sId) {
                      return secondLeafPointsToFirstLeaderDist[fId] < secondLeafPointsToFirstLeaderDist[sId];
                  });

        std::vector<Edge> leavesBridges;
        std::vector<double> bridgesEuclideanDistance;
        if (H[firstLeafLeaderP].find(secondLeafLeaderP) != H[firstLeafLeaderP].end()) {
            leavesBridges.emplace_back(firstLeafLeaderLocalPId, secondLeafLeaderLocalPId);
            bridgesEuclideanDistance.emplace_back(
                    std::sqrt(squared_distance(P[firstLeafLeaderP], P[secondLeafLeaderP])));
            for (unsigned i = 0; i < firstLeafInfoPIds.size(); i++) {
                firstLeafPointsIntraDist[i][firstLeafLeaderLocalPId] = firstLeafPointsIntraDist[firstLeafLeaderLocalPId][i] =
                        std::sqrt(CGAL::squared_distance(P[firstLeafInfoPIds[i]], P[firstLeafLeaderP]));
            }
            for (unsigned i = 0; i < secondLeafInfoPIds.size(); i++) {
                secondLeafPointsIntraDist[i][secondLeafLeaderLocalPId] = secondLeafPointsIntraDist[secondLeafLeaderLocalPId][i] =
                        std::sqrt(CGAL::squared_distance(P[secondLeafInfoPIds[i]], P[secondLeafLeaderP]));
            }
        }

        // Main processing inter Leaf pairs for merging of two Leafs
        for (unsigned firstLeafLocalPId: firstLeafLocalPIds) {
            unsigned u = firstLeafLocalPId, uGlobal = firstLeafLocalToGlobalPId[u];
            for (unsigned secondLeafLocalPId: secondLeafLocalPIds) {
                unsigned v = secondLeafLocalPId, vGlobal = secondLeafLocalToGlobal[v];
                double euclideanDistOfUV = std::sqrt(CGAL::squared_distance(P[uGlobal], P[vGlobal]));
                bool isNewBridgeNeeded = true;

                for (int b = (int) leavesBridges.size() - 1; b >= 0; b--) {
                    unsigned leavesBridgesFirstPLocalId = leavesBridges[b].first,
                            leavesBridgesSecondPLocalId = leavesBridges[b].second;
                    if (CGAL::is_zero(firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId])) {
                        firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId] = firstLeafPointsIntraDist[leavesBridgesFirstPLocalId][u] =
                                std::sqrt(CGAL::squared_distance(P[uGlobal],
                                                                 P[firstLeafLocalToGlobalPId[leavesBridgesFirstPLocalId]]));
                    }

                    if (CGAL::is_zero(secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v])) {
                        secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v] = secondLeafPointsIntraDist[v][leavesBridgesSecondPLocalId] =
                                std::sqrt(CGAL::squared_distance(P[vGlobal],
                                                                 P[secondLeafLocalToGlobal[leavesBridgesSecondPLocalId]]));
                    }
                    double leftSide = t * ((firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId] +
                                            secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v]) -
                                           euclideanDistOfUV);
                    if (leftSide <= -bridgesEuclideanDistance[b]) {
                        isNewBridgeNeeded = false;
                        break;
                    }
                }
                if (isNewBridgeNeeded) {
                    leavesBridges.emplace_back(u, v);
                    auto greedyHeuristicDist = DBL_MAX;
                    greedyHeuristicDist = greedyHeuristic(P, G, uGlobal, vGlobal);

                    if (greedyHeuristicDist <= t * euclideanDistOfUV) {
                        bridgesEuclideanDistance.emplace_back(greedyHeuristicDist);
                    } else {
                        double aStarShortestPathDist = aStar(P, G, uGlobal, vGlobal);
                        if (aStarShortestPathDist <= t * euclideanDistOfUV) {
                            bridgesEuclideanDistance.emplace_back(aStarShortestPathDist);
                        } else {
                            bridgesEuclideanDistance.emplace_back(euclideanDistOfUV);
                            G[uGlobal].insert(vGlobal);
                            G[vGlobal].insert(uGlobal);
                            edgesByThread.emplace_back(uGlobal, vGlobal);
//                            countE++;
                        }
                    }
                }
            }
        }
    }


    double computeStretchFactorOfTwoLeaves(
            const QuadTreeLeaf &firstLeaf,
            const QuadTreeLeaf &secondLeaf,
            std::vector<std::vector<double>> &firstLeafPointsIntraDist, bool isFirstLeafPointsIntraDistSaved = true) {
        double exactSFOfFSS = t;
        QuadTreeLeafInfo &firstLeafInfo = PointToLeafMap[keyOf(firstLeaf)], &secondLeafInfo = PointToLeafMap[keyOf(
                secondLeaf)];

        const std::vector<unsigned> &firstLeafInfoPIds = firstLeafInfo.internalPIDs,
                &secondLeafInfoPIds = secondLeafInfo.internalPIDs;

        std::unordered_map<unsigned, unsigned> firstLeafGlobalToLocalPId(firstLeafInfoPIds.size());
        std::vector<unsigned> firstLeafLocalToGlobalPId(firstLeafInfoPIds.size());

        std::unordered_map<unsigned, unsigned> secondLeafGlobalToLocal(secondLeafInfoPIds.size());
        std::vector<unsigned> secondLeafLocalToGlobal(secondLeafInfoPIds.size());

        // First Leaf points mapping
        for (unsigned i = 0; i < firstLeafInfoPIds.size(); ++i) {
            firstLeafGlobalToLocalPId.insert({firstLeafInfoPIds[i], i});
            firstLeafLocalToGlobalPId[i] = firstLeafInfoPIds[i];
        }
        // Second Leaf points mapping
        for (unsigned i = 0; i < secondLeafInfoPIds.size(); ++i) {
            secondLeafGlobalToLocal.insert({secondLeafInfoPIds[i], i});
            secondLeafLocalToGlobal[i] = secondLeafInfoPIds[i];
        }

        std::vector<std::vector<double>> secondLeafPointsIntraDist(secondLeafInfoPIds.size(),
                                                                   std::vector<double>(secondLeafInfoPIds.size(), 0.0));

        // First and second Leaf leader points info
        const unsigned &firstLeafLeaderP = firstLeafInfo.leaderPId, &secondLeafLeaderP = secondLeafInfo.leaderPId;
        unsigned firstLeafLeaderLocalPId = firstLeafGlobalToLocalPId[firstLeafLeaderP],
                secondLeafLeaderLocalPId = secondLeafGlobalToLocal[secondLeafLeaderP];

        std::vector<unsigned> firstLeafLocalPIds(firstLeafInfoPIds.size()), secondLeafLocalPIds(
                secondLeafInfoPIds.size());

        std::vector<double> firstLeafPointsToSecondLeaderDist(firstLeafInfoPIds.size(), DBL_MAX);

        for (unsigned i = 0; i < firstLeafInfoPIds.size(); i++) {
            firstLeafLocalPIds[i] = i;
            firstLeafPointsToSecondLeaderDist[i] = std::sqrt(
                    CGAL::squared_distance(P[firstLeafInfoPIds[i]], P[secondLeafLeaderP]));
        }

        std::vector<double> secondLeafPointsToFirstLeaderDist(secondLeafInfoPIds.size(), DBL_MAX);
        for (unsigned i = 0; i < secondLeafInfoPIds.size(); i++) {
            secondLeafLocalPIds[i] = i;
            secondLeafPointsToFirstLeaderDist[i] = std::sqrt(squared_distance(P[secondLeafInfoPIds[i]],
                                                                              P[firstLeafLeaderP]));
        }
        // Sort the first Leaf based on the second Leaf leader point
        std::sort(firstLeafLocalPIds.begin(), firstLeafLocalPIds.end(),
                  [&firstLeafPointsToSecondLeaderDist](const unsigned fId, const unsigned sId) {
                      return firstLeafPointsToSecondLeaderDist[fId] < firstLeafPointsToSecondLeaderDist[sId];
                  });
        // Sort the second Leaf based on the first Leaf leader point
        std::sort(secondLeafLocalPIds.begin(), secondLeafLocalPIds.end(),
                  [&secondLeafPointsToFirstLeaderDist](const unsigned fId, const unsigned sId) {
                      return secondLeafPointsToFirstLeaderDist[fId] < secondLeafPointsToFirstLeaderDist[sId];
                  });

        std::vector<Edge> leavesBridges;
        std::vector<double> bridgesED;
        if (H[firstLeafLeaderP].find(secondLeafLeaderP) != H[firstLeafLeaderP].end()) {
            leavesBridges.emplace_back(firstLeafLeaderLocalPId, secondLeafLeaderLocalPId);
            bridgesED.emplace_back(std::sqrt(squared_distance(P[firstLeafLeaderP], P[secondLeafLeaderP])));
            if (!isFirstLeafPointsIntraDistSaved) {
                for (unsigned i = 0; i < firstLeafInfoPIds.size(); i++) {
                    firstLeafPointsIntraDist[i][firstLeafLeaderLocalPId] = firstLeafPointsIntraDist[firstLeafLeaderLocalPId][i] =
                            std::sqrt(CGAL::squared_distance(P[firstLeafInfoPIds[i]], P[firstLeafLeaderP]));
                }
            }
            for (unsigned i = 0; i < secondLeafInfoPIds.size(); i++) {
                secondLeafPointsIntraDist[i][secondLeafLeaderLocalPId] = secondLeafPointsIntraDist[secondLeafLeaderLocalPId][i] =
                        std::sqrt(CGAL::squared_distance(P[secondLeafInfoPIds[i]], P[secondLeafLeaderP]));
            }
        }

        // Main processing inter Leaf pairs for merging of two Leafs
        for (unsigned firstLeafLocalPId: firstLeafLocalPIds) {
            unsigned u = firstLeafLocalPId, uGlobal = firstLeafLocalToGlobalPId[u];
            for (unsigned secondLeafLocalPId: secondLeafLocalPIds) {
                unsigned v = secondLeafLocalPId, vGlobal = secondLeafLocalToGlobal[v];
                double euclideanDistOfUV = std::sqrt(CGAL::squared_distance(P[uGlobal], P[vGlobal]));
                bool isNewBridgeNeeded = true;

                for (int b = (int) leavesBridges.size() - 1; b >= 0; b--) {
                    unsigned leavesBridgesFirstPLocalId = leavesBridges[b].first,
                            leavesBridgesSecondPLocalId = leavesBridges[b].second;
                    if (CGAL::is_zero(firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId])) {
                        firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId] = firstLeafPointsIntraDist[leavesBridgesFirstPLocalId][u] =
                                std::sqrt(CGAL::squared_distance(P[uGlobal],
                                                                 P[firstLeafLocalToGlobalPId[leavesBridgesFirstPLocalId]]));
                    }

                    if (CGAL::is_zero(secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v])) {
                        secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v] = secondLeafPointsIntraDist[v][leavesBridgesSecondPLocalId] =
                                std::sqrt(CGAL::squared_distance(P[vGlobal],
                                                                 P[secondLeafLocalToGlobal[leavesBridgesSecondPLocalId]]));
                    }
                    double leftSide = t * ((firstLeafPointsIntraDist[u][leavesBridgesFirstPLocalId] +
                                            secondLeafPointsIntraDist[leavesBridgesSecondPLocalId][v]) -
                                           euclideanDistOfUV);
                    if (leftSide <= -bridgesED[b]) {
                        isNewBridgeNeeded = false;
                        break;
                    }
                }

                if (isNewBridgeNeeded) {
                    leavesBridges.emplace_back(u, v);
                    auto greedyHeuristicDist = DBL_MAX;
                    greedyHeuristicDist = greedyHeuristic(P, H, uGlobal, vGlobal);

                    if (greedyHeuristicDist <= t * euclideanDistOfUV)
                        bridgesED.emplace_back(greedyHeuristicDist);
                    else {
                        double aStarShortestPathDist = aStar(P, H, uGlobal, vGlobal);
                        if (aStarShortestPathDist <= t * euclideanDistOfUV) {
                            bridgesED.emplace_back(aStarShortestPathDist);
                        } else
                            exactSFOfFSS = std::max(exactSFOfFSS, (double) aStarShortestPathDist / euclideanDistOfUV);
                    }
                }
            }
        }
        return exactSFOfFSS;
    }


    auto
    fastSparseSpanner(std::vector<Edge> &E, const double t, bool mayHaveDuplicates = true, unsigned threadCount = 1,
                      bool log = false) {
        struct FSSStats {
            double FSS_T;
            double FSS_M;
            double FSS_Avg_deg;
        };
        try {
            if (P.empty())
                throw std::runtime_error("Point set is empty! Aborting.");
            else if (P.size() < 2)
                throw std::runtime_error("Pointset must be greater than 1! Aborting.");

            if (t < 1)
                throw std::runtime_error("t must be > 1; Aborting!");
        }
        catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return FSSStats{0.0, 0.0, 0.0};
        }


        CGAL::Real_timer clock;
        this->t = t;
        clock.start();


        if (mayHaveDuplicates) {
            std::unordered_set<Point> S;
            for (auto p: P)
                S.insert(p);

            if (S.size() != P.size())
                throw std::runtime_error("Point set contains duplicates! Aborting.");
        }
        if (log)
            std::cout << timeStamp()
                      << "--------------------------------------FSS Start--------------------------------------------------------"
                      << std::endl;

        FSSStats fssStats{};
        double previousTime = 0.0, fssComputationTime;
        unsigned countE = 0, previousEdgeCount = 0, mergingCnt = 0;
        PointToIndexMap.reserve(P.size());
        H.resize(P.size());
        if (log)
            std::cout << timeStamp() << "Step 1: Quad tree construction has started..." << std::endl;

        T->refine(UINT_MAX, 2500);

        if (log) {
            std::cout << timeStamp() << "Done. Took " << clock.time() - previousTime << "s = "
                      << (clock.time() - previousTime) / 60.0 << "m = " << (clock.time() - previousTime) / 60.0 / 60.0
                      << "h " << std::endl << std::endl;
            previousTime = clock.time();
        }

        if (log)
            std::cout << timeStamp() << "Step 2: Greedy construction and leader points computation has started..."
                      << std::endl;

        for (unsigned i = 0; i < P.size(); ++i)
            PointToIndexMap[P[i]] = i;

        unsigned leafIndexCounter = 0;
        std::vector<Point> leaderPoints;
        leaderPoints.reserve(P.size() / 2500);

        std::unordered_map<Point, QuadTreeLeaf> leaderPointToLeafLeaf;
        leaderPointToLeafLeaf.reserve(P.size() / 2500);

        std::vector<std::unordered_set<unsigned >> G_T;
        G_T.reserve(P.size() / 2500);

        leafIdToLeaf.reserve(P.size() / 2500);
        size_t sizeOfLargestLeaf = 0;

        for (auto &curLeaf: T->traverse<CGAL::Orthtrees::Leaves_traversal>()) {
            PointToLeafMap.insert({keyOf(curLeaf), QuadTreeLeafInfo()});
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap[keyOf(curLeaf)];
            if (!curLeaf.is_null() && !curLeaf.empty()) sizeOfLargestLeaf = std::max(sizeOfLargestLeaf, curLeaf.size());
            curLeafInfo.leafID = leafIndexCounter++;
            leafIdToLeaf.emplace_back(curLeaf);
            G_T.emplace_back(std::unordered_set<unsigned>());
        }

        std::vector<std::vector<std::vector<Edge>>> step1Edges(threadCount,
                                                               std::vector<std::vector<Edge>>(leafIdToLeaf.size(),
                                                                                              std::vector<Edge>()));
        std::vector<std::vector<Point>> leaderPointsPL(threadCount, std::vector<Point>());
#pragma omp parallel for default(none) shared(std::cout, PointToLeafMap, sizeOfSmallestLeaf, leaderPointToLeafLeaf, leaderPointsPL, step1Edges, leafIndexCounter, sizeOfLargestLeaf, t) reduction(+:countE) schedule(dynamic) num_threads(threadCount)
        for (auto &curLeaf: leafIdToLeaf) {
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap[keyOf(curLeaf)];
            curLeafInfo.bboxOfLeaf = T->bbox(curLeaf);
            sizeOfSmallestLeaf = std::min(curLeafInfo.bboxOfLeaf.xmax() - curLeafInfo.bboxOfLeaf.xmin(),
                                          sizeOfSmallestLeaf);
            if (!curLeaf.is_null() && !curLeaf.empty()) {
                for (Point const &p: curLeaf)
                    curLeafInfo.internalPIDs.emplace_back(PointToIndexMap[p]);
                leaderPointsPL[omp_get_thread_num()].emplace_back(computeLeader(curLeaf));
                step1Edges[omp_get_thread_num()].emplace_back(
                        constructFG_GreedySpanner(P, curLeaf, *T, PointToLeafMap, t, countE));
            }
        }
        for (auto &lPLItr: leaderPointsPL) {
            for (auto &lP: lPLItr) {
                leaderPoints.emplace_back(lP);
            }
        }
        for (auto &step1Edge: step1Edges) {
            for (auto &edgeV: step1Edge) {
                for (auto &edge: edgeV) {
                    H[edge.first].insert(edge.second);
                    H[edge.second].insert(edge.first);
                }
            }
        }

        totalNumOfLeaf = leafIndexCounter;

        double enlargeFactor = sizeOfSmallestLeaf * 0.25;
        for (QuadTreeLeaf const &curLeaf: T->traverse<CGAL::Orthtrees::Leaves_traversal>()) {
            if (!curLeaf.is_null()) {
                QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(curLeaf));
                CGAL::Bbox_2 curLeafBbox = enlargeBbox(curLeafInfo.bboxOfLeaf, enlargeFactor);
                curLeafInfo.extendedBboxOfLeaf = curLeafBbox;
                curLeafInfo.centerOfBbox = {
                        (curLeafInfo.extendedBboxOfLeaf.xmin() + curLeafInfo.extendedBboxOfLeaf.xmax()) / 2.0,
                        (curLeafInfo.extendedBboxOfLeaf.ymin() + curLeafInfo.extendedBboxOfLeaf.ymax()) / 2.0};
                curLeafInfo.bboxOfVertices.emplace_back(curLeafInfo.extendedBboxOfLeaf.xmin(),
                                                        curLeafInfo.extendedBboxOfLeaf.ymin());
                curLeafInfo.bboxOfVertices.emplace_back(curLeafInfo.extendedBboxOfLeaf.xmin(),
                                                        curLeafInfo.extendedBboxOfLeaf.ymax());
                curLeafInfo.bboxOfVertices.emplace_back(curLeafInfo.extendedBboxOfLeaf.xmax(),
                                                        curLeafInfo.extendedBboxOfLeaf.ymin());
                curLeafInfo.bboxOfVertices.emplace_back(curLeafInfo.extendedBboxOfLeaf.xmax(),
                                                        curLeafInfo.extendedBboxOfLeaf.ymax());
                std::vector<QuadTreeLeaf> allAdjacentNeighbors;
                intersectedLeafs(curLeafBbox, T->root(), allAdjacentNeighbors);
                for (const auto &neighbor: allAdjacentNeighbors) {
                    if (!neighbor.is_null() &&
                        std::find(leafIdToLeaf.begin(), leafIdToLeaf.end(), neighbor) != leafIdToLeaf.end()) {
                        QuadTreeLeafInfo &neighborLeafInfo = PointToLeafMap.at(keyOf(neighbor));
                        if (curLeafInfo.leafID != neighborLeafInfo.leafID)
                            G_T[curLeafInfo.leafID].insert(neighborLeafInfo.leafID);
                    }
                }
            }
        }

        previousEdgeCount = (countE - previousEdgeCount);
        if (log) {
            std::cout << timeStamp() << "Done. Took " << clock.time() - previousTime << "s = "
                      << (clock.time() - previousTime) / 60.0 << "m = " << (clock.time() - previousTime) / 60.0 / 60.0
                      << "h " << std::endl;
            std::cout << timeStamp() << "Placed " << previousEdgeCount << " Edges" << ", Added avg. deg : "
                      << (double) (2 * previousEdgeCount) / (double) P.size() << std::endl << std::endl;
        }

        previousTime = clock.time();
        previousEdgeCount = countE;

        std::vector<Edge> wspdEdges;
        double WSPD_SF = (t >= 1.1 && t <= 1.25) ? 1.25 : t;
        if (log) {
            std::cout << timeStamp() << "Step 3: WSPD spanner construction has started... " << std::endl;
            std::cout << timeStamp() << "tPrime: " << WSPD_SF << std::endl;
        }
        if (leafIndexCounter >= 2)
            constructWSPDSpanner(leaderPoints, H, countE, PointToIndexMap, WSPD_SF, wspdEdges);

        previousEdgeCount = (countE - previousEdgeCount);
        if (log) {
            std::cout << timeStamp() << "Done. Took " << clock.time() - previousTime << "s = "
                      << (clock.time() - previousTime) / 60.0 << "m = " << (clock.time() - previousTime) / 60.0 / 60.0
                      << "h " << std::endl;
            std::cout << timeStamp() << "Placed " << previousEdgeCount << " Edges" << ", Added avg. deg : "
                      << (double) (2 * previousEdgeCount) / (double) P.size() << std::endl << std::endl;
        }

        previousTime = clock.time();
        previousEdgeCount = countE;

        if (log)
            std::cout << timeStamp() << "Step 4: Merging neighbors in G_T... " << std::endl;

        // first layer of all Leafs
        double phaseOneEdgeCount = 0.0;
        std::vector<std::pair<QuadTreeLeaf, QuadTreeLeaf>> nonDiagonalPairs, diagonalPairs;
        for (QuadTreeLeaf const &curLeaf: T->traverse<CGAL::Orthtrees::Leaves_traversal>()) {
            if (!curLeaf.is_null() && !curLeaf.empty()) {
                QuadTreeLeafInfo &curLeafInfo = PointToLeafMap[keyOf(curLeaf)];
                std::unordered_set<unsigned> &neighbors = G_T[curLeafInfo.leafID];
                for (const auto &neighbor: neighbors) {
                    auto &neighborLeaf = leafIdToLeaf[neighbor];
                    if (!neighborLeaf.is_null() && !neighborLeaf.empty()) {
                        if (isDiagonal(curLeaf, neighborLeaf) && curLeafInfo.leafID < neighbor)
                            diagonalPairs.emplace_back(curLeaf, neighborLeaf);
                        else if (curLeafInfo.leafID < neighbor) nonDiagonalPairs.emplace_back(curLeaf, neighborLeaf);
                    }
                }
            }
        }

        // non-diagonal pairs
        std::vector<std::vector<std::unordered_set<unsigned>>> GOfPhaseOneND(threadCount, H);
        std::vector<std::vector<Edge>> phaseOneEdgesND(threadCount, std::vector<Edge>());
        std::vector<std::unordered_set<index_tPair, index_tPairHash >> phaseOneFSSMergedLeafPairsOfThreadsND(
                threadCount, FSSMergedLeafPairs);
#pragma omp parallel for default(none) shared(nonDiagonalPairs, phaseOneFSSMergedLeafPairsOfThreadsND, GOfPhaseOneND, phaseOneEdgesND) reduction(+:countE, mergingCnt) schedule(dynamic) num_threads(threadCount)
        for (auto &nonDiagonalPair: nonDiagonalPairs) {
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(nonDiagonalPair.first)),
                    &neighborInfo = PointToLeafMap[keyOf(nonDiagonalPair.second)];
            greedyMerge(nonDiagonalPair.first, nonDiagonalPair.second, countE, GOfPhaseOneND[omp_get_thread_num()],
                        phaseOneEdgesND[omp_get_thread_num()]);
            mergingCnt++;
        }

        for (auto &phaseOneEdge: phaseOneEdgesND) {
            for (auto &edge: phaseOneEdge) {
                H[edge.first].insert(edge.second);
                H[edge.second].insert(edge.first);
            }
        }
        for (auto &nonDiagonalPair: nonDiagonalPairs) {
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(nonDiagonalPair.first)),
                    &neighborInfo = PointToLeafMap[keyOf(nonDiagonalPair.second)];
            FSSMergedLeafPairs.insert({curLeafInfo.leafID, neighborInfo.leafID});
        }

        GOfPhaseOneND.clear();
        phaseOneEdgesND.clear();
        phaseOneFSSMergedLeafPairsOfThreadsND.clear();

        std::vector<std::vector<std::unordered_set<unsigned>>> GOfPhaseOneD(threadCount, H);
        std::vector<std::vector<Edge>> phaseOneEdgesD(threadCount, std::vector<Edge>());
        std::vector<std::unordered_set<index_tPair, index_tPairHash >> phaseOneFSSMergedLeafPairsOfThreadsD(threadCount,
                                                                                                            FSSMergedLeafPairs);
#pragma omp parallel for default(none) shared(diagonalPairs, phaseOneFSSMergedLeafPairsOfThreadsD, GOfPhaseOneD, phaseOneEdgesD) reduction(+:countE, mergingCnt) schedule(dynamic) num_threads(threadCount)
        for (auto &diagonalPair: diagonalPairs) {
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(diagonalPair.first)),
                    &neighborInfo = PointToLeafMap[keyOf(diagonalPair.second)];
            greedyMerge(diagonalPair.first, diagonalPair.second, countE, GOfPhaseOneD[omp_get_thread_num()],
                        phaseOneEdgesD[omp_get_thread_num()]);
            mergingCnt++;
        }
        for (auto &phaseOneEdge: phaseOneEdgesD) {
            for (auto &edge: phaseOneEdge) {
                H[edge.first].insert(edge.second);
                H[edge.second].insert(edge.first);
            }
        }
        for (auto &diagonalPair: diagonalPairs) {
            QuadTreeLeafInfo &curLeafInfo = PointToLeafMap.at(keyOf(diagonalPair.first)),
                    &neighborInfo = PointToLeafMap[keyOf(diagonalPair.second)];
            FSSMergedLeafPairs.insert({curLeafInfo.leafID, neighborInfo.leafID});
        }

        GOfPhaseOneD.clear();
        phaseOneEdgesD.clear();
        phaseOneFSSMergedLeafPairsOfThreadsD.clear();

        phaseOneEdgeCount = previousEdgeCount = (countE - previousEdgeCount);
        if (log) {
            std::cout << timeStamp() << "Done. Took " << clock.time() - previousTime << "s = "
                      << (clock.time() - previousTime) / 60.0 << "m = " << (clock.time() - previousTime) / 60.0 / 60.0
                      << "h " << std::endl;
            std::cout << timeStamp() << "Placed " << previousEdgeCount << " Edges" << ", Added avg. deg : "
                      << (double) (2 * previousEdgeCount) / (double) P.size() << std::endl << std::endl;
        }
        previousTime = clock.time();
        previousEdgeCount = countE;

        if (log)
            std::cout << timeStamp() << "Step 5: Light merging distant leaves... " << std::endl;

        // second layer to onward
        unsigned numberOfLayers;
        if (t >= 2) numberOfLayers = 1;
        else if (t > 1.25) numberOfLayers = 3;
        else if (t > 1.05) numberOfLayers = 5;
        else if (t > 1.04) numberOfLayers = 6;
        else if (t > 1.03) numberOfLayers = 7;
        else if (t > 1.02) numberOfLayers = 8;
        else if (t > 1.01) numberOfLayers = 9;
        else numberOfLayers = 10;

        double phaseTwoEdgeCount = 0.0;
        for (unsigned nextLayer = 2; nextLayer <= numberOfLayers; nextLayer++) {
            std::vector<std::vector<Edge>> phaseTwoEdges(threadCount, std::vector<Edge>());
            std::vector<std::unordered_set<index_tPair, index_tPairHash>> FSSMergedLeafPairsOfThreads(threadCount,
                                                                                                      FSSMergedLeafPairs);
            std::vector<std::vector<std::unordered_set<unsigned>>> GOfPhase2(threadCount, H);
#pragma omp parallel for default(none) shared(GOfPhase2, G_T, leafIndexCounter, nextLayer, threadCount, previousEdgeCount, std::cout, phaseTwoEdges, FSSMergedLeafPairsOfThreads) reduction(+:mergingCnt, countE) schedule(dynamic) num_threads(threadCount)
            for (QuadTreeLeaf const &curLeaf: leafIdToLeaf) {
                if (!curLeaf.is_null() && !curLeaf.empty()) {
                    QuadTreeLeafInfo &curLeafInfo = PointToLeafMap[keyOf(curLeaf)];
                    std::unordered_set<unsigned> &neighbors = G_T[curLeafInfo.leafID];
                    std::vector<QuadTreeLeaf> nextLayerNeighbors;
                    boost::dynamic_bitset<> processedLeafs(leafIndexCounter);
                    processedLeafs[curLeafInfo.leafID] = true;
                    for (auto &neighbor: neighbors) {
                        nextLayerNeighbors.emplace_back(leafIdToLeaf[neighbor]);
                        processedLeafs[neighbor] = true;
                    }
                    std::vector<QuadTreeLeaf> expectedLayerNeighbors;
                    unsigned layerNo = 1;
                    while (layerNo <= nextLayer - 1) {
                        layerNo++;
                        getNextLayerNeighbors(G_T, leafIdToLeaf, processedLeafs, nextLayerNeighbors,
                                              expectedLayerNeighbors, layerNo, nextLayer);

                    }
                    std::vector<std::vector<double>> firstLeafDistance(curLeafInfo.internalPIDs.size(),
                                                                       std::vector<double>(
                                                                               curLeafInfo.internalPIDs.size(), 0.0));
                    unsigned prevMergingCount = mergingCnt;
                    for (auto &nextLayerNeighbor: expectedLayerNeighbors) {
                        QuadTreeLeafInfo &neighborInfo = PointToLeafMap[keyOf(nextLayerNeighbor)];
                        if (curLeafInfo.leafID < neighborInfo.leafID && !leafIdToLeaf[neighborInfo.leafID].empty()) {
                            unsigned firstPart = curLeafInfo.leafID, secondPart = neighborInfo.leafID;
                            index_tPair checkPair(std::min(firstPart, secondPart), std::max(firstPart, secondPart));
                            FSSMergedLeafPairsOfThreads[omp_get_thread_num()].insert(checkPair);
                            greedyMergeLight(curLeaf, nextLayerNeighbor, GOfPhase2[omp_get_thread_num()],
                                             phaseTwoEdges[omp_get_thread_num()]);
                            mergingCnt++;
                        }
                    }
                }
            }
            for (auto &phaseTwoEdge: phaseTwoEdges) {
                for (auto &edge: phaseTwoEdge) {
                    H[edge.first].insert(edge.second);
                    H[edge.second].insert(edge.first);
                    countE++;
                }
            }
            for (auto &fssMergedLeafPairsPerThread: FSSMergedLeafPairsOfThreads)
                for (auto &mergedLeafPair: fssMergedLeafPairsPerThread)
                    FSSMergedLeafPairs.insert(mergedLeafPair);
        }

        phaseTwoEdgeCount = previousEdgeCount = (countE - previousEdgeCount);

        if (log) {
            std::cout << timeStamp() << "Done. Took " << clock.time() - previousTime << "s = "
                      << (clock.time() - previousTime) / 60.0 << "m = " << (clock.time() - previousTime) / 60.0 / 60.0
                      << "h " << std::endl;
            std::cout << timeStamp() << "Placed " << previousEdgeCount << " Edges" << ", Added avg. deg : "
                      << (double) (2 * previousEdgeCount) / (double) P.size() << std::endl;
        }

        previousEdgeCount = countE;
        previousTime = clock.time();

        leaderPoints.clear();
        leaderPointToLeafLeaf.clear();
        E.clear();

        for (unsigned i = 0; i < H.size(); i++)
            for (unsigned const j: H[i])
                if (i < j)
                    E.emplace_back(i, j);
        expectedStretchFactor = t;
        isFSSCompleted = true;

        clock.stop();
        fssComputationTime = clock.time();
        fssStats.FSS_T = fssComputationTime;
        fssStats.FSS_Avg_deg = double(2 * E.size()) / (double) P.size();

        if (log)
            std::cout << timeStamp()
                      << "--------------------------------------FSS End--------------------------------------------------------"
                      << std::endl;
return fssStats;
    }

    double fastStretchFactor(unsigned threadCount = 1) {

        if (fastStretchFactor_t != DBL_MIN)
            return 0.0;

        if (!isFSSCompleted) {
            std::cout << "FSS is not completed yet. Need to construct FSS first." << std::endl;
            return 0.0;
        }
        if (isFSFCompleted) return 0.0;

        CGAL::Real_timer fssExactSFComputation;
        fssExactSFComputation.start();

        std::vector<double> exactSFOfFSS(threadCount, expectedStretchFactor);
        unsigned mergingCnt = 0;
        fastStretchFactor_t = expectedStretchFactor;
#pragma omp parallel for default(none) shared(H, expectedStretchFactor, T, PointToLeafMap, exactSFOfFSS, std::cout, FSSMergedLeafPairs, mergingCnt, leafIdToLeaf, totalNumOfLeaf) schedule(dynamic) num_threads(threadCount)
        for (unsigned i = 0; i < totalNumOfLeaf; i++) {
            QuadTreeLeaf &firstLeaf = leafIdToLeaf[i];
            if (firstLeaf.is_null()) continue;
            std::vector<std::vector<double>> firstLeafPointsIntraDist(firstLeaf.size(),
                                                                      std::vector<double>(firstLeaf.size(), 0.0));
            unsigned prevMergingCnt = mergingCnt;
            for (unsigned j = i + 1; j < totalNumOfLeaf; j++) {
                index_tPair checkPair(std::min(i, j), std::max(i, j));
                if (FSSMergedLeafPairs.find(checkPair) == FSSMergedLeafPairs.end()) {
                    QuadTreeLeaf &secondLeaf = leafIdToLeaf[j];
                    if (secondLeaf.is_null()) continue;
                    if (!firstLeaf.empty() && !secondLeaf.empty()) {
                        exactSFOfFSS[omp_get_thread_num()] = std::max(exactSFOfFSS[omp_get_thread_num()],
                                                                      computeStretchFactorOfTwoLeaves(
                                                                              firstLeaf,
                                                                              secondLeaf,
                                                                              firstLeafPointsIntraDist,
                                                                              (mergingCnt - prevMergingCnt) > 0));
                        mergingCnt++;
                    }
                }
            }
        }

        fastStretchFactor_t = expectedStretchFactor;
        for (double sf_t: exactSFOfFSS)
            fastStretchFactor_t = std::max(sf_t, fastStretchFactor_t);

        isFSFCompleted = true;
        fssExactSFComputation.stop();

        return fastStretchFactor_t;
    }
};

#endif //FASTSPARSESPANNER_FSS_H
