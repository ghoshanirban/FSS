#ifndef FASTSPARSESPANNER_WSPDSPANNER_H
#define FASTSPARSESPANNER_WSPDSPANNER_H

#include "../Utilities.h"

typedef K::Point_2 Point;

class linearithmicWSPD {
    struct FastSplitTreeNode;
    unsigned dimension;

    struct DllNode {
        unsigned long pointID;
        std::vector<std::list<DllNode>::iterator> crossPtr, ptrBuffer;
        std::list<DllNode>::iterator crossPtrToTheSameNodeInTheCopy = std::list<DllNode>::iterator{};
        FastSplitTreeNode *ptrToLeaf = nullptr;

        DllNode(unsigned long id, unsigned dimension) {
            pointID = id;
            crossPtr.resize(dimension, std::list<DllNode>::iterator{});
        };
    };

    struct FastSplitTreeNode {
        CGAL::Bbox_2 R; // CGAL::Bbox_2(0,0,0,0);
        bool isBoundingBoxSet = false;
        FastSplitTreeNode *left = nullptr, *right = nullptr;
        unsigned long pointIDofRep = ULONG_MAX; // a representative point for the node (can be used for other spanner algorithms)
        std::vector<std::list<DllNode> > LS;

        explicit FastSplitTreeNode(const unsigned dimension) {
            LS = std::vector<std::list<DllNode> >{dimension};
        }
    };

    struct coordinateComparator {
        unsigned i;
        std::vector<Point> *pointSet;

        explicit coordinateComparator(unsigned whichCoordinate, std::vector<Point> *P) {
            i = whichCoordinate; pointSet = P;
        }

        bool operator()(const DllNode &first, const DllNode &second) const {
            return ( (*pointSet)[first.pointID].cartesian((int)i) < (*pointSet)[second.pointID].cartesian((int) i));
        }
    };

    std::vector<Point> P;
    double sep;

    inline static Point centerOfBox(const CGAL::Bbox_2 &B) {
        return {B.xmin() + (B.xmax() - B.xmin()) / 2, B.ymin() + (B.ymax() - B.ymin()) / 2};
    }

    inline static double radiusOfBox(const CGAL::Bbox_2 &B) {
        return L2Distance(centerOfBox(B), Point(B.xmax(), B.ymax()));
    }

    inline static double distanceBetweenBoxes(const CGAL::Bbox_2 &B1, const CGAL::Bbox_2 &B2) {
        double commonRadius = std::max(radiusOfBox(B1), radiusOfBox(B2));
        return L2Distance(centerOfBox(B1), centerOfBox(B2)) - 2 * commonRadius ;
    }

    bool static isWellSeparated(const CGAL::Bbox_2 &B1, const CGAL::Bbox_2 &B2, const double s) {
        return (distanceBetweenBoxes(B1, B2) >= s * std::max(radiusOfBox(B1), radiusOfBox(B2)));
    }

    inline static double lMax(const CGAL::Bbox_2 &B) {
        return std::max((B.xmax() - B.xmin()), (B.ymax() - B.ymin()));
    }

    void findPairs(const FastSplitTreeNode *v, const FastSplitTreeNode *w, std::vector<Edge> &pairs) {
        if (isWellSeparated(v->R, w->R, sep))
            pairs.emplace_back(std::make_pair<unsigned, unsigned>(std::min(v->pointIDofRep, w->pointIDofRep),
                                                                  std::max(v->pointIDofRep, w->pointIDofRep)));
        else if (lMax(v->R) <= lMax(w->R)) {
            if (w->left  != nullptr) findPairs(v, w->left, pairs);
            if (w->right != nullptr) findPairs(v, w->right, pairs);
        } else {
            if (v->left  != nullptr) findPairs(v->left, w, pairs);
            if (v->right != nullptr) findPairs(v->right, w, pairs);
        }
    }

    inline unsigned findTheLongestSide (const CGAL::Bbox_2 &B) const {
        unsigned longestSide = 0;
        double lengthOflongestSide = B.max(0) - B.min(0);
        for (int i = 1; i < (int) dimension; ++i)
            if (B.max(i) - B.min(i) > lengthOflongestSide) {
                longestSide = i;
                lengthOflongestSide = B.max(i) - B.min(i);
            }
        return longestSide;
    }

    inline void copyLists(std::vector<std::list<DllNode>> &LS, std::vector<std::list<DllNode>> &CLS) {

        // setup space for cross-pointing in CLS; for this we are using extra space in the nodes of LS[0]
        for (DllNode &node: LS[0])
            node.ptrBuffer.resize(dimension); // to be used to set up the crosspointers in CLS

        CLS.resize(dimension);

        for (unsigned i = 0; i < dimension; ++i)
            for (DllNode &node: LS[i]) {
                CLS[i].emplace_back(node.pointID, dimension);
                node.crossPtrToTheSameNodeInTheCopy = prev(CLS[i].end());
                node.crossPtr[0]->ptrBuffer[i] = prev(CLS[i].end());
            }

        // set up the crosspointers in CLS
        for (DllNode &node: LS[0]) {
            for (unsigned i = 0; i < dimension; ++i)
                for (unsigned j = 0; j < dimension; ++j)
                    node.ptrBuffer[i]->crossPtr[j] = node.ptrBuffer[j];
        }
    }

    void constructPartialSplitTreeRecursive(FastSplitTreeNode &u, const CGAL::Bbox_2 &R,
                                            std::vector<std::list<DllNode>> &LS, std::vector<std::list<DllNode>> &CLS,
                                            std::list<FastSplitTreeNode *> &leavesInPartialSplitTree, unsigned size) {
        // step 2
        if ( size <= CLS[0].size() / 2 ) {
            for (unsigned i = 0; i < dimension; ++i)
                for (const DllNode &node: LS[i])
                    node.crossPtrToTheSameNodeInTheCopy->ptrToLeaf = &u;

            leavesInPartialSplitTree.push_back(&u);
            return ;
        }

        double xMin = P[LS[0].front().pointID].x(), xMax = P[LS[0].back().pointID].x();
        double yMin = P[LS[1].front().pointID].y(), yMax = P[LS[1].back().pointID].y();
        u.R = CGAL::Bbox_2(xMin, yMin, xMax, yMax);
        u.isBoundingBoxSet = true;

        int i = (int) findTheLongestSide(u.R);
        double H = u.R.min(i) + (u.R.max(i) - u.R.min(i)) / 2;

        auto p = LS[i].begin(),     pPrime = std::next(p),
                q = prev(LS[i].end()), qPrime = std::prev(q);

        unsigned sizePrime = 1;

        while (P[pPrime->pointID].cartesian(i) <= H && P[qPrime->pointID].cartesian(i) >= H) {
            p = pPrime; pPrime = std::next(p);
            q = qPrime; qPrime = std::prev(q);
            sizePrime++;
        }

        auto *v = new FastSplitTreeNode(dimension), *w = new FastSplitTreeNode(dimension);
        u.left = v; u.right = w;
        u.pointIDofRep = LS[0].front().pointID;
        v->isBoundingBoxSet = w->isBoundingBoxSet = true;

        if (i == 0) {
            v->R = CGAL::Bbox_2(R.xmin(), R.ymin(), H, R.ymax());
            w->R = CGAL::Bbox_2(H, R.ymin(), R.xmax(), R.ymax());
        } else if (i == 1) {
            v->R = CGAL::Bbox_2(R.xmin(), R.ymin(), R.xmax(), H);
            w->R = CGAL::Bbox_2(R.xmin(), H, R.xmax(), R.ymax());
        }

        if (P[pPrime->pointID].cartesian(i) > H) {
            std::vector<std::list<DllNode>::iterator> toBeErased;
            toBeErased.reserve(LS[i].size());

            for(auto it = LS[i].begin(); it != LS[i].end(); ++it) {
                DllNode &z = *it;
                for (unsigned j = 0; j < dimension; ++j) {
                    z.crossPtrToTheSameNodeInTheCopy->crossPtr[j]->ptrToLeaf = v;
                    if ((int) j != i)
                        LS[j].erase(z.crossPtr[j]);
                }
                toBeErased.push_back(it);
                if (z.pointID == p->pointID)
                    break;
            }
            for(auto it : toBeErased )
                LS[i].erase(it);
            leavesInPartialSplitTree.emplace_back(v);
            constructPartialSplitTreeRecursive(*w, w->R, LS, CLS, leavesInPartialSplitTree, size - sizePrime);

        } else {
            std::vector<std::list<DllNode>::iterator> toBeErased;
            toBeErased.reserve(LS[i].size());

            for (auto it = prev(LS[i].end()); ; --it) {
                DllNode &z = *it;
                for (unsigned j = 0; j < dimension; ++j) {
                    z.crossPtrToTheSameNodeInTheCopy->crossPtr[j]->ptrToLeaf = w;
                    if ((int) j != i)
                        LS[j].erase(z.crossPtr[j]);
                }
                toBeErased.push_back(it);
                if (z.pointID == q->pointID || it == LS[i].begin())
                    break;
            }

            for(auto it : toBeErased )
                LS[i].erase(it);
            leavesInPartialSplitTree.emplace_back(w);
            constructPartialSplitTreeRecursive(*v, v->R, LS, CLS, leavesInPartialSplitTree, size - sizePrime);
        }
    }

    void constructPartialSplitTree(FastSplitTreeNode &root, const CGAL::Bbox_2 &R, std::vector<std::list<DllNode>> &LS,
                                   std::list<FastSplitTreeNode *> &leaves) {

        if (LS[0].size() == 1) {
            root.pointIDofRep = LS[0].front().pointID;
            Point p = P[LS[0].front().pointID];
            root.R = CGAL::Bbox_2(p.x(), p.y(), p.x(), p.y());
            root.isBoundingBoxSet = true;
            return;
        }

        std::vector<std::list<DllNode>> CLS;
        copyLists(LS, CLS);

        constructPartialSplitTreeRecursive(root, R, LS, CLS, leaves, LS[0].size());

        // allocate space for setting up cross-pointers for nodes of LS
        for (DllNode &node: CLS[0])
            node.ptrBuffer.resize(dimension);

        for (unsigned i = 0; i < dimension; ++i)
            for (DllNode &node: CLS[i]) {
                FastSplitTreeNode *u = node.ptrToLeaf;
                DllNode newNode(node.pointID, dimension);
                u->LS[i].emplace_back(newNode);
                node.crossPtr[0]->ptrBuffer[i] = prev(u->LS[i].end());
            }

        // set up the crosspointers in LS using the buffers in the nodes of CLS[0]
        for (DllNode &node: CLS[0]) {
            for (unsigned i = 0; i < dimension; ++i)
                for (unsigned j = 0; j < dimension; ++j)
                    node.ptrBuffer[i]->crossPtr[j] = node.ptrBuffer[j];

            if( ! node.ptrToLeaf->isBoundingBoxSet ) {
                auto u = node.ptrToLeaf;
                double xMin = P[u->LS[0].front().pointID].x(), xMax = P[u->LS[0].back().pointID].x(),
                        yMin = P[u->LS[1].front().pointID].y(), yMax = P[u->LS[1].back().pointID].y();
                node.ptrToLeaf->R = CGAL::Bbox_2(xMin, yMin, xMax, yMax);
            }
        }
    }

    void constructFastSplitTree(FastSplitTreeNode &root,std::unordered_map<Point, unsigned> &PointToIndex) {
        // set up the LS[i]s and the cross-pointers
        std::vector<std::list<DllNode>> LS(dimension);
        for (unsigned i = 0; i < dimension; ++i)
            LS.emplace_back(std::list<DllNode>());

        for (Point &p: P) {
            for (unsigned i = 0; i < dimension; ++i)
                LS[i].emplace_back(DllNode(PointToIndex[p], dimension));

            for (unsigned i = 0; i < dimension; ++i)
                for (unsigned j = 0; j < dimension; ++j)
                    LS[i].back().crossPtr[j] = prev(LS[j].end());
        }


        // sort the LS[i]s
        for (unsigned i = 0; i < dimension; ++i)
            LS[i].sort(coordinateComparator(i, &P));

        double xMin = P[LS[0].front().pointID].x(), xMax = P[LS[0].back().pointID].x(),
               yMin = P[LS[1].front().pointID].y(), yMax = P[LS[1].back().pointID].y();
        root.R = CGAL::Bbox_2(xMin, yMin, xMax, yMax);
        root.isBoundingBoxSet = true;
        std::list<FastSplitTreeNode*> leaves;
        constructPartialSplitTree(root, root.R, LS, leaves);
        while( !leaves.empty() ) {
            auto leaf = leaves.front();
            leaves.pop_front();

            if ( !leaf->LS[0].empty() )
                constructPartialSplitTree(*leaf, leaf->R, leaf->LS, leaves);
        }
    }

public:
    linearithmicWSPD(const std::vector<Point> &pointSet, const double separationRatio, const unsigned dimension) {
        assert( pointSet.size() >= 2 );
        P = pointSet;
        sep = separationRatio;
        this->dimension = dimension;
    }

    void computePairsAndStoreIn(std::vector<Edge> &pairs, std::unordered_map<Point, unsigned> &PointToIndex) {

        pairs.reserve(50 * P.size());

        FastSplitTreeNode root(dimension);
        constructFastSplitTree(root,PointToIndex);

        std::queue<FastSplitTreeNode*> Q;
        Q.push(&root);

        while (!Q.empty()) {
            FastSplitTreeNode *current = Q.front();
            Q.pop();

            if (current->left != nullptr)  Q.push(current->left);
            if (current->right != nullptr) Q.push(current->right);
            if (current->left != nullptr && current->right != nullptr)
                findPairs(current->left, current->right, pairs);
        }
    }
};

void constructWSPDSpanner(std::vector<Point> &leaderPoints,
                                       std::vector<std::unordered_set<unsigned>> &G,
                                       unsigned &countE,
                                       std::unordered_map<Point, unsigned> &PointToIndex, double t, std::vector<Edge> &wspdEdges) {

    std::vector<Edge> wspdpairs;
    linearithmicWSPD wspd(leaderPoints, 4 * (t + 1) / (t - 1), 2);
    std::unordered_map<Point ,unsigned> leaderPLocalIds(leaderPoints.size());
    std::vector<unsigned > leaderLocalToGlobalIds(leaderPoints.size());
    unsigned lIds = 0;
    for( Point &leaderPoint: leaderPoints){
        leaderPLocalIds[leaderPoint] = lIds;
        leaderLocalToGlobalIds[lIds] = PointToIndex[leaderPoint];
        lIds++;
    }

    wspd.computePairsAndStoreIn(wspdpairs, leaderPLocalIds);

    for( auto &p : wspdpairs ) {
        G[leaderLocalToGlobalIds[p.first]].insert(leaderLocalToGlobalIds[p.second]);
        G[leaderLocalToGlobalIds[p.second]].insert(leaderLocalToGlobalIds[p.first]);
        wspdEdges.emplace_back(leaderLocalToGlobalIds[p.first],leaderLocalToGlobalIds[p.second] );
        countE++;
    }
}

void constructWSPDSpanner(std::vector<Point> &P, std::vector<Edge> &E, double t) {

    std::vector<Edge> wspdpairs;
    linearithmicWSPD wspd(P, 4 * (t + 1) / (t - 1), 2);
    std::unordered_map<Point ,unsigned> leaderPLocalIds(P.size());
    std::vector<unsigned > leaderLocalToGlobalIds(P.size());
    std::unordered_map<Point,unsigned> PointToIndex;

    for (unsigned i = 0; i < P.size(); ++i)
        PointToIndex[P[i]] = i;

    unsigned lIds = 0;
    for( Point &leaderPoint: P){
        leaderPLocalIds[leaderPoint] = lIds;
        leaderLocalToGlobalIds[lIds] = PointToIndex[leaderPoint];
        lIds++;
    }

    wspd.computePairsAndStoreIn(wspdpairs, leaderPLocalIds);

    for( auto &p : wspdpairs )
        E.emplace_back(leaderLocalToGlobalIds[p.first], leaderLocalToGlobalIds[p.second] );
}

#endif //FASTSPARSESPANNER_WSPDSPANNER_H
