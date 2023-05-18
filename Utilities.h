#ifndef FASTSPARSESPANNER_UTILITIES_H
#define FASTSPARSESPANNER_UTILITIES_H

#include <iostream>
#include <unordered_set>
#include <vector>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::FT number_t;
typedef unsigned index_t;

typedef std::pair<index_t, index_t> index_tPair;
typedef index_tPair Edge;
typedef boost::hash<index_tPair> index_tPairHash;
const size_t SIZE_T_MAXX = std::numeric_limits<size_t>::max();

template<typename Point_2>
inline number_t L2Distance(const Point_2 &p, const Point_2 &q) {
    return p == q ? 0 : CGAL::sqrt(CGAL::squared_distance(p, q));
}

void readPointsFromFile(const std::string &inputFileName, std::vector<Point> &points, const std::size_t n = SIZE_T_MAXX) {

    std::ifstream in(inputFileName);
    std::unordered_set<Point> S;

    unsigned numberOfLines = 0;
    if (in.is_open()) {
        double x, y;
        size_t i = 0;
        while (i < n && in >> x >> y) {
            numberOfLines++;
            S.insert(Point(x, y));
            ++i;
        }
        in.close();
    }
    for (auto p: S)
        points.emplace_back(p);

    if (numberOfLines - S.size() > 0)
        std::cout << "Oops! Duplicates found in the input file!" << std::endl;

}

static inline std::string timeStamp() {
    std::string str;
    std::time_t result = std::time(nullptr);
    str = std::asctime(std::localtime(&result));
    str.pop_back();
    return "[ " + str + " ] ";
}


#endif // FASTSPARSESPANNER_UTILITIES_H
