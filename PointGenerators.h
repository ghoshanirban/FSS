#ifndef FASTSPARSESPANNER_POINTGENERATORS_H
#define FASTSPARSESPANNER_POINTGENERATORS_H

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <unordered_set>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include "Utilities.h"


class RandomPointGenerator {

public:
    RandomPointGenerator() : m_randCgal(std::rand()) {}

    void uni_square(const index_t n, const double sizeOfSquare, std::vector<Point> &P) {
        typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<number_t, Point> > Point_generator;
        Point_generator g(sizeOfSquare / 2, m_randCgal);

        std::unordered_set<Point> P_unique;
        size_t remaining;
        while ((remaining = n - P_unique.size()) > 0)
            std::copy_n(g, remaining, inserter(P_unique));

        std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
    }    

    void normal_clustered(const unsigned n,
                          std::vector<Point> &P,
                          const number_t xStdDev = 2.0,
                          const number_t yStdDev = 2.0) {
        std::mt19937 rngX(seed());
        std::mt19937 rngY(seed());
        std::default_random_engine generatorX(rngX()), generatorY(rngY());
        std::normal_distribution<double> distributionX(0.0, xStdDev), distributionY(5.0, yStdDev);

        std::mt19937 rngShift(std::random_device{}());
        std::uniform_int_distribution<> shiftDistribution(0, INT32_MAX);

        index_t shiftX, shiftY;
        std::unordered_set<std::pair<index_t, index_t>, boost::hash<std::pair<index_t, index_t>>> S;

        std::unordered_set<Point> P_unique;
        unsigned totalPs = n;
        const unsigned numberOfClusters = (unsigned) std::sqrt(n);
        const unsigned pointsInACuster = (unsigned) std::sqrt(n);

        for (index_t c = 0; c < numberOfClusters; c++) {
            if (c != 0) {
                shiftX = shiftDistribution(rngShift) % (1000 * numberOfClusters);
                shiftY = shiftDistribution(rngShift) % (1000 * numberOfClusters);

                while (!(S.find(std::make_pair(shiftX, shiftY)) == S.end())) {
                    shiftX = shiftDistribution(rngShift) % (2000 * numberOfClusters); // set 740 for sink-spanner (prev 2000)
                    shiftY = shiftDistribution(rngShift) % (2000 * numberOfClusters); // set 740 for sink-spanner (prev 2000)
                }
            } else
                shiftX = shiftY = 0;

            S.insert(std::make_pair(shiftX, shiftY));

            for (index_t i = 0; i < pointsInACuster; i++) {
                double x = distributionX(generatorX) + (double) shiftX;
                double y = distributionY(generatorY) + (double) shiftY;
                P_unique.emplace(x, y);
            }
            if (c == numberOfClusters - 1) {
                const unsigned remainingPoints = totalPs - (numberOfClusters * pointsInACuster);
                for (index_t i = 0; i < remainingPoints; i++) {
                    double x = distributionX(generatorX) + (double) shiftX;
                    double y = distributionY(generatorY) + (double) shiftY;
                    P_unique.emplace(x, y);
                }
            }
        }

        std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
    }
    
    void grid_random(const index_t n, std::vector<Point> &P) {
        std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> S;
        std::unordered_set<Point> P_unique;

        std::mt19937 rngX(seed());
        std::mt19937 rngY(seed());
        std::uniform_int_distribution<> xDistribution(0, (int) ceil(0.7 * (double) n)), yDistribution(0, (int) ceil(
                0.7 * (double) n));

        index_t count = 0;

        while (count < n) {
            int x = xDistribution(rngX), y = yDistribution(rngY);

            if (S.find(std::make_pair(x, y)) == S.end()) {
                P_unique.emplace(x, y);
                S.insert(std::make_pair(x, y));
                count++;
            }
        }
        std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
    }
    void annulus(const index_t n, const double r2, const double r1, std::vector<Point> &P) {
        assert(r2 > r1);
        std::unordered_set<Point> P_unique;

        std::default_random_engine generator(seed());
        std::uniform_real_distribution<double> distributionR(r1, r2), distributionT(0, 1);

        for (index_t i = 0; i < n; i++) {
            double t = 2 * M_PI * distributionT(generator);
            double r = distributionR(generator);
            P_unique.emplace(r * cos(t), r * sin(t));
        }

        std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
    }    
        
    void galaxy(const index_t n, const unsigned numSpokes, std::vector<Point> &P) {
        // see https://itinerantgames.tumblr.com/post/78592276402/a-2d-procedural-galaxy-with-c
        //srand(seed());
        const double spokeAngle = 2 * M_PI / numSpokes,
                armOffsetMax = 0.5,
                rotationFactor = 5;

        std::unordered_set<Point> P_unique;

        while (P_unique.size() < n) {
            //for(index_t i=0; i<n; ++i) {
            double distance = randFloat();
            distance = pow(distance, 2);

            double angle = randFloat() * 2 * M_PI;
            double armOffset = randFloat() * armOffsetMax;
            armOffset -= armOffsetMax / 2;
            armOffset *= (1 / distance);

            double squaredArmOffset = pow(armOffset, 2);
            squaredArmOffset *= -1 * int(armOffset < 0);
            armOffset = squaredArmOffset;

            double rotation = distance * rotationFactor;

            angle = ((unsigned) (angle / spokeAngle)) * spokeAngle;
            angle += armOffset + rotation;

            P_unique.emplace(cos(angle) * distance * 500, sin(angle) * distance * 500);
        }

        std::copy(P_unique.begin(), P_unique.end(), std::back_inserter(P));
        //perturb(P, perturbationValue);
    }
    
    void convex(const index_t &n, const double &areaOfConvex, std::vector<Point> &P) {
        P.reserve(n);
        typedef CGAL::Random_points_in_square_2<Point, CGAL::Creator_uniform_2<double, Point> > Point_generator;
        Point_generator g(areaOfConvex, m_randCgal);
        CGAL::random_convex_set_2(n, std::back_inserter(P), g);
    }
    
    void spokes(const index_t n, const unsigned numSpokes, std::vector<Point> &P) {
        double spokeAngle = 2 * M_PI / numSpokes;

        std::unordered_set<Point> P_unique;

        while (P_unique.size() < n) {
            double distance = 500 * randFloat();
            double angle = randFloat() * 2 * M_PI;
            angle = ((unsigned) (angle / spokeAngle)) * spokeAngle;
            P_unique.emplace(cos(angle) * distance,
                             sin(angle) * distance);
        }

        std::copy(P_unique.begin(), P_unique.end(), back_inserter(P));
    }
private:
    CGAL::Random m_randCgal;

    size_t seed() {
        return std::rand();
    }

    double randFloat() {
        return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    }
};


#endif //FASTSPARSESPANNER_POINTGENERATORS_H
