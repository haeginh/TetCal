#include "SourceDistribution.hh"

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

void Require(bool condition, const std::string& message)
{
    if(!condition) throw std::runtime_error(message);
}

void RequireNear(double actual, double expected, const std::string& message)
{
    if(std::fabs(actual - expected) > 1.e-12)
        throw std::runtime_error(message + ": expected " + std::to_string(expected) +
                                 ", got " + std::to_string(actual));
}

std::map<int, double> Probabilities(const std::vector<tetcal::CDFEntry>& cdf)
{
    std::map<int, double> probabilities;
    double previous = 0.;
    for(std::vector<tetcal::CDFEntry>::const_iterator entry = cdf.begin();
        entry != cdf.end(); ++entry){
        probabilities[entry->second] = entry->first - previous;
        previous = entry->first;
    }
    return probabilities;
}

void RequireThrows(const std::function<void()>& operation, const std::string& message)
{
    try {
        operation();
    } catch(const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

}

int main()
{
    try {
        const std::vector<tetcal::SourceTet> tetrahedra = {
            {1., 10, 0}, {3., 10, 1}, {100., 20, 2}, {5., 30, 3}
        };
        const std::vector<int> sourceIDs = {10, 20};

        const std::map<int, double> volumeProbabilities = Probabilities(
            tetcal::BuildVolumeWeightedCDF(tetrahedra, sourceIDs));
        RequireNear(volumeProbabilities.at(0), 1. / 104., "volume weight for copy 0");
        RequireNear(volumeProbabilities.at(1), 3. / 104., "volume weight for copy 1");
        RequireNear(volumeProbabilities.at(2), 100. / 104., "volume weight for copy 2");
        Require(volumeProbabilities.find(3) == volumeProbabilities.end(),
                "an unrequested source was included");

        const std::map<int, double> fractionProbabilities = Probabilities(
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {1., 3.}));
        RequireNear(fractionProbabilities.at(0), 0.0625, "fraction weight for copy 0");
        RequireNear(fractionProbabilities.at(1), 0.1875, "fraction weight for copy 1");
        RequireNear(fractionProbabilities.at(2), 0.75, "fraction weight for copy 2");

        const std::map<int, double> scaledProbabilities = Probabilities(
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {10., 30.}));
        RequireNear(scaledProbabilities.at(0), fractionProbabilities.at(0),
                    "scaled fractions changed copy 0 probability");
        RequireNear(scaledProbabilities.at(1), fractionProbabilities.at(1),
                    "scaled fractions changed copy 1 probability");
        RequireNear(scaledProbabilities.at(2), fractionProbabilities.at(2),
                    "scaled fractions changed copy 2 probability");

        const std::map<int, double> zeroProbabilities = Probabilities(
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {0., 2.}));
        Require(zeroProbabilities.size() == 1 && zeroProbabilities.at(2) == 1.,
                "zero-weight source was not excluded");

        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {1.});
        }, "fraction-count mismatch was accepted");
        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {});
        }, "empty fraction list was accepted");
        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {1., -1.});
        }, "negative fraction was accepted");
        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, {10, 10}, {1., 1.});
        }, "duplicate source ID was accepted");
        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, {10, 99}, {1., 1.});
        }, "missing positive-weight source was accepted");
        const std::map<int, double> missingZeroProbabilities = Probabilities(
            tetcal::BuildFractionWeightedCDF(tetrahedra, {10, 99}, {1., 0.}));
        RequireNear(missingZeroProbabilities.at(0), 0.25,
                    "missing zero-weight source changed copy 0 probability");
        RequireNear(missingZeroProbabilities.at(1), 0.75,
                    "missing zero-weight source changed copy 1 probability");
        RequireThrows([&]() {
            tetcal::BuildFractionWeightedCDF(tetrahedra, sourceIDs, {0., 0.});
        }, "all-zero fractions were accepted");

        std::cout << "SourceDistribution tests passed" << std::endl;
        return EXIT_SUCCESS;
    } catch(const std::exception& error) {
        std::cerr << "SourceDistribution test failed: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
}
