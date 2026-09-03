#include "SourceDistribution.hh"

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace
{

void ValidateSourceIDs(const std::vector<int>& sourceIDs)
{
    if(sourceIDs.empty())
        throw std::invalid_argument("No source organ ID was defined");

    std::unordered_set<int> uniqueIDs;
    for(std::vector<int>::const_iterator id = sourceIDs.begin(); id != sourceIDs.end(); ++id){
        if(!uniqueIDs.insert(*id).second)
            throw std::invalid_argument("Duplicate source organ ID: " + std::to_string(*id));
    }
}

void ValidateTet(const tetcal::SourceTet& tet)
{
    if(!std::isfinite(tet.volume) || tet.volume <= 0.)
        throw std::invalid_argument("A source tetrahedron has a non-positive or non-finite volume");
}

std::vector<tetcal::CDFEntry> MakeCDF(std::vector<tetcal::CDFEntry> weightedCopies)
{
    if(weightedCopies.empty())
        throw std::invalid_argument("No tetrahedra were found for the requested source organs");

    double totalWeight = 0.;
    for(std::vector<tetcal::CDFEntry>::const_iterator entry = weightedCopies.begin();
        entry != weightedCopies.end(); ++entry)
        totalWeight += entry->first;

    if(!std::isfinite(totalWeight) || totalWeight <= 0.)
        throw std::invalid_argument("The total source weight must be finite and greater than zero");

    double cumulativeWeight = 0.;
    for(std::vector<tetcal::CDFEntry>::iterator entry = weightedCopies.begin();
        entry != weightedCopies.end(); ++entry){
        cumulativeWeight += entry->first;
        entry->first = cumulativeWeight / totalWeight;
    }
    weightedCopies.back().first = 1.;
    return weightedCopies;
}

}

namespace tetcal
{

std::vector<CDFEntry> BuildVolumeWeightedCDF(
    const std::vector<SourceTet>& tetrahedra,
    const std::vector<int>& sourceIDs)
{
    ValidateSourceIDs(sourceIDs);

    std::unordered_set<int> requested(sourceIDs.begin(), sourceIDs.end());
    std::unordered_set<int> found;
    std::vector<CDFEntry> weightedCopies;
    weightedCopies.reserve(tetrahedra.size());
    for(std::vector<SourceTet>::const_iterator tet = tetrahedra.begin(); tet != tetrahedra.end(); ++tet){
        if(requested.find(tet->sourceID) == requested.end()) continue;
        ValidateTet(*tet);
        found.insert(tet->sourceID);
        weightedCopies.push_back(CDFEntry(tet->volume, tet->copyNumber));
    }

    for(std::vector<int>::const_iterator id = sourceIDs.begin(); id != sourceIDs.end(); ++id){
        if(found.find(*id) == found.end())
            throw std::invalid_argument("No tetrahedra were found for source organ ID " +
                                        std::to_string(*id));
    }
    return MakeCDF(weightedCopies);
}

std::vector<CDFEntry> BuildFractionWeightedCDF(
    const std::vector<SourceTet>& tetrahedra,
    const std::vector<int>& sourceIDs,
    const std::vector<double>& sourceWeights)
{
    ValidateSourceIDs(sourceIDs);
    if(sourceWeights.size() != sourceIDs.size())
        throw std::invalid_argument("The number of source fractions must match the number of source organ IDs");

    double totalSourceWeight = 0.;
    std::unordered_map<int, double> weights;
    std::unordered_map<int, double> totalVolumes;
    for(std::size_t i = 0; i < sourceIDs.size(); ++i){
        if(!std::isfinite(sourceWeights[i]) || sourceWeights[i] < 0.)
            throw std::invalid_argument("Source fractions must be finite and non-negative");
        weights[sourceIDs[i]] = sourceWeights[i];
        totalSourceWeight += sourceWeights[i];
    }
    if(!std::isfinite(totalSourceWeight) || totalSourceWeight <= 0.)
        throw std::invalid_argument("At least one source fraction must be greater than zero");

    std::unordered_set<int> requested(sourceIDs.begin(), sourceIDs.end());
    for(std::vector<SourceTet>::const_iterator tet = tetrahedra.begin(); tet != tetrahedra.end(); ++tet){
        if(requested.find(tet->sourceID) == requested.end()) continue;
        ValidateTet(*tet);
        totalVolumes[tet->sourceID] += tet->volume;
    }

    for(std::size_t i = 0; i < sourceIDs.size(); ++i){
        if(sourceWeights[i] > 0. && totalVolumes.find(sourceIDs[i]) == totalVolumes.end())
            throw std::invalid_argument("No tetrahedra were found for source organ ID " +
                                        std::to_string(sourceIDs[i]));
    }

    std::vector<CDFEntry> weightedCopies;
    weightedCopies.reserve(tetrahedra.size());
    for(std::vector<SourceTet>::const_iterator tet = tetrahedra.begin(); tet != tetrahedra.end(); ++tet){
        if(requested.find(tet->sourceID) == requested.end()) continue;
        const double sourceWeight = weights[tet->sourceID];
        if(sourceWeight == 0.) continue;
        const double normalizedSourceWeight = sourceWeight / totalSourceWeight;
        const double tetWeight = normalizedSourceWeight * tet->volume / totalVolumes[tet->sourceID];
        weightedCopies.push_back(CDFEntry(tetWeight, tet->copyNumber));
    }
    return MakeCDF(weightedCopies);
}

}
