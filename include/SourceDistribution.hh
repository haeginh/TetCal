#ifndef SOURCE_DISTRIBUTION_HH
#define SOURCE_DISTRIBUTION_HH

#include <utility>
#include <vector>

namespace tetcal
{

struct SourceTet
{
    double volume;
    int sourceID;
    int copyNumber;
};

typedef std::pair<double, int> CDFEntry;

// Selects all requested tetrahedra in proportion to their volume.
std::vector<CDFEntry> BuildVolumeWeightedCDF(
    const std::vector<SourceTet>& tetrahedra,
    const std::vector<int>& sourceIDs);

// Selects a source organ according to sourceWeights, then selects a
// tetrahedron inside that organ in proportion to its volume.
std::vector<CDFEntry> BuildFractionWeightedCDF(
    const std::vector<SourceTet>& tetrahedra,
    const std::vector<int>& sourceIDs,
    const std::vector<double>& sourceWeights);

}

#endif
