#include <stdio.h>
#include "overlaps.h"
#include "fileutilities.h"
#include "utilities.h"

using namespace std;
using namespace MST;

// ========= CAResidueHasher ===========

template <typename hash_type, typename distance_metric>
CAResidueHasher<hash_type, distance_metric>::CAResidueHasher(vector<mstreal> bbox, mstreal increment): _bbox(bbox), _increment(increment) {
    _xBinSize = (hash_type)ceil((bbox[1] - bbox[0]) / increment);
    _yBinSize = (hash_type)ceil((bbox[3] - bbox[2]) / increment);
    _zBinSize = (hash_type)ceil((bbox[5] - bbox[4]) / increment);
}

template <typename hash_type, typename distance_metric>
hash_type CAResidueHasher<hash_type, distance_metric>::hash(Residue *res) const {
    Atom *ca = res->findAtom("CA");
    if (!ca)
        return 0;
    
    // Find the bin position in terms of x, y, z
    hash_type xBin = (hash_type)floor((ca->getX() - _bbox[0]) / _increment);
    hash_type yBin = (hash_type)floor((ca->getY() - _bbox[2]) / _increment);
    hash_type zBin = (hash_type)floor((ca->getZ() - _bbox[4]) / _increment);
    return xBin * _yBinSize * _zBinSize + yBin * _zBinSize + zBin;
}

template <typename hash_type, typename distance_metric>
vector<hash_type> CAResidueHasher<hash_type, distance_metric>::region(Residue *res, distance_metric cutoff) const {
    vector<hash_type> result;
    
    // Get CA position
    Atom *ca = res->findAtom("CA");
    if (!ca) {
        result.push_back(0);
        return result;
    }
    mstreal cax = ca->getX();
    mstreal cay = ca->getY();
    mstreal caz = ca->getZ();
    
    // Compute bin index bounds in each direction
    hash_type minX = (hash_type) max(0.0,
                                     floor((cax - cutoff - _bbox[0]) / _increment));
    hash_type maxX = (hash_type) min((double)_xBinSize,
                                     ceil((cax + cutoff - _bbox[0]) / _increment) + 1);
    hash_type minY = (hash_type) max(0.0,
                                     floor((cay - cutoff - _bbox[2]) / _increment));
    hash_type maxY = (hash_type) min((double)_yBinSize,
                                     ceil((cay + cutoff - _bbox[2]) / _increment) + 1);
    hash_type minZ = (hash_type) max(0.0,
                                     floor((caz - cutoff - _bbox[4]) / _increment));
    hash_type maxZ = (hash_type) min((double)_zBinSize,
                                     ceil((caz + cutoff - _bbox[4]) / _increment) + 1);
    
    // Add all hash values
    for (hash_type x = minX; x < maxX; x++) {
        for (hash_type y = minY; y < maxY; y++) {
            for (hash_type z = minZ; z < maxZ; z++) {
                result.push_back(x * _yBinSize * _zBinSize + y * _zBinSize + z);
            }
        }
    }
    
    return result;
}

// This only allows CAResidueHasher to be instantiated with its default template
// parameters. Otherwise, the implementations will not be accessible.
template class CAResidueHasher<unsigned long, mstreal>;

// ========= MaxDeviationVerifier ===========


bool MaxDeviationVerifier::verify(const vector<Residue *> &segment1, const vector<Residue *> &segment2) const {
    MstUtils::assert(segment1.size() == segment2.size(), "Segments must be same size");
    
    // Check all pairs of alpha carbons
    mstreal maxDistance = 0.0;
    for (int i = 0; i < segment1.size(); i++) {
        Atom *ca1 = segment1[i]->findAtom("CA");
        Atom *ca2 = segment2[i]->findAtom("CA");
        maxDistance = max(maxDistance, ca1->distance(*ca2));
    }
    
    return maxDistance <= _cutoff;
}
