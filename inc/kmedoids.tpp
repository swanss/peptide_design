/**
 A helper file defining implementations for the KmedoidsCluster class template.
 The implementations must be directly included in the header file (clustertree.h) in order
 for the compiler to synthesize class definitions, so they are provided here to avoid 
 implementing the clustering algorithm directly into the header.
*/
#include "ranklist.h"
#include <float.h>

template <typename Element>
void KMedoidsClusterer<Element>::computeMedoidDistances() {
    // Pairwise compute first- and second-closest selected objects
    closestDistances.clear();
    secondClosestDistances.clear();
    
    for (Element *el: allElements) {
        RankList<Element *> rankings(2, false);
        for (Element *med: selected) {
            rankings.insert(med, distanceFn(el, med));
        }
        vector<pair<Element *, double>> rankedItems = rankings.values();
        if (rankedItems.size() >= 1) {
            closestDistances[el] = rankedItems[0];
        } else {
            closestDistances[el] = make_pair(nullptr, DBL_MAX);
        }

        if (rankedItems.size() >= 2) {
            secondClosestDistances[el] = rankedItems[1];
        } else {
            secondClosestDistances[el] = make_pair(nullptr, DBL_MAX);
        }
    }
}

template <typename Element>
void KMedoidsClusterer<Element>::cluster() {
    unselected.insert(allElements.begin(), allElements.end());

    //selected.clear();
    //selected.insert(selected.end(), elements.begin(), numMedoids < elements.size() ? elements.begin() + numMedoids : elements.end());

    //int medoidIndex = 0;
    //for (int i = numMedoids; i < elements.size(); i++) {
    //    clusters[selected[medoidIndex % selected.size()]].push_back(elements[i]);
    //    medoidIndex++;
    //}

    // PAM algorithm
    buildInitialClusters();
    int numIterations = 0;
    while (performSwap()) {
        numIterations++;
        if (numIterations > 10)
            break;
    }
    // cout << "Converged in " << numIterations << " iterations" << endl;

    // Set the final return values
    centroids.insert(centroids.end(), selected.begin(), selected.end());
    for (auto it: closestDistances) {
        Element *child = it.first;
        pair<Element *, double> closest = it.second;
        clusters[closest.first].push_back(child);
    } 
}

template <typename Element>
void KMedoidsClusterer<Element>::buildInitialClusters() {
    // Initialize by finding the object whose distance to all others is minimal
    Element *bestElement = nullptr;
    double bestDistanceSum = DBL_MAX;
    for (Element *el: allElements) {
        double distanceSum = 0.0;
        for (Element *el2: allElements) {
            if (el != el2)
                distanceSum += distanceFn(el, el2);
        }
        if (distanceSum < bestDistanceSum) {
            bestElement = el;
            bestDistanceSum = distanceSum;
        }
    }   
    if (!bestElement) {
        return;
    }
    select(bestElement);

    while (selected.size() < min(numMedoids, (int)allElements.size())) {
        // Define gain as the amount by which clustering would benefit if
        // the object were selected; choose object to maximize gain
        Element *nextSelection = nullptr;
        double bestGain = -DBL_MAX;
        for (Element *el: unselected) {
            double gain = 0.0;
            for (Element *el2: unselected) {
                if (el == el2) continue;
                double dj = closestDistances[el2].second;
                gain += max(dj - distanceFn(el, el2), 0.0);
            }
            if (gain > bestGain) {
                nextSelection = el;
                bestGain = gain;
            }
        }
        select(nextSelection);
    }
}

template <typename Element>
bool KMedoidsClusterer<Element>::performSwap() {
    // Find an element in selected and an element in unselected, such that
    // swapping the two elements reduces the total dissimilarity between
    // each cluster element and its centroid.
    Element *bestSelected = nullptr;
    Element *bestUnselected = nullptr;
    double bestDelta = 0.0;

    for (Element *i: selected) {
        for (Element *h: unselected) {
            double delta = 0.0;
            for (Element *j: unselected) {
                if (j == h) continue;
                double dij = distanceFn(i, j);
                double dj = closestDistances[j].second;
                double djh = distanceFn(j, h);
                if (fabs(dij - dj) < 1e-4) {
                    // i is element j's centroid. Swapping i would move its
                    // centroid to its second closest option.
                    double ej = secondClosestDistances[j].second;
                    delta += min(djh, ej) - dj;
                } else {
                    // j isn't in i's cluster. It might move to the new cluster
                    // h, or not.
                    delta += min(djh - dj, 0.0);
                }
            }
            
            if (delta < bestDelta) {
                bestSelected = i;
                bestUnselected = h;
                bestDelta = delta;
            }
        }
    }

    if (!bestSelected || !bestUnselected)
        return false;
    
    swap(bestSelected, bestUnselected);
    return true;
}

