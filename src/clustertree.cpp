#include "clustertree.h"
#include "Util.h"
#include "Fragments.h"
#include <unordered_set>
#include "findpaths.h"
#include <iostream>

using Node = ClusterNode<FragmentInfo>;

// ClusterSearchResults

Structure *ClusterSearchResults::getFullStructure(int i) {
    return _fetcher->getFullStructure(results[i].info);
}

Structure ClusterSearchResults::getResultStructure(int i) {
    return _fetcher->getResultStructure(results[i].info);
}

vector<Residue *> ClusterSearchResults::getResultResidues(int i) {
    AtomPointerVector apv = _fetcher->getAPV(results[i].info);
    vector<Residue *> result;
    for (int i = 0; i < apv.size(); i++) {
        Residue *res = apv[i]->getResidue();
        if (res != nullptr && (result.empty() || res != result.back()))
            result.push_back(res);
    }
    return result;
}

FragmentInfo *ClusterSearchResults::getFragmentInfo(int i) {
    return results[i].info->copy();
}

AtomPointerVector ClusterSearchResults::getAPV(int i) {
    return _fetcher->getAPV(results[i].info);
}

Sequence ClusterSearchResults::getResultSequence(int i) {
    Structure s = _fetcher->getResultStructure(results[i].info);
    const vector<Residue*> &residues = s.getResidues();
	vector<res_t> aas;
	for (int i = 0; i < residues.size(); i++) {
		aas.push_back(SeqTools::aaToIdx(residues[i]->getName()));
	}
	return Sequence(aas, "result_" + to_string(i) + ": " + results[i].info->toString());

}

double ClusterSearchResults::getRMSD(int i) {
    if (results[i].rmsd == DBL_MAX)
        results[i].rmsd = _rmsdCalc->bestRMSD(results[i].info, query);
    return results[i].rmsd;
}

// ClusterTree

ClusterTree::~ClusterTree() {
    if (_root != nullptr) {
        delete _root;    
    }
}

void ClusterTree::cluster() {
    _fragmentFetcher->reset();

    long numInserted = 0;
    while (_fragmentFetcher->hasNext()) {    
        pair<AtomPointerVector, FragmentInfo *> item = _fragmentFetcher->next();
        insert(item.first, item.second); 
        numInserted++; 
        double treeDepth = log(numInserted) / log(_refineLogBase);
        if (numInserted > 1 && (double)(int)treeDepth == treeDepth) {
            cout << "Inserted " << numInserted << endl;
            cout << "Average nest depth: " << (depthTotal / (float)depthCount) << endl;
            clusteringOps = 0;
            calculateRMSDProps({ _root });
            refineClustering(_root);
            cout << "Performed " << clusteringOps << " clustering ops" << endl;
            depthCount = 0;
            depthTotal = 0;
        }
        if (numInserted % 10000 == 0) 
            cout << "Inserted " << numInserted << " APVs" << endl;
        if (sizeLimit > 0 && numInserted >= sizeLimit) 
            break;
    }
    cout << "Final refinements" << endl;
    calculateRMSDProps({ _root });
    clusteringOps = 0;
    refineClustering(_root);
    calculateRMSDProps({ _root }, true);    
}

void ClusterTree::insert(const AtomPointerVector &item, FragmentInfo *idInfo) {
    Node *newNode = new Node(idInfo, _childLimit); 
    if (!_root) {
        _root = newNode;
        return;
    }
    // Find the top k best candidate clusters
    vector<Node *> candidates = { _root };

    RankList<Node *> bestParents(_numInsertOptions, false);
    
    // One iteration per level of the tree
    while (!candidates.empty()) {
        RankList<Node *> newCandidateRanking(_numInsertOptions, false);
        for (auto node: candidates) {
            // Add the node's best children to the running for the next tree level
            for (auto child: node->getChildren()) {
                double rmsd = rmsdCalc.bestRMSD(child->getItem(), item);
                newCandidateRanking.insert(child, rmsd);
            }

            // If all of the node's children have children, then don't add to this
            // node directly. Otherwise it can be a candidate for item's parent
            vector<Node *> children = node->getChildren();
            if (children.size() >= _childLimit) {
                continue;
            }
            bool allWithChildren = all_of(children.begin(), children.end(), [](Node *item) { return item->hasChildren(); });
            if (children.empty() || !allWithChildren) {
                bestParents.insert(node, rmsdCalc.bestRMSD(node->getItem(), item));
            }
        }
        vector<Node *> newCandidates = newCandidateRanking.items();
        candidates.clear();
        candidates.insert(candidates.end(), newCandidates.begin(), newCandidates.end());
        depthTotal++;
    }
    depthCount++;

    // Take the best parent and add this node as a child
    auto bestParentItem = bestParents.values()[0];
    bestParentItem.first->addChild(newNode);
}

void ClusterTree::addClusteringCandidates(Node *parent, vector<Node *> &elements) {
    // We want to assemble a list of candidates for the new children of parent.
    // The goal is to get as close to _childLimit + _childLimit ^ 2 candidates
    // as possible, since this will allow for a roughly even breadth/depth
    // balance in the tree. If the tree is currently imbalanced, we may need to
    // go more than 2 levels deep to get enough nodes.
    int maxSize = _childLimit + _childLimit * _childLimit;
    if (elements.size() >= maxSize)
        return;

    vector<Node *> children = parent->getChildren();
    elements.insert(elements.end(), children.begin(), children.end());
    int levelsTraversed = 1;
    while (elements.size() < maxSize && !children.empty()) {
        // Since we have space, go through the children in order of decreasing
        // subtree radius and add their children.
        sort(children.begin(), children.end(), [](const Node *c1, const Node *c2) { return c1->subtreeRadius > c2->subtreeRadius; });

        vector<Node *> newChildren;
        for (Node *child: children) {
            vector<Node *> grandchildren = child->getChildren();
            elements.insert(elements.end(), grandchildren.begin(), grandchildren.end());
            newChildren.insert(newChildren.end(), grandchildren.begin(), grandchildren.end());
            if (elements.size() >= maxSize) break;
        }
        // Now the children are replaced with grandchildren, and we can go 
        // through another generation if there is still space.
        children = newChildren;
        levelsTraversed++;
    }
    // cout << "Traversed " << levelsTraversed << " to get " << elements.size() << " clustering candidates" << endl;
}

void ClusterTree::refineClustering(Node *parent, bool bottomUp, int levelLimit) {
    if (levelLimit == 0)
        return;

    clusteringOps++;
    if (clusteringOps % 1000 == 0)
        cout << clusteringOps << " clustering ops so far" << endl;

    if (bottomUp) {
        // Recurse downward first
        for (Node *child: parent->getChildren()) {
            refineClustering(child, bottomUp, levelLimit - 1);
        }
    } else {
        parent = reassignParent(parent);
    }

    //cout << "Clustering from root " << parent->getItem()->toString() << ", " << bottomUp << endl;
    vector<Node *> elements;
    addClusteringCandidates(parent, elements);
    // cout << "Elements to cluster: ";
    // for (auto el: elements) {
    //     cout << el->getItem()->toString() << ", ";
    // }
    // cout << endl;
    if (elements.empty()) return;
    
    // Perform the clustering
    KMedoidsClusterer<Node> clusterer(elements, _childLimit, [this](Node *c1, Node *c2) {
        return rmsdCalc.bestRMSD(c1->getItem(), c2->getItem());
    });

    // Reorganize the parent's children to reflect the new clusters
    for (Node *centroid: clusterer.getCentroids()) {
        centroid->moveToParent(parent);
        for (Node *element: clusterer.getCluster(centroid)) {
            if (element == centroid) continue;
            element->moveToParent(centroid);
        }
    }

    if (!bottomUp) {
        for (Node *centroid: clusterer.getCentroids()) {
            refineClustering(centroid, bottomUp, levelLimit - 1);
        }
    }
}

Node * ClusterTree::reassignParent(Node *parent) {
    // Assemble a list of candidates
    vector<Node *> elements; 
    elements.push_back(parent);
    addClusteringCandidates(parent, elements);
    // Find the best parent
    return assignBestParent(parent, elements);
}

Node * ClusterTree::assignBestParent(Node *parent, const vector<Node *> &elements) {
    // Find the element that minimizes the max distance to all of the elements
    Node *newParent = nullptr;
    double bestMaxDistance = DBL_MAX;

    for (Node *el: elements) {
        double maxDistance = 0.0;
        for (Node *el2: elements) {
            double dist = rmsdCalc.bestRMSD(el->getItem(), el2->getItem());
            if (dist > maxDistance) {
                maxDistance = dist;
            }
        }
        if (maxDistance < bestMaxDistance) {
            newParent = el;
            bestMaxDistance = maxDistance;
        }
    }

    if (newParent != nullptr && newParent != parent) {
        // We need to replace parent with newParent (just leave the old
        // parent as a child of newParent for now)
        if (parent == _root) {
            // Special case for root node
            _root = newParent;
            newParent->moveToParent(nullptr);
        } else {
            newParent->moveToParent(parent->getParent());
        }
        parent->moveToParent(newParent);
        // Also move all of the old parent's children to the new parent
        for (auto child: parent->getChildren())
            child->moveToParent(newParent);
    }
    return newParent;
}

void ClusterTree::calculateRMSDProps(vector<Node *> nodes, bool fullAncestry) {
    if (nodes.size() == 1) {
        // Root node
        nodes[0]->parentRMSD = 0.0;
    }
    
    Node *subtreeRoot = nodes.back();
    AtomPointerVector parentAPV = _fragmentFetcher->getAPV(subtreeRoot->getItem());
    if (fullAncestry && nodes.size() > 1) {
        // Calculate RMSD with respect to all ancestors
        if (subtreeRoot->ancestorRMSDs != nullptr)
            delete[] subtreeRoot->ancestorRMSDs;
        subtreeRoot->ancestorRMSDs = new double[nodes.size() - 1];
        int ancestorIndex = 0;
        for (int i = nodes.size() - 2; i >= 0; i--) {
            subtreeRoot->ancestorRMSDs[ancestorIndex++] = rmsdCalc.bestRMSD(nodes[i]->getItem(), subtreeRoot->getItem());
        }
    }

    // Recurse through the hierarchy
    for (Node *child: subtreeRoot->getChildren()) {
        // Set parentRMSD for children
        child->parentRMSD = rmsdCalc.bestRMSD(child->getItem(), parentAPV);

        // Reset subtree radius, since it will be determined by future calls
        child->subtreeRadius = 0.0;

        // Traverse the path to the subtree root and update subtreeRadius
        for (Node *node: nodes) {
            node->subtreeRadius = max(node->subtreeRadius, rmsdCalc.bestRMSD(node->getItem(), child->getItem()));
        }

        vector<Node *> newPath = nodes;
        newPath.push_back(child);
        calculateRMSDProps(newPath, fullAncestry);
    }
}

void ClusterTree::sortSubtreesByRadius(Node *root) {
    Node *parent = root ? root : _root;

    vector<Node *> children = parent->getChildren();
    sort(children.begin(), children.end(), [](const Node *c1, const Node *c2) { return c1->subtreeRadius > c2->subtreeRadius; });
    parent->setChildren(children);

    // Sort subtrees of children
    for (Node *child: children) {
        sortSubtreesByRadius(child);
    }
}

// Merge constructor
ClusterTree::ClusterTree(const vector<ClusterTree> &trees, int mergeLimit): rmsdCalc(trees[0]._fragmentFetcher) {
    // Make sure all trees have the same properties
    FragmentFetcher *fetcher = nullptr;
    int childLim = -1;
    for (const ClusterTree &tree: trees) {
        if (!fetcher) fetcher = tree._fragmentFetcher;
        if (childLim < 0) childLim = tree._childLimit;
        MstUtils::assert(tree._fragmentFetcher == fetcher, "All trees to be merged must have the same fragment fetcher reference");
        MstUtils::assert(tree._childLimit == childLim, "All trees to be merged must have the same child limit");
    }
    _fragmentFetcher = fetcher;
    _childLimit = childLim;
    
    // Make a dummy root node for now
    Node *dummyRoot = new Node(nullptr, _childLimit);
    _root = dummyRoot;

    // Find the best root node
    vector<Node *> rootCandidates;
    for (const ClusterTree &tree: trees) {
        Node *treeRoot = copySubtree(tree._root);
        _root->addChild(treeRoot);
    }
    addClusteringCandidates(_root, rootCandidates);    
    // This changes the value of _root
    Node *newRoot = assignBestParent(_root, rootCandidates);
    dummyRoot->moveToParent(nullptr);

    // Now we have a 'dirty' tree rooted at newRoot - clean it up
    calculateRMSDProps({ _root });
    refineClustering(_root, false, mergeLimit);
    if (mergeLimit >= 0) {
        vector<Node *> parents;
        vector<Node *> children;
        getNodesAtLevel(mergeLimit, parents);
        getNodesAtLevel(mergeLimit + 1, children);
        greedyCluster(parents, children);
    }
    calculateRMSDProps({ _root }, true);    
}

Node *ClusterTree::copySubtree(Node *parent) {
    Node *newParent = new Node(parent->getItem()->copy(), _childLimit);
    newParent->parentRMSD = parent->parentRMSD;
    newParent->subtreeRadius = parent->subtreeRadius;
    for (auto child: parent->getChildren()) {
        newParent->addChild(copySubtree(child));
    }
    return newParent;
}

void ClusterTree::getNodesAtLevel(int level, vector<Node *> &results, Node *parent) const {
    if (parent == nullptr)
        parent = _root;

    if (level == 0) {
        results.push_back(parent);
        return;
    }

    for (Node *child: parent->getChildren())
        getNodesAtLevel(level - 1, results, child);
}

vector<Node *> ClusterTree::getNodesAtLevel(int level) const {
    vector<Node *> results;
    getNodesAtLevel(level, results);
    return results;
}

void ClusterTree::greedyCluster(vector<Node *> parents, vector<Node *> children) {
    cout << "Clustering " << children.size() << " nodes into " << parents.size() << " parents" << endl;
    bool restrictChildCounts = children.size() <= parents.size() * _childLimit;
    if (restrictChildCounts)
        cout << "Restricting child counts to " << _childLimit << " per parent" << endl;
    else
        cout << "Too many children to cluster, so removing child count restriction" << endl;

    sort(children.begin(), children.end(), [](const Node *c1, const Node *c2) { return c1->subtreeRadius > c2->subtreeRadius; });
    for (Node *child: children) {
        child->moveToParent(nullptr);
    }

    for (int i = 0; i < children.size(); i++) {
        if (i % 1000 == 0)
            cout << "Assigning node " << i << endl;
        Node *child = children[i];
        // Assign the child to the parent to which it is closest that has
        // open slots
        double bestRMSD = DBL_MAX;
        Node *bestParent = nullptr;
        for (Node *parent: parents) {
            if (restrictChildCounts && parent->childSize() >= _childLimit)
                continue;

            double rmsd = rmsdCalc.bestRMSD(child->getItem(), parent->getItem());
            if (rmsd < bestRMSD) {
                bestRMSD = rmsd;
                bestParent = parent;
            }
        }

        MstUtils::assert(bestParent != nullptr, "No best parent found!");
        child->moveToParent(bestParent);
    }
}

string ClusterTree::toString() const {
    stringstream ss;
    ss << std::fixed << std::setprecision(3);
    vector<pair<Node *, int>> printStack = { make_pair(_root, 0) };
    unordered_set<Node *> visited;
    while (!printStack.empty()) {
        auto nextItem = printStack.back();
        Node *item = nextItem.first;
        printStack.pop_back();
        
        if (!visited.count(item)) {
            for (int i = 0; i < nextItem.second; i++) ss << " ";
            ss << item->getItem()->toString() << "," << item->parentRMSD << "," << item->subtreeRadius;
            if (item->ancestorRMSDs != nullptr) {
                for (int i = 1; i < nextItem.second; i++)
                    ss << "," << item->ancestorRMSDs[i];
            }
            ss << endl; 
            visited.insert(item);
        }

        for (auto child: item->getChildren()) {
            printStack.push_back(make_pair(child, nextItem.second + 1));
        }
    }
    return ss.str();
}

void ClusterTree::read(string clusterFilePath) {
    ifstream readstream(clusterFilePath);
    if (!readstream.is_open()) {
        cerr << "couldn't open in stream" << endl;
        return;
    }
    
    string line;
    Node *currentParent = nullptr;
    int previousNumSpaces = 0;
    while (getline(readstream, line)) {
        auto textStart = find_if(line.begin(), line.end(), [](int c) { return !isspace(c); });
        if (textStart == line.end()) continue;
        int numSpaces = distance(line.begin(), textStart);
        line.erase(line.begin(), textStart);
        vector<string> comps = splitString(line, ",");
        MstUtils::assert(comps.size() >= 3, "Need 3 comma-separated components, got " + to_string(comps.size()));

        FragmentInfo *info = _fragmentFetcher->makeFragmentInfo(comps[0]);
        Node *node = new Node(info, _childLimit);
        node->parentRMSD = stod(comps[1]);
        node->subtreeRadius = stod(comps[2]);
        if (currentParent == nullptr) {
            // Root node
            currentParent = node;
            _root = node;
        } else if (numSpaces > previousNumSpaces) {
            currentParent->addChild(node);
            currentParent = node;
        } else {
            while (numSpaces < previousNumSpaces) {
                currentParent = currentParent->getParent();
                previousNumSpaces--;
            }
            currentParent->getParent()->addChild(node);
            currentParent = node;
        }

        // Read ancestor RMSDs
        if (comps.size() > 3) {
            node->ancestorRMSDs = new double[int(comps.size()) - 2];
            int ancestorIndex = 0;
            node->ancestorRMSDs[ancestorIndex++] = node->parentRMSD;
            for (int i = 3; i < comps.size(); i++) {
                node->ancestorRMSDs[ancestorIndex++] = stod(comps[i]);
            }
        }
        previousNumSpaces = numSpaces;
    } 
}

void ClusterTree::write(string clusterFilePath) const {
    ofstream outFile(clusterFilePath, ofstream::out);
    outFile << toString() << endl;
    outFile.close();
}

ClusterSearchResults ClusterTree::search(const AtomPointerVector &query, double rmsdCutoff) {
    numNodesTotal = 0;
    numNodesSearched = 0;
    numNodesEliminated = 0;
    numNodesAutoAdded = 0;
    numNodesEliminatedByAncestor = 0;
    numNodesEliminatedByParent = 0;

    vector<_SearchResult> matches;
    rmsdCalc.setQuery(&query);
    _search(query, rmsdCutoff, _root, matches);
    ClusterSearchResults results(_fragmentFetcher, &rmsdCalc, query, matches);
    rmsdCalc.setQuery(nullptr);

    if (verbose) {
        cout << numNodesTotal << " nodes searched, " << numNodesSearched << " touched. " << numNodesEliminated << " eliminated, " << numNodesAutoAdded << " auto added." << endl;
        cout << numNodesEliminatedByAncestor << " eliminated by ancestor, vs " << numNodesEliminatedByParent << " by parent." << endl;
        cout << matches.size() << " results" << endl;
    }
    return results;
}

void ClusterTree::addSubtreeToResults(ClusterNode<FragmentInfo> *parent, vector<_SearchResult> &results, double parentRMSD) {
    // Add this node and all of its children
    results.emplace_back(parent->getItem(), parentRMSD);

    for (auto child: parent->getChildren())
       addSubtreeToResults(child, results); 
}

void ClusterTree::_search(const AtomPointerVector &query, double rmsdCutoff, Node *parent, vector<_SearchResult> &results) {
    // Check the RMSD with the parent
    if (sufficientMatches > 0 && results.size() >= sufficientMatches)
        return;
    numNodesSearched++;
    double parentRMSD = rmsdCalc.bestQueryRMSD(parent->getItem()); // , query);

    if (parentRMSD - parent->subtreeRadius > rmsdCutoff) {
        // The entire cluster is too far away to contain any matches!
        numNodesEliminated += parent->subtreeSize() - 1;
        numNodesTotal += parent->subtreeSize();
        return;
    } else if (parentRMSD + parent->subtreeRadius <= rmsdCutoff) {
        // The entire cluster is within the rmsdCutoff!
        addSubtreeToResults(parent, results, parentRMSD);
        numNodesAutoAdded += parent->subtreeSize();
        numNodesTotal += parent->subtreeSize();
    } else {
        // Some of the cluster children might be in the results.
        if (parentRMSD <= rmsdCutoff)
            results.emplace_back(parent->getItem(), parentRMSD);

        numNodesTotal++;
        // Enumerate the children and recurse
        for (auto child: parent->getChildren()) {
            if (child->ancestorRMSDs != nullptr) {
                double maxRMSD = 0.0;
                Node *currentParent = parent;
                int ancestorIndex = 0;
                while (currentParent != nullptr) {
                    double ancestorChildRMSD = child->ancestorRMSDs[ancestorIndex++]; // rmsdCalc.bestRMSD(child->getItem(), currentParent->getItem());
                    double ancestorQueryRMSD = rmsdCalc.bestQueryRMSD(currentParent->getItem());
                    maxRMSD = max(maxRMSD, ancestorQueryRMSD - ancestorChildRMSD);
                    currentParent = currentParent->getParent();
                    if (maxRMSD - child->subtreeRadius > rmsdCutoff)
                        break;
                    if (ancestorIndex == 4)
                        break;
                }
                if (maxRMSD - child->subtreeRadius > rmsdCutoff) {
                //if (parentRMSD - child->parentRMSD - child->subtreeRadius > rmsdCutoff) {
                    // The child and its children are all too far away
                    if (parentRMSD - child->parentRMSD - child->subtreeRadius > rmsdCutoff)
                        numNodesEliminatedByParent += child->subtreeSize();
                    else
                        numNodesEliminatedByAncestor += child->subtreeSize();
                    numNodesEliminated += child->subtreeSize();
                    numNodesTotal += child->subtreeSize();
                    continue;
                }
            }
            _search(query, rmsdCutoff, child, results);
        }
    }
}

// ClusterIterator

bool ClusterIterator::hasNext() {
    if (!_nextResult)
        makeNextResult();
    return _nextResult != nullptr;
}

Node * ClusterIterator::next() {
    if (!_nextResult)
        makeNextResult();
    Node *res = _nextResult;
    destroyNextResult();
    return res;
}

void ClusterIterator::makeNextResult() {
    if (_nextResult != nullptr)
        return;

    if (_stack.empty()) {
        // Nothing left to do
        return;
    }

    Node *topNode = _stack.back();
    _nextResult = topNode;

    if (topNode->hasChildren()) {
        _stack.push_back(topNode->getChildren()[0]);
        _indexes.push_back(0);
    } else {
        int index = _indexes.back();
        _stack.pop_back();
        _indexes.pop_back();
        while (!_stack.empty()) {
            Node *parent = _stack.back();
            if (index < parent->childSize() - 1) {
                // Add the next child of parent
                _stack.push_back(parent->getChildren()[index + 1]);
                _indexes.push_back(index + 1);
                break;
            }
            index = _indexes.back();
            _stack.pop_back();
            _indexes.pop_back();
        }
    }
}

void ClusterIterator::destroyNextResult() {
    _nextResult = nullptr;
}

// ClusterPairIterator

ClusterPairIterator::~ClusterPairIterator() {
    if (_nextResult != nullptr) {
        delete _nextResult;
        _nextResult = nullptr;
    }
}

bool ClusterPairIterator::hasNext() {
    if (!_nextResult)
        makeNextResult();
    return _nextResult != nullptr;
}

pair<Node *, Node *> ClusterPairIterator::next() {
    if (!_nextResult)
        makeNextResult();
    pair<Node *, Node *> res = *_nextResult;
    destroyNextResult();
    return res;
}

bool ClusterPairIterator::updateStacks(vector<ClusterNode<FragmentInfo> *> &nodeStack, vector<int> &indexStack, Node *markerNode, bool skipSubtree, bool *pushed) {
    Node *topNode = nodeStack.back();
    bool visitedMarker = false;

    // If this node has children, we will go to them next;
    // otherwise, if it has siblings, we will visit them;
    // and if not, we will go up the tree and find cousins to visit.
    // We don't need to directly visit ancestors.
    if (topNode->hasChildren() && !skipSubtree) {
        nodeStack.push_back(topNode->getChildren()[0]);
        indexStack.push_back(0);
        if (pushed != nullptr)
            *pushed = true;
    } else {
        int index = indexStack.back();
        if (topNode == markerNode)
            visitedMarker = true;
        nodeStack.pop_back();
        indexStack.pop_back();
        while (!nodeStack.empty()) {
            Node *parent = nodeStack.back();
            if (index < parent->childSize() - 1) {
                // Add the next child of parent
                nodeStack.push_back(parent->getChildren()[index + 1]);
                indexStack.push_back(index + 1);
                break;
            }
            index = indexStack.back();
            if (nodeStack.back() == markerNode)
                visitedMarker = true;
            nodeStack.pop_back();
            indexStack.pop_back();
        }
        if (pushed != nullptr)
            *pushed = false;
    }
    return visitedMarker;
}

void ClusterPairIterator::makeNextResult() {
    if (_nextResult != nullptr)
        return;

    _currentIsAncestor = _nextIsAncestor;

    Node *nextFirst = nullptr;
    if (!_firstStack.empty()) {
        nextFirst = _firstStack.back();
    } else {
        // There is nothing left to do
        return;
    }
    
    Node *nextSecond = nullptr;
    if (!_secondStack.empty()) {
        nextSecond = _secondStack.back();

        // We've established what the next item will be;
        // now update the stack for the next iteration.
        bool crossedMarker = updateStacks(_secondStack, _secondIndexes, nextFirst, _secondSubtreeSkipped, &_secondSubtreePushed);
        _secondSubtreeSkipped = false;
        _nextIsAncestor = _currentIsAncestor && !crossedMarker;
    }

    // Update first stacks if second is not available
    if (nextSecond == nullptr) {
        updateStacks(_firstStack, _firstIndexes);
        
        // The second stack should start at the point of the first stack
        _secondStack.insert(_secondStack.end(), _firstStack.begin(), _firstStack.end());
        _secondIndexes.insert(_secondIndexes.end(), _firstIndexes.begin(), _firstIndexes.end());
        if (_firstStack.empty()) {
            // Nothing left to do
            return;
        } else {
            nextFirst = _firstStack.back();
            nextSecond = _secondStack.back();
            _currentIsAncestor = true;
            // Need to update second stack again
            _nextIsAncestor = !updateStacks(_secondStack, _secondIndexes, nextFirst, false, &_secondSubtreePushed);
            _secondSubtreeSkipped = false;
        }
    }

    // At this point both first and second values should be present
    MstUtils::assert(nextFirst != nullptr, "Missing nextFirst");
    MstUtils::assert(nextSecond != nullptr, "Missing nextSecond");

    _nextResult = new pair<Node *, Node *>(nextFirst, nextSecond);
}

void ClusterPairIterator::destroyNextResult() {
    if (_nextResult == nullptr)
        return;
    delete _nextResult;
    _nextResult = nullptr;
}

void ClusterPairIterator::skipSecondSubtree() {
    if (_secondStack.empty())
        return;
    if (_secondSubtreePushed) {
        _secondStack.pop_back();
        _secondIndexes.pop_back();
        if (!_firstStack.empty())
            _nextIsAncestor = !updateStacks(_secondStack, _secondIndexes, _firstStack.back(), true, &_secondSubtreePushed);
    }
}
