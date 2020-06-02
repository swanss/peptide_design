#ifndef ClusterTree_H
#define ClusterTree_H

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <numeric>

#include "msttypes.h"
#include "mstsequence.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "seedutils.h"
#include "clusterutils.h"

using namespace MST;

/**
 * ClusterNode represents a node in a ClusterTree. ClusterNodes maintain 
 * references to their children, the RMSD with respect to their parent, and 
 * the RMSD radius of their cluster (i.e. the maximum RMSD of any child with
 * respect to this node). ClusterNodes are a template class and can be
 * defined to store any type of structural information.
*/
template <typename T>
class ClusterNode {

public:
    ClusterNode(T *item, int childLimit = 2): item(item), childLimit(childLimit) {};
    
    ~ClusterNode() {
        if (item != nullptr) {
            delete item;
        }
        if (ancestorRMSDs != nullptr)
            delete[] ancestorRMSDs;
        for (ClusterNode<T> *child: children) {
            delete child;
        }
    };

    T *getItem() const { return item; };
    ClusterNode<T> *getParent() const { return parent; }
    vector<ClusterNode<T> *> getChildren() const { return children; };
    size_t childSize() const { return children.size(); }
    bool hasChildren() const { return children.size() > 0; } 

    void addChild(ClusterNode<T> *child) {
        children.push_back(child);
        child->parent = this;
    }
    
    void removeChild(ClusterNode<T> *child) {
        auto it = find(children.begin(), children.end(), child);
        if (it != children.end())
            children.erase(it);
    } 

    void moveToParent(ClusterNode<T> *newParent) {
        if (parent != nullptr)
            parent->removeChild(this);
        if (newParent != nullptr)
            newParent->addChild(this);
    }

    void setChildren(vector<ClusterNode<T> *> &newChildren) {
        children.clear();
        children.insert(children.end(), newChildren.begin(), newChildren.end());
    }
    
    double parentRMSD = 0.0;
    double subtreeRadius = 0.0;
    double *ancestorRMSDs = nullptr;

    /**
     * @return the number of nodes in the subtree rooted by this node
     */
    int subtreeSize() const { return 1 + accumulate(children.begin(), children.end(), 0, [](int val, ClusterNode<T> *child) { return val + child->subtreeSize(); }); }

private:
    vector<ClusterNode<T> *> children;
    ClusterNode<T> *parent = nullptr;
    T *item = nullptr;
    int childLimit = 2; 
};

/**
 * A helper class that performs k-medoids clustering on an arbitrary data
 * type Element. The types supported for Element are declared in the
 * implementation file and currently support ClusterNode<FragmentInfo>.
 */
template <typename Element>
class KMedoidsClusterer {
public:
    KMedoidsClusterer(vector<Element *> elements, int numMedoids, std::function<double(Element *, Element *)> distanceFn): numMedoids(numMedoids), distanceFn(distanceFn), allElements(elements) {
        cluster();
    };

    // Retrieving clusters
    vector<Element *> getCentroids() { return centroids; }
    vector<Element *> getCluster(Element *centroid) { return clusters[centroid]; } 
private:
    int numMedoids;
    std::function<double(Element *, Element *)> distanceFn;

    vector<Element *> allElements;
    unordered_set<Element *> selected;
    unordered_set<Element *> unselected;
    unordered_map<Element *, pair<Element *, double>> closestDistances;
    unordered_map<Element *, pair<Element *, double>> secondClosestDistances;
    
    // Performing clustering
    void cluster();
    
    // Used to store final results
    vector<Element *> centroids;
    unordered_map<Element *, vector<Element *>> clusters;

    // Helpers for clustering
    void computeMedoidDistances();
    void select(Element *el) {
        unselected.erase(el);
        selected.insert(el);
        computeMedoidDistances();
    }
    void unselect(Element *el) {
        selected.erase(el);
        unselected.insert(el);
        computeMedoidDistances();
    }
    void swap(Element *sel, Element *unsel) {
        selected.erase(sel);
        selected.insert(unsel);
        unselected.erase(unsel);
        unselected.insert(sel);
        computeMedoidDistances();
    }

    // Cluster steps
    void buildInitialClusters();
    bool performSwap();
};

#include "kmedoids.tpp"

struct _SearchResult {
    _SearchResult(FragmentInfo *info, double rmsd = DBL_MAX): info(info), rmsd(rmsd) {};

    FragmentInfo *info;
    double rmsd;
};

class ClusterTree;

/**
 * Stores search results from a query on a ClusterTree. The search results
 * are stored in a memory-efficient way so that only the desired properties
 * will be loaded.
 */
class ClusterSearchResults {
    friend class ClusterTree;
public:
    /**
     * @return the number of search results
     */
    long size() { return (long)results.size(); };
    /**
     * @param i the index of the search result to retrieve
     * @return the full structure of the protein from which search result i
     *  originated
     */
    Structure *getFullStructure(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return the structure of the fragment for search result i
     */
    Structure getResultStructure(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return the residues that comprise the search result i
     */
    vector<Residue *> getResultResidues(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return a fragment info object representing search result i, which can
     *  be cast to a concrete subclass to retrieve specific information or
     *  used to query a FragmentFetcher
     */
    FragmentInfo *getFragmentInfo(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return an AtomPointerVector representing the structure of result i
     */
    AtomPointerVector getAPV(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return the RMSD between the original query and result i. This may have
     *  to be calculated during this call.
     */
    double getRMSD(int i);
    /**
     * @param i the index of the search result to retrieve
     * @return the Sequence of the given search result, with the residues mapping to
     *  the order of the original query
     */
    Sequence getResultSequence(int i);
private:
    ClusterSearchResults(FragmentFetcher *fetcher, RMSDCache *rmsdCalc, AtomPointerVector query, vector<_SearchResult> results): _fetcher(fetcher), _rmsdCalc(rmsdCalc), results(results), query(query) {}

    FragmentFetcher *_fetcher;
    RMSDCache *_rmsdCalc;
    vector<_SearchResult> results;
    AtomPointerVector query;
};

/**
 * A ClusterTree stores fragment structures in a tree format for optimized
 * lookup.
 */
class ClusterTree {
public:
    ClusterTree(FragmentFetcher *fragmentFetcher, int childLimit = 2, bool sharedCoordinates = false): _fragmentFetcher(fragmentFetcher), rmsdCalc(fragmentFetcher, sharedCoordinates), _childLimit(childLimit), _sharedCoordinates(sharedCoordinates) {};
    // Move constructor
    ClusterTree(ClusterTree &&other): _fragmentFetcher(other._fragmentFetcher), rmsdCalc(other.rmsdCalc), _childLimit(other._childLimit), _root(other._root) {
        // Remove ownership of the root node
        other._root = nullptr;
    }

    /**
     * Initializes a cluster tree by merging the given set of trees. All the
     * trees must have the same fragmentFetcher reference and child limit.
     *
     * @param mergeLimit if positive, indicates the number of levels to merge
     *  before performing a greedy assignment of the remaining children
     */
    ClusterTree(const vector<ClusterTree> &trees, int mergeLimit = -1);
    ~ClusterTree();

    /**
     * Reads the cluster tree from a file.
     */
    void read(string clusterFilePath);
    /**
     * Writes the cluster tree to file using the toString method.
     */
    void write(string clusterFilePath) const;
    /**
     * Performs the clustering.
     */
    void cluster();
    /**
     * Inserts a node into the tree with the given structure and info. The 
     * structure is used to determine where to place the node, but the info
     * object is what is actually stored.
     */
    void insert(const AtomPointerVector &item, FragmentInfo *idInfo);
    /**
     * Generates a tree-like string representation of the cluster tree. The
     * format of each line is [space]FRAGMENT,RMSD,RADIUS where [space] 
     * includes one space for each nest level, FRAGMENT is the toString value
     * for the node's fragment info, RMSD is the distance from this node to its
     * parent, and RADIUS is the maximum RMSD from this node to any of its
     * descendants.
     */
    string toString() const;

    /**
     * Sorts all children of each node in the tree by their subtree radius. This may
     * improve performance on some tree traversal tasks.
     */
    void sortSubtreesByRadius(ClusterNode<FragmentInfo> *root = nullptr);

    /**
     * Searches the database for a given query.
     *
     * @param query the structure to look for
     * @param rmsdCutoff the maximum RMSD to allow in results
     * 
     * @return a ClusterSearchResults object containing all of the fragments
     *  in the database with RMSD at most rmsdCutoff away from the query
     */
    ClusterSearchResults search(const AtomPointerVector &query, double rmsdCutoff);

    /**
     * Provide a positive size limit to stop the cluster() method after this
     * many fragments have been added.
     */
    int sizeLimit = -1;

    /**
     * If positive, represents the number of matches after which the search algorithm
     * terminates.
     */
    int sufficientMatches = -1;

    /**
     * If true, print statistics on the results of each query to stdout.
     */
    bool verbose = false;

    ClusterNode<FragmentInfo> *getRoot() const { return _root; }

    vector<ClusterNode<FragmentInfo> *> getNodesAtLevel(int level) const;

    FragmentFetcher *getFragmentFetcher() const { return _fragmentFetcher; };

private:
    FragmentFetcher *_fragmentFetcher;
    ClusterNode<FragmentInfo> *_root = nullptr;
    RMSDCache rmsdCalc;

    /**
     * If the size of the tree is an integer power of _refineLogBase, the 
     * tree's clustering will be refined.
     */
    int _refineLogBase = 4;
    /**
     * The maximum number of children any node in the tree (post-refinement) 
     * can have.
     */
    int _childLimit;
    /**
     * When inserting a new node into the tree, this provides the number of
     * nodes to consider for further recursion at each level of the tree.
     */
    int _numInsertOptions = 3;

    /**
     * Indicates whether all elements in the tree are
     * positioned in the same coordinate system. If so, an
     * in-place RMSD function will be used instead of the 
     * default optimal superposition.
     */
    bool _sharedCoordinates = false;

    // Stats for clustering
    int depthCount = 0;
    int depthTotal = 0;
    int clusteringOps = 0;

    // Stats for searching
    int numNodesTotal = 0;
    int numNodesSearched = 0;
    int numNodesEliminated = 0;
    int numNodesAutoAdded = 0;
    int numNodesEliminatedByAncestor = 0;
    int numNodesEliminatedByParent = 0;

    /**
     * Copies all the nodes and fragment info objects in the given parent's
     * subtree.
     */
    ClusterNode<FragmentInfo> *copySubtree(ClusterNode<FragmentInfo> *parent);

    /**
     * Builds a list of clustering candidates by adding the given node's
     * children, then adding the descendants one level at a time until there
     * are approximately _childLimit + _childLimit ^ 2 elements to cluster.
     *
     * @param parent the node whose children are to be clustered
     * @param elements the vector in which to place results
     */
    void addClusteringCandidates(ClusterNode<FragmentInfo> *parent, vector<ClusterNode<FragmentInfo> *> &elements);
    
    /**
     * Refines the clustering of the given node's children by performing k-
     * medoids clustering on its children and grandchildren, which determines
     * the new list of children nodes. This method also restructures the
     * tree so that the node closest to a "centroid" takes the place of a 
     * parent.
     *
     * This method supports both top-down and bottom-up recursive approaches,
     * which will lead to different clustering results. The top-down approach
     * guarantees that parent (or the node in its place) will have no more
     * than _childLimit children.
     *
     * @param parent the node whose children to re-cluster.
     * @param bottomUp if true, omit parent reassignment and cluster
     *  children first, then handle this parent - effectively a bottom-up
     *  recursive approach. If false, perform a top-down recursion.
     * @param levelLimit if nonnegative, indicates the level at which to
     *  stop refining the clustering. This may lead to a dirty tree (with
     *  more than _childLimit children per node).
     */
    void refineClustering(ClusterNode<FragmentInfo> *parent, bool bottomUp = false, int levelLimit = -1);

    /**
     * Assembles a list of candidates that could potentially replace the 
     * given node as a centroid, including the parent node itself, then calls
     * assignBestParent.
     *
     * @param parent the node to potentially reassign
     */
    ClusterNode<FragmentInfo>* reassignParent(ClusterNode<FragmentInfo> *parent);

    /**
     * Determines which of parent's clustering candidates has the smallest
     * maximum distance from the other candidates, and replaces parent in the
     * tree with that node.
     *
     * @param parent the node to potentially reassign
     * @param candidates the list of candidates that could replace the node
     */
    ClusterNode<FragmentInfo>* assignBestParent(ClusterNode<FragmentInfo> *parent, const vector<ClusterNode<FragmentInfo> *> &candidates);

    /**
     * Fills the given results vector with all nodes from the given level
     * of the tree, where level = 0 corresponds to the root.
     */
    void getNodesAtLevel(int level, vector<ClusterNode<FragmentInfo> *> &results, ClusterNode<FragmentInfo> *parent = nullptr) const;

    /**
     * Greedily assigns the children to the given set of parents such that
     * the subtree radii of the parents are roughly minimized.
     */
    void greedyCluster(vector<ClusterNode<FragmentInfo> *> parents, vector<ClusterNode<FragmentInfo> *> children);

    /**
     * Helper that sets parentRMSD and subtreeRadius properties for all
     * nodes. The nodes parameter indicates the path to the subtree to set,
     * where the last element is the root of the subtree.
     *
     * @param fullAncestry if true, calculate RMSDs to every node in the
     *  ancestor hierarchy for the last element
     */
    void calculateRMSDProps(vector<ClusterNode<FragmentInfo> *> nodes, bool fullAncestry = false);

    /**
     * Search helper function that adds an entire subtree recursively to the
     * search results.
     *
     * @param parent the root node of the subtree to add
     * @param results the vector into which to add results
     * @param parentRMSD the distance from the query to the parent, which will
     *  be stored in the result for the parent node
     */
    void addSubtreeToResults(ClusterNode<FragmentInfo> *parent, vector<_SearchResult> &results, double parentRMSD = DBL_MAX);

    /**
     * Recursive function that implements the search() method.
     *
     * @param query the structure to look for
     * @param rmsdCutoff the maximum distance from query to allow in
     *  returned search results
     * @param parent the root node of the subtree to search
     * @param results the vector in which to place results
     */
    void _search(const AtomPointerVector &query, double rmsdCutoff, ClusterNode<FragmentInfo> *parent, vector<_SearchResult> &results);
}; 

/**
 * A simple utility for iterating over all nodes in a cluster tree.
 * Iterates in DFS order.
 */
class ClusterIterator {
public:
    ClusterIterator(ClusterNode<FragmentInfo> *root): _root(root) {
        _stack.push_back(root);
        _indexes.push_back(0);
    };

    bool hasNext();
    ClusterNode<FragmentInfo> * next();

    void makeNextResult();
    void destroyNextResult();

private:
    ClusterNode<FragmentInfo> *_root;
    ClusterNode<FragmentInfo> *_nextResult = nullptr;

    vector<ClusterNode<FragmentInfo> *> _stack;
    vector<int> _indexes;
};

/**
 * Iterates over all pairs of nodes in a cluster tree.
 * Every pair of nodes will be visited exactly once
 * (regardless of order).
 */
class ClusterPairIterator {
public:
    ClusterPairIterator(const ClusterTree &tree): _tree(tree) {
        _firstStack.push_back(_tree.getRoot());
        _firstIndexes.push_back(0);
        _secondStack.push_back(_tree.getRoot());
        _secondIndexes.push_back(0);
        _nextIsAncestor = true;
    };
    ~ClusterPairIterator();

    bool hasNext();
    pair<ClusterNode<FragmentInfo> *, ClusterNode<FragmentInfo> *> next();
    bool secondDescendsFromFirst() { return _currentIsAncestor; };

    void makeNextResult();
    void destroyNextResult();

    void skipSecondSubtree();

private:
    const ClusterTree &_tree;

    pair<ClusterNode<FragmentInfo> *, ClusterNode<FragmentInfo> *> *_nextResult = nullptr;
    // Indicates whether in the next result, the first node
    // is an ancestor of the second node
    bool _currentIsAncestor = false;
    bool _nextIsAncestor = false;

    // If asked to skip second subtree when the recent stack involved a
    // push on top of the previous second node, remove the back node.
    // If the last stack update involved traversing upward at all, then
    // we just need to tell future updates not to push into the current
    // node (using _secondSubtreeSkipped).
    bool _secondSubtreePushed = false;
    bool _secondSubtreeSkipped = false;

    vector<ClusterNode<FragmentInfo> *> _firstStack;
    vector<int> _firstIndexes;
    vector<ClusterNode<FragmentInfo> *> _secondStack;
    vector<int> _secondIndexes;

    // Increments the tree traversal using a DFS iteration order.
    // If markerNode is removed from the stack, returns true. If markerNode
    // is not provided or not removed from the stack, returns false.
    // If skipSubtree is true, don't go into the current top node's children.
    bool updateStacks(vector<ClusterNode<FragmentInfo> *> &nodeStack, vector<int> &indexStack, ClusterNode<FragmentInfo> *markerNode = nullptr, bool skipSubtree = false, bool *pushed = nullptr);
};

#endif
