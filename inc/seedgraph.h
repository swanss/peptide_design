//
//  seedgraph.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 1/17/19.
//

#ifndef seedgraph_h
#define seedgraph_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <fstream>

#include "msttypes.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "fileutilities.h"
#include "structure_iter.h"
#include "clustertree.h"
#include "clusterutils.h"
#include "findfuseable.h"

using namespace std;
using namespace MST;

/**
 Represents a graph in which each node is a residue, and edges between nodes
 represent the possibility of a bond between two residues. Peptides can thus be
 thought of as paths through this graph.
 */
class SeedGraph {
public:
    /**
     Initialize an empty seed graph.
     
     @param adjSameResidues whether the edges in this graph connect mutually-
            interchangeable residues (true) or residues that can be bonded
            together (false)
     @param sCache a structure cache to use as a source of residues. If this is
            non-null, the new seed graph will NOT manage its memory.
     @param centerOnly edges should only be defined between the two most center residues of the overlap
     */
    SeedGraph(bool adjSameResidues = false, StructureCache *sCache = nullptr, bool centerOnly = false);
    
    SeedGraph(string adjacencyPath, bool adjSameResidues, StructureCache *sCache = nullptr, string pdbPrefix = "", bool centerOnly = false);
    
    /// Copy constructor
    SeedGraph(const SeedGraph& other): adjSameResidues(other.adjSameResidues), centerOnly(other.centerOnly), structures(other.structures), ownsStructureCache(false), adjacencies(other.adjacencies), reverseAdjacencies(other.reverseAdjacencies) {}

    /**
     Initialize a seed graph by reading a CSV file containing fuse candidates.
     
     @param neighborsFile the CSV file to read
     @param pdbPrefix a path prefix to use before each path found in the
            fuse candidate metadata
     @param adjSameResidues whether to form the graph by connecting residues
            with good overlap to each other (true), or by connecting residues
            that could be joined with a bond (false)
     */
    SeedGraph(FuseCandidateFile neighborsFile, string pdbPrefix = "", bool adjSameResidues = false, bool centerOnly = false): adjSameResidues(adjSameResidues), centerOnly(centerOnly) {
        if (adjSameResidues == true and centerOnly == true) MstUtils::error("centerOnly cannot be applied while adjSameResidues is true");
        structures = new StructureCache(pdbPrefix); load(neighborsFile); }
    ~SeedGraph();
    /** Print the graph's adjacency list to cout */
    void print();

    StructureCache *getStructures() { return structures; }

    /**
     Adds to the seed graph from the given CSV file, prepending pdbPrefix before
     each path. If nearStructure is not null, it defines a structure with which all
     entries loaded into the graph must overlap (by bounding box with inset
     `inset`). If usedCandidatesFile is not null, the candidates that are
     loaded in will be written out to that file.
     */
    void load(FuseCandidateFile neighborsFile, string pdbPrefix = "", Structure *nearStructure = nullptr, float inset = 0.0, FuseCandidateFile *usedCandidatesFile = nullptr);

    /**
     Adds to the seed graph from the given list of structures, in the same
     manner as the FuseCandidateFile version but without any connections between
     seeds.
     */
    void load(SeedListFile structuresFile, string pdbPrefix = "", Structure *nearStructure = nullptr, float inset = 0.0);

    /**
     * Adds to the seed graph from the given binary file without any connections between seeds.
     */
    void load(StructuresBinaryFile *binaryFile, Structure *nearStructure = nullptr, float inset = 0.0);
    
    void loadCache();

    /**
     * Adds nodes and edges to the seed graph using a cluster tree. Pairs of nodes in
     * the tree are searched to check for RMSD less than the given cutoff, as well as
     * an optional cosine angle criterion.
     *
     * @param clusterTree a tree of overlaps in which to search. The tree should be created
     *  with shared coordinates = true for an appropriate definition of overlaps.
     * @param numResOverlap the number of residues that need to overlap. This is the same
     *  parameter used to create the cluster tree (i.e. the number of residues in each
     *  node).
     * @param rmsdCutoff the maximum RMSD allowed for an overlap
     * @param minCosAngle the minimum cos(angle) between the first and last alpha carbons
     *  in two k-mers needed to count as an overlap, in addition to the rmsd criterion
     */
    void load(ClusterTree &clusterTree, int numResOverlap, mstreal rmsdCutoff, mstreal minCosAngle = -1.0);

    /**
     Reads the graph from the given adjacency path, such as written out from
     the write() method.
     
     @param adjacencyPath the path to the file to read
     @param pdbPrefix a path to prepend paths found in the file
     @param convertToAdjSameResidue if true, assume the file was created with
            C-to-N adjacencies, and convert them to the adjSameResidue mode
     */
    void read(string adjacencyPath, string pdbPrefix = "", SeedGraph *convertToAdjSameResidue = nullptr);
    
    /**
     Writes the seed graph to file at the given path.
     
     @param path the path to write to
     */
    void write(string path);
    
    /**
     Gets a writable identifier for the given residue.
     
     @param res the residue to write
     @return a string to write representing the residue
     */
    string writeCodeForResidue(Residue *res);
    
    /**
     Generates and returns a list of seed graphs that are isolated from
     each other (i.e. no edges between them in the original graph).
     */
    vector<SeedGraph> subgraphs();
    /**
     Generates the graph consisting of the neighborhood of the given
     residue.
     */
    SeedGraph neighborhood(Residue *residue);
    
    /**
     Returns the result of combining this graph with the given one, without
     modifying either one.
     
     @param graph the graph to take the union with
     @return the union of this graph and the given one
     */
    SeedGraph unionWith(SeedGraph &graph);
    
    /**
     Returns the result of removing all residues and adjacencies associated with
     the given subgraph. Does not modify either original graph.
     
     @param subgraph the subgraph to remove
     @return the result of removing the given subgraph from this one
     */
    SeedGraph removing(SeedGraph &subgraph);
    
    /** @return the number of residues in this graph */
    size_t residueSize() { return adjacencies.size(); };
    
    /** @return the set of residues in this graph */
    unordered_set<Residue *> getResidues();

    /** @return whether or not the graph contains the given residue */
    bool contains(Residue *res) { return adjacencies.count(res) != 0; }

    /** @return the number of seeds that contribute residues to this graph */
    size_t seedSize();
    
    /** @return the adjacencies and reverse adjacencies for the given residue */
    unordered_set<Residue *> bidirectionalNeighbors(Residue *res);
    
    /**
     @return the residues to which this residue is connected from C to N
     @throws assertion error if the graph is adjSameResidue
     */
    unordered_set<Residue *> forwardNeighbors(Residue *res);
    
    /**
     @return the residues to which this residue is connected from N to C
     @throws assertion error if the graph is adjSameResidue
     */
    unordered_set<Residue *> backwardNeighbors(Residue *res);

    /**
     Returns a graph with transformed edges; requires adjSameResidues to be false
     for this graph.
     
     @return the equivalent graph where edges between residues indicate direct
        overlap, not the possibility of a bond
     */
    SeedGraph withAdjSameResidue();
    
    /**
     Finds residues whose flanking regions overlap those of the given residue
     down to the given RMSD cutoff.
     
     @param residue the residue to compare to
     @param flankLength the number of residues on each side of residue to compare
     @param rmsdCutoff the maximum RMSD between the flanking regions of residue
            and those of a candidate in the graph
     @param rmsds a vector that upon return will contain each result segment's
            RMSD to the target residue segment
     
     @return a list of residues whose flanking regions match those of the given
             residue
     */
    vector<Residue *> similarContactResidues(Residue *residue, int flankLength, float rmsdCutoff, vector<float> *rmsds = nullptr);
    
    /**
     Gets a residue from a file-written identifier of the form PDB:index. If
     loadIfNeeded is true, loads the structure if not already loaded. Otherwise,
     returns nullptr.
     
     @param residueID the file-written identifier
     @param loadIfNeeded whether to load the structure if not already present
     @return a residue for that identifier
     */
    Residue * getResidueFromFile(string residueID, bool loadIfNeeded = true);
    
    /*
     added by sebastian on 08/29/20
     
     Gets a structure by its name in the structureCache. If loadIfNeeded is true, loads the
     structure if not already loaded. Otherwise, returns nullptr.
     */
    Structure * getStructureFromFile(string structureName, bool loadIfNeeded = true);
    
    /**
     Assembles the _representativeResidues map so that every residue is mapped
     to a single residue in its overlap neighborhood.
     */
    void computeRepresentativeResidues();

    /**
     Gets a single residue that represents this residue and its overlap cluster.
     */
    Residue * representativeResidue(Residue *res) {
        return _representativeResidues[res];
    }
    
    /**
     @return the number of edges in the graph
     */
    int edgeSize() {
        return numEdges;
    }
    
    /**
     Computes a topological order for the graph, assuming it is acyclic.
     
     @return the vertices in such an order that no vertex has an edge to a
             vertex that appears after it
     */
    vector<Residue *> topologicalOrder();
    
private:
    typedef unordered_map<Residue *, unordered_set<Residue *>> ResidueAdjacencyList;
    
    /// Stores a single residue that represents every residue and its associated cluster
    unordered_map<Residue *, Residue *> _representativeResidues;

    /**
     Cache of Structure objects loaded during the building of the graph. May
     be read in order to use the same structure objects in later steps.
     */
    StructureCache *structures = nullptr;
    bool ownsStructureCache = true;

    ResidueAdjacencyList adjacencies;
    ResidueAdjacencyList reverseAdjacencies;
    bool adjSameResidues = false; // If true, edges connect analogous residues (with low RMSD between themselves); otherwise, edges connect sequential residues
    
    
    /*
     The following parameter is only applicable if adjSameResidues = false. When centerOnly is set to
     true, edges are only defined in the center of the overlap.
     
     e.g.
     
     overlap1 = 3
     overlap2 = 2
     overlapSize = 4
     
     seed 1: 1->2->3->4->5->6
                       X
     seed 2:    1->2->3->4->5->6
     
     two edges are defined: seed 1 residue 4 to seed 2 residue 4, and seed 2 residue 3 to seed 1 residue 5.
     
     note: this is only defined when overlapSize is even.
     */
    
    bool centerOnly = false;
    
    int numEdges = 0;
    
    string seed_chain_ID = "0";
    
    /**
     Adds adjacencies from chain1, including edges that join chain1 to chain2.
     */
    void addAdjacencies(Chain *chain1, Chain *chain2, int overlap1, int overlap2, int overlapSize);
    
    /**
     Adds an edge from res1 to res2.
     
     @param res1 the source node
     @param res2 the target node
     */
    void addEdge(Residue *res1, Residue *res2);
    
    /** Adds the given residue to the graph. */
    void addResidue(Residue *res);
};

template <class T>
class SeedGraphMap: public SeedGraph {
public:
    SeedGraphMap(StructureCache *sCache = nullptr): SeedGraph(true, sCache) {};
    
    void setValue(Residue *res, T value) {
        MstUtils::assert(contains(res), "adding value to residue not contained in graph");
        _storage[res] = value;
    }
    
    T* value(Residue *res) {
        unordered_set<Residue *> visited;
        return recursiveValue(res, visited);
    }
    
    /**
     Allows clients to iterate over the scores available in this map (in no particular order).
     */
    class iterator {
    public:
        typedef iterator self_type;
        typedef pair<Residue *, T> value_type;
        typedef int difference_type;
        typedef forward_iterator_tag iterator_category;
        iterator(unordered_map<string, Structure *>::iterator it) : it_(it) { }
        iterator(const iterator& other): it_(other.it_) { }
        self_type operator++() { it_++; return *this; }
        self_type operator++(int junk) { self_type i = *this; it_++; return i; }
        Structure * operator*() { return it_->second; }
        Structure * const * operator->() { return &(it_->second); }
        bool operator==(const self_type& rhs) { return it_ == rhs.it_; }
        bool operator!=(const self_type& rhs) { return it_ != rhs.it_; }
    private:
        typename unordered_map<Residue *, T>::iterator it_;
    };
    
    /**
     Iterator pointing to the beginning of the structure cache.
     */
    iterator begin() {
        iterator it(_storage.begin());
        return it;
    }
    
    /**
     Iterator pointing to the end of the structure cache.
     */
    iterator end() {
        iterator it(_storage.end());
        return it;
    }


private:
    unordered_map<Residue *, T> _storage;

    T* recursiveValue(Residue *res, unordered_set<Residue *> &visitedResidues) {
        visitedResidues.insert(res);
        if (_storage.count(res) != 0) {
            return &_storage[res];
        } else {
            for (Residue *res: bidirectionalNeighbors(res)) {
                if (visitedResidues.count(res) != 0) continue;
                
                T* subValue = recursiveValue(res, visitedResidues);
                if (subValue != nullptr)
                    return subValue;
            }
        }
        return nullptr;
    }
};

#endif /* seedgraph_h */
