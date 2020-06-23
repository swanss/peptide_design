/*
 * This inclusion should be put at the beginning.  It will include <Python.h>.
 */
#include <boost/python.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <unordered_set>
#include "structure_iter.h"
#include "fileutilities.h"
#include "clustertree.h"
#include "clusterutils.h"
#include "seedgraph.h"

/// Wrapper that allows FragmentInfo to be polymorphic
struct FragmentInfoWrapper : FragmentInfo, boost::python::wrapper<FragmentInfo> {
    FragmentInfoWrapper(): boost::python::wrapper<FragmentInfo>() {}
    FragmentInfoWrapper(const FragmentInfo& rhs): FragmentInfo(rhs), boost::python::wrapper<FragmentInfo>() {}
    string toString() {
        return this->get_override("toString")();
    }
};

/// Code to allow converting an std::pair to a tuple
template <typename T1, typename T2>
struct std_pair_to_tuple {
    static PyObject* convert(std::pair<T1, T2> const& p) {
        return boost::python::incref(boost::python::make_tuple(p.first, p.second).ptr());
    }
    static PyTypeObject const *get_pytype () {return &PyTuple_Type; }
};

template <typename T1, typename T2>
struct std_pair_to_python_converter {
    std_pair_to_python_converter() {
        boost::python::to_python_converter<std::pair<T1, T2>, std_pair_to_tuple<T1, T2>, true>();
    }
};

// class ResidueSet : public unordered_set<Residue *>, boost::python::wrapper<unordered_set<Residue *>> {
//     ResidueSet(): boost::python::wrapper<unordered_set<Residue *>>() {}
//     ResidueSet(const unordered_set<Residue *> &rhs): unordered_set<Residue *>(rhs), boost::python::wrapper<unordered_set<Residue *>>() {}
// 
//     bool contains(Residue *res) { return count(res) != 0; }
// };

/*
 * This is a macro Boost.Python provides to signify a Python extension module.
 */
BOOST_PYTHON_MODULE(peptide_design) {
    // An established convention for using boost.python.
    using namespace boost::python;

    std_pair_to_python_converter<AtomPointerVector, FragmentInfo *>();

    class_<unordered_set<Residue *>>("ResidueSet", no_init)
        .def("__contains__", +[](const unordered_set<Residue *> &s, Residue *res) { return s.count(res) != 0; })
        .def("__len__", &unordered_set<Residue *>::size)
        .def("__iter__", range<return_value_policy<reference_existing_object>>(
                    static_cast<unordered_set<Residue *>::iterator (unordered_set<Residue *>::*)()>(&unordered_set<Residue *>::begin), 
                    static_cast<unordered_set<Residue *>::iterator (unordered_set<Residue *>::*)()>(&unordered_set<Residue *>::end)))
    ;

    class_<StructureCache>("StructureCache", init<StructuresBinaryFile *, long>())
        .def("__getitem__", +[](StructureCache &c, string name) { return c.getStructure(name, ""); }, return_value_policy<reference_existing_object>())
        .def("__contains__", &StructureCache::hasStructure)
        .def("__delitem__", &StructureCache::removeStructure)
    ;

    class_<SeedGraph>("SeedGraph", init<string, bool, StructureCache *, string>())
        .def("write", &SeedGraph::write)
        .def("writeCodeForResidue", &SeedGraph::writeCodeForResidue)
        .def("neighborhood", &SeedGraph::neighborhood)
        .def("unionWith", &SeedGraph::unionWith)
        .def("removing", &SeedGraph::removing)
        .def("__len__", &SeedGraph::residueSize)
        .def("__contains__", &SeedGraph::contains)
        .def("numSeeds", &SeedGraph::seedSize)
        .def("bidirectionalNeighbors", &SeedGraph::bidirectionalNeighbors)
        .def("forwardNeighbors", &SeedGraph::forwardNeighbors)
        .def("backwardNeighbors", &SeedGraph::backwardNeighbors)
        .def("computeRepresentativeResidues", &SeedGraph::computeRepresentativeResidues)
        .def("representativeResidue", &SeedGraph::representativeResidue, return_value_policy<reference_existing_object>())
        .def("numEdges", &SeedGraph::edgeSize)
    ;
    /*class_<vector<MST::Structure>>("StructureList")
        .def(vector_indexing_suite<vector<Structure>>());

    class_<Fragmenter>("Fragmenter", init<Structure&, FragmentParams&>())
            .def("fragment", &Fragmenter::fragment)
            .def("getFragmentStructures", &Fragmenter::getFragmentStructures)
    ;*/
    class_<StructuresBinaryFile, boost::noncopyable>("StructuresBinaryFile", init<string>())
        .def("hasNext", &StructuresBinaryFile::hasNext)
        .def("next", &StructuresBinaryFile::next, return_value_policy<reference_existing_object>())
        .def("reset", &StructuresBinaryFile::reset)
        .def("__getitem__", &StructuresBinaryFile::getStructureNamed, return_value_policy<reference_existing_object>())
        .def("skip", &StructuresBinaryFile::skip)
        .def("scanFilePositions", &StructuresBinaryFile::scanFilePositions)
        .def("__len__", &StructuresBinaryFile::structureCount)
    ;

    class_<ClusterNode<FragmentInfo>, boost::noncopyable>("ClusterNode", no_init)
        .add_property("item", make_function(&ClusterNode<FragmentInfo>::getItem, return_value_policy<reference_existing_object>()))
        .add_property("children", &ClusterNode<FragmentInfo>::getChildren)
        .add_property("parent", make_function(&ClusterNode<FragmentInfo>::getParent, return_value_policy<reference_existing_object>()))
        .def_readonly("subtreeRadius", &ClusterNode<FragmentInfo>::subtreeRadius)
        .def_readonly("parentRMSD", &ClusterNode<FragmentInfo>::parentRMSD)
        .def("subtreeSize", &ClusterNode<FragmentInfo>::subtreeSize)
    ;

    class_<vector<ClusterNode<FragmentInfo> *>>("ClusterNodeList")
        .def(vector_indexing_suite<vector<ClusterNode<FragmentInfo> *>>());

    class_<ClusterTree, boost::noncopyable>("ClusterTree", init<FragmentFetcher *, int, bool>())
        .def("cluster", &ClusterTree::cluster)
        .def("read", &ClusterTree::read)
        .def("write", &ClusterTree::write)
        .def_readwrite("sufficientMatches", &ClusterTree::sufficientMatches)
        .def("__str__", &ClusterTree::toString)
        .def("search", &ClusterTree::search)
        .add_property("root", make_function(&ClusterTree::getRoot, return_value_policy<reference_existing_object>()))
        .def("getNodesAtLevel", static_cast<vector<ClusterNode<FragmentInfo> *> (ClusterTree::*) (int) const>(&ClusterTree::getNodesAtLevel))
    ;

    class_<FragmentFetcher, boost::noncopyable>("FragmentFetcher", no_init);
    class_<FragmentInfoWrapper, boost::noncopyable>("FragmentInfo", no_init)
        .def("__str__", &FragmentInfoWrapper::toString);

    class_<PairFragmentFetcher, bases<FragmentFetcher>>("PairFragmentFetcher", init<string, string, int>())
        .def("getAPV", &PairFragmentFetcher::getAPV)
        .def("next", &PairFragmentFetcher::next)
        .def("hasNext", &PairFragmentFetcher::hasNext)
        .def("skip", &PairFragmentFetcher::skip)
        .def("reset", &PairFragmentFetcher::reset)
        .def("numAtoms", &PairFragmentFetcher::numAtoms)
    ;
    class_<SingleFragmentFetcher, bases<FragmentFetcher>>("SingleFragmentFetcher", init<StructuresBinaryFile *, int, string>())
        .def("getAPV", &SingleFragmentFetcher::getAPV)
        .def("next", &SingleFragmentFetcher::next)
        .def("hasNext", &SingleFragmentFetcher::hasNext)
        .def("skip", &SingleFragmentFetcher::skip)
        .def("reset", &SingleFragmentFetcher::reset)
        .def("numAtoms", &SingleFragmentFetcher::numAtoms)
        .def("infoToString", +[](const SingleFragmentFetcher &fetcher, FragmentInfo *info) { return info->toString(); })
    ;
    class_<BatchWorkerFragmentFetcher, bases<FragmentFetcher>>("BatchWorkerFragmentFetcher", init<FragmentFetcher *, int, int>())
        .def("getAPV", &BatchWorkerFragmentFetcher::getAPV)
        .def("next", &BatchWorkerFragmentFetcher::next)
        .def("hasNext", &BatchWorkerFragmentFetcher::hasNext)
        .def("skip", &BatchWorkerFragmentFetcher::skip)
        .def("reset", &BatchWorkerFragmentFetcher::reset)
        .def("numAtoms", &BatchWorkerFragmentFetcher::numAtoms)
    ;

    class_<ClusterSearchResults>("ClusterSearchResults", no_init)
        .def("__len__", &ClusterSearchResults::size)
        .def("__getitem__", &ClusterSearchResults::getFragmentInfo, return_value_policy<manage_new_object>())
        .def("getFullStructure", &ClusterSearchResults::getFullStructure, return_value_policy<reference_existing_object>())
        .def("getResultStructure", &ClusterSearchResults::getResultStructure)
        .def("getResultSequence", &ClusterSearchResults::getResultSequence)
        .def("getAPV", &ClusterSearchResults::getAPV)
        .def("getRMSD", &ClusterSearchResults::getRMSD)
    ;
}
