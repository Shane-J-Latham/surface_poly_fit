
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <surface_poly_fit/polyhedral_surface_py.h>
#include <surface_poly_fit/spf_cgal.h>

namespace spf
{

namespace py = pybind11;

template <typename TStlVec>
std::shared_ptr<TStlVec>
array_to_stlvec(py::object threeDObj)
{
    typedef TStlVec StlVec;
    typedef std::shared_ptr<StlVec> StlVecPtr;
    typedef typename StlVec::value_type ThreeD;

    StlVecPtr stlVecPtr(new StlVec());
    py::object j_objs[] = {py::cast(long(0)), py::cast(long(1)), py::cast(long(2))};
    for (long i = 0; i < py::len(threeDObj); ++i)
    {
        ThreeD p;
        py::object i_obj(py::cast(i));
        for (long j = 0; j < 3; ++j)
        {
            p[j] = py::cast<typename ThreeD::value_type>(threeDObj[i_obj][j_objs[j]]);
        }
        stlVecPtr->push_back(p);
    }
    return stlVecPtr;
}

template <typename TStlVec>
py::object
stlvec_to_array(const TStlVec & stlVec)
{
    typedef TStlVec StlVec;
    typedef typename StlVec::value_type ThreeD;

    py::dtype dtyp(py::dtype::of<typename ThreeD::value_type>());

    size_t shape[2]{stlVec.size(), 3};
    // py::array ary(dtyp, shape);
    auto ary = py::array_t<typename ThreeD::value_type>(shape);

    for (std::size_t i = 0; i < stlVec.size(); i++)
    {
        for (std::size_t j = 0; j < 3; j++)
        {
            ary.mutable_at(i, j) = stlVec[i][j];
        }
    }
    return py::object(ary);
}

PointStlVecPtr array_to_points(py::object pointsObj)
{
    return array_to_stlvec<PointStlVec>(pointsObj);
}

VectorStlVecPtr array_to_vectors(py::object vectorsObj)
{
    return array_to_stlvec<VectorStlVec>(vectorsObj);
}

py::object
trimeshpair_to_tuple(const TriMeshPair & meshPair)
{
    PointStlVecPtr points;
    TriFaceStlVecPtr faces;
    points = meshPair.first;
    faces = meshPair.second;
    py::object pointsObj(stlvec_to_array(*points));
    py::object facesObj(stlvec_to_array(*faces));

    return py::make_tuple(pointsObj, facesObj);
}

py::object
points_normal_pair_to_tuple(const PointVectorStlVecPair & pointNormalPair)
{
    PointStlVecPtr points;
    VectorStlVecPtr normals;
    points = pointNormalPair.first;
    normals = pointNormalPair.second;
    py::object pointsObj(stlvec_to_array(*points));
    py::object normalsObj(stlvec_to_array(*normals));

    return py::make_tuple(pointsObj, normalsObj);
}

}

PYBIND11_MODULE(_spf_cgal, m)
{
    spf::export_polyhedral_surface(m);
}
