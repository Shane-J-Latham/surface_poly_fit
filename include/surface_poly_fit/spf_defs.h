
#ifndef SURFACE_POLY_FIT_SPF_DEFS_H
#define SURFACE_POLY_FIT_SPF_DEFS_H

#include <cstdint>
#include <memory>
#include <vector>
#include <array>

namespace spf
{
    typedef std::array<double, 3> Point;
    typedef std::array<double, 3> Vector;
    typedef std::array<std::int64_t, 3> TriFace;

    typedef std::vector<Point> PointStlVec;
    typedef std::shared_ptr<PointStlVec> PointStlVecPtr;
    typedef std::vector<Vector> VectorStlVec;
    typedef std::shared_ptr<VectorStlVec> VectorStlVecPtr;
    typedef std::vector<TriFace> TriFaceStlVec;
    typedef std::shared_ptr<TriFaceStlVec> TriFaceStlVecPtr;

    typedef std::pair<PointStlVecPtr,VectorStlVecPtr> PointVectorStlVecPair;
    typedef std::pair<PointStlVecPtr,TriFaceStlVecPtr> TriMeshPair;
}

#endif
