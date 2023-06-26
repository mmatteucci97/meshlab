#ifndef CUT_AND_FILL_H
#define CUT_AND_FILL_H

#include "filter_cut_and_fill.h"

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/bitquad_support.h>
#include <vcg/complex/algorithms/bitquad_creation.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>
#include <vcg/complex/algorithms/isotropic_remeshing.h>
#include <vcg/space/fitting3.h>
#include <wrap/gl/glu_tessellator_cap.h>

#include <array>
#include <map>

using namespace std;
using namespace vcg;
using namespace vcg::tri;

float FilterCutAndFillPlugin::DetermineAverageEdgeLength(const MeshModel &mesh)
{
    CMeshO m = mesh.cm;
    float averageEdgeLength = 0;
    int count = 0;

    for(auto fi = m.face.begin(); fi!=m.face.end(); ++fi)
    {
       for(int i=0; i<3; i++)
       {
           float dx = fi->P1(i).X() - fi->P0(i).X();
           float dy = fi->P1(i).Y() - fi->P0(i).Y();
           float dz = fi->P1(i).Z() - fi->P0(i).Z();

           float dist = sqrt(dx*dx + dy*dy + dz*dz);
           averageEdgeLength += dist;

           count++;
       }
    }

    averageEdgeLength /= count / 2;

    return averageEdgeLength;
}

/**
 *
*/
void CheckMeshRequirement(CMeshO m)
{
    tri::UpdateTopology<CMeshO>::FaceFace(m);

    for(auto fi = m.face.begin(); fi!=m.face.end(); ++fi)
    {
        for(int i=0; i<3; i++)
        {
            CMeshO::VertexType *v0 = (*fi).V0(i);
            CMeshO::VertexType *v1 = (*fi).V1(i);
            CMeshO::VertexType *v2 = (*fi).V2(i);

            if((v0->Q()==0 && v1->Q()==0) || (v1->Q()==0 && v2->Q()==0) || (v2->Q()==0 && v0->Q()==0))
            {
                if(face::IsBorder(*fi, i))
                {
                    throw MLException("Cannot apply filter. Mesh has a border on the plane to slice");
                }

                else if(!face::IsManifold(*fi, i))
                {
                    throw MLException("Cannot apply filter. Mesh is not manifold on the plane to slice");
                }
            }
        }
    }


}

void UpdateMesh(CMeshO &m)
{
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(m);
    tri::UpdateTopology<CMeshO>::FaceFace(m);
    tri::UpdateTopology<CMeshO>::VertexFace(m);
    tri::Clean<CMeshO>::RemoveDuplicateVertex(m);
    tri::UpdateBounding<CMeshO>::Box(m);
}

void SetMeshRequirements(CMeshO &m)
{
    m.vert.EnableNormal();
    m.vert.EnableVFAdjacency();
    m.vert.EnableNormal();
    m.vert.EnableMark();
    m.vert.EnableQuality();

    m.face.EnableNormal();
    m.face.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableQuality();
    m.face.EnableMark();
}

void SetMeshRequirements(MeshModel &m)
{
    m.updateBoxAndNormals();

    m.updateDataMask
            (
                MeshModel::MM_VERTFACETOPO |
                MeshModel::MM_FACEFACETOPO |
                MeshModel::MM_FACEQUALITY |
                MeshModel::MM_VERTQUALITY |
                MeshModel::MM_FACEMARK |
                MeshModel::MM_VERTMARK |
                MeshModel::MM_VERTNORMAL |
                MeshModel::MM_VERTCOLOR |
                MeshModel::MM_WEDGTEXCOORD |
                MeshModel::MM_FACENUMBER |
                MeshModel::MM_GEOMETRY_AND_TOPOLOGY_CHANGE
            );
    m.cm.face.EnableMark();
}

static void CreateSection(CMeshO &section, CMeshO &m)
{
    int newVertices=0;
    int newEdges=0;
    for(auto vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
    {
        if((*vi).Q() == 0)
        {
            CMeshO::VertexType newV = *tri::Allocator<CMeshO>::AddVertices(section, 1);
            newV = *vi;
            newVertices++;
        }
    }

    std::set<std::pair<CMeshO::VertexType *, CMeshO::VertexType *> > addedEdges;

    for(auto fi = m.face.begin(); fi!=m.face.end(); ++fi)
    {
        if(!(*fi).IsD())
        {
            for(int j=0; j<3; j++)
            {
                CMeshO::VertexType *v0 = (*fi).V0(j);
                CMeshO::VertexType *v1 = (*fi).V1(j);

                if(v0->Q() == 0 && v1->Q()==0)
                {
                    if(addedEdges.find(std::make_pair(v0, v1)) == addedEdges.end() && addedEdges.find(std::make_pair(v1, v0)) == addedEdges.end())
                    {
                        CMeshO::EdgeType &newE = *vcg::tri::Allocator<CMeshO>::AddEdges(section, 1);
                        newE.V(0) = v0;
                        newE.V(1) = v1;
                        addedEdges.insert(std::make_pair(v0, v1));
                        newEdges++;
                    }
                }
            }
        }
    }
    tri::UpdateNormal<CMeshO>::PerVertex(section);
}

static void FillMesh(CMeshO &overMesh, CMeshO &underMesh, CMeshO &section, Point3m planeCenter)
{
    tri::Append<CMeshO, CMeshO>::Mesh(overMesh, section);
    tri::Clean<CMeshO>::FlipMesh(section);
    tri::Append<CMeshO, CMeshO>::Mesh(underMesh, section);

    tri::Clean<CMeshO>::RemoveDuplicateVertex(overMesh);
    tri::Clean<CMeshO>::RemoveDuplicateVertex(overMesh);

    tri::Clean<CMeshO>::RemoveDuplicateFace(overMesh);

    tri::Clean<CMeshO>::RemoveDuplicateVertex(underMesh);
    tri::Clean<CMeshO>::RemoveDuplicateVertex(underMesh);

    tri::Clean<CMeshO>::RemoveDuplicateFace(underMesh);
}

void BoundaryExpand(CMeshO &m)
{
    tri::UpdateTopology<CMeshO>::FaceFace(m);
    tri::UpdateSelection<CMeshO>::FaceAll(m);

    std::vector<std::tuple<size_t, size_t, size_t>> newTrianglesVector;
    std::vector<Point3f> newTriangleCoordinates;
    int count = 0;

    for(auto fi = m.face.begin(); fi!=m.face.end(); ++fi)
    {
        for(int i=0; i<3; i++)
        {
            if(face::IsBorder(*fi, i))
            {
                newTrianglesVector.emplace_back(tri::Index(m, fi->V0(i)), tri::Index(m, fi->V1(i)), m.vert.size()+ newTriangleCoordinates.size());
                newTriangleCoordinates.push_back(fi->P0(i) + (fi->P1(i) - fi->P2(i)));
            }
        }
    }

    for(int i=0; i<newTriangleCoordinates.size(); i++)
    {
        tri::Allocator<CMeshO>::AddVertex(m, newTriangleCoordinates[i]);
    }

    for(int i=0; i<newTrianglesVector.size(); i++)
    {
        tri:Allocator<CMeshO>::AddFace(m, std::get<0>(newTrianglesVector[i]), std::get<2>(newTrianglesVector[i]), std::get<1>(newTrianglesVector[i]));
    }
}
#endif // CUT_AND_FILL_H
