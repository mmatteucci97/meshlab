/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

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

using namespace std;
using namespace vcg;
using namespace vcg::tri;

/**
 * @brief Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions. If you want to add icons to
 *  your filtering actions you can do here by construction the QActions accordingly
 */
FilterCutAndFillPlugin::FilterCutAndFillPlugin()
{ 
    typeList = {
        FP_CUT_AND_FILL,
        FP_CUT_FILL_AND_REMESH
    };

	for(ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QString FilterCutAndFillPlugin::pluginName() const
{
    return "FilterCutAndFill";
}

/**
* @brief The FilterClass describes in which generic class of filters it fits.
* This choice affect the submenu in which each filter will be placed
* More than a single class can be chosen.
* @param a: the action of the filter
* @return the class od the filter
*/
FilterCutAndFillPlugin::FilterClass FilterCutAndFillPlugin::getClass(const QAction *a) const
{
    switch(ID(a))
    {
        case FP_CUT_AND_FILL :return FilterPlugin::Smoothing;
        case FP_CUT_FILL_AND_REMESH: return FilterPlugin::Smoothing;
        default :
            assert(0);
            return FilterPlugin::Generic;
    }
    return FilterPlugin::Generic;
}

/**
 * @brief FilterSamplePlugin::getPreConditions
 * @return
 */
int FilterCutAndFillPlugin::getPreConditions(const QAction *filter) const
{
    switch(ID(filter))
    {
        case FP_CUT_AND_FILL: return MeshModel::MM_VERTFACETOPO;
        case FP_CUT_FILL_AND_REMESH: return MeshModel::MM_FACENUMBER;
    }

    return MeshModel::MM_NONE;
}

int FilterCutAndFillPlugin::getRequirements(const QAction* filter) const
{
    switch (ID(filter)){
        case FP_CUT_AND_FILL: return MeshModel::MM_NONE;
        case FP_CUT_FILL_AND_REMESH: return MeshModel::MM_FACEQUALITY | MeshModel::MM_VERTQUALITY;
        default:
            return MeshModel::MM_NONE;
    }
}

/**
 * @brief FilterSamplePlugin::pythonFilterName if you want that your filter should have a different
 * name on pymeshlab, use this function to return its python name.
 * @param f
 * @return
 */
QString FilterCutAndFillPlugin::pythonFilterName(ActionIDType f) const
{
    switch(f) {
    case FP_CUT_AND_FILL :
        return tr("generate_two_parts_of_the_mesh");
    case FP_CUT_FILL_AND_REMESH:
        return tr("generate_two_parts_of_the_mesh_than_remesh_the_surface");
    default :
        assert(0);
        return QString();
    }
}

/**
 * @brief ST() must return the very short string describing each filtering action
 * (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterCutAndFillPlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
    case FP_CUT_AND_FILL :
        return "Cut And Fill";
    case FP_CUT_FILL_AND_REMESH:
        return "Cut Fill and Remesh";
    default :
		assert(0);
		return QString();
	}
}

/**
 * @brief // Info() must return the longer string describing each filtering action
 * (this string is used in the About plugin dialog)
 * @param filterId: the id of the filter
 * @return an info string of the filter
 */
 QString FilterCutAndFillPlugin::filterInfo(ActionIDType filterId) const
{
	switch(filterId) {
    case FP_CUT_AND_FILL :
        return "Cut the mesh along a plane and fill both the parts";
    case FP_CUT_FILL_AND_REMESH:
        return tr("Cut the mesh along a plane and fill both of the parts; teìhen remesh the surface");
	default :
		assert(0);
		return "Unknown Filter";
	}
}

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
* @brief This function define the needed parameters for each filter. Return true if the filter has some parameters
* it is called every time, so you can set the default value of parameters according to the mesh
* For each parameter you need to define,
* - the name of the parameter,
* - the default value
* - the string shown in the dialog
* - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
* @param action
* @param m
* @param parlst
*/
RichParameterList FilterCutAndFillPlugin::initParameterList(const QAction *action,const MeshModel &m)
{
    RichParameterList parlst;
    float maxVal;
    QStringList curvCalcMethods;
    QStringList curvColorMethods;
    QStringList loopWeightLst;
    float averageLen;
    float averageArea;

    QStringList axis = QStringList() <<"X Axis"<<"Y Axis"<<"Z Axis"<<"Custom Axis";
    parlst.addParam(RichEnum   ("planeAxis", 1, axis, tr("Plane perpendicular to"), tr("The Slicing plane will be done perpendicular to the axis")));
    parlst.addParam(RichDirection("customAxis",Point3f(0,1,0),"Custom axis","Specify a custom axis, this is only valid if the above parameter is set to Custom"));

    parlst.addParam(RichFloat  ("planeOffset", 0.0, "Cross plane offset", "Specify an offset of the cross-plane. The offset corresponds to the distance from the point specified in the plane reference parameter. By default (Cross plane offset == 0)"));
    parlst.addParam(RichEnum   ("relativeTo",0,QStringList()<<"Bounding box center"<<"Bounding box min"<<"Origin","plane reference","Specify the reference from which the planes are shifted"));

    parlst.addParam(RichBool("createSectionSurface",false,"Create also section surface","If selected, in addition to a layer with the section polyline, it will be created also a layer with a triangulated version of the section polyline. This only works if the section polyline is closed"));
    parlst.addParam(RichBool("createOverMesh", true, "Create also over mesh", "If selected, it will create another layer with the portion of the mesh over the section plane. It requires manifoldness of the mesh "));
    parlst.addParam(RichBool("createUnderMesh", true, "Create also under mesh", "If selected, it will create another layer with the portion of the mesh under the section plane. It requires manifoldness of the mesh "));

    switch(ID(action)) {
    case FP_CUT_AND_FILL :
    {
    }
        break;
    case FP_CUT_FILL_AND_REMESH:
    {
        maxVal = m.cm.bbox.Diag();
        averageLen = DetermineAverageEdgeLength(m);

        parlst.addParam(RichBool ("TargetLenFromAverageLen", true, "Calculate the remeshing maximum edge length from the average edge length of the original mesh"));
        parlst.addParam(RichPercentage("TargetLenAverage", averageLen, 0.01*averageLen, averageLen, "Target Edge Len Average", "Use this value of the average edge length"));
    }
        break;
    default :
        assert(0);
    }
    return parlst;
}

/*
 * todo: check only the vertices and the faces marked with quality=0 after determining the items that will be cut
 * check if they are non manifold and if they are of border
*/
void CheckMeshRequirement(MeshModel *m)
{
    if (tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m->cm)>0 || (tri::Clean<CMeshO>::CountNonManifoldVertexFF(m->cm,false) != 0))
    {
        throw MLException("Mesh is not two manifold, cannot apply filter");
    }
}

void UpdateMesh(CMeshO &m)
{
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(m);
    tri::UpdateTopology<CMeshO>::FaceFace(m);
    tri::UpdateTopology<CMeshO>::VertexFace(m);
    tri::UpdateBounding<CMeshO>::Box(m);
}

void SetMeshRequirements(CMeshO &m)
{
    tri::UpdateBounding<CMeshO>::Box(m);
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(m);

    m.vert.EnableVFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<CMeshO>::VertexFace(m);

    m.face.EnableFFAdjacency();
    tri::UpdateTopology<CMeshO>::FaceFace(m);

    m.face.EnableQuality();
    m.vert.EnableQuality();

    m.face.EnableMark();
    m.vert.EnableMark();

    m.vert.EnableNormal();
    m.face.EnableWedgeTexCoord();
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
                MeshModel::MM_WEDGTEXCOORD |
                MeshModel::MM_FACENUMBER |
                MeshModel::MM_GEOMETRY_AND_TOPOLOGY_CHANGE
            );
    m.cm.face.EnableMark();
}

static void CreateSection(CMeshO &section, CMeshO &m)
{
    // Add in section all vertices with Q==0
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

    // Add in section all edges composed by 2 vertices with Q==0
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
}

static void FillMesh(CMeshO &outputMesh, CMeshO &section, CMeshO &inputMesh, bool flipSection=false)
{
    CMeshO sectionCopy;
    tri::Append<CMeshO, CMeshO>::Mesh(sectionCopy, section);

    if(flipSection)
    {
        tri::Clean<CMeshO>::FlipMesh(sectionCopy);
    }

    tri::Append<CMeshO, CMeshO>::Mesh(outputMesh, sectionCopy);
    tri::Append<CMeshO, CMeshO>::Mesh(outputMesh, inputMesh);

    tri::Clean<CMeshO>::RemoveDuplicateVertex(outputMesh);
    tri::Clean<CMeshO>::RemoveDuplicateVertex(outputMesh);
    tri::Clean<CMeshO>::RemoveDuplicateFace(outputMesh);
    UpdateMesh(outputMesh);
}

static void BoundaryExpand(CMeshO &m)
{
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
                //AGGIUNGI UNA FACCIA PER OGNI EDGE DI BOUNDARY

                newTrianglesVector.emplace_back(tri::Index(m, fi->V0(i)), tri::Index(m, fi->V1(i)), m.vert.size()+ newTriangleCoordinates.size());


                newTriangleCoordinates.push_back(fi->P0(i) + (fi->P1(i) - fi->P2(i)));
            }
        }
    }

    //AGGIUNGI I TRIANGOLI A M
    for(int i=0; i<newTriangleCoordinates.size(); i++)
    {
        tri::Allocator<CMeshO>::AddVertex(m, newTriangleCoordinates[i]);
    }

    for(int i=0; i<newTrianglesVector.size(); i++)
    {
        tri:Allocator<CMeshO>::AddFace(m, std::get<0>(newTrianglesVector[i]), std::get<2>(newTrianglesVector[i]), std::get<1>(newTrianglesVector[i]));
    }
}


/**
* @brief The Real Core Function doing the actual mesh processing.
* @param action
* @param md: an object containing all the meshes and rasters of MeshLab
* @param par: the set of parameters of each filter
* @param cb: callback object to tell MeshLab the percentage of execution of the filter
* @return true if the filter has been applied correctly, false otherwise
*/
std::map<std::string, QVariant> FilterCutAndFillPlugin::applyFilter(const QAction * filter, const RichParameterList & par, MeshDocument &md, unsigned int& /*postConditionMask*/, vcg::CallBackPos *cb)
{
    std::map<std::string, QVariant> outputValues;
    MeshModel & m = *md.mm();
    // the mesh has to return to its original position
    if (m.cm.Tr != Matrix44m::Identity())
        tri::UpdatePosition<CMeshO>::Matrix(m.cm, Inverse(m.cm.Tr), true);

    // the mesh has to be correctly transformed
    if (m.cm.Tr != Matrix44m::Identity())
        tri::UpdatePosition<CMeshO>::Matrix(m.cm, m.cm.Tr, true);


    Point3m planeAxis(0,0,0);
    Scalarm planeOffset;
    Point3m planeCenter;
    Plane3m slicingPlane;

    Box3m bbox=m.cm.bbox;
    MeshModel* base=&m;
    MeshModel* orig=&m;

    // making up new layer name
    QString sectionName = QFileInfo(base->shortName()).baseName() + "_section_surface";

    CMeshO sectionPolyline;
    CMeshO sectionSurface;
    CMeshO underM;
    CMeshO overM;
    CMeshO underFM;
    CMeshO overFM;

    //MeshModel* sectionPolyline= md.addNewMesh("",sectionName,true);
    //MeshModel* sectionSurface= md.addNewMesh("",sectionName+"_filled");
    //MeshModel* underM= md.addNewMesh("",sectionName+"_under");
    //MeshModel* underFM= md.addNewMesh("",sectionName+"_under_filed");
    //MeshModel* overM= md.addNewMesh("",sectionName+"_over");
    //MeshModel* overFM= md.addNewMesh("",sectionName+"_over_filled");

    SetMeshRequirements(m);

    SetMeshRequirements(sectionPolyline);
    SetMeshRequirements(sectionSurface);

    SetMeshRequirements(underM);
    SetMeshRequirements(overM);

    SetMeshRequirements(underFM);
    SetMeshRequirements(overFM);

    tri::QualityMidPointFunctor<CMeshO> slicingfunc(0.0);
    tri::QualityEdgePredicate<CMeshO> slicingpred(0.0,0.0);

    int ind;

    ind = par.getEnum("planeAxis");
    if(ind>=0 && ind<3)
        planeAxis[ind] = 1.0f;
    else
        planeAxis=par.getPoint3m("customAxis");

    planeAxis.Normalize();
    planeOffset = par.getFloat("planeOffset");

    // Check if the mesh has the correct topology to perform the algorithm
    // TODO: check only if the mesh has boundary or if is it two-manifold only on the part of the mesh around the plane
    CheckMeshRequirement(base);

    switch(RefPlane(par.getEnum("relativeTo")))
    {
        case REF_CENTER:  planeCenter = bbox.Center()+ planeAxis*planeOffset*(bbox.Diag()/2.0);      break;
        case REF_MIN:     planeCenter = bbox.min+planeAxis*planeOffset*(bbox.Diag()/2.0);    break;
        case REF_ORIG:    planeCenter = planeAxis*planeOffset;  break;
    }



    //planeCenter+=planeAxis*planeDist ;
    slicingPlane.Init(planeCenter,planeAxis);

    tri::Append<CMeshO,CMeshO>::Mesh(underM,orig->cm);
    tri::UpdateQuality<CMeshO>::VertexFromPlane(underM, slicingPlane);

    tri::UpdateTopology<CMeshO>::FaceFace(underM);
    tri::RefineE<CMeshO, tri::QualityMidPointFunctor<CMeshO>, tri::QualityEdgePredicate<CMeshO> > (underM, slicingfunc, slicingpred, false);

    tri::UpdateSelection<CMeshO>::VertexFromQualityRange(underM,0,std::numeric_limits<float>::max());
    tri::UpdateSelection<CMeshO>::FaceFromVertexStrict(underM);

    tri::UpdateSelection<CMeshO>::FaceInvert(underM);

    CreateSection(sectionPolyline, underM);

    // Calculate the section surface
    tri::CapEdgeMesh(sectionPolyline, sectionSurface);
    tri::UpdateBounding<CMeshO>::Box(sectionSurface);
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(sectionSurface);

    switch(ID(filter)) {
        case FP_CUT_AND_FILL :
        {

        }
        break;
    case FP_CUT_FILL_AND_REMESH:
    {
        CMeshO toProjectCopy;

        //Variables used to remesh the section surface
        float maxVal;
        float targetLen;

        if(par.getBool("TargetLenFromAverageLen"))
        {
            targetLen = par.getAbsPerc("TargetLenAverage");
        }
        else
        {
            targetLen = par.getAbsPerc("TargetLen");
        }

        tri::IsotropicRemeshing<CMeshO>::Params params;

        params.SetTargetLen(targetLen);
        params.SetFeatureAngleDeg(creaseAngle);

        params.iter         = iterations;
        params.adapt        = adaptive;
        params.selectedOnly = selectedOnly;
        params.splitFlag    = splitFlag;
        params.collapseFlag = collapseFlag;
        params.swapFlag     = swapFlag;
        params.smoothFlag   = smoothFlag;
        params.projectFlag  = reprojectFlag;
        params.surfDistCheck= false;

        BoundaryExpand(sectionSurface);

        tri::Append<CMeshO, CMeshO>::Mesh(toProjectCopy, sectionSurface);

        SetMeshRequirements(toProjectCopy);

        toProjectCopy.face.EnableMark();

        try
        {
            tri::IsotropicRemeshing<CMeshO>::Do(sectionSurface, toProjectCopy, params, cb);
        }
        catch(vcg::MissingPreconditionException& excp)
        {
            log(excp.what());
            throw MLException(excp.what());
        }
        tri::UpdateSelection<CMeshO>::FaceInvert(sectionSurface);

        for(auto fi = sectionSurface.face.begin(); fi!=sectionSurface.face.end(); ++fi)
        {
            if((*fi).IsS())
            {
                tri::Allocator<CMeshO>::DeleteFace(sectionSurface, *fi);
            }
        }
        tri::Clean<CMeshO>::RemoveUnreferencedVertex(sectionSurface);

        tri::Allocator<CMeshO>::CompactEveryVector(sectionSurface);

        tri::UpdateBounding<CMeshO>::Box(sectionSurface);
    }
        break;
    default :
        wrongActionCalled(filter);
    }

    //copy the selected faces of underMesh on overMesh
    tri::Append<CMeshO, CMeshO>::Mesh(overM, underM, true);

    //clear the selected faces on underMesh
    tri::UpdateSelection<CMeshO>::VertexClear(underM);
    tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(underM);

    for(auto fi = underM.face.begin(); fi!=underM.face.end(); ++fi)
    {
        if(!(*fi).IsD() && (*fi).IsS())
        {
            tri::Allocator<CMeshO>::DeleteFace(underM, *fi);
        }
    }

    for(auto vi = underM.vert.begin(); vi != underM.vert.end(); ++vi)
    {
        if(!(*vi).IsD() && (*vi).IsS())
        {
            tri::Allocator<CMeshO>::DeleteVertex(underM, *vi);
        }
    }

    for(auto vii = overM.vert.begin(); vii != overM.vert.end(); ++vii)
    {
        if((*vii).IsD())
        {
            vcg::tri::Allocator<CMeshO>::DeleteVertex(overM,*vii);
        }
    }

    vcg::tri::Clean<CMeshO>::RemoveUnreferencedVertex(overM);
    tri::UpdateBounding<CMeshO>::Box(underM);
    if(underM.fn >0)
    {
        tri::UpdateNormal<CMeshO>::PerFaceNormalized(underM);
        tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(underM);
    }

    tri::UpdateBounding<CMeshO>::Box(overM);
     if(overM.fn >0)
    {
        tri::UpdateNormal<CMeshO>::PerFaceNormalized(overM);
        tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(overM);
    }

    FillMesh(underFM, sectionSurface, underM);
    FillMesh(overFM, sectionSurface, overM, true);

    tri::UpdateBounding<CMeshO>::Box(underFM);
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(underFM);

    tri::UpdateBounding<CMeshO>::Box(overFM);
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(overFM);

    tri::UpdateSelection<CMeshO>::Clear(underM);
    tri::UpdateSelection<CMeshO>::Clear(underFM);
    tri::UpdateSelection<CMeshO>::Clear(sectionSurface);
    tri::UpdateSelection<CMeshO>::Clear(overM);
    tri::UpdateSelection<CMeshO>::Clear(overFM);


    if(par.getBool("createSectionSurface"))
    {
        // making up new layer name
        QString sectionName = QFileInfo(base->shortName()).baseName() + "_section_surface";
        MeshModel* sectionSurfaceModel= md.addNewMesh("",sectionName,true);
        sectionSurfaceModel->cm = sectionSurface;
    }

    if(par.getBool("createUnderMesh"))
    {
        // making up new layer name
        QString underMName = QFileInfo(base->shortName()).baseName() + "_under_part";
        MeshModel* underMeshModel= md.addNewMesh("",underMName);
        underMeshModel->cm = underFM;
    }

    if(par.getBool("createOverMesh"))
    {
        // making up new layer name
        QString overMName = QFileInfo(base->shortName()).baseName() + "_under_part";
        MeshModel* overMeshModel= md.addNewMesh("",overMName);
        overMeshModel->cm = overFM;
    }

    return std::map<std::string, QVariant>();
}

 /**
 * @brief FilterSamplePlugin::postCondition
 * @return
 */
int FilterCutAndFillPlugin::postCondition(const QAction *filter) const
{
    switch(ID(filter))
    {
    case FP_CUT_AND_FILL: 	return MeshModel::MM_VERTCOORD | MeshModel::MM_FACENORMAL | MeshModel::MM_VERTNORMAL;
    }

	return MeshModel::MM_VERTCOORD | MeshModel::MM_FACENORMAL | MeshModel::MM_VERTNORMAL;
}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterCutAndFillPlugin)
