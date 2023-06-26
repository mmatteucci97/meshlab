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
#include "cut_and_fill.h"

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

typedef vcg::tri::QualityEdgePredicate<CMeshO> EdgePredicate;
typedef vcg::tri::QualityMidPointFunctor<CMeshO> MidPointFunctor;
typedef vcg::tri::IsotropicRemeshing<CMeshO>::Params RemeshingParams;

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

RichParameterList FilterCutAndFillPlugin::initParameterList(const QAction *action,const MeshModel &m)
{
    RichParameterList parlst;
    float maxVal;

    QStringList axis = QStringList() <<"X Axis"<<"Y Axis"<<"Z Axis"<<"Custom Axis";
    parlst.addParam(RichEnum   ("planeAxis", 0, axis, tr("Plane perpendicular to"), tr("The Slicing plane will be done perpendicular to the axis")));
    parlst.addParam(RichDirection("customAxis",Point3f(0,1,0),"Custom axis","Specify a custom axis, this is only valid if the above parameter is set to Custom"));

    parlst.addParam(RichFloat  ("planeOffset", 0.0, "Cross plane offset", "Specify an offset of the cross-plane. The offset corresponds to the distance from the point specified in the plane reference parameter. By default (Cross plane offset == 0)"));
    parlst.addParam(RichEnum   ("relativeTo",2,QStringList()<<"Bounding box center"<<"Bounding box min"<<"Origin","plane reference","Specify the reference from which the planes are shifted"));

    parlst.addParam(RichBool("createSectionSurface",false,"Create also section surface","If selected, in addition to a layer with the section polyline, it will be created also a layer with a triangulated version of the section polyline. This only works if the section polyline is closed"));
    parlst.addParam(RichBool("createOverMesh", true, "Create also over mesh", "If selected, it will create another layer with the portion of the mesh over the section plane. It requires manifoldness of the mesh "));
    parlst.addParam(RichBool("createUnderMesh", true, "Create also under mesh", "If selected, it will create another layer with the portion of the mesh under the section plane. It requires manifoldness of the mesh "));

    parlst.addParam(RichBool("allLayers", false, "Apply to all visible layers", "if selected, the filter will be applied to all visible mesh Layers."));

    switch(ID(action)) {
    case FP_CUT_AND_FILL :
    {
    }
        break;
    case FP_CUT_FILL_AND_REMESH:
    {
        float averageLen;
        averageLen = DetermineAverageEdgeLength(m);
        parlst.addParam(RichPercentage("TargetLenAverage", averageLen, 0.01*averageLen, 100*averageLen, "Target Edge Len Average", "Use this value of the average edge length"));
    }
        break;
    default :
        assert(0);
    }
    return parlst;
}

int FilterCutAndFillPlugin::postCondition(const QAction *filter) const
{
    switch(ID(filter))
    {
    case FP_CUT_AND_FILL: 	return MeshModel::MM_VERTCOORD | MeshModel::MM_FACENORMAL | MeshModel::MM_VERTNORMAL;
    }

    return MeshModel::MM_VERTCOORD | MeshModel::MM_FACENORMAL | MeshModel::MM_VERTNORMAL;
}

std::set<std::pair<CMeshO *, const char *>> FilterCutAndFillPlugin::SliceMesh(MeshModel &model, Point3m planeNormal, Point3m planeCenter, const RichParameterList &par, vcg::CallBackPos *cb, bool remesh)
{
    //------- parameters
    CMeshO m = model.cm;

    std::set<std::pair<CMeshO *, const char *>> retValues;

    CMeshO *sectionSurface = new CMeshO();
    CMeshO *underMesh = new CMeshO();
    CMeshO *overMesh = new CMeshO();
    CMeshO section;

    MidPointFunctor slicingFunc(0.0);
    EdgePredicate slicingPred(0.0, 0.0);

    Plane3m slicingPlane;

    //-----------------requirements
    SetMeshRequirements(model);

    SetMeshRequirements(m);
    SetMeshRequirements(*sectionSurface);
    SetMeshRequirements(*underMesh);
    SetMeshRequirements(*overMesh);
    SetMeshRequirements(section);

    //--------------


    planeNormal.Normalize();

    slicingPlane.Init(planeCenter, planeNormal);

    tri::UpdateQuality<CMeshO>::VertexFromPlane(m, slicingPlane);
    tri::UpdateTopology<CMeshO>::FaceFace(m);
    tri::RefineE<CMeshO, MidPointFunctor, EdgePredicate>(m, slicingFunc, slicingPred);

    CheckMeshRequirement(m);

    CreateSection(section, m);

    CapEdgeMesh(section, *sectionSurface);

    UpdateNormal<CMeshO>::PerFaceNormalized(*sectionSurface);
    UpdateNormal<CMeshO>::PerVertexAngleWeighted(*sectionSurface);

    Point3f sectionAverageNormal = Point3f(0.0, 0.0, 0.0);
    for(auto fi = sectionSurface->face.begin(); fi != sectionSurface->face.end(); ++fi)
    {
        if(fi->IsD()) continue;
        sectionAverageNormal += fi->N();
    }
    sectionAverageNormal /= sectionSurface->FN();

    if(sectionAverageNormal.dot(planeNormal) < 0.0)
    {
        tri::Clean<CMeshO>::FlipMesh(*sectionSurface);
    }

    if(remesh)
    {
        tri::IsotropicRemeshing<CMeshO>::Params remeshParams;
        float averageEdgeLength;
        CMeshO toProjectCopy;

        SetMeshRequirements(toProjectCopy);

        averageEdgeLength = par.getAbsPerc("TargetLenFromAverageLen");

        remeshParams.SetTargetLen(averageEdgeLength);
        remeshParams.SetFeatureAngleDeg(creaseAngle);
        remeshParams.selectedOnly = true;
        remeshParams.iter = iterations;
        remeshParams.surfDistCheck = false;
        remeshParams.adapt = true;
        remeshParams.splitFlag = splitFlag;
        remeshParams.smoothFlag = smoothFlag;
        remeshParams.collapseFlag = collapseFlag;
        remeshParams.swapFlag = swapFlag;
        remeshParams.projectFlag = reprojectFlag;

        BoundaryExpand(*sectionSurface);
        tri::Append<CMeshO, CMeshO>::Mesh(toProjectCopy, *sectionSurface);

        try
        {
            tri::IsotropicRemeshing<CMeshO>::Do(*sectionSurface, toProjectCopy, remeshParams);
        }
        catch(vcg::MissingPreconditionException& e)
        {
            cout << e.what() << endl;
        }
        tri::UpdateSelection<CMeshO>::FaceInvert(*sectionSurface);

        for(auto fi = sectionSurface->face.begin(); fi != sectionSurface->face.end(); ++fi)
        {
            if(fi->IsD()) continue;
            if(fi->IsS())
            {
                tri::Allocator<CMeshO>::DeleteFace(*sectionSurface, *fi);
            }
        }
        tri::Clean<CMeshO>::RemoveUnreferencedVertex(*sectionSurface);
        tri::Allocator<CMeshO>::CompactEveryVector(*sectionSurface);
        tri::UpdateBounding<CMeshO>::Box(*sectionSurface);


    }


    tri::UpdateSelection<CMeshO>::VertexFromQualityRange(m, 0.0, std::numeric_limits<float>::max());
    tri::UpdateSelection<CMeshO>::FaceFromVertexStrict(m);

    tri::Append<CMeshO, CMeshO>::Mesh(*underMesh, m, true);
    tri::UpdateSelection<CMeshO>::FaceInvert(m);
    tri::Append<CMeshO, CMeshO>::Mesh(*overMesh, m, true);

    FillMesh(*overMesh, *underMesh, *sectionSurface, m.bbox.Center());

    tri::UpdateSelection<CMeshO>::Clear(*underMesh);
    tri::UpdateSelection<CMeshO>::Clear(*overMesh);
    tri::UpdateSelection<CMeshO>::Clear(*sectionSurface);



    if(par.getBool("createSectionSurface"))
    {
        const char * sectionName = "_section_surface";
        tri::UpdateColor<CMeshO>::PerVertexQualityRamp(*sectionSurface);
        retValues.insert(std::make_pair(sectionSurface, sectionName));
    }

    if(par.getBool("createUnderMesh"))
    {
        const char * underMName = "_under_part";
        retValues.insert(std::make_pair(underMesh, underMName));
    }

    if(par.getBool("createOverMesh"))
    {
        const char * overMName = "_over_part";
        retValues.insert(std::make_pair(overMesh, overMName));
    }
    return retValues;
}


std::map<std::string, QVariant> FilterCutAndFillPlugin::applyFilter(const QAction * filter, const RichParameterList & par, MeshDocument &md, unsigned int& /*postConditionMask*/, vcg::CallBackPos *cb)
{
    std::map<std::string, QVariant> outputValues;

    MeshModel & m = *md.mm();

    if (m.cm.Tr != Matrix44m::Identity())
        tri::UpdatePosition<CMeshO>::Matrix(m.cm, m.cm.Tr, true);

    SetMeshRequirements(m);

    Box3m bbox=m.cm.bbox;

    Point3m planeNormal(0,0,0);
    Scalarm planeOffset;
    Point3m planeCenter;
    int ind;

    ind = par.getEnum("planeAxis");
    if(ind>=0 && ind<3)
        planeNormal[ind] = 1.0f;
    else
        planeNormal=par.getPoint3m("customAxis");

    planeNormal.Normalize();

    planeOffset = par.getFloat("planeOffset");

    switch(RefPlane(par.getEnum("relativeTo")))
    {
        case REF_CENTER:
            planeCenter = bbox.Center()+ planeNormal*planeOffset;
            break;
        case REF_MIN:
            planeCenter = bbox.min + planeNormal*planeOffset;
            break;
        case REF_ORIG:
        planeCenter = planeNormal*planeOffset;
        break;
    }
    bool applyToAllVisibleLayers = par.getBool("allLayers");
    int cnt=0;

    bool remesh=false;
    if(ID(filter) == FP_CUT_FILL_AND_REMESH)
    {
        remesh=true;
    }

    std::set<std::pair<CMeshO *, const char *>> newMeshes;

    if(!applyToAllVisibleLayers)
    {
        newMeshes = SliceMesh(m, planeNormal, planeCenter, par, cb, remesh);
    }
    else
    {
        for(auto mmp = md.meshBegin(); mmp!=md.meshEnd(); ++mmp)
        {
            if(mmp->isVisible())
            {
                std::set<std::pair<CMeshO *, const char *>> ret = SliceMesh(*mmp, planeNormal, planeCenter, par, cb, remesh);
                newMeshes.insert(ret.begin(), ret.end());
            }
        }
    }

    if(!newMeshes.empty())
    {
        for(auto mm : newMeshes)
        {
            QString name = QFileInfo(m.shortName()).baseName() + mm.second;
            MeshModel * mM = md.addNewMesh("", name);
            mM->cm = *mm.first;
            mM->updateBoxAndNormals();
        }
    }

    return std::map<std::string, QVariant>();
}


MESHLAB_PLUGIN_NAME_EXPORTER(FilterCutAndFillPlugin)
