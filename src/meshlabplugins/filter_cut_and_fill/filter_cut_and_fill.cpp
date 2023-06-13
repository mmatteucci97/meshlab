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

    parlst.addParam(RichBool("allLayers", false, "Apply to all visible layers", "if selected, the filter will be applied to all visible mesh Layers."));

    switch(ID(action)) {
    case FP_CUT_AND_FILL :
    {
    }
        break;
    case FP_CUT_FILL_AND_REMESH:
    {
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

std::set<std::pair<CMeshO *, const char *>> FilterCutAndFillPlugin::SliceMesh(MeshModel &m, Plane3m slicingPlane, const RichParameterList &par, vcg::CallBackPos *cb, bool remesh)
{
    cout << "Slicing " << m.shortName().toStdString() << ". VN=" << m.cm.VN() << "; FN=" << m.cm.FN() << endl;
    MeshModel* base=&m;
    MeshModel* orig=&m;

    cout << 1 << endl;

    // making up new layer name
    CMeshO sectionPolyline;
    CMeshO *sectionSurface = new CMeshO();
    CMeshO underM;
    CMeshO overM;
    CMeshO *underFM = new CMeshO();
    CMeshO *overFM = new CMeshO();

    cout << 2 << endl;


    SetMeshRequirements(sectionPolyline);
    SetMeshRequirements(*sectionSurface);

    SetMeshRequirements(underM);
    SetMeshRequirements(overM);

    SetMeshRequirements(*underFM);
    SetMeshRequirements(*overFM);

    cout << 3 << endl;

    tri::QualityMidPointFunctor<CMeshO> slicingfunc(0.0);
    tri::QualityEdgePredicate<CMeshO> slicingpred(0.0,0.0);

    cout << 4 << endl;


    // Check if the mesh has the correct topology to perform the algorithm
    // TODO: check only if the mesh has boundary or if is it two-manifold only on the part of the mesh around the plane
    CheckMeshRequirement(base);

    cout << 5 << endl;

    tri::Append<CMeshO,CMeshO>::Mesh(underM,orig->cm);
    tri::UpdateQuality<CMeshO>::VertexFromPlane(underM, slicingPlane);

    cout << 6 << endl;

    tri::UpdateTopology<CMeshO>::FaceFace(underM);
    tri::RefineE<CMeshO, tri::QualityMidPointFunctor<CMeshO>, tri::QualityEdgePredicate<CMeshO> > (underM, slicingfunc, slicingpred, false);

    cout << 7 << endl;

    tri::UpdateSelection<CMeshO>::VertexFromQualityRange(underM,0,std::numeric_limits<float>::max());
    tri::UpdateSelection<CMeshO>::FaceFromVertexStrict(underM);

    tri::UpdateSelection<CMeshO>::FaceInvert(underM);

    CreateSection(sectionPolyline, underM);

    // Calculate the section surface
    tri::CapEdgeMesh(sectionPolyline, *sectionSurface);
    tri::UpdateBounding<CMeshO>::Box(*sectionSurface);
    tri::UpdateNormal<CMeshO>::PerVertexNormalized(*sectionSurface);
    tri::UpdateTopology<CMeshO>::FaceFace(*sectionSurface);

    if(remesh)
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

        cout << "section surface vertices: " << sectionSurface->VN() << endl;
        cout << "section surface faces: " << sectionSurface->FN() << endl;

        BoundaryExpand(*sectionSurface);

        cout << "section surface vertices: " << sectionSurface->VN() << endl;
        cout << "section surface faces: " << sectionSurface->FN() << endl;

        tri::Append<CMeshO, CMeshO>::Mesh(toProjectCopy, *sectionSurface);

        SetMeshRequirements(toProjectCopy);

        toProjectCopy.face.EnableMark();

        try
        {
            tri::IsotropicRemeshing<CMeshO>::Do(*sectionSurface, toProjectCopy, params, cb);
        }
        catch(vcg::MissingPreconditionException& excp)
        {
            log(excp.what());
            throw MLException(excp.what());
        }
        tri::UpdateSelection<CMeshO>::FaceInvert(*sectionSurface);

        for(auto fi = sectionSurface->face.begin(); fi!=sectionSurface->face.end(); ++fi)
        {
            if((*fi).IsS())
            {
                tri::Allocator<CMeshO>::DeleteFace(*sectionSurface, *fi);
            }
        }
        tri::Clean<CMeshO>::RemoveUnreferencedVertex(*sectionSurface);

        tri::Allocator<CMeshO>::CompactEveryVector(*sectionSurface);

        tri::UpdateBounding<CMeshO>::Box(*sectionSurface);
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

    FillMesh(*underFM, *sectionSurface, underM);
    FillMesh(*overFM, *sectionSurface, overM, true);

    tri::UpdateBounding<CMeshO>::Box(*underFM);
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(*underFM);

    tri::UpdateBounding<CMeshO>::Box(*overFM);
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(*overFM);

    tri::UpdateSelection<CMeshO>::Clear(underM);
    tri::UpdateSelection<CMeshO>::Clear(*underFM);
    tri::UpdateSelection<CMeshO>::Clear(*sectionSurface);
    tri::UpdateSelection<CMeshO>::Clear(overM);
    tri::UpdateSelection<CMeshO>::Clear(*overFM);

    std::set<std::pair<CMeshO *, const char *>> retValues;


    if(par.getBool("createSectionSurface"))
    {
        // making up new layer name
        const char * sectionName = "_section_surface";
        retValues.insert(std::make_pair(sectionSurface, sectionName));
    }

    if(par.getBool("createUnderMesh"))
    {
        // making up new layer name
        const char * underMName = "_under_part";
        retValues.insert(std::make_pair(underFM, underMName));
    }

    if(par.getBool("createOverMesh"))
    {
        // making up new layer name
        const char * overMName = "_over_part";
        retValues.insert(std::make_pair(overFM, overMName));
    }

    for(auto mm : retValues)
    {

        cout << "Sliced " << mm.second << " vertices: " << mm.first->VN() << endl;
        cout << "Sliced " << mm.second << " faces: " << mm.first->FN() << endl;

    }

    return retValues;
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
    //CALCULATE THE PLANE USING THE SELECTED LAYER
    MeshModel & m = *md.mm();
    Box3m bbox=m.cm.bbox;
    // the mesh has to return to its original position
    if (m.cm.Tr != Matrix44m::Identity())
        tri::UpdatePosition<CMeshO>::Matrix(m.cm, Inverse(m.cm.Tr), true);

    // the mesh has to be correctly transformed
    if (m.cm.Tr != Matrix44m::Identity())
        tri::UpdatePosition<CMeshO>::Matrix(m.cm, m.cm.Tr, true);

    SetMeshRequirements(m);

    Point3m planeAxis(0,0,0);
    Scalarm planeOffset;
    Point3m planeCenter;
    Plane3m slicingPlane;
    int ind;

    ind = par.getEnum("planeAxis");
    if(ind>=0 && ind<3)
        planeAxis[ind] = 1.0f;
    else
        planeAxis=par.getPoint3m("customAxis");

    planeAxis.Normalize();
    planeOffset = par.getFloat("planeOffset");

    switch(RefPlane(par.getEnum("relativeTo")))
    {
        case REF_CENTER:  planeCenter = bbox.Center()+ planeAxis*planeOffset*(bbox.Diag()/2.0);      break;
        case REF_MIN:     planeCenter = bbox.min+planeAxis*planeOffset*(bbox.Diag()/2.0);    break;
        case REF_ORIG:    planeCenter = planeAxis*planeOffset;  break;
    }
    slicingPlane.Init(planeCenter,planeAxis);

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
        newMeshes = SliceMesh(m, slicingPlane, par, cb, remesh);
    }
    else
    {
        for(auto mmp = md.meshBegin(); mmp!=md.meshEnd(); ++mmp)
        {
            if(mmp->isVisible())
            {
                std::set<std::pair<CMeshO *, const char *>> ret = SliceMesh(*mmp, slicingPlane, par, cb, remesh);
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

            cout << "created mesh " << mM->shortName().toStdString() << ". VN=" << mM->cm.VN() << "; FN=" << mM->cm.FN() << endl;
        }
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
