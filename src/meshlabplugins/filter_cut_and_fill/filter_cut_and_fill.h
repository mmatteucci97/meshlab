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

/****************************************************************************
History
$Log: sampleplugins.h,v $

Revision 1,3 2020/05/20
Reorganization of the filter, comments in doxygen format

Revision 1.2  2006/11/29 00:59:21  cignoni
Cleaned plugins interface; changed useless help class into a plain string

Revision 1.1  2006/09/25 09:24:39  e_cerisoli
add sampleplugins

****************************************************************************/

#ifndef FILTER_CUT_AND_FILL_PLUGIN_H
#define FILTER_CUT_AND_FILL_PLUGIN_H

#include <common/plugins/interfaces/filter_plugin.h>

class FilterCutAndFillPlugin : public QObject, public FilterPlugin
{
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
    Q_INTERFACES(FilterPlugin)

    enum RefPlane {REF_CENTER, REF_MIN, REF_ORIG};

public:
    enum
    {
        FP_CUT_AND_FILL,
        FP_CUT_FILL_AND_REMESH
    } ;

    FilterCutAndFillPlugin();
    virtual ~FilterCutAndFillPlugin() {}

	QString pluginName() const;
    QString pythonFilterName(ActionIDType f) const;
    QString filterName(ActionIDType filter) const;
    QString filterInfo(ActionIDType filter) const;

    FilterClass getClass(const QAction* a) const;
    RichParameterList initParameterList(const QAction*, const MeshModel &/*m*/);
    std::map<std::string, QVariant> applyFilter(
            const QAction* action,
            const RichParameterList & parameters,
            MeshDocument &md,
            unsigned int& postConditionMask,
            vcg::CallBackPos * cb);

    int postCondition(const QAction *filter) const;
    int getPreConditions(const QAction *filter) const;
    int getRequirements(const QAction *filter) const;

    static float DetermineAverageEdgeLength(const MeshModel &mesh);

    FilterArity filterArity(const QAction*) const {return SINGLE_MESH;}

protected:

    int iterations = 10;
    bool adaptive= true;
    bool selectedOnly = true;
    float creaseAngle = 30;

    bool splitFlag = true;
    bool collapseFlag = true;
    bool swapFlag = true;
    bool smoothFlag = true;
    bool reprojectFlag = true;

    Scalarm lastisor_FeatureDeg;
};

#endif
