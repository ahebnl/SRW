/************************************************************************//**
 * File: sroptzp.h
 * Description: Optical element: Zone Plate (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTZPD_H
#define __SROPTZPD_H

#include "sroptzp.h"

//*************************************************************************

class srTZonePlateD : public srTZonePlate {
	double dftLen, pdcenter, pdedge;
	int nxdiv, nzdiv;

public:
	srTZonePlateD(srTStringVect* pElemInfo) : srTZonePlate(pElemInfo) {}
	srTZonePlateD(int _nZones, double _rn, double _thick, double _atLen1, double _atLen2, double _delta1, double _delta2, double _x = 0, double _y = 0, double _e0 = 0, double _dftLen = 0, int _nxdiv = 1, int _nzdiv = 1, double _pdcenter = 1, double _pdedge = 1)
		:srTZonePlate(_nZones, _rn, _thick, _atLen1, _atLen2, _delta1, _delta2, _x, _y, _e0) {
		dftLen = _dftLen;
		nxdiv = _nxdiv;
		nzdiv = _nzdiv;
		pdcenter = _pdcenter;
		pdedge = _pdedge;
	}


	virtual int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);
	
};

//*************************************************************************

#endif

