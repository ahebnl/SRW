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
#include <utility>

//*************************************************************************

class srTConnectDrift : public srTZonePlate {
	double dftLen;
	int nxdiv, nzdiv;
	double rdivs[32]; // 16 layers

public:
	srTConnectDrift(srTStringVect* pElemInfo) : srTZonePlate(pElemInfo) {}
	srTConnectDrift(int _nZones, double _rn, double _thick, double _atLen1, double _atLen2, double _delta1, double _delta2, double _x = 0, double _y = 0, double _e0 = 0, double _dftLen = 0,
		int _nxdiv = 0, int _nzdiv = 0, double *_rdivs = nullptr)
		:srTZonePlate(_nZones, _rn, _thick, _atLen1, _atLen2, _delta1, _delta2, _x, _y, _e0) {
		dftLen = _dftLen;
		nxdiv = _nxdiv;
		nzdiv = _nzdiv;
		for (int i = 0; i < 32; ++i) rdivs[i] = _rdivs[i];
	}

	int PropagateRad1(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);
	int PropagateRad2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);

	int test_kick(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);

	virtual int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);
	
private:
	double max_div_width() const;
	// div 1, 2, ..., 2*ndiv-1. find its coordinate range
	void fill_div_range(vector<int>& idx, int nx) const;
	// pair<double, double> get_div_range(int ix, double xStart, double xStep, int nx) const;

	double maxpd() const;
	//void resize_dest_rad(srTSRWRadStructAccessData& rad); // const;
	void init_dest_rad(srTSRWRadStructAccessData& rad, vector<int>& xidxdiv, vector<int>& zidxdiv) const;
};

// select subsection, accumulate and pad a srTSRWRadStructAccessData
struct CDRadStructHelper {
public:
	static void add(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src);
	static void assign(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src);
};

//*************************************************************************

#endif

