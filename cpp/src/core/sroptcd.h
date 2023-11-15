/************************************************************************//**
 * File: sroptcd.h
 * Description: Optical element: Element Combined with Drift Space (header)
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

#include "sroptelm.h"
#include "sroptshp.h"
#include <utility>
#include <cassert>

//*************************************************************************

class srTCombinedDrift : public srTGenOptElem {
	double dftLen;
	int nxdiv, nzdiv;
	double xdivs[32], zdivs[32]; // maximum 32 divides in each dimension
	double crsz[5 * 32 * 32]; // each cell has (method, px_range, px_density, pz_range, pz_density)
	double obsgrid[4]; // x half-range, x stepsize, then z

	srTShapedOptElem* elem;

public:
	srTCombinedDrift(srTStringVect* pElemInfo) {}
	srTCombinedDrift(srTShapedOptElem* _elem, double _dftLen = 0,
		int _nxdiv = 0, double *_xdivs = nullptr,
		int _nzdiv = 0, double *_zdivs = nullptr, double *_crsz=nullptr, double*_obsgrid=nullptr) {
		elem = _elem;
		assert(elem);
		dftLen = _dftLen;
		nxdiv = _nxdiv;
		memcpy(xdivs, _xdivs, nxdiv * sizeof (double));
		nzdiv = _nzdiv;
		memcpy(zdivs, _zdivs, nzdiv * sizeof(double));
		memcpy(crsz, _crsz, nxdiv * nzdiv * 5 * sizeof(double));
		memcpy(obsgrid, _obsgrid, 4 * sizeof(double));
	}

  int PropagateRad2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);

	virtual int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect);
	
private:
	double max_div_width() const;
	// div 1, 2, ..., 2*ndiv-1. find its coordinate range
	void fill_div_range(vector<int>& idx, int nx) const;
	// pair<double, double> get_div_range(int ix, double xStart, double xStep, int nx) const;

	double maxpd() const;
	//void resize_dest_rad(srTSRWRadStructAccessData& rad); // const;
	void init_dest_rad(srTSRWRadStructAccessData& rad, const srTSRWRadStructAccessData* pRadAccessData, int wg = 0) const;
	void init_dest_rad2(srTSRWRadStructAccessData& rad, const srTSRWRadStructAccessData* pRadAccessData) const;
};

// select subsection, accumulate and pad a srTSRWRadStructAccessData
struct CDRadStructHelper {
public:
	static void add(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src);
	static void assign(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src);
  static void sample(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src,
	  double xStart, double xEnd, double xStep, int nxpad,
	  double zStart, double zEnd, double zStep, int nzpad);
};

//*************************************************************************

#endif

