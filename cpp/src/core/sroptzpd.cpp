/************************************************************************//**
 * File: sroptzp.cpp
 * Description: Optical element: Zone Plate
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptzpd.h"
#include "sroptdrf.h"

#include <cassert>
#include <algorithm>

void do_efield_2d_sample(float *peout, float z, float x, float *pe, int ne, int nz, int nx, float z0, float dz, float x0, float dx,
	bool additive = false)
{
	// peout has len == ne*2
	// pe has dimension (nz, nx, ne*2), C-mem layout
	// x0, z0 are the starting coordinate. dx, dz are step size
	// z, x are sampling point
	int ixL = floor((x - x0) / dx); // left index in pe
	int izL = floor((z - z0) / dz); // left index
	long long SZ = nx * ne * 2;
	// Bilinear interpolation, see https://en.wikipedia.org/wiki/Bilinear_interpolation
	float x1 = x0 + ixL * dx, x2 = x1 + dx;
	float z1 = z0 + izL * dz, z2 = z1 + dz;
	if (!additive) memset(peout, 0, ne * 2 * sizeof(float));
	for (int ie = 0; ie < ne*2; ++ie) {
		// treat as 0 if out of boundary
		if (ixL < 0 || ixL + 1 >= nx || izL < 0 || izL + 1 >= nz) {
			peout[ie] += 0.0; continue;
		}
		float f11 = pe[izL * SZ + ixL * ne * 2 + ie];
		float f12 = pe[izL * SZ + (ixL+1) * ne * 2 + ie];
		float f21 = pe[(izL+1) * SZ + ixL * ne * 2 + ie];
		float f22 = pe[(izL+1) * SZ + (ixL + 1) * ne * 2 + ie];

		peout[ie] += (f11 * (z2 - z) * (x2 - x) + f21 * (z - z1) * (x2 - x) + f12 * (z2 - z) * (x - x1) + f22*(z - z1) * (x - x1)) / (z2 - z1) / (x2 - x1);
	}
}

void accumulate_rad(float* rdestx, float *rdestz, int nz, int nx, float x0, float dx, float z0, float dz, const srTSRWRadStructAccessData& r /*, const string& component */)
{
	// float* psrc = nullptr;
	// if (component == "X") psrc = r.pBaseRadX;
	// for each point rdest, we sample a result from srTSRWRadStructAccessData r and accumulate
	for (int iz = 0; iz < nz; ++iz) {
		float z = z0 + iz * dz;
		for (int ix = 0; ix < nx; ++ix) {
			float x = x0 + ix * dx;
			do_efield_2d_sample(rdestx, z, x, r.pBaseRadX, r.ne, r.nz, r.nx, r.zStart, r.zStep, r.xStart, r.xStep, true /* additive */);
			do_efield_2d_sample(rdestz, z, x, r.pBaseRadZ, r.ne, r.nz, r.nx, r.zStart, r.zStep, r.xStart, r.xStep, true /* additive */);
			rdestx += r.ne * 2;
			rdestz += r.ne * 2;
		}
		//
	}
}

//*************************************************************************
int srTZonePlateD::PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!!] The ZP with aperture and drift L= " << dftLen << endl;

	float xStep = pRadAccessData->xStep;
	float xStart = pRadAccessData->xStart;
	float zStep = pRadAccessData->zStep;
	float zStart = pRadAccessData->zStart;
	int nx = pRadAccessData->nx;
	int nz = pRadAccessData->nz;

	const int NXDIV = 3, NZDIV = 3; // we divide x and z into nx, nz subblocks
	double xwidth = (nx - 1) * xStep / NXDIV;
	double zwidth = (nz - 1) * zStep / NZDIV;
	//double pdedge = 10, pdcenter = 0.1; // the outter edge side becomes denser (x2)
	double pdedge = 5, pdcenter = 1;

	//fprintf(stderr, "center: xc= %g start+n/2*L= %g diff= %g\n", pRadAccessData->xc,
	//	(xStart + (nx - 1) / 2.0 * xStep), std::fabs(pRadAccessData->xc - (xStart + (nx - 1) / 2.0 * xStep)));
	//fprintf(stderr, "center: zc= %g start+n/2*L= %g diff= %g\n", pRadAccessData->zc,
	//	(zStart + (nz - 1) / 2.0 * zStep), std::fabs(pRadAccessData->zc - (zStart + (nz - 1) / 2.0 * zStep)));

	const long long RADSZ = pRadAccessData->nz * pRadAccessData->nx * (pRadAccessData->ne * 2);
	float * radx = (float*) calloc(RADSZ, sizeof(float));
	float * radz = (float*) calloc(RADSZ, sizeof(float));

	srTDriftSpace internal_drift(dftLen);

	srTRadResize shifter;
	shifter.ShiftTypeBeforeRes |= 3;
	shifter.pxm = 1.0 / NXDIV;
	shifter.pzm = 1.0 / NZDIV;
	for (int iz = 0; iz < NZDIV; ++iz) {
		shifter.zCenShift = -(pRadAccessData->zc - zStart - zwidth * (iz + 0.5));
		for (int ix = 0; ix < NXDIV; ++ix) {
			shifter.xCenShift = -(pRadAccessData->xc - xStart - xwidth * (ix + 0.5)); 
			double r = std::hypot(shifter.xCenShift, shifter.zCenShift);
			double Rfull = std::hypot(pRadAccessData->xStep * nx / 2.0, pRadAccessData->zStep * nz / 2.0);
			double pd = r / Rfull * (pdedge - pdcenter) + pdcenter;
			shifter.pxd = shifter.pzd = pd;

			srTSRWRadStructAccessData newRad(pRadAccessData);
			RadResizeGen(newRad, shifter);
			if (int result = srTZonePlate::PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}
			if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}

			accumulate_rad(radx, radz, nz, nx, xStart, xStep, zStart, zStep, newRad);
		}
	}
	
	memcpy(pRadAccessData->pBaseRadX, radx, RADSZ * sizeof(float));
	memcpy(pRadAccessData->pBaseRadZ, radz, RADSZ * sizeof(float));

	free(radx);
	free(radz);

	return 0;
}
//*************************************************************************
