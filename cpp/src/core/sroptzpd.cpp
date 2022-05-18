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
#include "gmfft.h"
#include "sroptang.h"
#include <cassert>
#include <algorithm>
#include <fstream>

using namespace std;

#define DEBUG_ZPD 3

template<class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
	return (v < lo ? lo : (v > hi ? hi : v));
}

static long next_fft_num(long n) {
	CGenMathFFT::NextCorrectNumberForFFT(n); return n;
}

void do_efield_2d_sample(float *peout, double z, double x, float *pe, int ne, int nz, int nx, double z0, double dz, double x0, double dx,
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
	double x1 = x0 + ixL * dx, x2 = x1 + dx;
	double z1 = z0 + izL * dz, z2 = z1 + dz;
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

void accumulate_rad(float* rdestx, float *rdestz, int nz, int nx, double x0, double dx, double z0, double dz, const srTSRWRadStructAccessData& r /*, const string& component */,
	double zshift = 0.0, double xshift = 0.0)
{
	// float* psrc = nullptr;
	// if (component == "X") psrc = r.pBaseRadX;
	// for each point (z, x) rdest, we sample a result from srTSRWRadStructAccessData r and accumulate
	// we are sampling (z+zshift, x+xshift) at r
	for (int iz = 0; iz < nz; ++iz) {
		double z = z0 + iz * dz;
		for (int ix = 0; ix < nx; ++ix) {
			double x = x0 + ix * dx;
			do_efield_2d_sample(rdestx, z+zshift, x+xshift, r.pBaseRadX, r.ne, r.nz, r.nx, r.zStart, r.zStep, r.xStart, r.xStep, true /* additive */);
			do_efield_2d_sample(rdestz, z+zshift, x+xshift, r.pBaseRadZ, r.ne, r.nz, r.nx, r.zStart, r.zStep, r.xStart, r.xStep, true /* additive */);
			rdestx += r.ne * 2;
			rdestz += r.ne * 2;
		}
		//
	}
}

void select_cell(srTSRWRadStructAccessData& wfr, double z0, double z1, double x0, double x1)
{
	int nreset = 0, NTOT = wfr.nz * wfr.nx;
	for (int iz = 0; iz < wfr.nz; ++iz) {
		double z = wfr.zStart + iz * wfr.zStep;
		for (int ix = 0; ix < wfr.nx; ++ix) {
			double x = wfr.xStart + ix * wfr.xStep;
			if (z < z0 || z >= z1 || x < x0 || x >= x1) {
				int offset = iz * wfr.nx * wfr.ne * 2 + ix * wfr.ne * 2;
				memset(wfr.pBaseRadX + offset, 0.0, wfr.ne * 2 * sizeof(float));
				memset(wfr.pBaseRadZ + offset, 0.0, wfr.ne * 2 * sizeof(float));
				++nreset;
			}
		}
	}
#if DEBUG_ZPD > 2
	fprintf(stderr, "reset %d (%.2f%%) points out side of z [%f %f], x [%f %f] N = %d %d start %f %f step= %g %g\n",
		nreset, nreset * 100.0 / NTOT, z0, z1, x0, x1, wfr.nz, wfr.nx, wfr.zStart, wfr.xStart, wfr.zStep, wfr.xStep);
#endif
}

// psrc - full (nz, nx, ne) 3D array, e has Ex, Ez
// pdst - destination 
void copy_sub_cell(float* pdst, const float* psrc, int nz, int nx, int ne, int iz0, int iz1, int ix0, int ix1)
{
#if DEBUG_ZPD > 2
	{
		const auto itm = std::minmax_element(psrc, psrc + nz*nx*ne*2);
		fprintf(stderr, "src N (%d,%d,%d,2) %ld min= %g max %g\n", nz, nx, ne, nz*nx*ne*2, *itm.first, *itm.second);
	}
#endif

	iz0 = clamp(iz0, 0, nz);
	ix0 = clamp(ix0, 0, nx);
	iz1 = clamp(iz1, 0, nz);
	ix1 = clamp(ix1, 0, nx);
	const long long SZZFULL = nx * ne * 2;
	const long long SZ = (ix1 - ix0) * ne * 2;
	for (int iz = iz0; iz < iz1; ++iz) {
		memcpy(pdst, psrc+iz*SZZFULL+ix0*ne*2, SZ * (sizeof psrc[0]));
		pdst += SZ;
	}

#if DEBUG_ZPD > 2
	{
		const auto itm = std::minmax_element(pdst - SZ * (iz1 - iz0), pdst);
		fprintf(stderr, "dst N (%d,%d,%d,2) %ld min= %g max %g\n", iz1-iz0, ix1-ix0, ne, SZ * (iz1 - iz0), *itm.first, *itm.second);
	}
#endif
}

// psrc - full (nz, nx, ne) 3D array, e has Ex, Ez
// pdst - destination, is a bigger grid nz1, nx1, ne1, place the src "cubic" inside dst and starts from (jz, jx, je)
// in this way, pdst has padding 
void copy_sub_cell_gen(float* pdst, int nz1, int nx1, int jz, int jx,
	const float* psrc, int nz, int nx, int ne, int iz0, int iz1, int ix0, int ix1)
{
#if DEBUG_ZPD > 2
	{
		fprintf(stderr, "copy subset (%d %d) of (%d %d %d,2) to (%d %d ne 2) left padding z= %d x= %d\n",
			iz1 - iz0, ix1 - ix0, nz, nx, ne, nz1, nx1, jz, jx);
		const auto itm = std::minmax_element(psrc, psrc + nz * nx * ne * 2);
		fprintf(stderr, "src N (%d,%d,%d,2) %ld jz= %d jx= %d min= %g max %g\n", nz, nx, ne, nz * nx * ne * 2, jz, jx, *itm.first, *itm.second);
	}
#endif

	iz0 = clamp(iz0, 0, nz);
	ix0 = clamp(ix0, 0, nx);
	iz1 = clamp(iz1, 0, nz);
	ix1 = clamp(ix1, 0, nx);
	// full dimension (nz, nx, ne*2) -> (nz1, nx1, ne*2) 
	// data block ([iz0,iz1), [ix0,ix1), ...) -> ([jz,jz+iz1-iz0), [jx,jx+ix1-ix0), ...)
	// each z, block length=(ix1-ix0)*ne*2, from psrc + iz*(nx*ne*2) + ix0*(ne*2) 
	// destination: pdst + (iz-iz0+jz)*(nx1*ne*2) + (ix-ix0+jx)*ne*2
	const int ne1 = ne;
	const long long DZSRC = nx * ne * 2;
	const long long DXSRC = ne * 2;
	const long long ZSZSRC = (ix1 - ix0) * ne * 2; // each z copy data
	const long long DZDST = nx1 * ne * 2; // dest step size
	const long long DXDST = ne * 2;
	float *pdstz = pdst + jz * DZDST + jx * ne * 2;
	for (int iz = iz0; iz < iz1; ++iz) {
		memcpy(pdstz, psrc + iz * DZSRC + ix0 * (ne * 2), ZSZSRC * (sizeof psrc[0]));
		pdstz += DZDST;
	}
	// loop method
	/*
	for (int iz = iz0; iz < iz1; ++iz) {
		for (int ix = ix0; ix < ix1; ++ix) {
			for (int ie = 0; ie < ne * 2; ++ie) {
				pdst[(jz + iz - iz0) * DZDST + (jx + ix - ix0) * DXDST + ie] = psrc[iz * DZSRC + ix * DXSRC + ie];
			}
		}
	}
	*/
#if DEBUG_ZPD > 2
	{
		const auto itm = std::minmax_element(pdst, pdst + nz1 * nx1 * ne * 2);
		fprintf(stderr, "dst N (%d,%d,%d) %ld min= %g max %g\n", iz1 - iz0, ix1 - ix0, ne, DZDST * (iz1 - iz0+jz), *itm.first, *itm.second);
	}
#endif
}

void sel_sub_cell(srTSRWRadStructAccessData * dst, /* const */ srTSRWRadStructAccessData& wfr, int iz0, int iz1, int ix0, int ix1,
	int npadz0=-1, int npadz1 = 0, int npadx0 = -1, int npadx1 = 0)
{
	// const int npad = 4;
	// npad.. < 0 means does not count, use whatever left
	dst->nz = iz1 - iz0 + max(0, npadz0) + max(0, npadz1);
	dst->nx = ix1 - ix0 + max(0, npadx0) + max(0, npadx1);
	CGenMathFFT::NextCorrectNumberForFFT(dst->nz);
	CGenMathFFT::NextCorrectNumberForFFT(dst->nx);
	fprintf(stderr, "select sub cell iz= [%d, %d) or %d (next good fft num) tot src nz= %d\n", iz0, iz1, dst->nz, wfr.nz);
	fprintf(stderr, "select sub cell ix= [%d, %d) or %d (tot src nx= %d)\n", ix0, ix1, dst->nx, wfr.nx);

	long nz0 = iz1 - iz0;
	long nx0 = ix1 - ix0;
	int jz = npadz0 > 0 ? npadz0 : dst->nz - nz0 - max(0, npadz1); // actual padding
	int jx = npadx0 > 0 ? npadx0 : dst->nx - nx0 - max(0, npadx1); // actual padding
	copy_sub_cell_gen(dst->pBaseRadX, dst->nz, dst->nx, jz, jx, wfr.pBaseRadX, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);
	copy_sub_cell_gen(dst->pBaseRadZ, dst->nz, dst->nx, jz, jx, wfr.pBaseRadZ, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);

	dst->ne = wfr.ne;
	
	dst->xStep = wfr.xStep;
	dst->zStep = wfr.zStep; 
	dst->xStart = wfr.xStart + ix0 * wfr.xStep - jx*dst->xStep;
	dst->zStart = wfr.zStart + iz0 * wfr.zStep - jz*dst->zStep;
	
	dst->xWfrMin = dst->xStart;
	dst->zWfrMin = dst->zStart;
	dst->xWfrMax = dst->xStart + dst->nx * dst->xStep;
	dst->zWfrMax = dst->zStart + dst->nz * dst->zStep;

#if DEBUG_ZPD > 2
	fprintf(stderr, "src hash= 0x%016zx\n", wfr.hashcode());
	fprintf(stderr, "dst hash= 0x%016zx\n", dst->hashcode());
#endif

	// dst->ComputeRadMoments();

#if DEBUG_ZPD > 2
	fprintf(stderr, "dst hash= 0x%016zx\n", dst->hashcode());
#endif
}


int srTZonePlateD::PropagateRad1(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!!] The ZP (nzdiv=" << nzdiv << ",nxdiv=" << nxdiv
		<< ") with aperture and drift L= " << dftLen << endl;
#if DEBUG_ZPD > 1
	cout << "zp input hash=" << std::hex << pRadAccessData->hashcode() << endl;
#endif

	double xStep = pRadAccessData->xStep;
	double xStart = pRadAccessData->xStart;
	double zStep = pRadAccessData->zStep;
	double zStart = pRadAccessData->zStart;

	const long long RADSZ = pRadAccessData->nz * pRadAccessData->nx * (pRadAccessData->ne * 2);
	float* radx = (float*)calloc(RADSZ, sizeof(float));
	float* radz = (float*)calloc(RADSZ, sizeof(float));

	srTDriftSpace internal_drift(dftLen);

#if DEBUG_ZPD > 0
	srTSRWRadStructAccessData newRad0(pRadAccessData);
	fprintf(stderr, "ref inp hash= 0x%016zx\n", newRad0.hashcode());
	srTZonePlate::PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	fprintf(stderr, "ref zp  hash= 0x%016zx\n", newRad0.hashcode());
	internal_drift.PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	fprintf(stderr, "ref drf hash= 0x%016zx\n", newRad0.hashcode()); 
	
	newRad0.dumpBinData("junk.zpd.ref.bin", "zpd_reference");
#endif

#if DEBUG_ZPD > 1
	string fname = "junk.zpd0.bin";
	pRadAccessData->dumpBinData(fname, "zpd0");
#endif

	/*
	if (int result = srTZonePlate::PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
		fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
		return result;
	}
	*/

	

#if DEBUG_ZPD > 1
	// fname = "junk.zpd01.bin";
	// newRad.dumpBinData(fname, "zpd01");
	
	ofstream junkfdiv("junk.main.txt");
	junkfdiv << "#nzdiv,nxdiv " << nzdiv << " " << nxdiv << endl;
#endif

	//if (nzdiv > 1) nzdiv += (pRadAccessData->nz % nzdiv ? 1 : 0);
	//if (nxdiv > 1) nxdiv += (pRadAccessData->nx % nxdiv ? 1 : 0);

	// original center
	double xc0 = pRadAccessData->xStart + pRadAccessData->xStep * pRadAccessData->nx / 2.0;
	double zc0 = pRadAccessData->zStart + pRadAccessData->zStep * pRadAccessData->nz / 2.0;
	// average length of Width and Height
	double wdavg = (pRadAccessData->xStep * pRadAccessData->nx / 2.0 + pRadAccessData->zStep * pRadAccessData->nz / 2.0) / 2.0;
	
	const int SZZ = ceil(pRadAccessData->nz * 1.0 / nzdiv);
	const int SZX = ceil(pRadAccessData->nx * 1.0 / nxdiv);
	for (int iz = 0; iz < nzdiv; ++iz) {
		for (int ix = 0; ix < nxdiv; ++ix) {

#if DEBUG_ZPD > 1
			fprintf(stderr, "## split (%d,%d) of (%d,%d) sz %d sx %d\n", iz, ix, nzdiv, nxdiv, SZZ, SZX);
#endif
			srTSRWRadStructAccessData newRad(pRadAccessData);

			memset(newRad.pBaseRadX, 0.0, RADSZ * sizeof newRad.pBaseRadX[0]);
			memset(newRad.pBaseRadZ, 0.0, RADSZ * sizeof newRad.pBaseRadZ[0]);

			/*if (nzdiv != 1 || nxdiv != 1) {
				int npad = 2;
				if (ix == 0) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else if (ix == 2*SZX) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1, 1);
			}*/
			int iz1 = min((iz+1) * SZZ, pRadAccessData->nz);
			int ix1 = min((ix+1) * SZX, pRadAccessData->nx);
			sel_sub_cell(&newRad, *pRadAccessData, iz*SZZ, iz1, ix*SZX, ix1, 0, 0, 0, 0);

			if (nzdiv == 1 && nxdiv == 1) {
				newRad.xc = pRadAccessData->xc; // for hash calc
				newRad.zc = pRadAccessData->zc;
				cerr << "hash= " << std::hex << pRadAccessData->hashcode() << endl;
				cerr << "hash= " << std::hex << newRad.hashcode() << std::dec << " (iz= " << iz << "/" << nzdiv << ", ix= " << ix << "/" << nxdiv << ")\n";
			}

#if DEBUG_ZPD > 1
			fname = "junk.zpd00." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz*SZZ << " " << ix*SZX << " " << fname << endl;
#endif
			double xc = newRad.xStart + newRad.xStep * newRad.nx / 2;
			if (nxdiv > 1) newRad.xStart -= xc;
			double ang_x = -atan(xc / dftLen);
			double zc = newRad.zStart + newRad.zStep * newRad.nz / 2;
			if (nzdiv > 1) newRad.zStart -= zc;
			double ang_z = -atan(zc / dftLen);

			srTRadResize resz;
			
			double ri = hypot(xc - xc0, zc - zc0);

			resz.pxd = pdcenter + ri / wdavg * (pdedge - pdcenter);
			resz.pzd = pdcenter + ri / wdavg * (pdedge - pdcenter); 
			if (nzdiv > 1 || nxdiv > 1) {
				RadResizeGen(newRad, resz);

				fprintf(stderr, "scaling pxd= %g nx= %d pzd= %g nz= %d (ri=%g wdavg= %f, ratio= %g)\n",
					resz.pxd, newRad.nx, resz.pzd, newRad.nz, ri, wdavg, ri / wdavg);
			}

#if DEBUG_ZPD > 1
			junkfdiv << "#pdizix " << iz << " " << ix << " pzd " << resz.pzd << " pxd " << resz.pxd << endl;
#endif

#if DEBUG_ZPD > 1
			fprintf(stderr, "xc=%g zc=%g ang_x=%g ang_z=%g\n", xc, zc, ang_x, ang_z);
			fprintf(stderr, "Slice iz=%d+%d ix=%d+%d z=[%g %g], x=[%g %g]\n",
				iz, SZZ, ix, SZX, newRad.zStart, newRad.zStart+newRad.nz*newRad.zStep, newRad.xStart, newRad.xStart + newRad.nx * newRad.xStep);

			fname = "junk.zpd01." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz*SZZ << " " << ix*SZX << " " << fname << endl;
			if (nzdiv == 1 && nxdiv == 1) {
				newRad.xc = pRadAccessData->xc; // for hash calc
				newRad.zc = pRadAccessData->zc;
				cerr << "hash= " << std::hex << newRad.hashcode() << std::dec << " (iz= " << iz << "/" << nzdiv << ", ix= " << ix << "/" << nxdiv << ")\n";
			}

			fprintf(stderr, "zpd inp hash= 0x%016zx\n", newRad.hashcode());
#endif
//			if (nzdiv > 1 || nxdiv > 1) {
//				srTOptAngle inter_angle(-ang_x, -ang_z);
//				if (int result = inter_angle.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
//					fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
//					return result;
//				}
//			}
			
			if (int result = srTZonePlate::PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}

#if DEBUG_ZPD > 2

			fname = "junk.zpd02." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);

			fprintf(stderr, "zpd zp  hash= 0x%016zx\n", newRad.hashcode());
#endif

			if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}
			

			if (nzdiv > 1 || nxdiv > 1) {
				srTOptAngle inter_angle(ang_x, ang_z);
				if (int result = inter_angle.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
					fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
					return result;
				}
			}

#if DEBUG_ZPD > 2
			fprintf(stderr, "zpd dft  hash= 0x%016zx\n", newRad.hashcode());
#endif

			accumulate_rad(radx, radz, pRadAccessData->nz, pRadAccessData->nx, xStart, xStep, zStart, zStep, newRad);

#if DEBUG_ZPD > 1
			{
				fname = "junk.zpd1." + to_string(iz) + "_" + to_string(ix) + ".bin";
				newRad.dumpBinData(fname, fname);
				junkfdiv << "#fout " << iz * SZZ << " " << ix * SZX << " " << fname << endl;
				{
					const auto itm = std::minmax_element(radx, radx + RADSZ);
					fprintf(stderr, "min= %g max %g\n", *itm.first, *itm.second);
				}
				fprintf(stderr, "Slice (iz, ix) = %d %d Hash= 0x%016zx, Hash= 0x%016zx\n", iz, ix, pRadAccessData->hashcode(), newRad.hashcode());
			}
#endif

#if DEBUG_ZPD > 2
			{
				// dump the accumulated
				fname = "junk.zpd1.sum." + to_string(iz) + "_" + to_string(ix) + ".bin";
				ofstream out(fname, ios::out | ios::binary);
				out.write((char*)&pRadAccessData->nz, sizeof(long));
				out.write((char*)&pRadAccessData->nx, sizeof(long));
				out.write((char*)&pRadAccessData->ne, sizeof(long));
				out.write((char*)&RADSZ, sizeof(long long));
				fprintf(stderr, "write accum field %d %d %d RADSZ=%d\n", pRadAccessData->nz, pRadAccessData->nx, pRadAccessData->ne, RADSZ);
				out.write((char*)radz, RADSZ * sizeof(float));
				out.write((char*)radx, RADSZ * sizeof(float));
			}
#endif
		}
	}
	
	
	memcpy(pRadAccessData->pBaseRadX, radx, RADSZ * sizeof(float));
	memcpy(pRadAccessData->pBaseRadZ, radz, RADSZ * sizeof(float));

	free(radx);
	free(radz);

#if DEBUG_ZPD > 0
	fprintf(stderr, "DONE ZPD div=(%d %d) nz=%d nx= %d\n", nzdiv, nxdiv, pRadAccessData->nz, pRadAccessData->nx);
	fflush(stderr);
	pRadAccessData->dumpBinData("junk.zpd.bin", "junk.zpd.bin");
#endif

#if DEBUG_ZPD > 1
	junkfdiv.close(); 
#endif

	return 0;
}


int sel_sub_ring(srTSRWRadStructAccessData& dst, /* const */ srTSRWRadStructAccessData* wfr, double zx0, double wzx, int iring, int n, int npad = 0)
{
	// skip
	if (iring == 0 && zx0 > wfr->zStep / 2.0) return 0; // the center but zx0 is not origin
	if (iring > 0 && zx0 < wfr->zStep) return 0; // not the center but zx0 is too small

	// n is total grid, including npad. i.e. there are (n - 2*npad) non-zero points each dimension
	// const int npad = 4;
	// npad.. < 0 means does not count, use whatever left
	
	dst.ne = wfr->ne;
	dst.xStep = dst.zStep = wzx / n;
	double wz = 0.0, wx = 0.0;
	fprintf(stderr, "select ring %d outside of (%f, %f) width= %f n= %d\n", iring, zx0, zx0, wzx, n);
	switch (iring) {
	case 0: wz = wx = wzx;  break; // fix later
	case 1: // right block
		dst.zStart = -zx0; wz = 2 * zx0 + wzx; 
		dst.xStart = zx0 + dst.xStep; wx = wzx;
		break;
	case 2: // bottom block
		dst.zStart = -zx0 - wzx; wz = wzx;
		dst.xStart = -zx0; wx = 2 * zx0 + wzx;
		break;
	case 3: // left
		dst.zStart = -zx0 - wzx; wz = 2 * zx0 + wzx;
		dst.xStart = -zx0 - wzx; wx = wzx;
		break;
	case 4: // top
		dst.xStart = -zx0 - wzx; wz = wzx;
		dst.zStart = zx0 + dst.zStep; wx = 2 * zx0 + wzx;
		break;
	default: exit(-1);
	}
	if (iring == 0) {
		dst.zStep = dst.xStep = wzx / (n - 1);
		dst.zStart = dst.xStart = -wzx / 2.0;
		dst.nz = dst.nx = next_fft_num(n+2*npad);
	}
	else {
		dst.zStart -= npad * dst.zStep;
		dst.xStart -= npad * dst.xStep;
		long nx0 = round(wx / dst.xStep) + 1 + npad * 2;
		long nz0 = round(wz / dst.zStep) + 1 + npad * 2;
		dst.nz = next_fft_num(nz0);
		dst.nx = next_fft_num(nx0);
		fprintf(stderr, "FFT num iring %d, nx %d +%d, nz %d +%d\n", iring, nx0, dst.nx - nx0, nz0, dst.nz - nz0);
	}
	
	fprintf(stderr, "ring size: [%g %g], [%g %g]\n",
		dst.zStart + npad * dst.zStep, dst.zStart + npad * dst.zStep + wz,
		dst.xStart + npad * dst.xStep, dst.xStart + npad * dst.xStep + wx);
	
	
	dst.xWfrMin = dst.xStart;
	dst.zWfrMin = dst.zStart;
	dst.xWfrMax = dst.xStart + dst.nx * dst.xStep;
	dst.zWfrMax = dst.zStart + dst.nz * dst.zStep;

	// sampling
	float* rdestz = dst.pBaseRadZ, *rdestx = dst.pBaseRadX;
	for (int iz = 0; iz < dst.nz; ++iz) {
		double z = dst.zStart + iz * dst.zStep;
		if (z - dst.zStart  < npad * dst.zStep || z - dst.zStart > wz) continue;
		for (int ix = 0; ix < dst.nx; ++ix) {
			double x = dst.xStart + ix * dst.xStep;
			if (x - dst.xStart  < npad * dst.xStep || x - dst.xStart > wx) continue;
			do_efield_2d_sample(rdestx + iz * dst.nx * dst.ne * 2 + ix * dst.ne*2, z, x, wfr->pBaseRadX, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
			do_efield_2d_sample(rdestz + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadZ, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
			
			//rdestx += wfr->ne * 2;
			//rdestz += wfr->ne * 2;
		}
		//
	}
	return dst.nz * dst.nx;
}

// zxll0 - lower left z/x coordinate, if <0 it is the center region
// wzx - width (smaller one for rectangular)
int sel_sub_ring(srTSRWRadStructAccessData& dst, /* const */ srTSRWRadStructAccessData* wfr, int iring, double zxll0, double wzx, double minstepsz)
{	
	// the center block or other 4 blocks
	assert((iring == 0 && zxll0 < 0) || (iring > 0 && iring <= 4 && zxll0 > 0));
	dst.ne = wfr->ne;
	const double GAP = minstepsz;
	const double WDSHORT = (iring == 0 ? 2*fabs(zxll0) : wzx - GAP);
	const double WDLONG = (iring == 0 ? 2*fabs(zxll0) : zxll0 * 2 + wzx);

	double wz = 0.0, wx = 0.0;
	// fprintf(stderr, "select ring %d outside of (%f, %f) width= %f n= %d\n", iring, zx0, zx0, wzx, n);
	switch (iring) {
	case 0:
		dst.xStart = dst.zStart = zxll0;  wz = wx = -2*zxll0;
		break;
	case 1: // right block
		dst.zStart = -zxll0; wz = WDLONG;
		dst.xStart = zxll0 + GAP; wx = WDSHORT;
		break;
	case 2: // bottom block
		dst.zStart = -zxll0 - wzx; wz = WDSHORT;
		dst.xStart = -zxll0; wx = WDLONG;
		break;
	case 3: // left
		dst.zStart = -zxll0 - wzx; wz = WDLONG;
		dst.xStart = -zxll0 - wzx; wx = WDSHORT;
		break;
	case 4: // top
		dst.zStart = zxll0 + GAP; wz = WDSHORT;
		dst.xStart = -zxll0 - wzx; wx = WDLONG;
		break;
	default: exit(-1);
	}
	
	long nx0 = ceil(wx / minstepsz);
	long nz0 = ceil(wz / minstepsz);
	dst.nz = next_fft_num(nz0);
	dst.nx = next_fft_num(nx0);
	dst.xStep = wx / (dst.nx - 1);
	dst.zStep = wz / (dst.nz - 1);
	fprintf(stderr, "FFT num iring %d, nx %d +%d, nz %d +%d\n", iring, nx0, dst.nx - nx0, nz0, dst.nz - nz0);
	

	fprintf(stderr, "ring size: [%g %g], [%g %g]\n",
		dst.zStart, dst.zStart + wz, dst.xStart, dst.xStart + wx);


	dst.xWfrMin = dst.xStart;
	dst.zWfrMin = dst.zStart;
	dst.xWfrMax = dst.xStart + (dst.nx-1) * dst.xStep;
	dst.zWfrMax = dst.zStart + (dst.nz-1) * dst.zStep;

	// sampling
	int npad = 0;
	float* rdestz = dst.pBaseRadZ, * rdestx = dst.pBaseRadX;
	for (int iz = 0; iz < dst.nz; ++iz) {
		double z = dst.zStart + iz * dst.zStep;
		if (z - dst.zStart  < npad * dst.zStep || z - dst.zStart > wz) continue;
		for (int ix = 0; ix < dst.nx; ++ix) {
			double x = dst.xStart + ix * dst.xStep;
			if (x - dst.xStart  < npad * dst.xStep || x - dst.xStart > wx) continue;
			do_efield_2d_sample(rdestx + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadX, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
			do_efield_2d_sample(rdestz + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadZ, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);

			//rdestx += wfr->ne * 2;
			//rdestz += wfr->ne * 2;
		}
		//
	}
	return dst.nz * dst.nx;
}

// pdedge = 1, keep original stepsize
// pdcenter < 1, larger stepsize
// iring = 0, 1, 2, ..., numring - 1
double sub_ring_stepsize(double stepsz, double pdcenter, double pdedge, int numring, int iring)
{
	assert(numring > 0);
	if (numring <= 1) return stepsz;
	double pd = pdcenter + iring * (pdedge - pdcenter) / (numring - 1);
	return stepsz / pd;
}

int srTZonePlateD::PropagateRad2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!! Method=Rad2] The ZP (nzdiv=" << nzdiv << ",nxdiv=" << nxdiv
		<< ") with aperture and drift L= " << dftLen << endl;
#if DEBUG_ZPD > 1
	fprintf(stderr, "orig field z [%g %g] nz= %d\n", pRadAccessData->zStart, pRadAccessData->zStart + pRadAccessData->zStep * pRadAccessData->nz, pRadAccessData->nz);
	fprintf(stderr, "orig field x [%g %g] nx= %d\n", pRadAccessData->xStart, pRadAccessData->xStart + pRadAccessData->xStep * pRadAccessData->nx, pRadAccessData->nx);
	// cout << "zp input hash=" << std::hex << pRadAccessData->hashcode() << endl;
#endif

	double xStep = pRadAccessData->xStep;
	double xStart = pRadAccessData->xStart;
	double zStep = pRadAccessData->zStep;
	double zStart = pRadAccessData->zStart;

	const long long RADSZ = pRadAccessData->nz * pRadAccessData->nx * (pRadAccessData->ne * 2);
	float* radx = (float*)calloc(RADSZ, sizeof(float));
	float* radz = (float*)calloc(RADSZ, sizeof(float));

	srTDriftSpace internal_drift(dftLen);

#if DEBUG_ZPD > 0
	srTSRWRadStructAccessData newRad0(pRadAccessData);
	//fprintf(stderr, "ref inp hash= 0x%016zx\n", newRad0.hashcode());
	srTZonePlate::PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref zp  hash= 0x%016zx\n", newRad0.hashcode());
	internal_drift.PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref drf hash= 0x%016zx\n", newRad0.hashcode());

	newRad0.dumpBinData("junk.zpd.ref.bin", "zpd_reference");

#endif

#if DEBUG_ZPD > 1
	string fname = "junk.zpd0.bin";
	pRadAccessData->dumpBinData(fname, "zpd0");
#endif

	/*
	if (int result = srTZonePlate::PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
		fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
		return result;
	}
	*/



#if DEBUG_ZPD > 1
	// fname = "junk.zpd01.bin";
	// newRad.dumpBinData(fname, "zpd01");

	ofstream junkfdiv("junk.main.txt");
	junkfdiv << "#nzdiv,nxdiv " << nzdiv << " " << nxdiv << endl;
#endif

	//if (nzdiv > 1) nzdiv += (pRadAccessData->nz % nzdiv ? 1 : 0);
	//if (nxdiv > 1) nxdiv += (pRadAccessData->nx % nxdiv ? 1 : 0);

	// original center
	double xc0 = pRadAccessData->xStart + pRadAccessData->xStep * pRadAccessData->nx / 2.0;
	double zc0 = pRadAccessData->zStart + pRadAccessData->zStep * pRadAccessData->nz / 2.0;
	// average length of Width and Height
	// double wdavg = (pRadAccessData->xStep * pRadAccessData->nx / 2.0 + pRadAccessData->zStep * pRadAccessData->nz / 2.0) / 2.0;

	const int SZZ = ceil(pRadAccessData->nz * 1.0 / nzdiv);
	const int SZX = ceil(pRadAccessData->nx * 1.0 / nxdiv);

	const double fac = 1.0 / 1.001; // sqrt(2.0); // scaling factor of ring width
	const double FW = (pRadAccessData->nz - 1) * pRadAccessData->zStep;
	const double w0 = (nzdiv <= 1 ? FW : FW * (1-fac) / (1-pow(fac, nzdiv))) / 2;
	fprintf(stderr, "# center ring r = %g fac= %g, outer width= %g tot NZ= %d NX= %d\n",
		w0, fac, w0 * pow(fac, nzdiv-1), pRadAccessData->nz, pRadAccessData->nx);

	vector<double> wi{ 2 * w0, w0 * fac };
	vector<double> zxll0{ -w0, w0 }; // starting corner of ring i
	for (int i = 2; i < nzdiv; ++i) { wi.push_back(w0 * pow(fac, i)); zxll0.push_back(zxll0[i - 1] + wi[i - 1]); }

	for (int iz = 0; iz < nzdiv; ++iz) {

		double minstepsz = sub_ring_stepsize(pRadAccessData->zStep, pdcenter, pdedge, nzdiv, iz);
		// double wi = (iz == 0 ? 2*w0 : w0 * pow(fac, iz));
		// zxi = (iz == 0 ? -w0 : zxi + w0 * pow(fac, iz - 1));

#if DEBUG_ZPD > 1
		fprintf(stderr, "## split %d / %d xzll0= %g wd= %g\n", iz, nzdiv, zxll0[iz], wi[iz]);
#endif

		for (int iring = (iz == 0 ? 0 : 1); iring <= (iz == 0 ? 0 : 4); ++iring) {
			srTSRWRadStructAccessData newRad(pRadAccessData);
			memset(newRad.pBaseRadX, 0.0, RADSZ * sizeof newRad.pBaseRadX[0]);
			memset(newRad.pBaseRadZ, 0.0, RADSZ * sizeof newRad.pBaseRadZ[0]);

			sel_sub_ring(newRad, pRadAccessData, iring, zxll0[iz], wi[iz], minstepsz);
			
			if (nzdiv == 1) {
				newRad.xc = pRadAccessData->xc; // for hash calc
				newRad.zc = pRadAccessData->zc;
				cerr << "hash= " << std::hex << pRadAccessData->hashcode() << endl;
			}

			fprintf(stderr, "ring %d %d, z [%g %g] nz= %d\n", iz, iring, newRad.zStart, newRad.zStart + newRad.zStep * (newRad.nz-1), newRad.nz);
			fprintf(stderr, "ring %d %d, x [%g %g] nx= %d\n", iz, iring, newRad.xStart, newRad.xStart + newRad.xStep * (newRad.nx-1), newRad.nx);


#if DEBUG_ZPD > 1
			fname = "junk.zpd00." + to_string(iz) + "_" + to_string(iring) + ".bin";
			newRad.dumpBinData(fname, "after slice");
			junkfdiv << "#fin " << iz << " " << iring << " " << fname << endl;
#endif
			double xc = newRad.xStart + newRad.xStep * (newRad.nx-1.0) / 2;
			double ang_x = -atan(xc / dftLen);
			double zc = newRad.zStart + newRad.zStep * (newRad.nz-1.0) / 2;
			double ang_z = -atan(zc / dftLen);
			/*
			if (iring > 0 && iring % 2 == 1) {
				newRad.xStart -= xc; ang_z = 0;
			}
			if (iring > 0 && iring % 2 == 0) {
				newRad.zStart -= zc; ang_x = 0;
			}
			*/
			if (iring > 0) { newRad.xStart -= xc; newRad.zStart -= zc; }

			srTRadResize resz;

			double ri = hypot(xc - xc0, zc - zc0);

			// resz.pzd = resz.pxd = 1 / pow(fac, iz);
			/*
			resz.pzd = resz.pxd = (nzdiv <= 1 ? pdedge : pdcenter + (pdedge - pdcenter) * iz / (nzdiv-1.0));

			RadResizeGen(newRad, resz);
			fprintf(stderr, "zStep= %g scaling pxd= %g nx= %d pzd= %g nz= %d\n",
					newRad.zStep, resz.pxd, newRad.nx, resz.pzd, newRad.nz);
			*/

#if DEBUG_ZPD > 1
			// junkfdiv << "#pdizix " << iz << " " << ix << " pzd " << resz.pzd << " pxd " << resz.pxd << endl;
#endif

#if DEBUG_ZPD > 1
			fprintf(stderr, "xc=%g zc=%g ang_x=%g ang_z=%g\n", xc, zc, ang_x, ang_z);
			fprintf(stderr, "--- Slice iz=%d iring=%d z=[%g %g] %d, x=[%g %g] %d\n",
				iz, iring, newRad.zStart, newRad.zStart + newRad.nz * newRad.zStep, newRad.nz,
				newRad.xStart, newRad.xStart + newRad.nx * newRad.xStep, newRad.nx);

			fname = "junk.zpd01." + to_string(iz) + "_" + to_string(iring) + ".bin";
			newRad.dumpBinData(fname, "after resize and shift");
			junkfdiv << "#fin " << iz << " " << iring << " " << fname << endl;
			if (nzdiv == 1 && nxdiv == 1) {
				newRad.xc = pRadAccessData->xc; // for hash calc
				newRad.zc = pRadAccessData->zc;
				cerr << "hash= " << std::hex << newRad.hashcode() << std::dec << " iz= " << iz << "/" << nzdiv << "\n";
			}

			// fprintf(stderr, "zpd inp hash= 0x%016zx\n", newRad.hashcode());
#endif
			//			if (nzdiv > 1 || nxdiv > 1) {
			//				srTOptAngle inter_angle(-ang_x, -ang_z);
			//				if (int result = inter_angle.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
			//					fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
			//					return result;
			//				}
			//			}

			if (int result = srTZonePlate::PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}

#if DEBUG_ZPD > 2

			fname = "junk.zpd02." + to_string(iz) + "_" + to_string(iring) + ".bin";
			newRad.dumpBinData(fname, "after ZP, before drift");

			fprintf(stderr, "zpd zp  hash= 0x%016zx\n", newRad.hashcode());
#endif

			if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}


			if (iring > 0) {
				srTOptAngle inter_angle(ang_x, ang_z);
				if (int result = inter_angle.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
					fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
					return result;
				}
			}

#if DEBUG_ZPD > 2
			fprintf(stderr, "zpd dft  hash= 0x%016zx\n", newRad.hashcode());
			// fname = "junk.zpd03." + to_string(iz) + "_" + to_string(iring) + ".bin";
			//newRad.dumpBinData(fname, "after drift and phase/angle");
#endif

			accumulate_rad(radx, radz, pRadAccessData->nz, pRadAccessData->nx, xStart, xStep, zStart, zStep, newRad);

#if DEBUG_ZPD > 1
			{
				fname = "junk.zpd1." + to_string(iz) + "_" + to_string(iring) + ".bin";
				newRad.dumpBinData(fname, "after drift and phase/angle");
				junkfdiv << "#fout " << iz << " " << iring << " " << fname << endl;
				{
					const auto itm = std::minmax_element(radx, radx + RADSZ);
					fprintf(stderr, "min= %g max %g\n", *itm.first, *itm.second);
				}
				fprintf(stderr, "Slice (iz, iring) = %d %d Hash= 0x%016zx, Hash= 0x%016zx\n", iz, iring, pRadAccessData->hashcode(), newRad.hashcode());
			}
#endif

#if DEBUG_ZPD > 2
			{
				// dump the accumulated
				fname = "junk.zpd1.sum." + to_string(iz) + "_" + to_string(iring) + ".bin";
				ofstream out(fname, ios::out | ios::binary);
				out.write((char*)&pRadAccessData->nz, sizeof(long));
				out.write((char*)&pRadAccessData->nx, sizeof(long));
				out.write((char*)&pRadAccessData->ne, sizeof(long));
				out.write((char*)&RADSZ, sizeof(long long));
				fprintf(stderr, "write accum field %d %d %d RADSZ=%d\n", pRadAccessData->nz, pRadAccessData->nx, pRadAccessData->ne, RADSZ);
				out.write((char*)radz, RADSZ * sizeof(float));
				out.write((char*)radx, RADSZ * sizeof(float));
			}
#endif
		}
	}


	memcpy(pRadAccessData->pBaseRadX, radx, RADSZ * sizeof(float));
	memcpy(pRadAccessData->pBaseRadZ, radz, RADSZ * sizeof(float));

	free(radx);
	free(radz);

#if DEBUG_ZPD > 0
	fprintf(stderr, "DONE ZPD div=(%d %d) nz=%d nx= %d\n", nzdiv, nxdiv, pRadAccessData->nz, pRadAccessData->nx);
	fflush(stderr);
	pRadAccessData->dumpBinData("junk.zpd.bin", "junk.zpd.bin");
#endif

#if DEBUG_ZPD > 1
	junkfdiv.close();
#endif

	return 0;
}


int srTZonePlateD::PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	return PropagateRad2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
}


