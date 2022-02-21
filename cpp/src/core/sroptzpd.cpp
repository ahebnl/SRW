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

#include <cassert>
#include <algorithm>
#include <fstream>

using namespace std;

#define DEBUG_ZPD

template<class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
	return (v < lo ? lo : (v > hi ? hi : v));
}

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

void accumulate_rad(float* rdestx, float *rdestz, int nz, int nx, float x0, float dx, float z0, float dz, const srTSRWRadStructAccessData& r /*, const string& component */,
	float zshift = 0.0, float xshift = 0.0)
{
	// float* psrc = nullptr;
	// if (component == "X") psrc = r.pBaseRadX;
	// for each point (z, x) rdest, we sample a result from srTSRWRadStructAccessData r and accumulate
	// we are sampling (z+zshift, x+xshift) at r
	for (int iz = 0; iz < nz; ++iz) {
		float z = z0 + iz * dz;
		for (int ix = 0; ix < nx; ++ix) {
			float x = x0 + ix * dx;
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
		float z = wfr.zStart + iz * wfr.zStep;
		for (int ix = 0; ix < wfr.nx; ++ix) {
			float x = wfr.xStart + ix * wfr.xStep;
			if (z < z0 || z >= z1 || x < x0 || x >= x1) {
				int offset = iz * wfr.nx * wfr.ne * 2 + ix * wfr.ne * 2;
				memset(wfr.pBaseRadX + offset, 0.0, wfr.ne * 2 * sizeof(float));
				memset(wfr.pBaseRadZ + offset, 0.0, wfr.ne * 2 * sizeof(float));
				++nreset;
			}
		}
	}
#ifdef DEBUG_ZPD
	fprintf(stderr, "reset %d (%.2f%%) points out side of z [%f %f], x [%f %f] N = %d %d start %f %f step= %g %g\n",
		nreset, nreset * 100.0 / NTOT, z0, z1, x0, x1, wfr.nz, wfr.nx, wfr.zStart, wfr.xStart, wfr.zStep, wfr.xStep);
#endif
}

// psrc - full (nz, nx, ne) 3D array, e has Ex, Ez
// pdst - destination 
void copy_sub_cell(float* pdst, const float* psrc, int nz, int nx, int ne, int iz0, int iz1, int ix0, int ix1)
{
#ifdef DEBUG_ZPD
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

#ifdef DEBUG_ZPD
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
#ifdef DEBUG_ZPD
	{
		const auto itm = std::minmax_element(psrc, psrc + nz * nx * ne * 2);
		fprintf(stderr, "src N (%d,%d,%d,2) %ld min= %g max %g\n", nz, nx, ne, nz * nx * ne * 2, *itm.first, *itm.second);
	}
#endif

	iz0 = clamp(iz0, 0, nz);
	ix0 = clamp(ix0, 0, nx);
	iz1 = clamp(iz1, 0, nz);
	ix1 = clamp(ix1, 0, nx);
	const int ne1 = ne;
	const long long SZZFULL = nx * ne * 2;
	const long long SZ = (ix1 - ix0) * ne * 2;
	const long long SZ1 = nx1 * ne1 * 2; // dest step size
	pdst += jz * SZ1 + jx * ne1 * 2;
	for (int iz = iz0; iz < iz1; ++iz) {
		memcpy(pdst, psrc + iz * SZZFULL + ix0 * ne * 2, SZ * (sizeof psrc[0]));
		pdst += SZ1;
	}

#ifdef DEBUG_ZPD
	{
		const auto itm = std::minmax_element(pdst - SZ * (iz1 - iz0), pdst);
		fprintf(stderr, "dst N (%d,%d,%d,2) %ld min= %g max %g\n", iz1 - iz0, ix1 - ix0, ne, SZ * (iz1 - iz0), *itm.first, *itm.second);
	}
#endif
}

void sel_sub_cell(srTSRWRadStructAccessData * dst, /* const */ srTSRWRadStructAccessData& wfr, int iz0, int iz1, int ix0, int ix1, int npad=4)
{
	// const int npad = 4;
	dst->nz = iz1 - iz0 + npad * 2;
	dst->nx = ix1 - ix0 + npad * 2;
	CGenMathFFT::NextCorrectNumberForFFT(dst->nz);
	CGenMathFFT::NextCorrectNumberForFFT(dst->nx);
	fprintf(stderr, "nz = %d or %d (next good fft num)\n", iz1 - iz0, dst->nz);
	fprintf(stderr, "nx = %d or %d\n", ix1 - ix0, dst->nx);

	long nz0 = iz1 - iz0;
	long nx0 = ix1 - ix0;
	int jz = (dst->nz - nz0) / 2;
	int jx = (dst->nx - nx0) / 2;
	copy_sub_cell_gen(dst->pBaseRadX, dst->nz, dst->nx, jz, jx, wfr.pBaseRadX, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);
	copy_sub_cell_gen(dst->pBaseRadZ, dst->nz, dst->nx, jz, jx, wfr.pBaseRadZ, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);

	dst->ne = wfr.ne;
	
	dst->xStep = wfr.xStep;
	dst->zStep = wfr.zStep; 
	dst->xWfrMin = dst->xStart = wfr.xStart + ix0 * wfr.xStep - jx*dst->xStep;
	dst->zWfrMin = dst->zStart = wfr.zStart + iz0 * wfr.zStep - jz*dst->zStep;
	
	dst->xWfrMax = dst->xStart + dst->nx * dst->xStep;
	dst->zWfrMax = dst->zStart + dst->nz * dst->zStep;

	dst->ComputeRadMoments();
}

void dump_raddata_txt(const srTSRWRadStructAccessData& wfr, const char *fname, const char *title)
{
	fprintf(stderr, "dumping file: %s\n", fname);
	ofstream out(fname);
	out << "# " << wfr.nz << ' ' << wfr.nx << ' ' << wfr.ne << ' ' << wfr.hashcode() << endl;
	out << "# " << title << endl;
	out << "# zStart= " << wfr.zStart << " zStep= " << wfr.zStep << endl;
	out << "# xStart= " << wfr.xStart << " xStep= " << wfr.xStep << endl;
	out << "# zWfrMin= " << wfr.zWfrMin << " zWfrMax= " << wfr.zWfrMax << endl; 
	out << "# xWfrMin= " << wfr.xWfrMin << " xWfrMax= " << wfr.xWfrMax << endl;
	{
		const auto itm1 = std::minmax_element(wfr.pBaseRadZ, wfr.pBaseRadZ + wfr.nz * wfr.nx * wfr.ne * 2);
		out << "# EzMinMax " << *itm1.first << " " << *itm1.second << endl;
		const auto itm2 = std::minmax_element(wfr.pBaseRadX, wfr.pBaseRadX + wfr.nz * wfr.nx * wfr.ne * 2);
		out << "# ExMinMax " << *itm2.first << " " << *itm2.second << endl;
	}

	const long long N = wfr.nz * wfr.nx * wfr.ne * 2;
	out << "# Ez " << N << endl;
	for (long long i = 0; i < N; i += 2) {
		out << wfr.pBaseRadZ[i] << ' ' << wfr.pBaseRadZ[i+1] << '\n';
	}
	out << "# Ex\n";
	for (long long i = 0; i < N; i += 2) {
		out << wfr.pBaseRadX[i] << ' ' << wfr.pBaseRadX[i+1] << '\n';
	}
	out << endl;

	fprintf(stderr, "file dumped: %s %s\n", fname, title);
}

void dump_raddata(const srTSRWRadStructAccessData& wfr, const char* fname, const char* title)
{
	fprintf(stderr, "dumping bin file: %s\n", fname);
	ofstream out(fname, ios::out | ios::binary);
	const int HDSZ = 64;
	char buf[HDSZ];
	memset(buf, 0, HDSZ);
	strncpy(buf, "SRWB1", 5);
	strncpy(buf + 5, title, HDSZ-5);
	out.write(buf, HDSZ);
	out.write((char*)&wfr.nz, sizeof(long));
	out.write((char*)&wfr.nx, sizeof(long));
	out.write((char*)&wfr.ne, sizeof(long));
	
	out.write((char*)&wfr.zStart, sizeof(double));
	out.write((char*)&wfr.zStep, sizeof(double));
	out.write((char*)&wfr.xStart, sizeof(double));
	out.write((char*)&wfr.xStep, sizeof(double));
	
	out.write((char*)&wfr.zWfrMin, sizeof(double));
	out.write((char*)&wfr.zWfrMax, sizeof(double));
	out.write((char*)&wfr.xWfrMin, sizeof(double));
	out.write((char*)&wfr.xWfrMax, sizeof(double));
	
	const long long N = wfr.nz * wfr.nx * wfr.ne * 2;
	out.write((char*)wfr.pBaseRadZ, N * sizeof(float));
	out.write((char*)wfr.pBaseRadX, N * sizeof(float));

	out.close();

	fprintf(stderr, "file dumped: %s %s (bin)\n", fname, title);
}

int srTZonePlateD::PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!!] The ZP (nz=" << nzdiv << ",nx=" << nxdiv
		 << ") with aperture and drift L= " << dftLen << " hash=" << pRadAccessData->hashcode() << endl;
	
	float xStep = pRadAccessData->xStep;
	float xStart = pRadAccessData->xStart;
	float zStep = pRadAccessData->zStep;
	float zStart = pRadAccessData->zStart;

	const long long RADSZ = pRadAccessData->nz * pRadAccessData->nx * (pRadAccessData->ne * 2);
	float* radx = (float*)calloc(RADSZ, sizeof(float));
	float* radz = (float*)calloc(RADSZ, sizeof(float));

	srTDriftSpace internal_drift(dftLen);

#ifdef DEBUG_ZPD
	srTSRWRadStructAccessData newRad0(pRadAccessData);
	srTZonePlate::PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	internal_drift.PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	dump_raddata(newRad0, "junk.zpd.ref.bin", "zpd_reference");
#endif

	srTSRWRadStructAccessData newRad(pRadAccessData);

#ifdef DEBUG_ZPD
	string fname = "junk.zpd0.bin";
	dump_raddata(newRad, fname.c_str(), "zpd0");

	ofstream junkfdiv("junk.main.txt");
	junkfdiv << "#nzdiv,nxdiv " << nzdiv << " " << nxdiv << endl;
#endif

	if (nzdiv > 1) nzdiv += (pRadAccessData->nz % nzdiv ? 1 : 0);
	if (nxdiv > 1) nxdiv += (pRadAccessData->nx % nxdiv ? 1 : 0);
	for (int iz = 0, SZZ = pRadAccessData->nz / nzdiv; iz < pRadAccessData->nz; iz += SZZ) {
		for (int ix = 0, SZX = pRadAccessData->nx / nxdiv; ix < pRadAccessData->nx; ix += SZX) {
			// if (iz != 0 || ix != SZX) continue;

#ifndef DEBUG_ZPD
			fprintf(stderr, "split %d sz %d, %d sz %d\n", iz, SZZ, ix, SZX);
#endif

			memset(newRad.pBaseRadX, 0.0, RADSZ * sizeof newRad.pBaseRadX[0]);
			memset(newRad.pBaseRadZ, 0.0, RADSZ * sizeof newRad.pBaseRadZ[0]);
			int npad = 4;
			if (nzdiv == 1 && nxdiv == 1) npad = 0;
			sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, npad);

#ifdef DEBUG_ZPD
			fname = "junk.zpd0." + to_string(iz/SZZ) + "_" + to_string(ix/SZX) + ".bin";
			dump_raddata(newRad, fname.c_str(), fname.c_str());
			junkfdiv << "#fin " << iz << " " << ix << " " << fname << endl;
#endif
			if (int result = srTZonePlate::PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}
			if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}
			accumulate_rad(radx, radz, pRadAccessData->nz, pRadAccessData->nx, xStart, xStep, zStart, zStep, newRad);

#ifdef DEBUG_ZPD
			{
				fname = "junk.zpd1." + to_string(iz / SZZ) + "_" + to_string(ix / SZX) + ".bin";
				dump_raddata(newRad, fname.c_str(), fname.c_str());
				junkfdiv << "#fout " << iz << " " << ix << " " << fname << endl;
				{
					const auto itm = std::minmax_element(radx, radx + RADSZ);
					fprintf(stderr, "min= %g max %g\n", *itm.first, *itm.second);
				}
				fprintf(stderr, "Slice (iz, ix) = %d %d Hash= %ld, Hash= %ld\n", iz, ix, pRadAccessData->hashcode(), newRad.hashcode());

				// dump the accumulated
				fname = "junk.zpd1.sum." + to_string(iz / SZZ) + "_" + to_string(ix / SZX) + ".bin";
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

#ifdef DEBUG_ZPD
	fprintf(stderr, "DONE ZPD div=(%d %d) nz=%d nx= %d\n", nzdiv, nxdiv, pRadAccessData->nz, pRadAccessData->nx);
	fflush(stderr);
	
	dump_raddata(*pRadAccessData, "junk.zpd.bin", "junk.zpd.bin");
	junkfdiv.close(); 
#endif

	return 0;
}



