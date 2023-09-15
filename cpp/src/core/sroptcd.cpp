/************************************************************************//**
 * File: sroptzp.cpp
 * Description: Optical element: Element Combined with Drift Space (CD)
 * Author: AN HE, Brookhaven National Laboratory, 2022
 *
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 * 
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptcd.h"
#include "sroptdrf.h"
#include "gmfft.h"
#include "sroptang.h"
#include <cassert>
#include <algorithm>
#include <fstream>

#if ! _MSC_VER
#include <unistd.h>
#endif

using namespace std;

#define DEBUG_ZPD 9  // No output of the debug information when DEBUG_ZPD=0 

template<class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
	return (v < lo ? lo : (v > hi ? hi : v));
}

static long next_fft_num(long n) {
	CGenMathFFT::NextCorrectNumberForFFT(n); return n;
}


#if DEBUG_ZPD > 0 && !_MSC_VER
#include <sys/resource.h>
#include <sys/time.h>
void print_usage(const char *tag)
{
  struct rusage u;
  getrusage(RUSAGE_SELF, &u);
  fprintf(stdout, "CPU time [%s]: %ld sec user mem %ld\n", tag, u.ru_utime.tv_sec, u.ru_maxrss);
}
#else
void print_usage(const char*) {}
#endif


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

void accumulate_rad(float* rdestx0, float *rdestz0, int nx, int nz, double x0, double dx, double z0, double dz, const srTSRWRadStructAccessData& r /*, const string& component */,
	double zshift = 0.0, double xshift = 0.0)
{
	// float* psrc = nullptr;
	// if (component == "X") psrc = r.pBaseRadX;
	// for each point (z, x) rdest, we sample a result from srTSRWRadStructAccessData r and accumulate
	// we are sampling (z+zshift, x+xshift) at r
	float *rdestx = rdestx0, *rdestz = rdestz0;
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

void CDRadStructHelper::add(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src)
{
	accumulate_rad(dest->pBaseRadX, dest->pBaseRadZ, dest->nx, dest->nz, dest->xStart, dest->xStep, dest->zStart, dest->zStep, *src);
}

// dest will keep its mem size (nz, nx, ne), but new xStart, xStep
void CDRadStructHelper::assign(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src)
{
	// fprintf(stderr, "BaseRadWasEmulated=%d\n", dest->BaseRadWasEmulated);	
	dest->zStart = src->zStart;
	dest->zStep = (src->nz - 1) * src->zStep / (dest->nz - 1);
	dest->xStart = src->xStart;
	dest->xStep = (src->nx - 1) * src->xStep / (dest->nx - 1);

	long long SZ = dest->nx * dest->nz * dest->ne * 2;
	memset(dest->pBaseRadX, 0, SZ * sizeof dest->pBaseRadX[0]);
	memset(dest->pBaseRadZ, 0, SZ * sizeof dest->pBaseRadZ[0]);

	accumulate_rad(dest->pBaseRadX, dest->pBaseRadZ, dest->nx, dest->nz, dest->xStart, dest->xStep, dest->zStart, dest->zStep, *src);
}

// dest samples from src's [xStart, xEnd], then add extra padding
// after the padding, dest->xStart <= xStart, ...
void CDRadStructHelper::sample(srTSRWRadStructAccessData* dest, const srTSRWRadStructAccessData* src,
	double xStart, double xEnd, double xStep, int nxpad,
	double zStart, double zEnd, double zStep, int nzpad)
{
	int ncx = (xEnd - xStart) / xStep + 1;
	int ncz = (zEnd - zStart) / zStep + 1;
	dest->nx = next_fft_num(ncx + 2 * nxpad);
	dest->nz = next_fft_num(ncz + 2 * nzpad);
	
	dest->xStep = xStep;
	dest->zStep = zStep;
	dest->xStart = xStart - xStep * nxpad;
	dest->zStart = zStart - zStep * nzpad;

	dest->ReAllocBaseRadAccordingToNeNxNz();

	long long SZ = dest->nx * dest->nz * dest->ne * 2; 
	memset(dest->pBaseRadX, 0, SZ * sizeof dest->pBaseRadX[0]);
	memset(dest->pBaseRadZ, 0, SZ * sizeof dest->pBaseRadZ[0]);
	
	
	accumulate_rad(dest->pBaseRadX, dest->pBaseRadZ, dest->nx, dest->nz, dest->xStart, dest->xStep, dest->zStart, dest->zStep, *src);
  
	// set padding 0
	int nxpad2 = dest->nx - ncx; // right side padding
	int nzpad2 = dest->nz - ncz;
	for (int iz = 0; iz < dest->nz; ++iz) {
	  for (int ix = 0; ix < dest->nx; ++ix) {
		if (iz < nzpad || iz + nzpad2 >= dest->nz ||
            ix < nxpad || ix + nxpad2 >= dest->nx) { // skip the center region
			int c = iz * dest->nx * dest->ne * 2 + ix * dest->ne * 2;
			memset(dest->pBaseRadX + c, 0, dest->ne * 2 * sizeof(float));
			memset(dest->pBaseRadZ + c, 0, dest->ne * 2 * sizeof(float));
        }
      }
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
				memset(wfr.pBaseRadX + offset, 0, wfr.ne * 2 * sizeof(float));
				memset(wfr.pBaseRadZ + offset, 0, wfr.ne * 2 * sizeof(float));
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
// pdst - destination, is a bigger grid nz1, nx1, ne1==ne, place the src "cubic" inside dst and starts from (jz, jx, je)
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
		fprintf(stderr, "dst N (%d,%d,%d) srcnz= %d, srcnx=%d min= %g max %g\n", nz1, nx1, ne, iz1 - iz0, ix1 - ix0, *itm.first, *itm.second);
		assert(nz1 == next_fft_num(nz1));
		assert(nx1 == next_fft_num(nx1));
	}
#endif
}

// dst must have large enough 
void sel_sub_cell(srTSRWRadStructAccessData * dst, /* const */ srTSRWRadStructAccessData& wfr, int iz0, int iz1, int ix0, int ix1,
	int npadz0=-1, int npadz1 = 0, int npadx0 = -1, int npadx1 = 0)
{
	// npad.. < 0 means does not count, use whatever left
	dst->nz = iz1 - iz0 + max(0, npadz0) + max(0, npadz1);
	dst->nx = ix1 - ix0 + max(0, npadx0) + max(0, npadx1);
	CGenMathFFT::NextCorrectNumberForFFT(dst->nz);
	CGenMathFFT::NextCorrectNumberForFFT(dst->nx);
	fprintf(stderr, "select sub cell iz= [%d, %d) or %d (next good fft num) tot src nz= %d\n", iz0, iz1, dst->nz, wfr.nz);
	fprintf(stderr, "select sub cell ix= [%d, %d) or %d (tot src nx= %d)\n", ix0, ix1, dst->nx, wfr.nx);

	long nz0 = iz1 - iz0;
	long nx0 = ix1 - ix0;
	// int jz = npadz0 > 0 ? npadz0 : dst->nz - nz0 - max(0, npadz1); // actual left padding
	// int jx = npadx0 > 0 ? npadx0 : dst->nx - nx0 - max(0, npadx1); // actual left padding
	int jz = max(0, npadz0);
	int jx = max(0, npadx0);
	copy_sub_cell_gen(dst->pBaseRadX, dst->nz, dst->nx, jz, jx, wfr.pBaseRadX, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);
	copy_sub_cell_gen(dst->pBaseRadZ, dst->nz, dst->nx, jz, jx, wfr.pBaseRadZ, wfr.nz, wfr.nx, wfr.ne, iz0, iz1, ix0, ix1);

	dst->ne = wfr.ne;
	
	dst->xStep = wfr.xStep;
	dst->zStep = wfr.zStep; 
	dst->xStart = wfr.xStart + ix0 * wfr.xStep - jx*dst->xStep;
	dst->zStart = wfr.zStart + iz0 * wfr.zStep - jz*dst->zStep;
	
	dst->xWfrMin = dst->xStart;
	dst->zWfrMin = dst->zStart;
	dst->xWfrMax = dst->xStart + (dst->nx - 1) * dst->xStep;
	dst->zWfrMax = dst->zStart + (dst->nz - 1) * dst->zStep;

#if DEBUG_ZPD > 2
	fprintf(stderr, "src hash= 0x%016zx\n", wfr.hashcode());
	fprintf(stderr, "dst hash= 0x%016zx nx= %d nz= %d\n", dst->hashcode(), dst->nx, dst->nz);
#endif

	// dst->ComputeRadMoments();

#if DEBUG_ZPD > 2
	fprintf(stderr, "dst hash= 0x%016zx\n", dst->hashcode());
#endif
}

//void sel_sub_cell_zx(srTSRWRadStructAccessData* dst, /* const */ srTSRWRadStructAccessData& wfr, double z0, double z1, double x0, double x1,
//	int npadz0 = -1, int npadz1 = 0, int npadx0 = -1, int npadx1 = 0)
//{
//	// (std::min) avoid macro expansion, a Visual Studio fix
//	int iz0 = (std::max)(int(ceil((z0 - wfr.zStart) / wfr.zStep)), 0);
//	int iz1 = (std::min)(long(floor((z1 - wfr.zStart) / wfr.zStep)), wfr.nz);
//
//	int ix0 = (std::max)(int(ceil((x0 - wfr.xStart) / wfr.xStep)), 0);
//	int ix1 = (std::min)(long(floor((x1 - wfr.xStart) / wfr.xStep)), wfr.nx);
//
//	sel_sub_cell(dst, wfr, iz0, iz1, ix0, ix1, npadz0, npadz1, npadx0, npadx1);
//} no longer to be used, need to be deleted 


void srTCombinedDrift::init_dest_rad(srTSRWRadStructAccessData& rad, const srTSRWRadStructAccessData* pinrad, int wg) const
{
	// get the largest x/z Range and smallest stepsize
	double zStep = 1.0, xStep = 1.0; //  pinrad->zStep, xStep = pinrad->xStep;
	double wx = 0.0; //  pinrad->xStep* (pinrad->nx - 1);
	double wz = 0.0; //  pinrad->zStep* (pinrad->nz - 1);
	
	for (int ix = 0; ix < nxdiv; ++ix) {
		const double xr = (ix == 0 ? xdivs[ix] : (xdivs[ix] - xdivs[ix - 1])); // ratio of cell
		for (int iz = 0; iz < nzdiv; ++iz) {
			const double zr = (iz == 0 ? zdivs[iz] : (zdivs[iz] - zdivs[iz - 1]));
			int j = (ix * nzdiv + iz) * 5;
			// fprintf(stderr, " xr= %g zr= %g  zoom %g %g\n", xr, zr, crsz[j+2], crsz[j+4]);
			wx = max(wx, xr * crsz[j+1]); // 
			wz = max(wz, zr * crsz[j+3]);
			xStep = min(xStep, 1/crsz[j+2]); // new stepsize (factor)
			zStep = min(zStep, 1/crsz[j+4]);
			
			//fprintf(stderr, "ix,iz= %d %d wz= %g wx= %g, xStep= %g zStep=%g\n",
			//	ix, iz, wz, wx, xStep, zStep);
		}
	}

  if (wg) {
    fprintf(stderr, "resize step size %g %g (larger means smaller stepsize)\n", obsgrid[1], obsgrid[3]);
    xStep /= obsgrid[1]; zStep /= obsgrid[3];
    wx *= obsgrid[0]; wz *= obsgrid[2];
  }
	
  rad.xStep = xStep * pinrad->xStep;
	rad.zStep = zStep * pinrad->zStep;

	rad.nx = wx * pinrad->xStep * (pinrad->nx - 1) / rad.xStep + 1;
	rad.nz = wz * pinrad->zStep * (pinrad->nz - 1) / rad.zStep + 1;

	rad.xStart = -rad.xStep * (rad.nx - 1) / 2.0;
	rad.zStart = -rad.zStep * (rad.nz - 1) / 2.0;

	rad.ReAllocBaseRadAccordingToNeNxNz();

	long long sz = rad.nx * rad.nz * rad.ne * 2;
	memset(rad.pBaseRadX, 0, sz * sizeof rad.pBaseRadX[0]);
	memset(rad.pBaseRadZ, 0, sz * sizeof rad.pBaseRadZ[0]);
	double xEnd = rad.xStart + rad.xStep * (rad.nx - 1);
	double zEnd = rad.zStart + rad.zStep * (rad.nz - 1);
	fprintf(stderr, "the dest will have\n  step x %e z %e,\n  start x %e %e z %e %e,\n  width nx %d nz %d\n  "
		"fac wx= %g wz= %g xStep= %g zStep= %g\n",
		rad.xStep, rad.zStep, rad.xStart, xEnd, rad.zStart, zEnd, rad.nx, rad.nz, wx, wz, xStep, zStep);
}

void srTCombinedDrift::init_dest_rad2(srTSRWRadStructAccessData& rad, const srTSRWRadStructAccessData* pinrad) const
{
	fprintf(stderr, "init from xStart= %g xStep= %g\n", pinrad->xStart, pinrad->xStep);

	rad.xStart = obsgrid[0] * pinrad->xStart;
	rad.zStart = obsgrid[2] * pinrad->zStart;

  for(int i = 0; i < 10 && rad.xStart < -1e-7; ++i) { rad.xStart /= 2.0; rad.zStart /= 2.0; fprintf(stderr, "WARNING: rescale receiving plane: %g %g\n", rad.xStart, rad.zStart); }
  
	rad.xStep = pinrad->xStep / obsgrid[1];
	rad.zStep = pinrad->zStep / obsgrid[3];
	
	rad.nx = next_fft_num(2*std::abs(rad.xStart / rad.xStep)+1);
	rad.nz = next_fft_num(2 * std::abs(rad.zStart / rad.zStep) + 1);

	rad.xStart = -(rad.nx - 1) * rad.xStep / 2.0;
	rad.zStart = -(rad.nz - 1) * rad.zStep / 2.0;

	rad.ReAllocBaseRadAccordingToNeNxNz();

	long long sz = rad.nx * rad.nz * rad.ne * 2;
	memset(rad.pBaseRadX, 0, sz * sizeof rad.pBaseRadX[0]);
	memset(rad.pBaseRadZ, 0, sz * sizeof rad.pBaseRadZ[0]);
	double xEnd = rad.xStart + rad.xStep * (rad.nx - 1);
	double zEnd = rad.zStart + rad.zStep * (rad.nz - 1);
	fprintf(stderr, "the dest will have\n  step x %e z %e,\n  start x %e %e z %e %e,\n  width nx %d nz %d\n",
		rad.xStep, rad.zStep, rad.xStart, xEnd, rad.zStart, zEnd, rad.nx, rad.nz);
}

void treat_phase_shift(srTSRWRadStructAccessData* pRadAccessData, double phase)
{
	float* pEx0 = pRadAccessData->pBaseRadX;
	float* pEz0 = pRadAccessData->pBaseRadZ;
	//long PerX = pRadAccessData->ne << 1;
	//long PerZ = PerX*pRadAccessData->nx;
	long long PerX = pRadAccessData->ne << 1;
	long long PerZ = PerX * pRadAccessData->nx;

	// srTGenOptElem::TreatPhaseShift(srTEFieldPtrs& EPtrs, double PhShift)
	for (int iz = 0; iz < pRadAccessData->nz; iz++)
	{
		srTEFieldPtrs EFieldPtrs;
		long long izPerZ = iz * PerZ;
		srTEXZ EXZ;

		EXZ.z = (pRadAccessData->zStart) + iz * (pRadAccessData->zStep);

		float* pEx_StartForX = pEx0 + izPerZ;
		float* pEz_StartForX = pEz0 + izPerZ;
		EXZ.x = pRadAccessData->xStart;
		long long ixPerX = 0;

		for (int ix = 0; ix < pRadAccessData->nx; ix++)	{
			float* pEx_StartForE = pEx_StartForX + ixPerX;
			float* pEz_StartForE = pEz_StartForX + ixPerX;
			EXZ.e = pRadAccessData->eStart;
			long long iePerE = 0;

			for (int ie = 0; ie < pRadAccessData->ne; ie++)	{
				EFieldPtrs.pExRe = pEx_StartForE + iePerE;
				EFieldPtrs.pExIm = EFieldPtrs.pExRe + 1;
				EFieldPtrs.pEzRe = pEz_StartForE + iePerE;
				EFieldPtrs.pEzIm = EFieldPtrs.pEzRe + 1;
				
				EXZ.aux_offset = izPerZ + ixPerX + iePerE;
				// RadPointModifier(EXZ, EFieldPtrs, pBufVars); //OC29082019
				//RadPointModifier(EXZ, EFieldPtrs);
				double twoPi_d_Lambda = (5.06773065e+06) * EXZ.e;
				srTGenOptElem::TreatPhaseShift(EFieldPtrs, phase * twoPi_d_Lambda);

				iePerE += 2;
				EXZ.e += pRadAccessData->eStep;
			}
			ixPerX += PerX;
			EXZ.x += pRadAccessData->xStep;
		}
	}
}

int srTCombinedDrift::PropagateRad1(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!!] The ZP (nzdiv=" << nzdiv << ",nxdiv=" << nxdiv
		<< ") with aperture and drift L= " << dftLen << endl
		<< "  inp xStart= " << pRadAccessData->xStart << " xStep= " << pRadAccessData->xStep
		<< "  obs grid half-width fac= " << obsgrid[0] << " step fac= " << obsgrid[1] << endl;

	assert(dftLen > 0); // NOTE: the shift, then kick (as a inverse of shift) depends on dftLen
						// a dftLen==0 will make this inverse impossible

#if DEBUG_ZPD > 1
	cout << "zp input hash= " << std::hex << pRadAccessData->hashcode() << endl;
#endif
	bool kick_then_shift = true; // do tilt method 
	if (nxdiv <= 2 && nzdiv <= 2 && crsz[0] == 5) { // 2x2 and method==5, do padding zero method
		kick_then_shift = false;
	}

	if (crsz[0] == 5 && (nxdiv != 2 || nzdiv != 2)) {
		cerr << "mode==5 but given grid not 2x2: " << nzdiv << "x" << nxdiv << endl;
		return -1;
	}
	
	

#if DEBUG_ZPD > 1
	string fname = "junk.zpd00.bin";
	pRadAccessData->dumpBinData(fname, "zpd0");
#endif

  // sleep(30);
  print_usage("init");

	double xStep = pRadAccessData->xStep;
	double xStart = pRadAccessData->xStart;
	double zStep = pRadAccessData->zStep;
	double zStart = pRadAccessData->zStart;
	//pRadAccessData->BaseRadWasEmulated == false !!

	srTDriftSpace internal_drift(dftLen);
	// internal_drift.LocalPropMode = 1;

#if DEBUG_ZPD > 0
	fprintf(stderr, "input xStep= %g zStep= %g\n", pRadAccessData->xStep, pRadAccessData->zStep);
	srTSRWRadStructAccessData newRad0(pRadAccessData);
	/*
	{
		srTRadResize resz;

		resz.pxm = resz.pzm = 1;
		resz.pxd = resz.pzd = 1.4;
		RadResizeGen(newRad0, resz);
	}
	*/
	elem->PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref zp  hash= 0x%016zx\n", newRad0.hashcode());
	//ParPrecWfrPropag.MethNo = 0; // crsz[0];
	//ParPrecWfrPropag.AnalTreatment = 4; // LocalPropMode -> 1
	internal_drift.PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref drf hash= 0x%016zx\n", newRad0.hashcode()); 

	newRad0.dumpBinData("junk.zpd.ref.bin", "zpd_reference");
	fprintf(stderr, "dumpped reference output of ZPD nz,nx,ne=(%d,%d,%d) zStep= %g xStep= %g meth= %d\n",
		newRad0.nz, newRad0.nx, newRad0.ne, newRad0.zStep, newRad0.xStep, ParPrecWfrPropag.MethNo);
#endif

	// elem->PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);

#if DEBUG_ZPD > 1
	ofstream junkfdiv("junk.main.txt");
	junkfdiv << "#nzdiv,nxdiv " << nzdiv << " " << nxdiv << endl;
#endif

	// original center
	double xc0 = pRadAccessData->xStart + pRadAccessData->xStep * pRadAccessData->nx / 2.0;
	double zc0 = pRadAccessData->zStart + pRadAccessData->zStep * pRadAccessData->nz / 2.0;
	// average length of Width and Height
	double wdavg = (pRadAccessData->xStep * pRadAccessData->nx / 2.0 + pRadAccessData->zStep * pRadAccessData->nz / 2.0) / 2.0;
	
  // sleep(30);
  print_usage("before_init_destRad");

	// WARNING: assuming destRad is large enough to hold each cell
	srTSRWRadStructAccessData destRad(pRadAccessData);
  init_dest_rad(destRad, pRadAccessData, 1);

  // sleep(30);
  print_usage("after_init_destRad");
  if (false) {
    srTRadResize resz;
    resz.pxm = obsgrid[0]; resz.pzm = obsgrid[2];
    resz.pxd = obsgrid[1]; resz.pzd = obsgrid[3];
    fprintf(stderr, "resize accumulated field: resz.pxm= %f, resz.pxd=%f\n", resz.pxm, resz.pxd);
    RadResizeGen(destRad, resz);
    fprintf(stderr, "accum xStart= %g xStep= %g\n", destRad.xStart, destRad.xStep);
  }
  print_usage("resized_destRad");
  // sleep(30);

#if DEBUG_ZPD > 2
	srTSRWRadStructAccessData radAfterZP(&destRad);
#endif

	const long long RADSZ2 = destRad.nz * destRad.nx * destRad.ne * 2;
	//memset(destRad.pBaseRadX, 0, RADSZ2 * sizeof destRad.pBaseRadX[0]);
	//memset(destRad.pBaseRadZ, 0, RADSZ2 * sizeof destRad.pBaseRadZ[0]);
	fprintf(stderr, "dest dim (%d %d) RADSZ2= %d kick_then_shift=%d\n", destRad.nz, destRad.nx, RADSZ2, kick_then_shift);

	fprintf(stderr, "z in %d [%g %g] step= %g\n", pRadAccessData->nz,
		pRadAccessData->zStart, pRadAccessData->zStart + pRadAccessData->zStep * (pRadAccessData->nz - 1),
		pRadAccessData->zStep);
	fprintf(stderr, "x in %d [%g %g] step= %g\n", pRadAccessData->nx,
		pRadAccessData->xStart, pRadAccessData->xStart + pRadAccessData->xStep * (pRadAccessData->nx - 1),
		pRadAccessData->xStep);

	int zidxdiv[33] = { 0 }, xidxdiv[33] = { 0 };
	for (int ix = 0; ix < nxdiv; ++ix) {
		xidxdiv[ix + 1] = ceil(xdivs[ix] * pRadAccessData->nx);
		fprintf(stderr, "%d (f=%.4f) ", xidxdiv[ix + 1], xdivs[ix]);
	}
	fprintf(stderr, "\n");
	
	for (int iz = 0; iz < nzdiv; ++iz) {
		zidxdiv[iz + 1] = ceil(zdivs[iz] * pRadAccessData->nz);
		fprintf(stderr, "%d ", zidxdiv[iz + 1]);
	}
	fprintf(stderr, "\n");

	for (int iz = 0; iz < nzdiv; ++iz) {
		const int iz0 = zidxdiv[iz], iz1 = zidxdiv[iz + 1];

		for (int ix = 0; ix < nxdiv; ++ix) {
			int ix0 = xidxdiv[ix], ix1 = xidxdiv[ix+1];
			const int k = ix * nzdiv + iz;
			double pmx = crsz[5 * k + 1], pmz = crsz[5 * k + 3];
			const double pdx = crsz[5*k+2];
			const double pdz = crsz[5*k+4];
		
#if DEBUG_ZPD > 1
			fprintf(stderr, "## split_%d_%d pd= %f %f z [%.4e,%.4e] x [%g %g]  (%d,%d)\n", ix, iz, pdx, pdz,
				pRadAccessData->zStart + iz0 * pRadAccessData->zStep,
				pRadAccessData->zStart + (iz1-1) * pRadAccessData->zStep, 
				pRadAccessData->xStart + ix0 * pRadAccessData->xStep, 
				pRadAccessData->xStart + (ix1-1) * pRadAccessData->xStep,
				iz1 - iz0, ix1 - ix0);
#endif
			
			srTSRWRadStructAccessData newRad(pRadAccessData);// &destRad);
			// newRad.BaseRadWasEmulated = false;
			const long long RADSZ = newRad.nz * newRad.nx * newRad.ne * 2;
			assert((iz1 - iz0) * (ix1 - ix0) * pRadAccessData->ne * 2 <= RADSZ);
			memset(newRad.pBaseRadX, 0, RADSZ * sizeof newRad.pBaseRadX[0]);
			memset(newRad.pBaseRadZ, 0, RADSZ * sizeof newRad.pBaseRadZ[0]);

      print_usage("new_subset");

			/*if (nzdiv != 1 || nxdiv != 1) {
				int npad = 2;
				if (ix == 0) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else if (ix == 2*SZX) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1, 1);
			}*/

			int xlpad = 8, xrpad = 8, zlpad = 8, zrpad = 8;
			if (!kick_then_shift) { // make sure (0.0, 0.0) is in this selection
				xlpad = max(10.0, (pRadAccessData->xStart + ix0 * pRadAccessData->xStep) / pRadAccessData->xStep + 8);
				zlpad = max(10.0, (pRadAccessData->zStart + iz0 * pRadAccessData->zStep) / pRadAccessData->zStep + 8);
				xrpad = max(10.0, -(pRadAccessData->xStart + (ix1-1) * pRadAccessData->xStep) / pRadAccessData->xStep + 8);
				zrpad = max(10.0, -(pRadAccessData->zStart + (iz1-1) * pRadAccessData->zStep) / pRadAccessData->zStep + 8);
				fprintf(stderr, "padding with z %d %d x %d %d\n", zlpad, zrpad, xlpad, xrpad);
			}
			sel_sub_cell(&newRad, *pRadAccessData, iz0, iz1, ix0, ix1, zlpad, zrpad, xlpad, xrpad);

      print_usage("sel_sub_cell");

			if (!kick_then_shift) {// make sure (0.0, 0.0) is in this selection
				assert(newRad.xStart < 0.0);
				assert(newRad.xStart + (newRad.nx - 1) * newRad.xStep > 0.0);
				assert(newRad.zStart < 0.0);
				assert(newRad.zStart + (newRad.nz - 1) * newRad.zStep > 0.0);
			}
#if DEBUG_ZPD > 10
			fname = "junk.zpd10." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz << " " << ix << " z " << iz0 << " " << iz1 << " x " << ix0 << " " << ix1 << " " << fname << endl;
#endif
			srTRadResize resz;

			resz.pxm = pmx; resz.pzm = pmz;
			resz.pxd = pdx; resz.pzd = pdz;
			fprintf(stderr, "resz.pmd= %f, resz.pmd=%f\n", pmx, pmz);
			RadResizeGen(newRad, resz);

      print_usage("after_resize_newRad");

#if DEBUG_ZPD > 1
			junkfdiv << "#pdizix " << iz << " " << ix << " pzd " << resz.pzd << " pxd " << resz.pxd << endl;
#endif

			if (int result = elem->PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}

#if DEBUG_ZPD > 10
			fname = "junk.zpd11." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz << " " << ix << " z " << iz0 << " " << iz1 << " x " << ix0 << " " << ix1 << " " << fname << endl;

			fprintf(stderr, "zpd11_%d_%d hash= 0x%016zx Par method= %d %d\n", iz, ix, newRad.hashcode(), ParPrecWfrPropag.MethNo, ParPrecWfrPropag.AnalTreatment);
#endif

			// the real center
			double xc = newRad.xStart + newRad.xStep * (newRad.nx - 1) / 2.0;
			double zc = newRad.zStart + newRad.zStep * (newRad.nz - 1) / 2.0;
			
			xc = newRad.xStart + 0.5 * newRad.xStep * (newRad.nx-1);
			zc = newRad.zStart + 0.5 * newRad.zStep * (newRad.nz-1);
			
			/*if (kick_then_shift && (nxdiv > 1 || nzdiv > 1)) {
				//newRad.xStart -= xc;
				//newRad.zStart -= zc;
			}*/

			//double ang_x = atan(xc / dftLen);
			//double ang_z = atan(zc / dftLen);

			
#if DEBUG_ZPD > 1

			fprintf(stderr, "xc=%g zc=%g ang_x=%g ang_z=%g\n", xc, zc, atan(xc/dftLen), atan(zc/dftLen));
			fprintf(stderr, "Slice iz=%d ix=%d z=[%g %g], x=[%g %g] nz= %d nx= %d\n",
				iz, ix, newRad.zStart, newRad.zStart + newRad.nz * newRad.zStep, newRad.xStart, newRad.xStart + newRad.nx * newRad.xStep,
				newRad.nz, newRad.nx);

#endif
			/*
			if (kick_then_shift && dftLen > 0 && (nzdiv > 1 || nxdiv > 1)) {
				srTOptAngle inter_angle(ang_x, ang_z);
				fprintf(stderr, "adjust angle %f %f\n", ang_x, ang_z);
				if (int result = inter_angle.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
					fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
					return result;
				}
			}
			*/

#if DEBUG_ZPD > 15
			// CDRadStructHelper::add(&radAfterZP, &newRad);

			fname = "junk.zpd12." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);

			fprintf(stderr, "zpd zp  hash= 0x%016zx\n", newRad.hashcode());
#endif

			srTDriftSpace internal_drift(dftLen);
			if (kick_then_shift) {
				internal_drift.shift_obsx = -xc;
				internal_drift.shift_obsz = -zc;
			}
			// int Ann = ParPrecWfrPropag.AnalTreatment;
			// ParPrecWfrPropag.AnalTreatment = 4;
			if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}
			// ParPrecWfrPropag.AnalTreatment = Ann; // restore


#if DEBUG_ZPD > 15
			fname = "junk.zpd13." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
#endif

			// treat_phase_shift(&newRad, (xc*xc + zc*zc) / dftLen);
			
#if DEBUG_ZPD > 2
			fprintf(stderr, "zpd dft  hash= 0x%016zx\n", newRad.hashcode());
#endif
			
			if (kick_then_shift) {
				newRad.xStart -= xc;
				newRad.zStart -= zc;
			}
			

#if DEBUG_ZPD > 1
			{
				//srTRadResize resz;

				//resz.pxm = 0.05; resz.pzm = 0.05;
				//resz.pxd = 10; resz.pzd = 10;
				//RadResizeGen(newRad, resz);

				fname = "junk.zpd2." + to_string(iz) + "_" + to_string(ix) + ".bin";
				newRad.dumpBinDataCore(fname, fname, 1e-7, 1e-7);
				junkfdiv << "#fout " << iz << " " << ix << " " << fname << endl;
				{
					const auto itm = std::minmax_element(destRad.pBaseRadX, destRad.pBaseRadX + RADSZ2);
					fprintf(stderr, "min= %g max %g\n", *itm.first, *itm.second);
				}
				fprintf(stderr, "Slice (iz, ix) = %d %d Hash= 0x%016zx, Hash= 0x%016zx\n", iz, ix, pRadAccessData->hashcode(), newRad.hashcode());
			}
#endif

			{
				//debug
				/*
				int k = 0;
				for (int i = 0; i < newRad.nz; ++i) {
					for (int j = 0; j < newRad.nx * newRad.ne * 2; ++j, ++k) {
						if (fabs(newRad.pBaseRadX[k]) > 1e6 || fabs(newRad.pBaseRadZ[k]) > 1e6)
							fprintf(stderr, "(%d %d k=%d)  Ux= %g Uz= %g ", i, j, k, newRad.pBaseRadX[k], newRad.pBaseRadZ[k]);
					}
				}
				fprintf(stderr, "\n");
				*/
			}

			

			// accumulate_rad(destRad.pBaseRadX, destRad.pBaseRadZ, pRadAccessData->nx, pRadAccessData->nz, xStart, xStep, zStart, zStep, newRad);
      if (newRad.xStep > destRad.xStep || newRad.zStep > destRad.zStep) {
        fprintf(stderr, "WARNING: propagated slice has larger grid %g %g than observation point (%g %g), the accumulated field may have flat top\n", newRad.xStep, newRad.zStep, destRad.xStep, destRad.zStep);
      }
			CDRadStructHelper::add(&destRad, &newRad);

      print_usage("merged_to_destRad");

#if DEBUG_ZPD > 2
			{
				// dump the accumulated
				fname = "junk.zpd2.sum." + to_string(iz) + "_" + to_string(ix) + ".bin";
				ofstream out(fname, ios::out | ios::binary);
				out.write((char*)&destRad.nz, sizeof destRad.nz);
				out.write((char*)&destRad.nx, sizeof destRad.nx);
				out.write((char*)&destRad.ne, sizeof destRad.ne);
				out.write((char*)&RADSZ2, sizeof RADSZ2);
				fprintf(stderr, "write accum field %d %d %d (%d) RADSZ=%d (%d)\n",
					destRad.nz, destRad.nx, destRad.ne, sizeof destRad.nz, RADSZ2, sizeof RADSZ2);
				
				out.write((char*)&destRad.zStart, sizeof destRad.zStart);
				out.write((char*)&destRad.zStep, sizeof destRad.zStep);
				out.write((char*)&destRad.xStart, sizeof destRad.xStart);
				out.write((char*)&destRad.xStep, sizeof destRad.xStep);

				out.write((char*)destRad.pBaseRadZ, RADSZ2 * sizeof(float));
				out.write((char*)destRad.pBaseRadX, RADSZ2 * sizeof(float));
			}
#endif
		}
	}

  fprintf(stderr, "destRad: x %g %g z %g %g nx= %d nz= %d\n",
          pRadAccessData->xStart, pRadAccessData->xStep, pRadAccessData->zStart, pRadAccessData->zStep,
          pRadAccessData->nx, pRadAccessData->nz);

  /*
  srTRadResize resz;
  resz.pxm = obsgrid[0]; resz.pzm = obsgrid[2];
  resz.pxd = obsgrid[1]; resz.pzd = obsgrid[3];
  fprintf(stderr, "resize accumulated field: resz.pxm= %f, resz.pxd=%f\n", resz.pxm, resz.pxd);
  RadResizeGen(destRad, resz);
  */

  print_usage("before_copy_to_pRadAcc");
	pRadAccessData->nx = destRad.nx;
	pRadAccessData->nz = destRad.nz;
	pRadAccessData->zStart = destRad.zStart;
	pRadAccessData->zStep = destRad.zStep;
	pRadAccessData->xStart = destRad.xStart;
	pRadAccessData->xStep = destRad.xStep;
	// ERROR: memory leak, but otherwise the delete[] in ReAllocBaseRadAccordingToNeNxNz gives me coredump.
	pRadAccessData->pBaseRadX = nullptr;
	pRadAccessData->pBaseRadZ = nullptr;
	// pRadAccessData->ReAllocBaseRadAccordingToNeNxNz();
	pRadAccessData->ModifyWfrNeNxNz();

  print_usage("pRadAcc_resized");


	CDRadStructHelper::assign(pRadAccessData, &destRad);
	fprintf(stderr, "%g %g %g %g\n", pRadAccessData->xStart, pRadAccessData->xStep, pRadAccessData->zStart, pRadAccessData->zStep);
	//free(radx);
	//free(radz);
#if DEBUG_ZPD > 2
	radAfterZP.dumpBinData("junk.zpd.afterzp.bin", "zpd.afterzp.bin");
#endif

  print_usage("after_destRad_to_pRadAcc");
	
#if DEBUG_ZPD > 0
	fprintf(stderr, "DONE ZPD div=(%d %d) nz=%d nx= %d\n", nzdiv, nxdiv, pRadAccessData->nz, pRadAccessData->nx);
	fflush(stderr);
	// pRadAccessData->dumpBinData("junk.zpd.bin", "junk.zpd.bin");
	destRad.dumpBinData("junk.zpd.bin", "final, but before copy to pRadAccessData");
#endif

#if DEBUG_ZPD > 1
	junkfdiv.close(); 
#endif

	return 0;
}

int srTCombinedDrift::PropagateRad2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	cout << "[WARNING!!] The ZP (nzdiv=" << nzdiv << ",nxdiv=" << nxdiv
		<< ") with aperture and drift L= " << dftLen << endl
		<< "  inp xStart= " << pRadAccessData->xStart << " xStep= " << pRadAccessData->xStep
		<< " nx= " << pRadAccessData->nx
		<< "  zStart = " << pRadAccessData->zStart << " zStep = " << pRadAccessData->zStep
		<< " nz= " << pRadAccessData->nz
		<< "  obs grid half-width fac= " << obsgrid[0] << " step fac= " << obsgrid[1] << endl;

	assert(dftLen > 0); // NOTE: the shift, then kick (as a inverse of shift) depends on dftLen
						// a dftLen==0 will make this inverse impossible

#if DEBUG_ZPD > 1
	cout << "zp input hash= " << std::hex << pRadAccessData->hashcode() << endl;
#endif
	bool kick_then_shift = true; // do tilt method 
	if (nxdiv <= 2 && nzdiv <= 2 && crsz[0] == 5) { // 2x2 and method==5, do padding zero method
		kick_then_shift = false;
	}

	if (crsz[0] == 5 && (nxdiv != 2 || nzdiv != 2)) {
		cerr << "mode==5 but given grid not 2x2: " << nzdiv << "x" << nxdiv << endl;
		return -1;
	}
	
	

#if DEBUG_ZPD > 1
	string fname = "junk.zpd00.bin";
	pRadAccessData->dumpBinData(fname, "zpd0");
#endif

  // sleep(30);
  print_usage("init");

	double xStep = pRadAccessData->xStep;
	double xStart = pRadAccessData->xStart;
	double zStep = pRadAccessData->zStep;
	double zStart = pRadAccessData->zStart;
	//pRadAccessData->BaseRadWasEmulated == false !!

	srTDriftSpace internal_drift(dftLen);
	// internal_drift.LocalPropMode = 1;

#if DEBUG_ZPD > 0
	fprintf(stderr, "input xStep= %g zStep= %g\n", pRadAccessData->xStep, pRadAccessData->zStep);
	srTSRWRadStructAccessData newRad0(pRadAccessData);
	
  elem->PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref zp  hash= 0x%016zx\n", newRad0.hashcode());
	//ParPrecWfrPropag.MethNo = 0; // crsz[0];
	//ParPrecWfrPropag.AnalTreatment = 4; // LocalPropMode -> 1
	internal_drift.PropagateRadiation(&newRad0, ParPrecWfrPropag, ResBeforeAndAfterVect);
	//fprintf(stderr, "ref drf hash= 0x%016zx\n", newRad0.hashcode()); 

	newRad0.dumpBinData("junk.zpd.ref.bin", "zpd_reference");
	fprintf(stderr, "dumpped reference output of ZPD nz,nx,ne=(%d,%d,%d) zStep= %g xStep= %g meth= %d\n",
		newRad0.nz, newRad0.nx, newRad0.ne, newRad0.zStep, newRad0.xStep, ParPrecWfrPropag.MethNo);
#endif

	// elem->PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);

#if DEBUG_ZPD > 1
	ofstream junkfdiv("junk.main.txt");
	junkfdiv << "#nzdiv,nxdiv " << nzdiv << " " << nxdiv << endl;
#endif

	// original center
	double xc0 = pRadAccessData->xStart + pRadAccessData->xStep * pRadAccessData->nx / 2.0;
	double zc0 = pRadAccessData->zStart + pRadAccessData->zStep * pRadAccessData->nz / 2.0;
	// average length of Width and Height
	double wdavg = (pRadAccessData->xStep * pRadAccessData->nx / 2.0 + pRadAccessData->zStep * pRadAccessData->nz / 2.0) / 2.0;
	
  // sleep(30);
  print_usage("before_init_destRad");

	// WARNING: assuming destRad is large enough to hold each cell
	srTSRWRadStructAccessData destRad(pRadAccessData);
  init_dest_rad(destRad, pRadAccessData, 1);

  // sleep(30);
  print_usage("after_init_destRad");
  if (false) {
    srTRadResize resz;
    resz.pxm = obsgrid[0]; resz.pzm = obsgrid[2];
    resz.pxd = obsgrid[1]; resz.pzd = obsgrid[3];
    fprintf(stderr, "resize accumulated field: resz.pxm= %f, resz.pxd=%f\n", resz.pxm, resz.pxd);
    RadResizeGen(destRad, resz);
    fprintf(stderr, "accum xStart= %g xStep= %g\n", destRad.xStart, destRad.xStep);
  }
  print_usage("resized_destRad");
  // sleep(30);

#if DEBUG_ZPD > 2
	srTSRWRadStructAccessData radAfterZP(&destRad);
#endif

	const long long RADSZ2 = destRad.nz * destRad.nx * destRad.ne * 2;
	//memset(destRad.pBaseRadX, 0, RADSZ2 * sizeof destRad.pBaseRadX[0]);
	//memset(destRad.pBaseRadZ, 0, RADSZ2 * sizeof destRad.pBaseRadZ[0]);
	fprintf(stderr, "dest dim (%d %d) RADSZ2= %d kick_then_shift=%d\n", destRad.nz, destRad.nx, RADSZ2, kick_then_shift);

	fprintf(stderr, "z in %d [%g %g] step= %g\n", pRadAccessData->nz,
		pRadAccessData->zStart, pRadAccessData->zStart + pRadAccessData->zStep * (pRadAccessData->nz - 1),
		pRadAccessData->zStep);
	fprintf(stderr, "x in %d [%g %g] step= %g\n", pRadAccessData->nx,
		pRadAccessData->xStart, pRadAccessData->xStart + pRadAccessData->xStep * (pRadAccessData->nx - 1),
		pRadAccessData->xStep);

	int zidxdiv[33] = { 0 }, xidxdiv[33] = { 0 };
	for (int ix = 0; ix < nxdiv; ++ix) {
		xidxdiv[ix + 1] = ceil(xdivs[ix] * pRadAccessData->nx);
		fprintf(stderr, "%d (f=%.4f) ", xidxdiv[ix + 1], xdivs[ix]);
	}
	fprintf(stderr, "\n");
	
	for (int iz = 0; iz < nzdiv; ++iz) {
		zidxdiv[iz + 1] = ceil(zdivs[iz] * pRadAccessData->nz);
		fprintf(stderr, "%d ", zidxdiv[iz + 1]);
	}
	fprintf(stderr, "\n");

	for (int iz = 0; iz < nzdiv; ++iz) {
		const int iz0 = zidxdiv[iz], iz1 = zidxdiv[iz + 1];
		double zStart = pRadAccessData->zStart + iz0 * pRadAccessData->zStep;
		double zEnd = zStart + (iz1 - iz0) * pRadAccessData->zStep;
		
		for (int ix = 0; ix < nxdiv; ++ix) {
			int ix0 = xidxdiv[ix], ix1 = xidxdiv[ix+1];
			const int k = ix * nzdiv + iz;
			double pmx = crsz[5 * k + 1], pmz = crsz[5 * k + 3];
			const double pdx = crsz[5*k+2];
			const double pdz = crsz[5*k+4];
		
#if DEBUG_ZPD > 1
			fprintf(stderr, "## split_%d_%d pd= %f %f z [%.4e,%.4e] x [%g %g]  (%d,%d)\n", ix, iz, pdx, pdz,
				pRadAccessData->zStart + iz0 * pRadAccessData->zStep,
				pRadAccessData->zStart + (iz1-1) * pRadAccessData->zStep, 
				pRadAccessData->xStart + ix0 * pRadAccessData->xStep, 
				pRadAccessData->xStart + (ix1-1) * pRadAccessData->xStep,
				iz1 - iz0, ix1 - ix0);
#endif
			
			double xStart = pRadAccessData->xStart + ix0 * pRadAccessData->xStep;
			double xEnd = xStart + (ix1 - ix0) * pRadAccessData->xStep;

			int nzpad = max(0, (zEnd - zStart) * (pmz - 1) / pRadAccessData->zStep * pdz / 2);
			int nxpad = max(0, (xEnd - xStart) * (pmx - 1) / pRadAccessData->xStep * pdx / 2);

			srTSRWRadStructAccessData newRad(pRadAccessData);

      CDRadStructHelper::sample(&newRad, pRadAccessData,
		  xStart, xEnd, pRadAccessData->xStep / pdx, nxpad,
		  zStart, zEnd, pRadAccessData->zStep / pdz, nzpad);

      print_usage("new_subset");

			/*if (nzdiv != 1 || nxdiv != 1) {
				int npad = 2;
				if (ix == 0) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else if (ix == 2*SZX) sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1750, 1750);
				else sel_sub_cell(&newRad, *pRadAccessData, iz, iz + SZZ, ix, ix + SZX, 1, 1, 1, 1);
			}*/

#if DEBUG_ZPD > 5
			fname = "junk.zpd10." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz << " " << ix << " z " << iz0 << " " << iz1 << " x " << ix0 << " " << ix1 << " " << fname << endl;
#endif


			if (int result = elem->PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}

#if DEBUG_ZPD > 5
			fname = "junk.zpd11." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
			junkfdiv << "#fin " << iz << " " << ix << " z " << iz0 << " " << iz1 << " x " << ix0 << " " << ix1 << " " << fname << endl;

			fprintf(stderr, "zpd11_%d_%d hash= 0x%016zx Par method= %d %d\n", iz, ix, newRad.hashcode(), ParPrecWfrPropag.MethNo, ParPrecWfrPropag.AnalTreatment);
#endif

			// the real center
			double xc = newRad.xStart + newRad.xStep * (newRad.nx - 1) / 2.0;
			double zc = newRad.zStart + newRad.zStep * (newRad.nz - 1) / 2.0;
			
			xc = newRad.xStart + 0.5 * newRad.xStep * (newRad.nx-1);
			zc = newRad.zStart + 0.5 * newRad.zStep * (newRad.nz-1);
			
			
#if DEBUG_ZPD > 1
			fprintf(stderr, "xc=%g zc=%g ang_x=%g ang_z=%g\n", xc, zc, atan(xc/dftLen), atan(zc/dftLen));
			fprintf(stderr, "Slice iz=%d ix=%d z=[%g %g], x=[%g %g] nz= %d nx= %d\n",
				iz, ix, newRad.zStart, newRad.zStart + newRad.nz * newRad.zStep, newRad.xStart, newRad.xStart + newRad.nx * newRad.xStep,
				newRad.nz, newRad.nx);
#endif

#if DEBUG_ZPD > 5
			fname = "junk.zpd12." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);

			fprintf(stderr, "zpd zp  hash= 0x%016zx\n", newRad.hashcode());
#endif

			srTDriftSpace internal_drift(dftLen);
			if (kick_then_shift) {
				internal_drift.shift_obsx = -xc;
				internal_drift.shift_obsz = -zc;
			}

      if (int result = internal_drift.PropagateRadiation(&newRad, ParPrecWfrPropag, ResBeforeAndAfterVect)) {
				fprintf(stderr, "ERROR %d: %s", result, __FUNCTION__);
				return result;
			}


#if DEBUG_ZPD > 5
			fname = "junk.zpd13." + to_string(iz) + "_" + to_string(ix) + ".bin";
			newRad.dumpBinData(fname, fname);
#endif

			// treat_phase_shift(&newRad, (xc*xc + zc*zc) / dftLen);
			
#if DEBUG_ZPD > 2
			fprintf(stderr, "zpd dft  hash= 0x%016zx\n", newRad.hashcode());
#endif
			
			if (kick_then_shift) {
				newRad.xStart -= xc;
				newRad.zStart -= zc;
			}
			

#if DEBUG_ZPD > 1
			{
				//srTRadResize resz;

				//resz.pxm = 0.05; resz.pzm = 0.05;
				//resz.pxd = 10; resz.pzd = 10;
				//RadResizeGen(newRad, resz);

				fname = "junk.zpd2.core." + to_string(iz) + "_" + to_string(ix) + ".bin";
				// newRad.dumpBinDataCore(fname, fname, newRad.xStep*50, newRad.zStep*50);
				fname = "junk.zpd2." + to_string(iz) + "_" + to_string(ix) + ".bin";
				newRad.dumpBinData(fname, fname);
				junkfdiv << "#fout " << iz << " " << ix << " " << fname << endl;
				{
					const auto itm = std::minmax_element(destRad.pBaseRadX, destRad.pBaseRadX + RADSZ2);
					fprintf(stderr, "min= %g max %g\n", *itm.first, *itm.second);
				}
				// fprintf(stderr, "Slice (iz, ix) = %d %d Hash= 0x%016zx, Hash= 0x%016zx\n", iz, ix, pRadAccessData->hashcode(), newRad.hashcode());
			}
#endif

			{
				//debug
				/*
				int k = 0;
				for (int i = 0; i < newRad.nz; ++i) {
					for (int j = 0; j < newRad.nx * newRad.ne * 2; ++j, ++k) {
						if (fabs(newRad.pBaseRadX[k]) > 1e6 || fabs(newRad.pBaseRadZ[k]) > 1e6)
							fprintf(stderr, "(%d %d k=%d)  Ux= %g Uz= %g ", i, j, k, newRad.pBaseRadX[k], newRad.pBaseRadZ[k]);
					}
				}
				fprintf(stderr, "\n");
				*/
			}

			

			// accumulate_rad(destRad.pBaseRadX, destRad.pBaseRadZ, pRadAccessData->nx, pRadAccessData->nz, xStart, xStep, zStart, zStep, newRad);
      if (newRad.xStep > destRad.xStep || newRad.zStep > destRad.zStep) {
        fprintf(stderr, "WARNING: propagated slice has larger grid %g %g than observation point (%g %g), the accumulated field may have flat top\n", newRad.xStep, newRad.zStep, destRad.xStep, destRad.zStep);
      }
			CDRadStructHelper::add(&destRad, &newRad);

      print_usage("merged_to_destRad");

#if DEBUG_ZPD > 2
			{
				// dump the accumulated
				fname = "junk.zpd2.sum." + to_string(iz) + "_" + to_string(ix) + ".bin";
				ofstream out(fname, ios::out | ios::binary);
				out.write((char*)&destRad.nz, sizeof destRad.nz);
				out.write((char*)&destRad.nx, sizeof destRad.nx);
				out.write((char*)&destRad.ne, sizeof destRad.ne);
				out.write((char*)&RADSZ2, sizeof RADSZ2);
				fprintf(stderr, "write accum field %d %d %d (%d) RADSZ=%d (%d)\n",
					destRad.nz, destRad.nx, destRad.ne, sizeof destRad.nz, RADSZ2, sizeof RADSZ2);
				
				out.write((char*)&destRad.zStart, sizeof destRad.zStart);
				out.write((char*)&destRad.zStep, sizeof destRad.zStep);
				out.write((char*)&destRad.xStart, sizeof destRad.xStart);
				out.write((char*)&destRad.xStep, sizeof destRad.xStep);

				out.write((char*)destRad.pBaseRadZ, RADSZ2 * sizeof(float));
				out.write((char*)destRad.pBaseRadX, RADSZ2 * sizeof(float));
				out.close();

				fname = "junk.zpd3.sum." + to_string(iz) + "_" + to_string(ix) + ".bin";
				destRad.dumpBinData(fname, fname);
			}
#endif
		}
	}

  fprintf(stderr, "destRad: x %g %g z %g %g nx= %d nz= %d\n",
          pRadAccessData->xStart, pRadAccessData->xStep, pRadAccessData->zStart, pRadAccessData->zStep,
          pRadAccessData->nx, pRadAccessData->nz);

  /*
  srTRadResize resz;
  resz.pxm = obsgrid[0]; resz.pzm = obsgrid[2];
  resz.pxd = obsgrid[1]; resz.pzd = obsgrid[3];
  fprintf(stderr, "resize accumulated field: resz.pxm= %f, resz.pxd=%f\n", resz.pxm, resz.pxd);
  RadResizeGen(destRad, resz);
  */

  print_usage("before_copy_to_pRadAcc");
	pRadAccessData->nx = destRad.nx;
	pRadAccessData->nz = destRad.nz;
	pRadAccessData->zStart = destRad.zStart;
	pRadAccessData->zStep = destRad.zStep;
	pRadAccessData->xStart = destRad.xStart;
	pRadAccessData->xStep = destRad.xStep;
	// ERROR: memory leak, but otherwise the delete[] in ReAllocBaseRadAccordingToNeNxNz gives me coredump.
	pRadAccessData->pBaseRadX = nullptr;
	pRadAccessData->pBaseRadZ = nullptr;
	// pRadAccessData->ReAllocBaseRadAccordingToNeNxNz();
	pRadAccessData->ModifyWfrNeNxNz();

  print_usage("pRadAcc_resized");

	CDRadStructHelper::assign(pRadAccessData, &destRad);
	fprintf(stderr, "assigned destRad back to pRadAccessData: %g %g %g %g\n", pRadAccessData->xStart, pRadAccessData->xStep, pRadAccessData->zStart, pRadAccessData->zStep);
	//free(radx);
	//free(radz);
#if DEBUG_ZPD > 2
	// destRad.dumpBinData("junk.zpd.afterzp.bin", "zpd.afterzp.bin");
#endif

  print_usage("after_destRad_to_pRadAcc");
	
#if DEBUG_ZPD > 0
	fprintf(stderr, "DONE ZPD div=(%d %d) nz=%d nx= %d\n", nzdiv, nxdiv, pRadAccessData->nz, pRadAccessData->nx);
	fflush(stderr);
	// pRadAccessData->dumpBinData("junk.zpd.bin", "junk.zpd.bin");
	destRad.dumpBinData("junk.zpd.bin", "final, but before copy to pRadAccessData");
#endif

#if DEBUG_ZPD > 1
	junkfdiv.close(); 
#endif

	return 0;
}

// the following is related to slice the wavefront according to ring
//int sel_sub_ring(srTSRWRadStructAccessData& dst, /* const */ srTSRWRadStructAccessData* wfr, double zx0, double wzx, int iring, int n, int npad = 0)
//{
//	// skip
//	if (iring == 0 && zx0 > wfr->zStep / 2.0) return 0; // the center but zx0 is not origin
//	if (iring > 0 && zx0 < wfr->zStep) return 0; // not the center but zx0 is too small
//
//	// n is total grid, including npad. i.e. there are (n - 2*npad) non-zero points each dimension
//	// const int npad = 4;
//	// npad.. < 0 means does not count, use whatever left
//	
//	dst.ne = wfr->ne;
//	dst.xStep = dst.zStep = wzx / n;
//	double wz = 0.0, wx = 0.0;
//	fprintf(stderr, "select ring %d outside of (%f, %f) width= %f n= %d\n", iring, zx0, zx0, wzx, n);
//	switch (iring) {
//	case 0: wz = wx = wzx;  break; // fix later
//	case 1: // right block
//		dst.zStart = -zx0; wz = 2 * zx0 + wzx; 
//		dst.xStart = zx0 + dst.xStep; wx = wzx;
//		break;
//	case 2: // bottom block
//		dst.zStart = -zx0 - wzx; wz = wzx;
//		dst.xStart = -zx0; wx = 2 * zx0 + wzx;
//		break;
//	case 3: // left
//		dst.zStart = -zx0 - wzx; wz = 2 * zx0 + wzx;
//		dst.xStart = -zx0 - wzx; wx = wzx;
//		break;
//	case 4: // top
//		dst.xStart = -zx0 - wzx; wz = wzx;
//		dst.zStart = zx0 + dst.zStep; wx = 2 * zx0 + wzx;
//		break;
//	default: exit(-1);
//	}
//	if (iring == 0) {
//		dst.zStep = dst.xStep = wzx / (n - 1);
//		dst.zStart = dst.xStart = -wzx / 2.0;
//		dst.nz = dst.nx = next_fft_num(n+2*npad);
//	}
//	else {
//		dst.zStart -= npad * dst.zStep;
//		dst.xStart -= npad * dst.xStep;
//		long nx0 = round(wx / dst.xStep) + 1 + npad * 2;
//		long nz0 = round(wz / dst.zStep) + 1 + npad * 2;
//		dst.nz = next_fft_num(nz0);
//		dst.nx = next_fft_num(nx0);
//		fprintf(stderr, "FFT num iring %d, nx %d +%d, nz %d +%d\n", iring, nx0, dst.nx - nx0, nz0, dst.nz - nz0);
//	}
//	
//	fprintf(stderr, "ring size: [%g %g], [%g %g]\n",
//		dst.zStart + npad * dst.zStep, dst.zStart + npad * dst.zStep + wz,
//		dst.xStart + npad * dst.xStep, dst.xStart + npad * dst.xStep + wx);
//	
//	
//	dst.xWfrMin = dst.xStart;
//	dst.zWfrMin = dst.zStart;
//	dst.xWfrMax = dst.xStart + dst.nx * dst.xStep;
//	dst.zWfrMax = dst.zStart + dst.nz * dst.zStep;
//
//	// sampling
//	float* rdestz = dst.pBaseRadZ, *rdestx = dst.pBaseRadX;
//	for (int iz = 0; iz < dst.nz; ++iz) {
//		double z = dst.zStart + iz * dst.zStep;
//		if (z - dst.zStart  < npad * dst.zStep || z - dst.zStart > wz) continue;
//		for (int ix = 0; ix < dst.nx; ++ix) {
//			double x = dst.xStart + ix * dst.xStep;
//			if (x - dst.xStart  < npad * dst.xStep || x - dst.xStart > wx) continue;
//			do_efield_2d_sample(rdestx + iz * dst.nx * dst.ne * 2 + ix * dst.ne*2, z, x, wfr->pBaseRadX, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
//			do_efield_2d_sample(rdestz + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadZ, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
//			
//			//rdestx += wfr->ne * 2;
//			//rdestz += wfr->ne * 2;
//		}
//		//
//	}
//	return dst.nz * dst.nx;
//}
//
//// zxll0 - lower left z/x coordinate, if <0 it is the center region
//// wzx - width (smaller one for rectangular)
//int sel_sub_ring(srTSRWRadStructAccessData& dst, /* const */ srTSRWRadStructAccessData* wfr, int iring, double zxll0, double wzx, double minstepsz)
//{	
//	// the center block or other 4 blocks
//	assert((iring == 0 && zxll0 < 0) || (iring > 0 && iring <= 4 && zxll0 > 0));
//	dst.ne = wfr->ne;
//	const double GAP = minstepsz;
//	const double WDSHORT = (iring == 0 ? 2*fabs(zxll0) : wzx - GAP);
//	const double WDLONG = (iring == 0 ? 2*fabs(zxll0) : zxll0 * 2 + wzx);
//
//	double wz = 0.0, wx = 0.0;
//	// fprintf(stderr, "select ring %d outside of (%f, %f) width= %f n= %d\n", iring, zx0, zx0, wzx, n);
//	switch (iring) {
//	case 0:
//		dst.xStart = dst.zStart = zxll0;  wz = wx = -2*zxll0;
//		break;
//	case 1: // right block
//		dst.zStart = -zxll0; wz = WDLONG;
//		dst.xStart = zxll0 + GAP; wx = WDSHORT;
//		break;
//	case 2: // bottom block
//		dst.zStart = -zxll0 - wzx; wz = WDSHORT;
//		dst.xStart = -zxll0; wx = WDLONG;
//		break;
//	case 3: // left
//		dst.zStart = -zxll0 - wzx; wz = WDLONG;
//		dst.xStart = -zxll0 - wzx; wx = WDSHORT;
//		break;
//	case 4: // top
//		dst.zStart = zxll0 + GAP; wz = WDSHORT;
//		dst.xStart = -zxll0 - wzx; wx = WDLONG;
//		break;
//	default: exit(-1);
//	}
//	
//	long nx0 = ceil(wx / minstepsz);
//	long nz0 = ceil(wz / minstepsz);
//	dst.nz = next_fft_num(nz0);
//	dst.nx = next_fft_num(nx0);
//	dst.xStep = wx / (dst.nx - 1);
//	dst.zStep = wz / (dst.nz - 1);
//	fprintf(stderr, "FFT num iring %d, nx %d +%d, nz %d +%d\n", iring, nx0, dst.nx - nx0, nz0, dst.nz - nz0);
//	
//
//	fprintf(stderr, "ring size: [%g %g], [%g %g]\n",
//		dst.zStart, dst.zStart + wz, dst.xStart, dst.xStart + wx);
//
//
//	dst.xWfrMin = dst.xStart;
//	dst.zWfrMin = dst.zStart;
//	dst.xWfrMax = dst.xStart + (dst.nx-1) * dst.xStep;
//	dst.zWfrMax = dst.zStart + (dst.nz-1) * dst.zStep;
//
//	// sampling
//	int npad = 0;
//	float* rdestz = dst.pBaseRadZ, * rdestx = dst.pBaseRadX;
//	for (int iz = 0; iz < dst.nz; ++iz) {
//		double z = dst.zStart + iz * dst.zStep;
//		if (z - dst.zStart  < npad * dst.zStep || z - dst.zStart > wz) continue;
//		for (int ix = 0; ix < dst.nx; ++ix) {
//			double x = dst.xStart + ix * dst.xStep;
//			if (x - dst.xStart  < npad * dst.xStep || x - dst.xStart > wx) continue;
//			do_efield_2d_sample(rdestx + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadX, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
//			do_efield_2d_sample(rdestz + iz * dst.nx * dst.ne * 2 + ix * dst.ne * 2, z, x, wfr->pBaseRadZ, wfr->ne, wfr->nz, wfr->nx, wfr->zStart, wfr->zStep, wfr->xStart, wfr->xStep, false /* additive */);
//
//			//rdestx += wfr->ne * 2;
//			//rdestz += wfr->ne * 2;
//		}
//		//
//	}
//	return dst.nz * dst.nx;
//}
//
//// pdedge = 1, keep original stepsize
//// pdcenter < 1, larger stepsize
//// iring = 0, 1, 2, ..., numring - 1
//double sub_ring_stepsize(double stepsz, double pdcenter, double pdedge, int numring, int iring)
//{
//	assert(numring > 0);
//	if (numring <= 1) return stepsz;
//	double pd = pdcenter + iring * (pdedge - pdcenter) / (numring - 1);
//	return stepsz / pd;
//}
// the above is related to slice the wavefront according to ring


int srTCombinedDrift::PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
{
	assert(nzdiv > 0 && nxdiv > 0);
	/*if (crsz[0] == 5) { // method = 2x2 each has padding to include (0.0, 0.0)
		return PropagateRad22(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
		assert(nxdiv == 2 && nzdiv == 2);
		nzdiv = nxdiv = 2; // it has to be 2x2, special case
	}
	*/

	// return PropagateRad1(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
	return PropagateRad2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
}


