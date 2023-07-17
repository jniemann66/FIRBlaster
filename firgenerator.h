#ifndef FIRGENERATOR_H
#define FIRGENERATOR_H


/*
* Copyright (C) 2016 - 2023 Judd Niemann - All Rights Reserved.
* You may use, distribute and modify this code under the
* terms of the GNU Lesser General Public License, version 2.1
*
* You should have received a copy of GNU Lesser General Public License v2.1
* with this file. If not, please refer to: https://github.com/jniemann66/ReSampler
*/

// FIRFilter.h : simple FIR filter implementation



//#include "alignedmalloc.h"

#include "factorial.h"

#include <typeinfo>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <vector>
#include <immintrin.h>

#include <fftw3.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -- Functions beyond this point are for manipulating filter taps, and not for actually performing filtering -- //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// makeLPF() : generate low pass filter coefficients, using sinc function
template<typename FloatType> bool makeLPF(FloatType* filter, int Length, FloatType ft)
{
#ifdef FIR_QUAD_PRECISION

	// use quads internally, regardless of FloatType
	__float128 ft = transitionFreq / sampleRate; // normalised transition frequency
	// assert(ft < 0.5Q);
	int halfLength = Length / 2;
	__float128 halfM = 0.5Q * (Length - 1);
	__float128 M_TWOPIq = 2.0Q * M_PIq;

	if (Length & 1)
		filter[halfLength] = 2.0Q * ft; // if length is odd, avoid divide-by-zero at centre-tap

	for (int n = 0; n<halfLength; ++n) {
		__float128 sinc = sinq(fmodq(M_TWOPIq * ft * (n - halfM), M_TWOPIq)) / (M_PIq * (n - halfM));	// sinc function
		filter[Length - n - 1] = filter[n] = sinc;	// exploit symmetry
	}

#else

	// assert(ft < 0.5);
	int halfLength = Length / 2;
	double halfM = 0.5 * (Length - 1);
	double M_TWOPI = 2.0 * M_PI;

	if (Length & 1)
		filter[halfLength] = 2.0 * ft; // if length is odd, avoid divide-by-zero at centre-tap

	for (int n = 0; n < halfLength; ++n) {
		// sinc function
		double sinc = sin(fmod(M_TWOPI * ft * (n - halfM), M_TWOPI)) / (M_PI * (n - halfM));
		filter[Length - n - 1] = filter[n] = sinc;	// exploit symmetry
	}
#endif

	return true;
}

// makeLPF() : generate high pass filter coefficients, using sinc function
template<typename FloatType> bool makeHPF(FloatType* filter, int Length, FloatType ft)
{
#ifdef FIR_QUAD_PRECISION

	// use quads internally, regardless of FloatType
	__float128 ft = transitionFreq / sampleRate; // normalised transition frequency
	// assert(ft < 0.5Q);
	int halfLength = Length / 2;
	__float128 halfM = 0.5Q * (Length - 1);
	__float128 M_TWOPIq = 2.0Q * M_PIq;

	if (Length & 1)
		filter[halfLength] = 1.0Q - 2.0Q * ft; // if length is odd, avoid divide-by-zero at centre-tap

	for (int n = 0; n<halfLength; ++n) {
		__float128 sinc = sinq(fmodq(M_TWOPIq * ft * (n - halfM), M_TWOPIq)) / (M_PIq * (n - halfM));	// sinc function
		filter[Length - n - 1] = filter[n] = -sinc;	// exploit symmetry
	}

#else

	// assert(ft < 0.5);
	int halfLength = Length / 2;
	double halfM = 0.5 * (Length - 1);
	double M_TWOPI = 2.0 * M_PI;

	if (Length & 1)
		filter[halfLength] = 1.0 - 2.0 * ft; // if length is odd, avoid divide-by-zero at centre-tap

	for (int n = 0; n < halfLength; ++n) {
		// sinc function
		double sinc = sin(fmod(M_TWOPI * ft * (n - halfM), M_TWOPI)) / (M_PI * (n - halfM));
		filter[Length - n - 1] = filter[n] = -sinc;	// exploit symmetry
	}
#endif

	return true;
}

// This function converts a requested sidelobe height (in dB) to a value for the Beta parameter used in a Kaiser window:
template<typename FloatType> FloatType calcKaiserBeta(FloatType dB)
{
	if(dB < 21.0)
	{
		return 0;
	}
	if ((dB >= 21.0) && (dB <= 50.0)) {
		return 0.5842 * pow((dB - 21), 0.4) + 0.07886 * (dB - 21);
	}
	if (dB > 50.0) {
		return 0.1102 * (dB - 8.7);
	}

	return 0;
}

// I0() : 0th-order Modified Bessel function of the first kind:
inline double I0(double z)
{
	double result = 0.0;
	for (int k = 0; k < 34; ++k) {
		double kfact = factorial[k];
		double x = pow(z * z / 4.0, k) / (kfact * kfact);
		result += x;
	}
	return result;
}

#ifdef FIR_QUAD_PRECISION
inline __float128 I0q(__float128 x)
{
	__float128 result = 0.0Q;
	__float128 kfact = 1.0Q;
	__float128 xx_4 = x * x / 4.0Q;
	for (int k = 0; k < 60; ++k){
		result += powq(xx_4, k) / factorialSquaredq[k];
	}
	return result;
}
#endif

// applyKaiserWindow() - This function applies a Kaiser Window to an array of filter coefficients ("textbook" version):
template<typename FloatType> bool applyKaiserWindow(FloatType* filter, int Length, double Beta)
{
	// Note: sometimes, the Kaiser Window formula is defined in terms of Alpha (instead of Beta),
	// in which case, Alpha def= Beta / pi

	if (Length < 1)
		return false;

#ifdef FIR_QUAD_PRECISION

	for (int n = 0; n < Length; ++n) {
		filter[n] *= I0q(Beta * sqrtq(1.0Q - powq((2.0Q * n / (Length - 1) - 1), 2.0Q)))
					 / I0q(Beta);
	}

#else

	for (int n = 0; n < Length; ++n) {
		filter[n] *= I0(Beta * sqrt(1.0 - pow((2.0 * n / (Length - 1) - 1), 2.0)))
			/ I0(Beta);
	}

#endif

	return true;
}

// applyKaiserWindow2() - applies a Kaiser Window to an array of filter coefficients (alternative formula):
template<typename FloatType> bool applyKaiserWindow2(FloatType* filter, int Length, double Beta)
 {
	 double A;	// use double internally, regardless of FloatType (speed not an issue here; no reason not to)
	 double maxA = 0; // for diagnostics
	 for (int n = 0; n < Length; ++n) {

		 // simplified Kaiser Window Equation:
		 A = (2.0 * Beta / Length) * sqrt(n*(Length - n - 1));
		 maxA = std::max(maxA, A);
		 filter[n] *= I0(A) / I0(Beta);
	 }

	return true;
}

inline std::vector<double> makeHilbert(int length)
{
	std::vector<double> coeffs;
	int length_ = std::max(length, 3) | 1;
	coeffs.resize(length_);
	int c = length_ / 2;
	coeffs[c] = 0.0;
	for(int n = 0; n < c; n++) {
		double s = std::sin((n - c) * M_PI / 2.0);
		double s2 = 2.0 * s * s;
		coeffs[n] = s2 / (M_PI * (n - c));
		coeffs[length_ - n - 1] = s2 / (M_PI * (c - n));
	}

	return coeffs;
}

// the following is a set of Complex-In, Complex-Out transforms used for constructing a minimum-Phase FIR:

// logV() : logarithm of a vector of Complex doubles
inline std::vector<std::complex<double>>
logV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return std::log(x); });
	return output;
}

// limitDynRangeV() : set a limit (-dB) on how quiet signal is allowed to be below the peak.
// Guaranteed to never return zero.
inline std::vector<std::complex<double>>
limitDynRangeV(const std::vector<std::complex<double>>& input, double dynRangeDB) {
	double dynRangeLinear = pow(10, std::abs(dynRangeDB) / 20.0); // will give same result for positive or negative dB values.

	// find peak:
	double peak=0.0;
	for (auto &c : input) {
		peak = std::max(peak, std::abs(c));
	}

	// determine low threshold
	double lowThresh = peak / dynRangeLinear;	// a level which is dynRangeDB dB below peak
	std::complex<double> lastX = lowThresh;		// variable for storing last output value

	std::vector<std::complex<double>> output(input.size(), 0);

	std::transform(input.begin(), input.end(), output.begin(),
		[lowThresh, &lastX](std::complex<double> x) -> std::complex<double> {

		double level = std::abs(x);
		if (level < lowThresh) {
			if (level == 0.0) {		// when input is zero, we must somehow make the modulus of the output equal to lowThresh
				x = lastX;			// sticky output; use last output value instead of zero
			}
			else {
				x = (x / level) * lowThresh; // scale x such that |x| == lowThresh
				lastX = x;
			}
		}
		return x;
	});

	return output;
}

// realV() : real parts of a vector of Complex doubles
inline std::vector<std::complex<double>>
realV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return x.real(); });
	return output;
}

// imagV() : imaginary parts of a vector of Complex doubles (answer placed in imaginary part of output):
inline std::vector<std::complex<double>>
imagV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return{ 0,x.imag() }; });
	return output;
}

// expV() : exp of a vector of Complex doubles
inline std::vector<std::complex<double>>
expV(const std::vector<std::complex<double>>& input) {
	std::vector<std::complex<double>> output(input.size(), 0);
	std::transform(input.begin(), input.end(), output.begin(),
		[](std::complex<double> x) -> std::complex<double> {return exp(x); });
	return output;
}

// fftV() : FFT of vector of Complex doubles
inline std::vector<std::complex<double>>
fftV(std::vector<std::complex<double>> input) {

	std::vector<std::complex<double>> output(input.size(), {0.0, 0.0}); // output vector

	// create, execute, destroy plan:
	fftw_plan p = fftw_plan_dft_1d(static_cast<int>(input.size()),
		reinterpret_cast<fftw_complex*>(&input[0]),
		reinterpret_cast<fftw_complex*>(&output[0]),
		FFTW_FORWARD,
		FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);

	return output;
}

// ifftV() : Inverse FFT of vector of Complex doubles
inline std::vector<std::complex<double>>
ifftV(std::vector<std::complex<double>> input) {

	std::vector<std::complex<double>> output(input.size(), 0); // output vector

	// create, execute, destroy plan:
	fftw_plan p = fftw_plan_dft_1d(static_cast<int>(input.size()),
		reinterpret_cast<fftw_complex*>(&input[0]),
		reinterpret_cast<fftw_complex*>(&output[0]),
		FFTW_BACKWARD,
		FFTW_ESTIMATE);

	fftw_execute(p);
	fftw_destroy_plan(p);

	// scale output:
	double reciprocalSize = 1.0 / input.size();
	for (auto &c : output){
		c *= reciprocalSize;
	}

	return output;
}

// AnalyticSignalV() : Analytic signal of vector of Complex doubles
// (Note: This function is referred to as "hilbert()" in Matlab / Octave, but it is not exactly a hilbert transform.
// The hilbert Transform is placed in the imaginary part, and the original input is in the real part.)
// See Footnote* below for more information on algorithm ...

inline std::vector<std::complex<double>>
AnalyticSignalV(const std::vector<std::complex<double>>& input) {

	std::vector<std::complex<double>> U = fftV(input);

	size_t N = input.size();
	size_t halfN = N / 2;

	// Note: U[0], U[halfN] unchanged:
	for (size_t n = 1; n < N; ++n) {
		if (n > halfN)
			U[n] = 0;
		if (n < halfN)
			U[n] *= 2.0;
	}

	std::vector<std::complex<double>> output = ifftV(U);
	return output;
}

// makeMinPhase() : transform linear-phase FIR filter coefficients into minimum-phase (in-place version)
template<typename FloatType>
void makeMinPhase(FloatType* pFIRcoeffs, size_t length)
{
	auto fftLength = static_cast<size_t>(pow(2, 2.0 + ceil(log2(length)))); // use FFT 4x larger than (length rounded-up to power-of-2)

	std::vector <std::complex<double>> complexInput;
	std::vector <std::complex<double>> complexOutput;

	// Pad zeros on either side of FIR:

	size_t frontPaddingLength = (fftLength - length) / 2;
	size_t backPaddingLength = fftLength - frontPaddingLength - length;

	for (size_t n = 0; n < frontPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	for (size_t n = 0; n < length; ++n) {
		complexInput.push_back({ pFIRcoeffs[n], 0 });
	}

	for (size_t n = 0; n < backPaddingLength; ++n) {
		complexInput.emplace_back(0, 0);
	}

	assert(complexInput.size() == fftLength); // make sure padding worked properly.

	// Formula is as follows:

	// take the reversed array of
	// the real parts of
	// the ifft of
	// e to the power of
	// the Analytic Signal of
	// the real parts of
	// the log of
	// the dynamic-ranged limited version of
	// the fft of
	// the original filter

	complexOutput = realV(ifftV(expV(AnalyticSignalV(realV(logV(limitDynRangeV(fftV(complexInput),-190)))))));
	std::reverse(complexOutput.begin(), complexOutput.end());

	// write all the real parts back to coeff array:
	size_t n = 0;
	for (auto &c : complexOutput) {
		if (n < length)
			pFIRcoeffs[n] = c.real();
		else
			break;
		++n;
	}
}

///////////////////////////////////////////////////////////////////////
// utility functions:

// dumpKaiserWindow() - utility function for displaying Kaiser Window:
inline void dumpKaiserWindow(size_t length, double Beta) {
	std::vector<double> f(length, 1);
	applyKaiserWindow<double>(f.data(), static_cast<int>(length), Beta);
	for (size_t i = 0; i < length; ++i) {
		std::cout << i << ": " << f[i] << std::endl;
	}

	std::vector<double> g(length, 1);
	applyKaiserWindow<double>(g.data(), static_cast<int>(length), Beta);
	for (size_t i = 0; i < length; ++i) {
		std::cout << i << ": " << g[i] << std::endl;
	}
}

// asserts that the two Kaiser Window formulas agree with each other (within a specified tolerance)
inline void assertKaiserWindow(size_t length, double Beta) {

	const double tolerance = 0.001;
	const double upper = 1.0 + tolerance;
	const double lower = 1.0 - tolerance;

	std::vector<double> f(length, 1);
	applyKaiserWindow2<double>(f.data(), static_cast<int>(length), Beta);

	std::vector<double> g(length, 1);
	applyKaiserWindow<double>(g.data(), static_cast<int>(length), Beta);

	for (size_t i = 0; i < length; ++i) {
		double ratio = f[i] / g[i];
		assert(ratio < upper && ratio > lower);
	}
}

// dumpFilter() - utility function for displaying filter coefficients:
template<typename FloatType> void dumpFilter(const FloatType* Filter, int Length) {
	const auto default_precision {std::cout.precision()};
	std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	for (int i = 0; i < Length; ++i) {
		std::cout << Filter[i] << std::endl;
	}
	std::cout << std::setprecision(default_precision);
}

inline void dumpComplexVector(const std::vector<std::complex<double>>& v)
{
	for (auto &c : v) {
		std::cout << c.real() << "+" << c.imag() << "i" << std::endl;
	}
}

template<typename FloatType>
void dumpFFT(FloatType* data, size_t length)
{
	auto pow2length = static_cast<size_t>(pow(2, 1.0 + floor(log2(length))));

	std::vector <std::complex<double>> complexInput;
	std::vector <std::complex<double>> complexOutput;

	for (int n = 0; n < pow2length; ++n) {
		if (n<length)
			complexInput.push_back({ data[n], 0 });
		else
			complexInput.emplace_back(0, 0); // pad remainder with zeros (to-do: does it mattter where the zeros are put ?)
	}

	complexOutput = fftV(complexInput);

	std::setprecision(17);
	std::cout << "real,imag,mag,phase" << std::endl;
	for (auto &c : complexOutput) {
		std::cout << c.real() << "," << c.imag() << "," << std::abs(c) << "," << arg(c) << std::endl;
	}
}

inline void testSinAccuracy() {

	const int numSteps = 10000000;
	const double inc = M_PI / numSteps;
	double t = M_PI / -2.0;
	double maxError = 0.0;
	double worstT = 0.0;

	for (int i = 0; i < numSteps; ++i ) {
		// calc relative error of
		// |(sin 2t - 2 * sint * cost) / sin 2t|
		// (double-angle identity)

		double e = std::abs((std::sin(2.0 * t) - 2.0 * std::sin(t) * std::cos(t)) / std::sin(2.0 * t));
		//double e = std::abs((sin(2.0 * t) - 2.0 * sin(t) * cos(t)) / sin(2.0 * t));
		if (e > maxError) {
			worstT = t;
			maxError = e;
		}
		t += inc;
	}
	std::cout << "maxError: " << std::setprecision(33) << maxError << std::endl;
	std::cout << "worstT: " << worstT << std::endl;
}


#endif // FIRGENERATOR_H
