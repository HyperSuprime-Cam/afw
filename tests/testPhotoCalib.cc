/*
 * LSST Data Management System
 * See COPYRIGHT file at the top of the source tree.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program. If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PhotoCalib

#include "boost/test/unit_test.hpp"

#include "lsst/afw/image/PhotoCalib.h"

namespace lsst {
namespace afw {
namespace image {

/// Tests whether the default constructor creates a zeroed calib
BOOST_AUTO_TEST_CASE(PhotoCalibDefault) {
    PhotoCalib photoCalib;
    BOOST_TEST(photoCalib.getFluxMag0() == 0);
    BOOST_TEST(photoCalib.getFluxMag0Sigma() == 0);
}

/// Tests the non-spatially-varying zeropoint constructor
BOOST_AUTO_TEST_CASE(PhotoCalibNonVarying) {
    PhotoCalib photoCalibNoError(10);
    BOOST_TEST(photoCalibNoError.getFluxMag0() == 10.0);
    BOOST_TEST(photoCalibNoError.getFluxMag0Sigma() == 0);
    PhotoCalib photoCalibwithError(0.5, 0.1);
    BOOST_TEST(photoCalibwithError.getFluxMag0() == 0.5);
    BOOST_TEST(photoCalibwithError.getFluxMag0Sigma() == 0.1);
}

/// Tests whether a spatially constant calib converts correctly
BOOST_AUTO_TEST_CASE(PhotoCalibNonVaryingConvertCounts) {
    PhotoCalib photoCalib(100, 1);
    BOOST_TEST(photoCalib.countsToMaggies(.5) == 0.005);
    BOOST_TEST(photoCalib.countsToMagnitude(.5) == -2.5*log10(0.5/100));
    BOOST_TEST(photoCalib.countsToMaggies(.5, 3.) == std::make_pair(0.005, 2.));
    BOOST_TEST(photoCalib.countsToMagnitude(.5, 3.) == std::make_pair(-2.5*log10(0.5/100), 2.));
}

}}}
