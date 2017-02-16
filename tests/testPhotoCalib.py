#
# LSST Data Management System
# Copyright 2008-2016 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
from __future__ import absolute_import, division, print_function
from builtins import zip
from builtins import range

import unittest
import math
import numpy as np

import lsst.utils.tests
import lsst.afw.image
import lsst.afw.image.utils
import lsst.daf.base
import lsst.pex.exceptions


class PhotoCalibTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.photoCalib = lsst.afw.image.PhotoCalib()

    def tearDown(self):
        del self.photoCalib

    def testPhotom(self):
        """Test the zero-point information"""

        flux, fluxErr = 1000.0, 10.0
        flux0, flux0Err = 1e12, 1e10
        self.photoCalib.setFluxMag0(flux0)

        self.assertEqual(flux0, self.photoCalib.getFluxMag0()[0])
        self.assertEqual(0.0, self.photoCalib.getFluxMag0()[1])
        self.assertEqual(22.5, self.photoCalib.getMagnitude(flux))
        # Error just in flux
        self.assertAlmostEqual(self.photoCalib.getMagnitude(flux, fluxErr)[1], 2.5/math.log(10)*fluxErr/flux)
        # Error just in flux0
        self.photoCalib.setFluxMag0(flux0, flux0Err)
        self.assertEqual(flux0Err, self.photoCalib.getFluxMag0()[1])
        self.assertAlmostEqual(self.photoCalib.getMagnitude(flux, 0)[1], 2.5/math.log(10)*flux0Err/flux0)

        self.assertAlmostEqual(flux0, self.photoCalib.getFlux(0))
        self.assertAlmostEqual(flux, self.photoCalib.getFlux(22.5))

        # I don't know how to test round-trip if fluxMag0 is significant compared to fluxErr
        self.photoCalib.setFluxMag0(flux0, flux0 / 1e6)
        for fluxErr in (flux / 1e2, flux / 1e4):
            mag, magErr = self.photoCalib.getMagnitude(flux, fluxErr)
            self.assertAlmostEqual(flux, self.photoCalib.getFlux(mag, magErr)[0])
            self.assertLess(abs(fluxErr - self.photoCalib.getFlux(mag, magErr)[1]), 1.0e-4)

        # Test context manager; shouldn't raise an exception within the block, should outside
        with lsst.afw.image.utils.CalibNoThrow():
            self.assertTrue(np.isnan(self.photoCalib.getMagnitude(-50.0)))
        with self.assertRaises(lsst.pex.exceptions.DomainError):
            self.photoCalib.getMagnitude(-50.0)

    def testPhotomMulti(self):
        self.photoCalib.setFluxMag0(1e12, 1e10)
        flux, fluxErr = 1000.0, 10.0
        num = 5

        mag, magErr = self.photoCalib.getMagnitude(flux, fluxErr)  # Result assumed to be true: tested elsewhere

        fluxList = np.array([flux for i in range(num)], dtype=float)
        fluxErrList = np.array([fluxErr for i in range(num)], dtype=float)

        magList = self.photoCalib.getMagnitude(fluxList)
        for m in magList:
            self.assertEqual(m, mag)

        mags, magErrs = self.photoCalib.getMagnitude(fluxList, fluxErrList)

        for m, dm in zip(mags, magErrs):
            self.assertEqual(m, mag)
            self.assertEqual(dm, magErr)

        flux, fluxErr = self.photoCalib.getFlux(mag, magErr)  # Result assumed to be true: tested elsewhere

        fluxList = self.photoCalib.getFlux(magList)
        for f in fluxList:
            self.assertEqual(f, flux)

        fluxes = self.photoCalib.getFlux(mags, magErrs)
        for f, df in zip(fluxes[0], fluxes[1]):
            self.assertAlmostEqual(f, flux)
            self.assertAlmostEqual(df, fluxErr)

    def testCtorFromMetadata(self):
        """Test building a Calib from metadata"""

        flux0, flux0Err = 1e12, 1e10
        flux, fluxErr = 1000.0, 10.0

        metadata = lsst.daf.base.PropertySet()
        metadata.add("FLUXMAG0", flux0)
        metadata.add("FLUXMAG0ERR", flux0Err)

        self.photoCalib = lsst.afw.image.Calib(metadata)

        self.assertEqual(flux0, self.photoCalib.getFluxMag0()[0])
        self.assertEqual(flux0Err, self.photoCalib.getFluxMag0()[1])
        self.assertEqual(22.5, self.photoCalib.getMagnitude(flux))
        # Error just in flux
        self.photoCalib.setFluxMag0(flux0, 0)

        self.assertAlmostEqual(self.photoCalib.getMagnitude(flux, fluxErr)[1], 2.5/math.log(10)*fluxErr/flux)

        # Check that we can clean up metadata
        lsst.afw.image.stripCalibKeywords(metadata)
        self.assertEqual(len(metadata.names()), 0)

    def testCalibEquality(self):
        self.assertEqual(self.photoCalib, self.photoCalib)
        self.assertFalse(self.photoCalib != self.photoCalib)  # using assertFalse to directly test != operator

        calib2 = lsst.afw.image.Calib()
        calib2.setFluxMag0(1200)

        self.assertNotEqual(calib2, self.photoCalib)

    def testCalibFromCalibs(self):
        """Test creating a Calib from an array of Calibs"""
        mag0, mag0Sigma = 1.0, 0.01

        calibs = lsst.afw.image.vectorCalib()
        ncalib = 3
        for i in range(ncalib):
            calib = lsst.afw.image.Calib()
            calib.setFluxMag0(mag0, mag0Sigma)

            calibs.append(calib)

        ocalib = lsst.afw.image.Calib(calibs)
        # the following is surely incorrect; see DM-7619
        self.assertEqual(ocalib.getFluxMag0(), (0.0, 0.0))

        # Check that we can only merge Calibs with the same fluxMag0 values
        calibs[0].setFluxMag0(1.001*mag0, mag0Sigma)
        with self.assertRaises(lsst.pex.exceptions.InvalidParameterError):
            lsst.afw.image.Calib(calibs)

    def testCalibNegativeFlux(self):
        """Check that we can control if negative fluxes raise exceptions"""
        self.photoCalib.setFluxMag0(1e12)

        with self.assertRaises(lsst.pex.exceptions.DomainError):
            self.photoCalib.getMagnitude(-10)
        with self.assertRaises(lsst.pex.exceptions.DomainError):
            self.photoCalib.getMagnitude(-10, 1)

        lsst.afw.image.Calib.setThrowOnNegativeFlux(False)
        mags = self.photoCalib.getMagnitude(-10)
        self.assertTrue(np.isnan(mags))
        mags = self.photoCalib.getMagnitude(-10, 1)
        self.assertTrue(np.isnan(mags[0]))
        self.assertTrue(np.isnan(mags[1]))

        lsst.afw.image.Calib.setThrowOnNegativeFlux(True)

        # Re-check that we raise after resetting ThrowOnNegativeFlux.
        with self.assertRaises(lsst.pex.exceptions.DomainError):
            self.photoCalib.getMagnitude(-10)
        with self.assertRaises(lsst.pex.exceptions.DomainError):
            self.photoCalib.getMagnitude(-10, 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
