#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

"""
Tests for Masks

Run with:
   python Mask.py
or
   python
   >>> import Mask; Mask.run()
"""

import os
import os.path

import sys
import unittest
import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as pexExcept
import lsst.daf.base
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import eups
import lsst.afw.display.ds9 as ds9

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def showMaskDict(d=None, msg=None):
    if not d:
        d = afwImage.MaskU(0,0)
        if not msg:
            msg = "default"

    try:
        d = d.getMaskPlaneDict()
    except AttributeError:
        pass
        
    if msg:
        print "%-15s" % msg,
    print sorted([(d[p], p) for p in d])

class MaskTestCase(unittest.TestCase):
    """A test case for Mask"""

    def setUp(self):
        self.Mask = afwImage.MaskU

        self.Mask.clearMaskPlaneDict() # reset so tests will be deterministic

        for p in ("BAD", "SAT", "INTRP", "CR", "EDGE"):
            self.Mask.addMaskPlane(p)

        self.BAD  = afwImage.MaskU_getPlaneBitMask("BAD")
        self.CR   = afwImage.MaskU_getPlaneBitMask("CR")
        self.EDGE = afwImage.MaskU_getPlaneBitMask("EDGE")

        self.val1 = self.BAD | self.CR
        self.val2 = self.val1 | self.EDGE

        self.mask1 = afwImage.MaskU(100, 200)
        self.mask1.set(self.val1)
        self.mask2 = afwImage.MaskU(self.mask1.getDimensions())
        self.mask2.set(self.val2)

        dataDir = os.path.join(eups.productDir("afwdata"), "data")
        if dataDir:
            if True:
                self.maskFile = os.path.join(dataDir, "small_MI_msk.fits")
            else:
                self.maskFile = os.path.join(dataDir, "871034p_1_MI_msk.fits")
        else:
            self.maskFile = None

    def tearDown(self):
        del self.mask1
        del self.mask2

    def testArrays(self):
        image1 = afwImage.MaskU(afwGeom.ExtentI(5,6)) # could use MaskU(5, 6) but check extent(5, 6) form too
        array1 = image1.getArray()
        self.assertEqual(array1.shape[0], image1.getHeight())
        self.assertEqual(array1.shape[1], image1.getWidth())
        image2 = afwImage.MaskU(array1, False)
        self.assertEqual(array1.shape[0], image2.getHeight())
        self.assertEqual(array1.shape[1], image2.getWidth())
        image3 = afwImage.makeMaskFromArray(array1)
        self.assertEqual(array1.shape[0], image2.getHeight())
        self.assertEqual(array1.shape[1], image2.getWidth())
        self.assertEqual(type(image3), afwImage.MaskU)
        array1[:,:] = numpy.random.uniform(low=0, high=10, size=array1.shape)
        for j in range(image1.getHeight()):
            for i in range(image1.getWidth()):
                self.assertEqual(image1.get(i, j), array1[j, i])
                self.assertEqual(image2.get(i, j), array1[j, i])

    def testInitializeMasks(self):
        val = 0x1234
        msk = afwImage.MaskU(afwGeom.ExtentI(10, 10), val)
        self.assertEqual(msk.get(0, 0), val)
        
    def testSetGetMasks(self):
        self.assertEqual(self.mask1.get(0, 0), self.val1)
    
    def testOrMasks(self):
        self.mask2 |= self.mask1
        self.mask1 |= self.val2
        
        self.assertEqual(self.mask1.get(0, 0), self.val1 | self.val2)
        self.assertEqual(self.mask2.get(0, 0), self.val1 | self.val2)
    
    def testAndMasks(self):
        self.mask2 &= self.mask1
        self.mask1 &= self.val2
        
        self.assertEqual(self.mask1.get(0, 0), self.val1 & self.val2)
        self.assertEqual(self.mask1.get(0, 0), self.BAD | self.CR)
        self.assertEqual(self.mask2.get(0, 0), self.val1 & self.val2)

    def testXorMasks(self):
        self.mask2 ^= self.mask1
        self.mask1 ^= self.val2
        
        self.assertEqual(self.mask1.get(0, 0), self.val1 ^ self.val2)
        self.assertEqual(self.mask2.get(0, 0), self.val1 ^ self.val2)

    def testLogicalMasksMismatch(self):
        "Test logical operations on Masks of different sizes"
        i1 = afwImage.MaskU(afwGeom.ExtentI(100, 100))
        i1.set(100)
        i2 = afwImage.MaskU(afwGeom.ExtentI(10, 10))
        i2.set(10)
        
        def tst(i1, i2): i1 |= i2
        utilsTests.assertRaisesLsstCpp(self, lsst.pex.exceptions.LengthErrorException, tst, i1, i2)

        def tst2(i1, i2): i1 &= i2
        utilsTests.assertRaisesLsstCpp(self, lsst.pex.exceptions.LengthErrorException, tst2, i1, i2)
    
        def tst2(i1, i2): i1 ^= i2
        utilsTests.assertRaisesLsstCpp(self, lsst.pex.exceptions.LengthErrorException, tst2, i1, i2)
    
    def testMaskPlanes(self):
        planes = self.Mask().getMaskPlaneDict()
        self.assertEqual(len(planes), self.Mask.getNumPlanesUsed())
        
        for k in sorted(planes.keys()):
            self.assertEqual(planes[k], self.Mask.getMaskPlane(k))
            
    def testCopyConstructors(self):
        dmask = afwImage.MaskU(self.mask1, True) # deep copy
        smask = afwImage.MaskU(self.mask1) # shallow copy
        
        self.mask1 |= 32767             # should only change dmask
        self.assertEqual(dmask.get(0, 0), self.val1)
        self.assertEqual(smask.get(0, 0), self.val1 | 32767)

    def testSubmasks(self):
        smask = afwImage.MaskU(self.mask1, 
                               afwGeom.Box2I(afwGeom.Point2I(1, 1), afwGeom.ExtentI(3, 2)),
                               afwImage.LOCAL)
        mask2 = afwImage.MaskU(smask.getDimensions())

        mask2.set(666)
        smask <<= mask2
        
        del smask
        del mask2
        
        self.assertEqual(self.mask1.get(0, 0), self.val1)
        self.assertEqual(self.mask1.get(1, 1), 666)
        self.assertEqual(self.mask1.get(4, 1), self.val1)
        self.assertEqual(self.mask1.get(1, 2), 666)
        self.assertEqual(self.mask1.get(4, 2), self.val1)
        self.assertEqual(self.mask1.get(1, 3), self.val1)

    def testReadFits(self):
        if not self.maskFile:
            print >> sys.stderr, "Warning: afwdata is not set up; not running the FITS I/O tests"
            return

        nMaskPlanes0 = afwImage.MaskU_getNumPlanesUsed()
        mask = afwImage.MaskU(self.maskFile) # will shift any unrecognised mask planes into unused slots

        if False:
            for (k, v) in afwImage.MaskU_getMaskPlaneDict().items():
                print k, v

        self.assertEqual(mask.get(32, 1), 0)
        self.assertEqual(mask.get(50, 50), 0)
        self.assertEqual(mask.get(0,  0), (1<<nMaskPlanes0))

    def testReadFitsConform(self):
        if not self.maskFile:
            print >> sys.stderr, "Warning: afwdata is not set up; not running the FITS I/O tests"
            return
        
        hdu = 0
        mask = afwImage.MaskU(self.maskFile, hdu, None, afwGeom.Box2I(), afwImage.LOCAL, True)

        if False:
            import lsst.afw.display.ds9 as ds9
            ds9.mtv(mask)

        if False:
            for (k, v) in afwImage.MaskU_getMaskPlaneDict().items():
                print k, v

        self.assertEqual(mask.get(32, 1), 0)
        self.assertEqual(mask.get(50, 50), 0)
        self.assertEqual(mask.get(0, 0), 1)

    def testWriteFits(self):
        if not self.maskFile:
            print >> sys.stderr, "Warning: afwdata is not set up; not running the FITS I/O tests"
            return

        nMaskPlanes0 = afwImage.MaskU_getNumPlanesUsed()
        mask = afwImage.MaskU(self.maskFile)

        self.assertEqual(mask.get(32, 1), 0)
        self.assertEqual(mask.get(50, 50), 0)
        self.assertEqual(mask.get(0, 0), (1<<nMaskPlanes0)) # as header had none of the canonical planes

        tmpFile = "foo.fits"
        mask.writeFits(tmpFile)
        #
        # Read it back
        #
        md = lsst.daf.base.PropertySet()
        rmask = self.Mask(tmpFile, 0, md)
        os.remove(tmpFile)
        
        self.assertEqual(mask.get(0, 0), rmask.get(0, 0))
        #
        # Check that we wrote (and read) the metadata successfully
        #
        mp_ = "MP_" if True else self.Mask.maskPlanePrefix() # currently private
        for (k, v) in self.Mask().getMaskPlaneDict().items():
            self.assertEqual(md.get(mp_ + k), v)

    def testReadWriteXY0(self):
        """Test that we read and write (X0, Y0) correctly"""
        mask = afwImage.MaskU(afwGeom.ExtentI(10, 20))

        x0, y0 = 1, 2
        mask.setXY0(x0, y0)
        tmpFile = "foo.fits"
        mask.writeFits(tmpFile)

        mask2 = mask.Factory(tmpFile)
        os.remove(tmpFile)

        self.assertEqual(mask2.getX0(), x0)
        self.assertEqual(mask2.getY0(), y0)

    def testMaskInitialisation(self):
        dims = self.mask1.getDimensions()
        factory = self.mask1.Factory

        self.mask1.set(666)

        del self.mask1                 # tempt C++ to reuse the memory
        self.mask1 = factory(dims)
        self.assertEqual(self.mask1.get(10, 10), 0)

        del self.mask1
        self.mask1 = factory(afwGeom.ExtentI(20, 20))
        self.assertEqual(self.mask1.get(10, 10), 0)

    def testBoundsChecking(self):
        """Check that pixel indexes are checked in python"""
        tsts = []
        def tst():
            self.mask1.get(-1, 0)
        tsts.append(tst)

        def tst():
            self.mask1.get(0, -1)
        tsts.append(tst)

        def tst():
            self.mask1.get(self.mask1.getWidth(), 0)
        tsts.append(tst)

        def tst():
            self.mask1.get(0, self.mask1.getHeight())
        tsts.append(tst)

        for tst in tsts:
            utilsTests.assertRaisesLsstCpp(self, lsst.pex.exceptions.LengthErrorException, tst)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class OldMaskTestCase(unittest.TestCase):
    """A test case for Mask (based on MaskU_1.cc); these are taken over from the DC2 fw tests
    and modified to run with the new (DC3) APIs"""

    def setUp(self):
        self.Mask = afwImage.MaskU           # the class

        self.testMask = self.Mask(afwGeom.Extent2I(300, 400), 0)

        self.testMask.clearMaskPlaneDict() # reset so tests will be deterministic

        for p in ("CR", "BP"):
            self.Mask.addMaskPlane(p)

        self.region = afwGeom.Box2I(afwGeom.Point2I(100, 300), afwGeom.Extent2I(10, 40))
        self.subTestMask = self.Mask(self.testMask, self.region, afwImage.LOCAL)

        if False:
            self.pixelList = afwImage.listPixelCoord()
            for x in range(0, 300):
                for y in range(300, 400, 20):
                    self.pixelList.push_back(afwImage.PixelCoord(x, y))
                    
    def tearDown(self):
        del self.testMask
        del self.subTestMask
        del self.region

    def testPlaneAddition(self):
        """Test mask plane addition"""

        nplane = self.testMask.getNumPlanesUsed()
        for p in ("XCR", "XBP"):
            self.assertEqual(self.Mask.addMaskPlane(p), nplane, "Assigning plane %s" % (p))
            nplane += 1

        def pname(p):
            return "P%d" % p

        nextra = 8
        for p in range(0, nextra):
            plane = self.Mask.addMaskPlane(pname(p))

        for p in range(0, nextra):
            self.testMask.removeAndClearMaskPlane(pname(p))
        
        self.assertEqual(nplane + nextra, self.Mask.getNumPlanesUsed(), "Adding and removing planes")

        for p in range(0, nextra):
            self.Mask.removeMaskPlane(pname(p))
        
        self.assertEqual(nplane, self.testMask.getNumPlanesUsed(), "Adding and removing planes")

    def testMetadata(self):
        """Test mask plane metadata interchange with MaskPlaneDict"""
        #
        # Demonstrate that we can extract a MaskPlaneDict into metadata
        #
        metadata = lsst.daf.base.PropertySet()

        self.Mask.addMaskPlanesToMetadata(metadata)
        for (k, v) in self.Mask().getMaskPlaneDict().items():
            self.assertEqual(metadata.getInt("MP_%s" % k), v)
        #
        # Now add another plane to metadata and make it appear in the mask Dict, albeit
        # in general at another location (hence the getNumPlanesUsed call)
        #
        metadata.addInt("MP_" + "Whatever", afwImage.MaskU_getNumPlanesUsed())

        self.testMask.conformMaskPlanes(afwImage.MaskU_parseMaskPlaneMetadata(metadata))
        for (k, v) in self.Mask().getMaskPlaneDict().items():
            self.assertEqual(metadata.getInt("MP_%s" % k), v)

    def testPlaneOperations(self):
        """Test mask plane operations"""

        planes = self.Mask().getMaskPlaneDict()
        self.testMask.clearMaskPlane(planes['CR'])

        if False:
            for p in planes.keys():
                self.testMask.setMaskPlaneValues(planes[p], self.pixelList)

        #printMaskPlane(self.testMask, planes['CR'])

        #print "\nClearing mask"
        self.testMask.clearMaskPlane(planes['CR'])

        #printMaskPlane(self.testMask, planes['CR'])

    def testPlaneRemoval(self):
        """Test mask plane removal"""

        def checkPlaneBP():
            self.Mask.getMaskPlane("BP")

        testMask2 = self.Mask(self.testMask.getDimensions())
        self.testMask = self.Mask(self.testMask.getDimensions())
        self.testMask.removeAndClearMaskPlane("BP")

        d = testMask2.getMaskPlaneDict()
        
        checkPlaneBP                    # still present in default mask
        self.assertTrue("BP" in testMask2.getMaskPlaneDict()) # should still be in testMask2

        self.Mask.removeMaskPlane("BP") # remove from default mask too

        utilsTests.assertRaisesLsstCpp(self, pexExcept.InvalidParameterException, checkPlaneBP)
        #
        # Check that removeAndClearMaskPlane can clear the default too
        #
        self.Mask.addMaskPlane("BP")
        self.testMask.removeAndClearMaskPlane("BP", True)

        utilsTests.assertRaisesLsstCpp(self, pexExcept.InvalidParameterException, checkPlaneBP)

    def testInvalidPlaneOperations(self):
        """Test mask plane operations invalidated by Mask changes"""

        testMask3 = self.Mask(self.testMask.getDimensions())
        
        name = "Great Timothy"
        self.Mask.addMaskPlane(name)
        testMask3.removeAndClearMaskPlane(name)

        self.Mask.getMaskPlane(name)    # should be fine
        self.assertRaises(IndexError, lambda : testMask3.getMaskPlaneDict()[name])
            
        def tst():
            self.testMask |= testMask3

        utilsTests.assertRaisesLsstCpp(self, pexExcept.RuntimeErrorException, tst)

        self.Mask.addMaskPlane(name)    # The dictionary should be back to the same state, so ...
        tst                             # ... assertion should not fail

        self.testMask.removeAndClearMaskPlane(name, True)
        self.Mask.addMaskPlane("Mario") # takes name's slot
        self.Mask.addMaskPlane(name)

        utilsTests.assertRaisesLsstCpp(self, pexExcept.RuntimeErrorException, tst)

    def testInvalidPlaneOperations2(self):
        """Test mask plane operations invalidated by Mask changes"""

        testMask3 = self.Mask(self.testMask.getDimensions())

        name = "Great Timothy"
        name2 = "Our Boss"
        self.Mask.addMaskPlane(name)
        self.Mask.addMaskPlane(name2)
        oldDict = testMask3.getMaskPlaneDict() # a description of the Mask's current dictionary

        for n in (name, name2):
            self.testMask.removeAndClearMaskPlane(n)
            self.Mask.removeMaskPlane(n)

        self.Mask.addMaskPlane(name2)        # added in opposite order to the planes in testMask3
        self.Mask.addMaskPlane(name)

        self.assertNotEqual(self.testMask.getMaskPlaneDict()[name], oldDict[name])

        def tst():
            self.testMask |= testMask3

        self.testMask.removeAndClearMaskPlane("BP")

        utilsTests.assertRaisesLsstCpp(self, pexExcept.RuntimeErrorException, tst)
        #
        # OK, that failed as it should.  Fixup the dictionaries and try again
        #
        self.Mask.addMaskPlane("BP")
        testMask3.conformMaskPlanes(oldDict) # convert testMask3 from oldDict to current default

        self.testMask |= testMask3      # shouldn't throw

    def testConformMaskPlanes(self):
        """Test conformMaskPlanes() when the two planes are actually the same"""

        testMask3 = self.Mask(self.testMask.getDimensions())

        name = "XXX"
        self.Mask.addMaskPlane(name)
        oldDict = testMask3.getMaskPlaneDict()
        testMask3.removeAndClearMaskPlane(name) # invalidates dictionary version

        testMask3.conformMaskPlanes(oldDict)

        self.testMask |= testMask3

    def testConformMaskPlanes2(self):
        """Test conformMaskPlanes() when the two planes are different"""

        testMask3 = afwImage.MaskU(self.testMask.getDimensions())
        
        name1 = "Great Timothy"
        name2 = "Our Boss"
        p1 = self.Mask.addMaskPlane(name1)
        p2 = self.Mask.addMaskPlane(name2)
        oldDict = self.testMask.getMaskPlaneDict()

        testMask3.setMaskPlaneValues(p1, 0, 5, 0)
        testMask3.setMaskPlaneValues(p2, 0, 5, 1)

        if display:
            im = afwImage.ImageF(testMask3.getDimensions())
            ds9.mtv(im)                 # bug in ds9's Mask display; needs an Image first
            ds9.mtv(testMask3)

        self.assertEqual(testMask3.get(0, 0), testMask3.getPlaneBitMask(name1))
        self.assertEqual(testMask3.get(0, 1), testMask3.getPlaneBitMask(name2))

        self.testMask.removeAndClearMaskPlane(name1, True)
        self.testMask.removeAndClearMaskPlane(name2, True)
        self.Mask.addMaskPlane(name2) # added in opposite order to testMask3
        self.Mask.addMaskPlane(name1)

        self.assertEqual(self.testMask.get(0, 0), 0)

        if display:
            ds9.mtv(im, frame=1)
            ds9.mtv(testMask3, frame=1)

        self.assertNotEqual(testMask3.get(0, 0), testMask3.getPlaneBitMask(name1))
        self.assertNotEqual(testMask3.get(0, 1), testMask3.getPlaneBitMask(name2))

        testMask3.conformMaskPlanes(oldDict)

        self.assertEqual(testMask3.get(0, 0), testMask3.getPlaneBitMask(name1))
        self.assertEqual(testMask3.get(0, 1), testMask3.getPlaneBitMask(name2))

        if display:
            ds9.mtv(im, frame=2)
            ds9.mtv(testMask3, frame=2)

        self.testMask |= testMask3

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def printMaskPlane(mask, plane,
                   xrange=range(250, 300, 10), yrange=range(300, 400, 20)):
    """Print parts of the specified plane of the mask"""
    
    if True:
        xrange = range(min(xrange), max(xrange), 25)
        yrange = range(min(yrange), max(yrange), 25)

    for x in xrange:
        for y in yrange:
            if False:                   # mask(x,y) confuses swig
                print x, y, mask(x, y), mask(x, y, plane)
            else:
                print x, y, mask(x, y, plane)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(OldMaskTestCase) # test suite from vw-based Masks
    suites += unittest.makeSuite(MaskTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
