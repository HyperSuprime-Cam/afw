# This file is part of afw.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import copy
import gc
import unittest

import lsst.utils.tests

from lsst.afw.typehandling import Storable
import testGenericMapLib as cppLib


class DemoStorable(Storable):
    """Test that we can inherit from Storable in Python.
    """

    def __init__(self, state):
        super().__init__()
        self._state = state

    def __str__(self):
        return "value = %s" % self._state

    def __repr__(self):
        return "DemoStorable(%r)" % self._state

    def __hash__(self):
        return hash(self._state)

    def __copy__(self):
        return DemoStorable(self._state)

    def __deepcopy__(self, memo=None):
        return DemoStorable(copy.deepcopy(self._state, memo))

    def __eq__(self, other):
        return self._state == other._state


class SpecializedStorable(cppLib.CppStorable):
    """Test that we can inherit from C++ subclasses of Storable that
    are not Storable itself.
    """

    def __repr__(self):
        return "Pythonic " + super().__repr__()


class PythonStorableTestSuite(lsst.utils.tests.TestCase):

    def setUp(self):
        self.aList = [42]
        self.testbed = DemoStorable(self.aList)

    def testCopy(self):
        shallow = copy.copy(self.testbed)
        self.assertIsNot(shallow, self.testbed)
        self.assertEqual(shallow, self.testbed)

        deep = copy.deepcopy(self.testbed)
        self.assertIsNot(deep, self.testbed)
        self.assertEqual(deep, self.testbed)

        self.aList.append(43)
        self.assertEqual(shallow, DemoStorable([42, 43]))
        self.assertEqual(deep, DemoStorable([42]))

    def testStr(self):
        self.assertEqual(str(self.testbed), "value = [42]")

    def testRepr(self):
        self.assertEqual(repr(self.testbed), "DemoStorable([42])")
        cppLib.assertPythonStorable(self.testbed, "DemoStorable([42])")

    def testHash(self):
        with self.assertRaises(TypeError):
            hash(self.testbed)

    def testEq(self):
        self.assertEqual(self.testbed, DemoStorable([42]))
        self.assertNotEqual(self.testbed, DemoStorable(0))

    def testGarbageCollection(self):
        cppLib.keepStaticStorable(DemoStorable(3))

        gc.collect()

        retrieved = cppLib.keepStaticStorable()
        self.assertIsInstance(retrieved, Storable)
        self.assertIsInstance(retrieved, DemoStorable)
        self.assertEqual(retrieved, DemoStorable(3))

    def testInheritedGarbageCollection(self):
        cppLib.keepStaticStorable(SpecializedStorable("Foo"))

        gc.collect()

        retrieved = cppLib.keepStaticStorable()
        self.assertIsInstance(retrieved, Storable)
        self.assertIsInstance(retrieved, cppLib.CppStorable)
        self.assertIsInstance(retrieved, SpecializedStorable)
        self.assertEqual(repr(retrieved), "Pythonic Foo")
        cppLib.assertPythonStorable(retrieved, "Pythonic Foo")


class CppStorableTestSuite(lsst.utils.tests.TestCase):

    def setUp(self):
        self.initstr = "Just a string"
        self.testbed = cppLib.CppStorable(self.initstr)

    def testNewValue(self):
        """Test a Python-side state change in both C++ and Python.
        """
        self.assertEqual(self.testbed.value, self.initstr)
        cppLib.assertCppValue(self.testbed, self.initstr)

        newstr = "Stringly typed"
        self.testbed.value = newstr

        self.assertEqual(self.testbed.value, newstr)
        cppLib.assertCppValue(self.testbed, newstr)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
