=======
SpanSet
=======
A ``SpanSet`` is a class that represents regions of pixels on an image.  As the
name implies, a ``SpanSet`` is a collection of ``Span`` objects which behaive as
a mathematical set. The class provides methods for doing set operations, as well
as convienience method for operating on data products at the locations defined
by the ``SpanSet``.

When working with ``SpanSet`` objects it is important to note that they are
immutable once they are created. As such, all method calls on a ``SpanSet``
object which operate on a ``SpanSet`` create modify and return a new
``SpanSet`` object. Another important note is that regions defined by a
``SpanSet`` are not guarenteed to be contiguous. There is however a method which
will return if a given ``SpanSet`` is contiguous, and a method which will split
an non-contiguous ``SpanSet`` into a standard vector (or list in python) of
contiguous ``SpanSet`` objects.

Utility Methods
===============
``SpanSets`` contain many untily functions for common operations, or accessing
internal attributes. These utility functions are grouped below, and more
information for each method can be found in the methods documentation

- getArea
- getBBox
- isContiguous
- shiftedBy
- clippedTo
- transformedBy
- overlaps
- contains
- computeCentroid
- computeShape
- findEdgePixels
- spanSetFromShape :sup:`†`
- applyFunctor :sup:`††`

:sup:`†` This is a static method on the ``SpanSet`` class

:sup:`††` applyFunctor is only available in c++

Set Methods
===========
As mentioned above, a ``SpanSet`` behaves as a mathematical set, and accordingly
has the various set methods summarized below.

- dilate
- erode
- intersect
- intersectNot
- union :sup:`†`

:sup:`†` In c++ this method is called union\_

Mask Methods
============
Masks and ``SpanSet``\s can often be used in similar fassions to represent areas
of interest in an image. As such SpanSets provide many utility functions to make
transporting this information between these two data types similar. As Masks are
similar to ``SpanSet``\s they can even participate in ``SpanSet`` mathematical
operations. The list below contains both convience functions, and the set
operations Masks can participate in.

- intersect
- intersectNot
- union :sup:`†`
- setMask
- clearMask
- maskToSpanSet :sup:`††`

:sup:`†` In c++ this method is called union\_

:sup:`††` This is actually a free function in afw geom and not a ``SpanSet``
method but, is included in this list as it is logically grouped together.
