// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
 * \file
 *
 * \brief Utilities to detect sets of Footprint%s
 *
 * Create and use an lsst::afw::detection::FootprintSet, a collection of pixels above (or below) a threshold
 * in an Image
 *
 * The "collections of pixels" are represented as lsst::afw::detection::Footprint%s, so an example application
 * would be:
 * \code
    namespace image = lsst::afw::image; namespace detection = lsst::afw::detection;

    image::MaskedImage<float> img(10,20);
    *img.getImage() = 100;

    detection::FootprintSet<float> sources(img, 10);
    cout << "Found " << sources.getFootprints()->size() << " sources" << std::endl;
 * \endcode
 */
#include <algorithm>
#include <cassert>
#include <set>
#include <string>
#include <typeinfo>
#include "boost/format.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/afw/detection/Peak.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/detection/FootprintSet.h"
#include "lsst/afw/detection/FootprintCtrl.h"
#include "lsst/afw/detection/HeavyFootprint.h"

namespace detection = lsst::afw::detection;
namespace image = lsst::afw::image;
namespace math = lsst::afw::math;
namespace pexLogging = lsst::pex::logging;
namespace geom = lsst::afw::geom;

/************************************************************************************************************/
namespace {
    /// Don't let doxygen see this block  \cond

    typedef boost::uint64_t IdPixelT;    // Type of temporary Images used in merging Footprints

    struct Threshold_traits {
    };
    struct ThresholdLevel_traits : public Threshold_traits { // Threshold is a single number
    };
    struct ThresholdPixelLevel_traits : public Threshold_traits { // Threshold varies from pixel to pixel
    };
    struct ThresholdBitmask_traits : public Threshold_traits { // Threshold ORs with a bitmask
    };

    //
    // Define our own functions to handle NaN tests;  this gives us the
    // option to define a value for e.g. image::MaskPixel or int
    //
    template<typename T>
    inline bool isBadPixel(T) {
        return false;
    }

    template<>
    inline bool isBadPixel(float val) {
        return std::isnan(val);
    }

    template<>
    inline bool isBadPixel(double val) {
        return std::isnan(val);
    }

    /*
     * Return the number of bits required to represent a unsigned long
     */
    int nbit(unsigned long i) {
        int n = 0;
        while (i > 0) {
            ++n;
            i >>= 1;
        }

        return n;
    }
    /*
     * Find the list of pixel values that lie in a Footprint
     *
     * Used when the Footprints are constructed from an Image containing Footprint indices
     */
    template <typename ImageT>
    class FindIdsInFootprint: public detection::FootprintFunctor<ImageT> {
    public:
        explicit FindIdsInFootprint(ImageT const& image ///< The image the source lives in
                                   ) : detection::FootprintFunctor<ImageT>(image), _ids(), _old(0) {}
        /// \brief Reset everything for a new Footprint
        void reset() {
            _ids.clear();
            _old = 0;
        }
        
        /// \brief method called for each pixel by apply()
        void operator()(typename ImageT::xy_locator loc, ///< locator pointing at the pixel
                        int x,                           ///< column-position of pixel
                        int y                            ///< row-position of pixel
                       ) {
            typename ImageT::Pixel val = loc(0, 0);

            if (val != _old) {
                _ids.insert(val);
                _old = val;
            }
        }

        std::set<typename ImageT::Pixel> const& getIds() const {
            return _ids;
        }

    private:
        std::set<typename ImageT::Pixel> _ids;
        typename ImageT::Pixel _old;
    };

    /********************************************************************************************************/
    /*
     * Sort peaks by decreasing pixel value.  N.b. -ve peaks are sorted the same way as +ve ones
     */
    struct SortPeaks {
	bool operator()(CONST_PTR(detection::PeakRecord) a, CONST_PTR(detection::PeakRecord) b) {
            if (a->getPeakValue() != b->getPeakValue()) {
                return (a->getPeakValue() > b->getPeakValue());
            }

            if (a->getIx() != b->getIx()) {
                return (a->getIx() < b->getIx());
            }

            return (a->getIy() < b->getIy());
        }
    };
    /********************************************************************************************************/
    /*
     * Worker routine for merging two FootprintSets, possibly growing them as we proceed
     */
    detection::FootprintSet
    mergeFootprintSets(
        detection::FootprintSet const &lhs, // the FootprintSet to be merged to
        int rLhs,                                         // Grow lhs Footprints by this many pixels
        detection::FootprintSet const &rhs, // the FootprintSet to be merged into lhs
        int rRhs,                                         // Grow rhs Footprints by this many pixels
        detection::FootprintControl const& ctrl           // Control how the grow is done
                      )
    {
        typedef detection::Footprint Footprint;
        typedef detection::FootprintSet::FootprintList FootprintList;
        // The isXXX routines return <isset, value>
        bool const circular = ctrl.isCircular().first && ctrl.isCircular().second;
        bool const isotropic = ctrl.isIsotropic().second; // isotropic grow as opposed to a Manhattan metric
                                        // n.b. Isotropic grows are significantly slower
        bool const left =  ctrl.isLeft().first  && ctrl.isLeft().second;
        bool const right = ctrl.isRight().first && ctrl.isRight().second;
        bool const up =    ctrl.isUp().first    && ctrl.isUp().second;
        bool const down =  ctrl.isDown().first  && ctrl.isDown().second;

        geom::Box2I const region = lhs.getRegion();
        if (region != rhs.getRegion()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              boost::format("The two FootprintSets must have the same region").str());
        }
        
        image::Image<IdPixelT>::Ptr idImage(new image::Image<IdPixelT>(region));
        idImage->setXY0(region.getMinX(), region.getMinY());
        *idImage = 0;

        FootprintList const& lhsFootprints = *lhs.getFootprints();
        FootprintList const& rhsFootprints = *rhs.getFootprints();
        int const nLhs = lhsFootprints.size();
        int const nRhs = rhsFootprints.size();
        /*
         * In general the lists of Footprints overlap, so we need to make sure that the IDs can be
         * uniquely recovered from the idImage.  We do this by allocating a range of bits to the lhs IDs
         */
        int const lhsIdNbit = nbit(nLhs);
        int const lhsIdMask = (lhsIdNbit == 0) ? 0x0 : (1 << lhsIdNbit) - 1;

        if (std::size_t(nRhs << lhsIdNbit) > std::numeric_limits<IdPixelT>::max() - 1) {
            throw LSST_EXCEPT(lsst::pex::exceptions::OverflowErrorException,
                              (boost::format("%d + %d footprints need too many bits; change IdPixelT typedef")
                               % nLhs % nRhs).str());
        }
        /*
         * When we insert grown Footprints into the idImage we can potentially overwrite an entire Footprint,
         * losing any peaks that it might contain.  We'll preserve the overwritten Ids in case we need to
         * get them back (n.b. Footprints that overlap, but both if which survive, will appear in this list)
         */
        typedef std::map<int, std::set<boost::uint64_t> > OldIdMap;
        OldIdMap overwrittenIds;        // here's a map from id -> overwritten IDs

        IdPixelT id = 1;                     // the ID inserted into the image
        for (FootprintList::const_iterator ptr = lhsFootprints.begin(), end = lhsFootprints.end();
             ptr != end; ++ptr, ++id) {
            CONST_PTR(Footprint) foot = *ptr;

            if (rLhs > 0) {
                foot = circular ?
                    growFootprint(*foot, rLhs, isotropic) : growFootprint(*foot, rLhs, left, right, up, down);
            }

            std::set<boost::uint64_t> overwritten;
            foot->insertIntoImage(*idImage, id, true, 0x0, &overwritten);

            if (!overwritten.empty()) {
                overwrittenIds.insert(overwrittenIds.end(), std::make_pair(id, overwritten));
            }
        }

        assert (id <= std::size_t(1 << lhsIdNbit));
        id = (1 << lhsIdNbit);
        for (FootprintList::const_iterator ptr = rhsFootprints.begin(), end = rhsFootprints.end();
             ptr != end; ++ptr, id += (1 << lhsIdNbit)) {
            CONST_PTR(Footprint) foot = *ptr;

            if (rRhs > 0) {
                foot = circular ?
                    growFootprint(*foot, rRhs, isotropic) : growFootprint(*foot, rRhs, left, right, up, down);
            }

            std::set<boost::uint64_t> overwritten;
            foot->insertIntoImage(*idImage, id, true, lhsIdMask, &overwritten);

            if (!overwritten.empty()) {
                overwrittenIds.insert(overwrittenIds.end(), std::make_pair(id, overwritten));
            }
        }

        detection::FootprintSet fs(*idImage, detection::Threshold(1),
                                   1, false); // detect all pixels in rhs + lhs
        /*
         * Now go through the new Footprints looking up and remembering their progenitor's IDs; we'll use
         * these IDs to merge the peaks in a moment
         *
         * We can't do this as we go through the idFinder as the IDs it returns are
         *   (lhsId + 1) | ((rhsId + 1) << nbit)
         * and, depending on the geometry, values of lhsId and/or rhsId can appear multiple times
         * (e.g. if nbit is 2, idFinder IDs 0x5 and 0x6 both contain lhsId = 0) so we get duplicates
         * of peaks.  This is not too bad, but it's a bit of a pain to make the lists unique again,
         * and we avoid this by this two-step process.
         */
        FindIdsInFootprint<image::Image<IdPixelT> > idFinder(*idImage);
        for (FootprintList::iterator ptr = fs.getFootprints()->begin(),
                 end = fs.getFootprints()->end(); ptr != end; ++ptr) {
            PTR(Footprint) foot = *ptr;

            idFinder.apply(*foot);      // find the (mangled) [lr]hsFootprint IDs that contribute to foot

            std::set<boost::uint64_t> lhsFootprintIndxs, rhsFootprintIndxs; // indexes into [lr]hsFootprints

            for (std::set<IdPixelT>::iterator idptr = idFinder.getIds().begin(),
                     idend = idFinder.getIds().end(); idptr != idend; ++idptr) {
                unsigned int indx = *idptr;
                if ((indx & lhsIdMask) > 0) {
                    boost::uint64_t i = (indx & lhsIdMask) - 1;
                    lhsFootprintIndxs.insert(i);
                    /*
                     * Now allow for Footprints that vanished beneath this one
                     */
                    OldIdMap::iterator mapPtr = overwrittenIds.find(indx);
                    if (mapPtr != overwrittenIds.end()) {
                        std::set<boost::uint64_t> &overwritten = mapPtr->second;

                        for (std::set<boost::uint64_t>::iterator ptr = overwritten.begin(),
                                 end = overwritten.end(); ptr != end; ++ptr){
                            lhsFootprintIndxs.insert((*ptr & lhsIdMask) - 1);
                        }
                    }
                }
                indx >>= lhsIdNbit;

                if (indx > 0) {
                    boost::uint64_t i = indx - 1;
                    rhsFootprintIndxs.insert(i);
                    /*
                     * Now allow for Footprints that vanished beneath this one
                     */
                    OldIdMap::iterator mapPtr = overwrittenIds.find(indx);
                    if (mapPtr != overwrittenIds.end()) {
                        std::set<boost::uint64_t> &overwritten = mapPtr->second;

                        for (std::set<boost::uint64_t>::iterator ptr = overwritten.begin(),
                                 end = overwritten.end(); ptr != end; ++ptr) {
                            rhsFootprintIndxs.insert(*ptr - 1);
                        }
                    }
                }
            }
            /*
             * We now have a complete set of Footprints that contributed to this one, so merge
             * all their Peaks into the new one
             */
            detection::PeakCatalog &peaks = foot->getPeaks();

            for (std::set<boost::uint64_t>::iterator ptr = lhsFootprintIndxs.begin(),
                     end = lhsFootprintIndxs.end(); ptr != end; ++ptr) {
                boost::uint64_t i = *ptr;
                assert (i < lhsFootprints.size());
                detection::PeakCatalog const& oldPeaks = lhsFootprints[i]->getPeaks();

                int const nold = peaks.size();
                peaks.insert(peaks.end(), oldPeaks.begin(), oldPeaks.end());
                // We use getInternal() here to get the vector of shared_ptr that Catalog uses internally,
                // which causes the STL algorithm to copy pointers instead of PeakRecords (which is what
                // it'd try to do if we passed Catalog's own iterators).
                std::inplace_merge(peaks.getInternal().begin(), peaks.getInternal().begin() + nold,
                                   peaks.getInternal().end(), SortPeaks());
            }

            for (std::set<boost::uint64_t>::iterator ptr = rhsFootprintIndxs.begin(),
                     end = rhsFootprintIndxs.end(); ptr != end; ++ptr) {
                boost::uint64_t i = *ptr;
                assert (i < rhsFootprints.size());
                detection::PeakCatalog const& oldPeaks = rhsFootprints[i]->getPeaks();

                int const nold = peaks.size();
                peaks.insert(peaks.end(), oldPeaks.begin(), oldPeaks.end());
                // See note above on why we're using getInternal() here.
                std::inplace_merge(peaks.getInternal().begin(), peaks.getInternal().begin() + nold,
                                   peaks.getInternal().end(), SortPeaks());
            }
        }

        return fs;
    }
/*
 * run-length code for part of object
 */
    class IdSpan {
    public:
        typedef boost::shared_ptr<IdSpan> Ptr;
        
        explicit IdSpan(int id, int y, int x0, int x1, double good) : 
            id(id), y(y), x0(x0), x1(x1), good(good) {}
        int id;                         /* ID for object */
        int y;                          /* Row wherein IdSpan dwells */
        int x0, x1;                     /* inclusive range of columns */
        bool good;                      /* includes a value over the desired threshold? */
    };
/*
 * comparison functor; sort by ID then row
 */
    struct IdSpanCompar : public std::binary_function<const IdSpan::Ptr, const IdSpan::Ptr, bool> {
        bool operator()(IdSpan::Ptr const a, IdSpan::Ptr const b) {
            if (a->id < b->id) {
                return true;
            } else if (a->id > b->id) {
                return false;
            } else {
                return (a->y < b->y) ? true : false;
            }
        }
    };
/*
 * Follow a chain of aliases, returning the final resolved value.
 */
    int resolve_alias(std::vector<int> const &aliases, /* list of aliases */
                      int id) {         /* alias to look up */
        int resolved = id;              /* resolved alias */
        
        while (id != aliases[id]) {
            resolved = id = aliases[id];
        }
        
        return(resolved);
    }
    /// \endcond
}

/************************************************************************************************************/

namespace {
    /*
     * Find all the Peaks within a Footprint
     */
    template <typename ImageT>
    class FindPeaksInFootprint: public detection::FootprintFunctor<ImageT> {
    public:
        explicit FindPeaksInFootprint(ImageT const& image, ///< The image the source lives in
                                      bool polarity,       ///< true if we're looking for -ve "peaks"
                                      detection::PeakCatalog &peaks
                                     ) : detection::FootprintFunctor<ImageT>(image),
                                         _polarity(polarity), _peaks(peaks) {}
        
        /// \brief method called for each pixel by apply()
        void operator()(typename ImageT::xy_locator loc, ///< locator pointing at the pixel
                        int x,                           ///< column-position of pixel
                        int y                            ///< row-position of pixel
                       ) {
            typename ImageT::Pixel val = loc(0, 0);

            if (_polarity) {            // look for +ve peaks
                if (loc(-1,  1) > val || loc( 0,  1) > val || loc( 1,  1) > val || 
                    loc(-1,  0) > val ||                      loc( 1,  0) > val || 
                    loc(-1, -1) > val || loc( 0, -1) > val || loc( 1, -1) > val) {
                    return;
                }
            } else {                    // look for -ve "peaks" (pits)
                if (loc(-1,  1) < val || loc( 0,  1) < val || loc( 1,  1) < val || 
                    loc(-1,  0) < val ||                      loc( 1,  0) < val || 
                    loc(-1, -1) < val || loc( 0, -1) < val || loc( 1, -1) < val) {
                    return;
                }
            }

            PTR(detection::PeakRecord) newPeak = _peaks.addNew();
            newPeak->setIx(x);
            newPeak->setIy(y);
            newPeak->setFx(x);
            newPeak->setFy(y);
            newPeak->setPeakValue(val);
        }
    private:
        bool _polarity;
        detection::PeakCatalog &_peaks;
    };

    /*
     * Find the maximum (or minimum, if polarity is false) pixel in a Footprint
     */
    template <typename ImageT>
    class FindMaxInFootprint : public detection::FootprintFunctor<ImageT> {
    public:
        explicit FindMaxInFootprint(ImageT const& image, ///< The image the source lives in
                                    bool polarity        ///< true if we're looking for -ve "peaks"
                                   ) : detection::FootprintFunctor<ImageT>(image),
                                       _polarity(polarity), _x(0), _y(0),
                                       _min( std::numeric_limits<double>::max()),
                                       _max(-std::numeric_limits<double>::max()) {}

        /// \brief Reset everything for a new Footprint
        void reset() {
            _x = _y = 0;
            _min =  std::numeric_limits<double>::max();
            _max = -std::numeric_limits<double>::max();
        }
        virtual void reset(detection::Footprint const&) {}

        /// \brief method called for each pixel by apply()
        void operator()(typename ImageT::xy_locator loc, ///< locator pointing at the pixel
                        int x,                           ///< column-position of pixel
                        int y                            ///< row-position of pixel
                       ) {
            typename ImageT::Pixel val = loc(0, 0);

            if (_polarity) {
                if (val > _max) {
                    _max = val;
                    _x = x;
                    _y = y;
                }
            } else {
                if (val < _min) {
                    _min = val;
                    _x = x;
                    _y = y;
                }
            }
        }

        // Add the Footprint's peak to the given PeakCatalog
        void addPeak(detection::PeakCatalog & peakCat) const {
            PTR(detection::PeakRecord) newPeak = peakCat.addNew();
            newPeak->setIx(_x);
            newPeak->setIy(_y);
            newPeak->setFx(_x);
            newPeak->setFy(_y);
            newPeak->setPeakValue(_polarity ? _max : _min);
        }
    private:
        bool _polarity;
        int _x, _y;
        double _min, _max;
    };
    
    template<typename ImageT, typename ThresholdT>
    void findPeaks(PTR(detection::Footprint) foot, ImageT const& img, bool polarity, ThresholdT)
    {
        FindPeaksInFootprint<ImageT> peakFinder(img, polarity, foot->getPeaks());
        peakFinder.apply(*foot, 1);

        // We use getInternal() here to get the vector of shared_ptr that Catalog uses internally,
        // which causes the STL algorithm to copy pointers instead of PeakRecords (which is what
        // it'd try to do if we passed Catalog's own iterators).
        std::stable_sort(foot->getPeaks().getInternal().begin(), foot->getPeaks().getInternal().end(),
                         SortPeaks());

        if (foot->getPeaks().empty()) {
            FindMaxInFootprint<ImageT> maxFinder(img, polarity);
            maxFinder.apply(*foot);
            maxFinder.addPeak(foot->getPeaks());
        }
    }

    // No need to search for peaks when processing a Mask
    template<typename ImageT>
    void findPeaks(PTR(detection::Footprint), ImageT const&, bool, ThresholdBitmask_traits)
    {
        ;
    }
}

/************************************************************************************************************/
/*
 * Functions to determine if a pixel's in a Footprint
 */
template<typename ImagePixelT, typename IterT>
static inline bool inFootprint(ImagePixelT pixVal, IterT,
                               bool polarity, double thresholdVal, ThresholdLevel_traits) {
    return (polarity ? pixVal : -pixVal) >= thresholdVal;
}

template<typename ImagePixelT, typename IterT>
static inline bool inFootprint(ImagePixelT pixVal, IterT var,
                               bool polarity, double thresholdVal, ThresholdPixelLevel_traits) {
    return (polarity ? pixVal : -pixVal) >= thresholdVal*::sqrt(*var);
}

template<typename ImagePixelT, typename IterT>
static inline bool inFootprint(ImagePixelT pixVal, IterT,
                               bool, double thresholdVal, ThresholdBitmask_traits) {
    return (pixVal & static_cast<long>(thresholdVal));
}

/*
 * Advance the x_iterator to the variance image, when relevant (it may be NULL otherwise)
 */
template<typename IterT>
static inline IterT
advancePtr(IterT varPtr, Threshold_traits) {
    return varPtr;
}

template<typename IterT>
static inline IterT
advancePtr(IterT varPtr, ThresholdPixelLevel_traits) {
    return varPtr + 1;
}

/*
 * Here's the working routine for the FootprintSet constructors; see documentation
 * of the constructors themselves
 */
template<typename ImagePixelT, typename MaskPixelT, typename VariancePixelT, typename ThresholdTraitT>
static void findFootprints(
        typename detection::FootprintSet::FootprintList *_footprints, // Footprints
        geom::Box2I const& _region,               // BBox of pixels that are being searched
        image::ImageBase<ImagePixelT> const &img, // Image to search for objects
        image::Image<VariancePixelT> const *var,  // img's variance
        double const footprintThreshold,  // threshold value for footprint
        double const includeThresholdMultiplier,  // threshold (relative to footprintThreshold) for inclusion
        bool const polarity,                      // if false, search _below_ thresholdVal
        int const npixMin,                        // minimum number of pixels in an object
        bool const setPeaks                       // should I set the Peaks list?
)
{
    int id;                             /* object ID */
    int in_span;                        /* object ID of current IdSpan */
    int nobj = 0;                       /* number of objects found */
    int x0 = 0;                         /* unpacked from a IdSpan */

    typedef typename image::Image<ImagePixelT> ImageT;
    double includeThreshold = footprintThreshold * includeThresholdMultiplier; // Threshold for inclusion
    
    int const row0 = img.getY0();
    int const col0 = img.getX0();
    int const height = img.getHeight();
    int const width = img.getWidth();
/*
 * Storage for arrays that identify objects by ID. We want to be able to
 * refer to idp[-1] and idp[width], hence the (width + 2)
 */
    std::vector<int> id1(width + 2);
    std::fill(id1.begin(), id1.end(), 0);
    std::vector<int> id2(width + 2);
    std::fill(id2.begin(), id2.end(), 0);
    std::vector<int>::iterator idc = id1.begin() + 1; // object IDs in current/
    std::vector<int>::iterator idp = id2.begin() + 1; //                       previous row

    std::vector<int> aliases;           // aliases for initially disjoint parts of Footprints
    aliases.reserve(1 + height/20);     // initial size of aliases

    std::vector<IdSpan::Ptr> spans;     // y:x0,x1 for objects
    spans.reserve(aliases.capacity());  // initial size of spans

    aliases.push_back(0);               // 0 --> 0
/*
 * Go through image identifying objects
 */
    typedef typename image::Image<ImagePixelT>::x_iterator x_iterator;
    typedef typename image::Image<VariancePixelT>::x_iterator x_var_iterator;

    in_span = 0;                        // not in a span
    for (int y = 0; y != height; ++y) {
        if (idc == id1.begin() + 1) {
            idc = id2.begin() + 1;
            idp = id1.begin() + 1;
        } else {
            idc = id1.begin() + 1;
            idp = id2.begin() + 1;
        }
        std::fill_n(idc - 1, width + 2, 0);
        
        in_span = 0;                    /* not in a span */
        bool good = (includeThresholdMultiplier == 1.0); /* Span exceeds the threshold? */

        x_iterator pixPtr = img.row_begin(y);
        x_var_iterator varPtr = (var == NULL) ? NULL : var->row_begin(y);
        for (int x = 0; x < width; ++x, ++pixPtr, varPtr = advancePtr(varPtr, ThresholdTraitT())) {
            ImagePixelT const pixVal = *pixPtr;

            if (isBadPixel(pixVal) ||
                !inFootprint(pixVal, varPtr, polarity, footprintThreshold, ThresholdTraitT())) {
                if (in_span) {
                    IdSpan::Ptr sp(new IdSpan(in_span, y, x0, x - 1, good));
                    spans.push_back(sp);

                    in_span = 0;
                    good = false;
                }
            } else {                    /* a pixel to fix */
                if (idc[x - 1] != 0) {
                    id = idc[x - 1];
                } else if (idp[x - 1] != 0) {
                    id = idp[x - 1];
                } else if (idp[x] != 0) {
                    id = idp[x];
                } else if (idp[x + 1] != 0) {
                    id = idp[x + 1];
                } else {
                    id = ++nobj;
                    aliases.push_back(id);
                }

                idc[x] = id;
                if (!in_span) {
                    x0 = x;
                    in_span = id;
                }
/*
 * Do we need to merge ID numbers? If so, make suitable entries in aliases[]
 */
                if (idp[x + 1] != 0 && idp[x + 1] != id) {
                    aliases[resolve_alias(aliases, idp[x + 1])] = resolve_alias(aliases, id);
               
                    idc[x] = id = idp[x + 1];
                }

                if (!good && inFootprint(pixVal, varPtr, polarity, includeThreshold, ThresholdTraitT())) {
                    good = true;
                }
            }
        }

        if (in_span) {
            IdSpan::Ptr sp(new IdSpan(in_span, y, x0, width - 1, good));
            spans.push_back(sp);
        }
    }
/*
 * Resolve aliases; first alias chains, then the IDs in the spans
 */
    for (unsigned int i = 0; i < spans.size(); i++) {
        spans[i]->id = resolve_alias(aliases, spans[i]->id);
    }
/*
 * Sort spans by ID, so we can sweep through them once
 */
    if (spans.size() > 0) {
        std::sort(spans.begin(), spans.end(), IdSpanCompar());
    }
/*
 * Build Footprints from spans
 */
    unsigned int i0;                    // initial value of i
    if (spans.size() > 0) {
        id = spans[0]->id;
        i0 = 0;
        for (unsigned int i = 0; i <= spans.size(); i++) { // <= size to catch the last object
            if (i == spans.size() || spans[i]->id != id) {
                PTR(detection::Footprint) fp(new detection::Footprint(i - i0, _region));
            
                bool good = false;      // Span includes pixel sufficient to include footprint in set?
                for (; i0 < i; i0++) {
                    good |= spans[i0]->good;
                    fp->addSpan(spans[i0]->y + row0, spans[i0]->x0 + col0, spans[i0]->x1 + col0);
                }

                if (good && !(fp->getNpix() < npixMin)) {
                    _footprints->push_back(fp);
                }
            }

            if (i < spans.size()) {
                id = spans[i]->id;
            }
        }
    }
/*
 * Find all peaks within those Footprints
 */
    if (setPeaks) {
        typedef detection::FootprintSet::FootprintList::iterator fiterator;
        for (fiterator ptr = _footprints->begin(), end = _footprints->end(); ptr != end; ++ptr) {
            findPeaks(*ptr, img, polarity, ThresholdTraitT());
        }
    }
}

/************************************************************************************************************/
/*
 * \brief Find a FootprintSet given an Image and a threshold
 */
template<typename ImagePixelT>
detection::FootprintSet::FootprintSet(
    image::Image<ImagePixelT> const &img, //!< Image to search for objects
    Threshold const &threshold,     //!< threshold to find objects
    int const npixMin,              //!< minimum number of pixels in an object
    bool const setPeaks            //!< should I set the Peaks list?
) : lsst::daf::base::Citizen(typeid(this)),
    _footprints(new FootprintList()),
    _region(img.getBBox(image::PARENT))
{
    typedef float VariancePixelT;
     
    findFootprints<ImagePixelT, afw::image::MaskPixel, VariancePixelT, ThresholdLevel_traits>(
        _footprints.get(), 
        _region, 
        img,
        NULL,
        threshold.getValue(img), threshold.getIncludeMultiplier(), threshold.getPolarity(),
        npixMin,
        setPeaks
    );
}

// NOTE: not a template to appease swig (see note by instantiations at bottom)

/*
 * \brief Find a FootprintSet given a Mask and a threshold
 */
template <typename MaskPixelT>
detection::FootprintSet::FootprintSet(
    image::Mask<MaskPixelT> const &msk, //!< Image to search for objects
    Threshold const &threshold,     //!< threshold to find objects
    int const npixMin               //!< minimum number of pixels in an object
) : lsst::daf::base::Citizen(typeid(this)),
    _footprints(new FootprintList()),
    _region(msk.getBBox(image::PARENT))
{
    switch (threshold.getType()) {
      case Threshold::BITMASK:
          findFootprints<MaskPixelT, MaskPixelT, float, ThresholdBitmask_traits>(
            _footprints.get(), _region, msk, NULL, threshold.getValue(), threshold.getIncludeMultiplier(),
            threshold.getPolarity(), npixMin, false);
        break;

      case Threshold::VALUE:
        findFootprints<MaskPixelT, MaskPixelT, float, ThresholdLevel_traits>(
            _footprints.get(), _region, msk, NULL, threshold.getValue(), threshold.getIncludeMultiplier(),
            threshold.getPolarity(), npixMin, false);
        break;

      default:
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "You must specify a numerical threshold value with a Mask");
    }
}


/**
 * \brief Find a FootprintSet given a MaskedImage and a threshold
 *
 * Go through an image, finding sets of connected pixels above threshold
 * and assembling them into Footprint%s;  the resulting set of objects
 * is returned as an \c array<Footprint::Ptr>
 *
 * If threshold.getPolarity() is true, pixels above the Threshold are
 * assembled into Footprints; if it's false, then pixels \e below Threshold
 * are processed (Threshold will probably have to be below the background level
 * for this to make sense, e.g. for difference imaging)
 */
template<typename ImagePixelT, typename MaskPixelT>
detection::FootprintSet::FootprintSet(
    const image::MaskedImage<ImagePixelT, MaskPixelT> &maskedImg, //!< MaskedImage to search for objects
    Threshold const &threshold,     //!< threshold for footprints (controls size)
    std::string const &planeName,   //!< mask plane to set (if != "")
    int const npixMin,              //!< minimum number of pixels in an object
    bool const setPeaks            //!< should I set the Peaks list?
) : lsst::daf::base::Citizen(typeid(this)),
    _footprints(new FootprintList()),
    _region(
        geom::Point2I(maskedImg.getX0(), maskedImg.getY0()),
        geom::Extent2I(maskedImg.getWidth(), maskedImg.getHeight())
    )
{
    typedef typename image::MaskedImage<ImagePixelT, MaskPixelT>::Variance::Pixel VariancePixelT;
    // Find the Footprints    
    switch (threshold.getType()) {
      case Threshold::PIXEL_STDEV:
        findFootprints<ImagePixelT, MaskPixelT, VariancePixelT, ThresholdPixelLevel_traits>(
            _footprints.get(), 
            _region,
            *maskedImg.getImage(), 
            maskedImg.getVariance().get(), 
            threshold.getValue(maskedImg),
            threshold.getIncludeMultiplier(),
            threshold.getPolarity(),
            npixMin,
            setPeaks
                                                                                  );
        break;
      default:
        findFootprints<ImagePixelT, MaskPixelT, VariancePixelT, ThresholdLevel_traits>(
            _footprints.get(), 
            _region,
            *maskedImg.getImage(), 
            maskedImg.getVariance().get(), 
            threshold.getValue(maskedImg),
            threshold.getIncludeMultiplier(),
            threshold.getPolarity(),
            npixMin,
            setPeaks
                                                                                  );
        break;
    }
    // Set Mask if requested    
    if (planeName == "") {
        return;
    }
    //
    // Define the maskPlane
    //
    const typename image::Mask<MaskPixelT>::Ptr mask = maskedImg.getMask();
    mask->addMaskPlane(planeName);

    MaskPixelT const bitPlane = mask->getPlaneBitMask(planeName);
    //
    // Set the bits where objects are detected
    //
    typedef image::Mask<MaskPixelT> MaskT;

    class MaskFootprint : public detection::FootprintFunctor<MaskT> {
    public:
        MaskFootprint(MaskT const& mimage,
                      MaskPixelT bit) : detection::FootprintFunctor<MaskT>(mimage), _bit(bit) {}

        void operator()(typename MaskT::xy_locator loc, int, int) {
            *loc |= _bit;
        }
    private:
        MaskPixelT _bit;
    };

    MaskFootprint maskit(*maskedImg.getMask(), bitPlane);
    for (FootprintList::const_iterator fiter = _footprints->begin();         
         fiter != _footprints->end(); ++fiter
    ) {
        Footprint::Ptr const foot = *fiter;

        maskit.apply(*foot);
    }
}

/************************************************************************************************************/
/**
 * Construct an empty FootprintSet given a region that its footprints would have lived in
 */
detection::FootprintSet::FootprintSet(geom::Box2I region ///< the desired region
) :
    lsst::daf::base::Citizen(typeid(this)),
    _footprints(PTR(FootprintList)(new FootprintList)), _region(region) {
}

/**
 * Copy constructor
 */
detection::FootprintSet::FootprintSet(
    FootprintSet const &rhs         //!< the input FootprintSet
) :
    lsst::daf::base::Citizen(typeid(this)),
    _footprints(new FootprintList), _region(rhs._region)
{
    _footprints->reserve(rhs._footprints->size());
    for (FootprintSet::FootprintList::const_iterator ptr = rhs._footprints->begin(),
             end = rhs._footprints->end(); ptr != end; ++ptr) {
        _footprints->push_back(PTR(Footprint)(new Footprint(**ptr)));
    }
}

/// Assignment operator.
detection::FootprintSet &
detection::FootprintSet::operator=(FootprintSet const& rhs) {
    FootprintSet tmp(rhs);
    swap(tmp);                          // See Meyers, Effective C++, Item 11    
    return *this;
}

/************************************************************************************************************/
/**
 * Merge a FootprintSet into *this
 */
void detection::FootprintSet::merge(
        detection::FootprintSet const& rhs, ///< the Footprints to merge
        int tGrow,                          ///< No. of pixels to grow this Footprints
        int rGrow,                          ///< No. of pixels to grow rhs Footprints
        bool isotropic                      ///< Use (expensive) isotropic grow
)
{
    detection::FootprintControl const ctrl(true, isotropic);
    detection::FootprintSet fs = mergeFootprintSets(*this, tGrow, rhs, rGrow, ctrl);
    swap(fs);                           // Swap the new FootprintSet into place
}

/// Set the corners of the FootprintSet's MaskedImage to region
///
/// N.b. updates all the Footprints' regions too
//
void detection::FootprintSet::setRegion(
    geom::Box2I const& region ///< desired region
) {
    _region = region;

    for (FootprintSet::FootprintList::iterator ptr = _footprints->begin(),
             end = _footprints->end(); ptr != end; ++ptr
    ) {
        (*ptr)->setRegion(region);
    }
}

/************************************************************************************************************/
/**
 * Grow all the Footprints in the input FootprintSet, returning a new FootprintSet
 *
 * The output FootprintSet may contain fewer Footprints, as some may well have been merged
 */
detection::FootprintSet::FootprintSet(
    FootprintSet const &rhs,        //!< the input FootprintSet
    int r,                          //!< Grow Footprints by r pixels
    bool isotropic                  //!< Grow isotropically (as opposed to a Manhattan metric)
    //!< @note Isotropic grows are significantly slower
)
    : lsst::daf::base::Citizen(typeid(this)), _footprints(new FootprintList), _region(rhs._region)
{
    if (r == 0) {
        FootprintSet fs = rhs;
        swap(fs);                       // Swap the new FootprintSet into place
        return;
    } else if (r < 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          (boost::format("I cannot grow by negative numbers: %d") % r).str());
    }

    detection::FootprintControl const ctrl(true, isotropic);
    detection::FootprintSet fs = mergeFootprintSets(FootprintSet(rhs.getRegion()), 0, rhs, r, ctrl);
    swap(fs);                           // Swap the new FootprintSet into place
}

/************************************************************************************************************/

detection::FootprintSet::FootprintSet(detection::FootprintSet const& rhs,
                                      int ngrow,
                                      detection::FootprintControl const& ctrl)
    : lsst::daf::base::Citizen(typeid(this)), _footprints(new FootprintList), _region(rhs._region)
{
    if (ngrow == 0) {
        FootprintSet fs = rhs;
        swap(fs);                       // Swap the new FootprintSet into place
        return;
    } else if (ngrow < 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          str(boost::format("I cannot grow by negative numbers: %d") % ngrow));
    }

    detection::FootprintSet fs = mergeFootprintSets(FootprintSet(rhs.getRegion()), 0, rhs, ngrow, ctrl);
    swap(fs);                           // Swap the new FootprintSet into place
}

/************************************************************************************************************/
/**
 * Return the FootprintSet corresponding to the merge of two input FootprintSets
 *
 * \todo Implement this.  There's RHL Pan-STARRS code to do it, but it isn't yet converted to LSST C++
 */
detection::FootprintSet::FootprintSet(
        FootprintSet const& fs1,
        FootprintSet const& fs2,
        bool const 
                                                              )
    : lsst::daf::base::Citizen(typeid(this)),
      _footprints(new FootprintList()),
      _region(fs1._region)
{
    _region.include(fs2._region);
    throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException, "NOT IMPLEMENTED");
}

/************************************************************************************************************/
/**
 * Return an Image with pixels set to the Footprint%s in the FootprintSet
 *
 * \returns an image::Image::Ptr
 */
PTR(image::Image<detection::FootprintIdPixel>)
detection::FootprintSet::insertIntoImage(
    bool const relativeIDs          ///< Use IDs starting at 0 (rather than the ones in the Footprint%s)
) const {
    PTR(image::Image<detection::FootprintIdPixel>) im(
        new image::Image<detection::FootprintIdPixel>(_region)
    );
    *im = 0;

    detection::FootprintIdPixel id = 0;
    for (FootprintList::const_iterator fiter = _footprints->begin(); 
         fiter != _footprints->end(); fiter++
    ) {
        Footprint::Ptr const foot = *fiter;
        
        if (relativeIDs) {
            id++;
        } else {
            id = foot->getId();
        }
        
        foot->insertIntoImage(*im.get(), id);
    }
    
    return im;
}

/************************************************************************************************************/
/**
 * Convert all the Footprints in the FootprintSet to be HeavyFootprint%s
 */
template<typename ImagePixelT, typename MaskPixelT>
void
detection::FootprintSet::makeHeavy(
    image::MaskedImage<ImagePixelT, MaskPixelT> const& mimg, ///< the image providing pixel values
    HeavyFootprintCtrl const *ctrl     ///< Control how we manipulate HeavyFootprints
)
{
    HeavyFootprintCtrl ctrl_s = HeavyFootprintCtrl();

    if (!ctrl) {
        ctrl = &ctrl_s;
    }

    for (FootprintList::iterator ptr = _footprints->begin(),
                                          end = _footprints->end(); ptr != end; ++ptr) {
        ptr->reset(new detection::HeavyFootprint<ImagePixelT, MaskPixelT>(**ptr, mimg, ctrl));
    }
}

void detection::FootprintSet::makeSources(
    lsst::afw::table::SourceCatalog & cat
) const {
    for (FootprintList::const_iterator i = _footprints->begin(); i != _footprints->end(); ++i) {
        PTR(afw::table::SourceRecord) r = cat.addNew();
        r->setFootprint(*i);
    }
}


/************************************************************************************************************/
//
// Explicit instantiations
//

#ifndef DOXYGEN

#define INSTANTIATE(PIXEL)                      \
    template detection::FootprintSet::FootprintSet(                     \
        image::Image<PIXEL> const &, Threshold const &, int const, bool const); \
    template detection::FootprintSet::FootprintSet(                     \
        image::MaskedImage<PIXEL,image::MaskPixel> const &, Threshold const &, \
        std::string const &, int const, bool const);\
    template void detection::FootprintSet::makeHeavy(image::MaskedImage<PIXEL,image::MaskPixel> const &, \
                                                     HeavyFootprintCtrl const *)

template detection::FootprintSet::FootprintSet(image::Mask<image::MaskPixel> const &,
                                               Threshold const &, int const);

template void detection::FootprintSet::setMask(image::Mask<image::MaskPixel> *, std::string const &);
template void detection::FootprintSet::setMask(PTR(image::Mask<image::MaskPixel>), std::string const &);

INSTANTIATE(boost::uint16_t);
INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);

#endif // !DOXYGEN
