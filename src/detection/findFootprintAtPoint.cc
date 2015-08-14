#include <vector>
#include <deque>

#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Mask.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Span.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst {
namespace afw {
namespace detection {


namespace {

// Indicates we hit a stop bit
class StopException {};

// Functor to apply threshold
template <typename PixelT>
struct Thresholder {
public:
    Thresholder(PixelT const threshold, bool polarity) : _threshold(threshold), _polarity(polarity) {}

    bool operator()(PixelT const pixel) const {
        return (_polarity && pixel >= _threshold) || (!_polarity && pixel < _threshold);
    }

private:
    PixelT _threshold;        // Threshold value
    bool _polarity;           // Polarity of threshold (true for >, false for <)
};


// Mask plane names for detecting
//
// We're not using the mask planes defined in afw::image::Mask because there's a
// finite number of those, and we don't want to run out. This provides our own
// namespace.
enum MaskPlanes {
    DETECTED = 0x0001, // Pixel is detected
    STOP     = 0x0002, // Stop searching
};

// A search for spans about a point
//
// As we search, we add starting points to a queue (to trace all the fingers), and
// work until the queue is empty.
//
// If template parameter 'throwOnStop' is true, this will throw StopException
// if a stopping point appears in the footprint.
//
// This is the common code behind findFootprintAtPoint and checkFootprintAtPoint.
template <bool throwOnStop, typename PixelT>
class SpanSearch {
public:

    // Starting point for a search
    struct StartingPoint {
        StartingPoint(geom::Span const& _span, bool const _dir) : span(_span), direction(_dir) {}
        geom::Span span;     // The span of interest
        bool direction;      // Direction to search: true for up, false for down
    };

    typedef typename std::deque<StartingPoint> Queue;

    SpanSearch(Footprint & footprint, // Footprint to which to add spans
               image::Image<PixelT> const& image, // Image we're searching
               geom::Point2I const& point,        // Point around which to search
               Thresholder<PixelT> const& threshold, // Threshold to apply
               std::vector<geom::Point2I> const& stops=std::vector<geom::Point2I>() // Stopping points
               );

    // Add a new StartingPoint to the queue
    void add(geom::Span const& span, bool direction) {
        _queue.push_back(StartingPoint(span, direction));
    }

    // Restart a search
    //
    // Commonly used after encountering an overhanging span, so we need to search
    // both up and down.
    void restart(geom::Span const& span) {
        add(span, true);
        add(span, false);
    }

    // Mark a span as found
    geom::Span const& found(int ySpan, int xSpanStart, int xSpanStop);

    // Search for spans from a starting point
    //
    // Additional starting points may be uncovered, which are appended
    // to the queue.
    void process(StartingPoint const& start);

    // Search for new spans
    void run();

private:
    Footprint & _footprint;             // Footprint we'll construct
    image::Image<PixelT> const& _image; // Image we're searching
    image::Mask<image::MaskPixel> _mask; // Mask that tells us what we've done already
    Thresholder<PixelT> const& _threshold; // Threshold to apply
    Queue _queue;                        // Queue of StartingPoints
};

template <bool throwOnStop, typename PixelT>
SpanSearch<throwOnStop, PixelT>::SpanSearch(
    Footprint & fp,
    image::Image<PixelT> const& image,
    geom::Point2I const& point,
    Thresholder<PixelT> const& threshold,
    std::vector<geom::Point2I> const& stops
    ) :
    _footprint(fp),
    _image(image),
    _mask(image.getBBox(image::PARENT)),
    _threshold(threshold),
    _queue()
{
    int const x0 = image.getX0(), y0 = image.getY0();
    int const width = image.getWidth();

    // We need a mask for two purposes: to indicate which pixels are already detected,
    // and to store the "stop" pixels --- those that, once reached, should stop us
    // looking for the rest of the Footprint.  These are generally set from peaks.
    _mask = 0;
    for (std::vector<geom::Point2I>::const_iterator ii = stops.begin(); ii != stops.end(); ++ii) {
        _mask(ii->getX() - x0, ii->getY() - y0) |= STOP;
    }

    // Find starting span passing through (row, col)
    int xOrigin = point.getX() - x0;
    int xStart = xOrigin, y = point.getY() - y0;
    typename image::Image<PixelT>::const_x_iterator imageRow = image.x_at(xOrigin, y);
    typename image::Mask<image::MaskPixel>::const_x_iterator maskRow = _mask.x_at(xOrigin, y);
    for (; xStart >= 0; --xStart, --imageRow, --maskRow) {
        assert(_mask.getBBox(image::LOCAL).contains(geom::Point2I(xStart, y)));
        if (*maskRow & DETECTED || !threshold(*imageRow)) {
            break;
        }
    }
    int xStop = point.getX() - x0;
    imageRow = image.x_at(xOrigin, y);
    maskRow = _mask.x_at(xOrigin, y);
    for (; xStop < width; ++xStop, ++imageRow, ++maskRow) {
        assert(_mask.getBBox(image::LOCAL).contains(geom::Point2I(xStop, y)));
        if (*maskRow & DETECTED || !threshold(*imageRow)) {
            break;
        }
    }

    restart(found(y + y0, xStart + 1 + x0, xStop - 1 + x0));
}

// Marking a span as found involves recording it in the Footprint, and on the mask
template <bool throwOnStop, typename PixelT>
geom::Span const& SpanSearch<throwOnStop, PixelT>::found(int ySpan, int xSpanStart, int xSpanStop)
{
    geom::Span const& span = _footprint.addSpan(ySpan, xSpanStart, xSpanStop);

    int const x0 = _mask.getX0(), y0 = _mask.getY0();
    assert(xSpanStart >= x0 && xSpanStop < x0 + _mask.getWidth() &&
           ySpan >= y0 && ySpan < y0 + _mask.getHeight());

    // Record that we've detected these pixels
    int x = xSpanStart - x0;
    int const y = ySpan - y0;
    for (image::Mask<image::MaskPixel>::x_iterator ii = _mask.x_at(x, y); x <= xSpanStop - x0; ++x, ++ii) {
        *ii |= DETECTED;
        if (throwOnStop && *ii & STOP) {
            throw StopException();
        }
    }

    return span;
}


// Search the image for pixels above threshold, starting at a single StartingPoint.
// More StartingPoints will be added to the queue as we go along, which we will
// process on subsequent calls of this method.
//
// This is an 8-way scanline flood fill algorithm (https://en.wikipedia.org/wiki/Flood_fill#Scanline_fill),
// except the usual flood fill concern about colors are replaced by a concern about the threshold.
//
// This is the guts of findFootprintAtPoint.
template <bool throwOnStop, typename PixelT>
void SpanSearch<throwOnStop, PixelT>::process(StartingPoint const& sp)
{
    int const x0 = _image.getX0(), y0 = _image.getY0();
    int const width = _image.getWidth(), height = _image.getHeight();
    int xStartNext = sp.span.getX0() - x0, xStopNext = sp.span.getX1() - x0; // Span we'll work on next

    // Work up/down as far as we can go, looking for spans over threshold.
    int const di = sp.direction ? 1 : -1;  // How do we get to the next row?
    bool foundNext = true;            // Have we found a span for the next line?
    for (int i = sp.span.getY() - y0 + di; foundNext && i < height && i >= 0; i += di) {
        typename image::Image<PixelT>::const_x_iterator imageRow = _image.row_begin(i);
        typename image::Mask<image::MaskPixel>::const_x_iterator maskRow = _mask.row_begin(i);
        foundNext = false;
        int xStart = xStartNext, xStop = xStopNext; // Span for this iteration
        xStartNext = -1; xStopNext = -1;
        int spanStart = -1, spanStop = -1; // Start and stop x indices for working span

        // Search left from the pixel diagonally to the left of (i - di, xStart).
        int const startLeft = xStart - 1;
        imageRow = _image.x_at(startLeft, i);
        maskRow = _mask.x_at(startLeft, i);
        for (int j = startLeft; j >= -1; --j, --imageRow, --maskRow) {
            assert(j == -1 || _mask.getBBox(image::LOCAL).contains(geom::Point2I(j, i)));
            if (j < 0 || *maskRow & DETECTED || !_threshold(*imageRow)) {
                if (j < xStart - 1) {   // we found the end of a span of pixels above threshold
                    spanStart = j + 1;
                    foundNext = true;
                    xStartNext = spanStart;
                }
                break;
            }
        }

        if (foundNext) {
            // There's a span running to the left. Since it overhangs the previous line,
            // we'll need to push this onto the queue to process later. Now find where it ends.
            int const start = xStart;
            imageRow = _image.x_at(start, i);
            maskRow = _mask.x_at(start, i);
            for (int j = start; j <= width; ++j, ++imageRow, ++maskRow) {
                assert(j == width || _mask.getBBox(image::LOCAL).contains(geom::Point2I(j, i)));
                if (j >= width || *maskRow & DETECTED || !_threshold(*imageRow)) {
                    spanStop = j - 1;
                    break;
                }
            }

            add(found(i + y0, spanStart + x0, spanStop + x0), !sp.direction);
            xStopNext = spanStop;
        } else {
            // we'll resume searching at spanStop + 1
            spanStop = xStart - 1;
        }

        // Now search to the right for spans that are connected to the previous span.
        // For the first span we find that is entirely contained within the previous
        // span, we will use it to continue iterating up/down the image (we don't have to
        // check both up and down). Otherwise, we'll have to push it onto the queue to
        // process later.
        //
        // Note that column 'width' exists virtually, and always ends the last span; this
        // is why we claim below that xStopSpan is always set
        //

        int const startRight = spanStop + 1;
        imageRow = _image.x_at(startRight, i);
        maskRow = _mask.x_at(startRight, i);
        for (int j = startRight; j <= xStop + 1; ++j, ++imageRow, ++maskRow) {
            if (!(*maskRow & DETECTED) && j != width && _threshold(*imageRow)) {
                // We've found a span; now find where it ends
                spanStart = j;
                ++j; ++imageRow; ++maskRow;
                for (; j <= width; ++j, ++imageRow, ++maskRow) {
                    if ((j >= width) || (*maskRow & DETECTED) || !_threshold(*imageRow)) {
                        spanStop = j - 1;
                        break;
                    }
                }
                assert(spanStop >= 0);

                if (foundNext) {
                    // We've already got something to process on the next line, so put this on the queue.
                    if (spanStop > xStop) {
                        // Overhangs, so need to look in both directions
                        restart(found(i + y0, spanStart + x0, spanStop + x0));
                    } else {
                        // Entirely contained within the length of the previous span (i.e., no overhang),
                        // so only need to search in one direction
                        add(found(i + y0, spanStart + x0, spanStop + x0), sp.direction);
                    }
                } else {
                    // We'll use this span as the basis for checking the next line
                    foundNext = true;
                    xStartNext = spanStart;
                    xStopNext = spanStop;
                    geom::Span const& span = found(i + y0, spanStart + x0, spanStop + x0);
                    if (spanStop > xStop) {
                        // There's an overhang we'll come back to later
                        add(span, !sp.direction);
                    }
                }
            }
        }
    }
}



template <bool throwOnStop, typename PixelT>
void SpanSearch<throwOnStop, PixelT>::run()
{
    while (!_queue.empty()) {
        process(_queue.front());
        _queue.pop_front();
    }
}


} // anonymous namespace


template <typename PixelT>
PTR(Footprint)
findFootprintAtPoint(
    image::Image<PixelT> const& image,      // image to search
    geom::Point2I const& point, // starting position
    PixelT const thresholdValue, // Threshold level
    bool const polarity          // Polarity of threshold
    )
{
    if (!image.getBBox(image::PARENT).contains(point)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "Point does not lie in image");
    }
    PTR(Footprint) footprint = boost::make_shared<Footprint>();
    Thresholder<PixelT> threshold(thresholdValue, polarity);
    if (!threshold(image.get0(point))) {
        return footprint;
    }

    SpanSearch<false, PixelT> search(*footprint, image, point, threshold);
    search.run();
    return footprint;
}

template <typename PixelT>
bool checkFootprintAtPoint(
    image::Image<PixelT> const& image,      // image to search
    geom::Point2I const& point, // starting position
    PixelT const thresholdValue, // Threshold level
    bool const polarity,        // Polarity of threshold
    std::vector<geom::Point2I> const& stops // stopping points
    )
{
    if (!image.getBBox(image::PARENT).contains(point)) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "Point does not lie in image");
    }
    Thresholder<PixelT> threshold(thresholdValue, polarity);
    PTR(Footprint) footprint = boost::make_shared<Footprint>();
    if (!threshold(image.get0(point))) {
        return false;
    }

    try {
        SpanSearch<true, PixelT> searches(*footprint, image, point, threshold, stops);
        searches.run();
    } catch (StopException const&) {
        // We hit a stop
        return true;
    }
    return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Implicit instantiations

#define INSTANTIATE(PIXEL) \
template PTR(Footprint) \
findFootprintAtPoint<PIXEL>(image::Image<PIXEL> const&, \
                            geom::Point2I const&, \
                            PIXEL const, \
                            bool const \
    ); \
template bool \
checkFootprintAtPoint<PIXEL>(image::Image<PIXEL> const&, \
                             geom::Point2I const&, \
                             PIXEL const, \
                             bool const, \
                             std::vector<geom::Point2I> const&  \
    );


INSTANTIATE(float);

}}} // lsst::afw::detection
