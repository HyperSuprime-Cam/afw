// -*- LSST-C++ -*-
/**
 * \file
 *
 * \brief Definition of templated functions declared in KernelFunctions.h
 *
 * This file is meant to be included by lsst/fw/KernelFunctions.h
 *
 * \author Russell Owen
 *
 * \ingroup fw
 */
#include <iostream>

#include <boost/format.hpp>
#include <vw/Image.h>

#include <lsst/mwi/utils/Trace.h>
#include <lsst/fw/ImageUtils.h>

namespace mwiu = lsst::mwi::utils;

/**
 * \brief Apply convolution kernel to a masked image at one point
 *
 * Note: this is a high performance routine; the user is expected to:
 * - handle edge extension
 * - figure out the kernel center and adjust the supplied pixel accessors accordingly
 * For an example of how to do this see the convolve function.
 *
 * \ingroup fw
 */
template <typename ImageT, typename MaskT, typename KernelT>
inline void lsst::fw::kernel::apply(
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> &outAccessor,    ///< accessor for output pixel
    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> const &imageAccessor,
        ///< accessor to for image pixel that overlaps (0,0) pixel of kernel(!)
    typename lsst::fw::Image<KernelT>::pixel_accessor const &kernelAccessor,
        ///< accessor for (0,0) pixel of kernel
    unsigned int cols,      ///< number of columns in kernel
    unsigned int rows,      ///< number of rows in kernel
    KernelT threshold   ///< if a kernel pixel > threshold then the corresponding image mask pixel
                        ///< is ORd into the output pixel, otherwise the mask pixel is ignored
) {
    typedef typename lsst::fw::Image<KernelT>::pixel_accessor kernelAccessorType;
    *outAccessor.image = 0;
    *outAccessor.variance = 0;
    *outAccessor.mask = 0;

    lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imRow = imageAccessor;
    kernelAccessorType kRow = kernelAccessor;
    for (unsigned int row = 0; row < rows; ++row) {
        MaskedPixelAccessor<ImageT, MaskT> imCol = imRow;
        kernelAccessorType kCol = kRow;
        for (unsigned int col = 0; col < cols; ++col) {
            *outAccessor.image += (*kCol) * (*imCol.image);
            *outAccessor.variance += (*kCol) * (*kCol) * (*imCol.variance);
            if ((*imCol.mask) && (*kCol > threshold)) {
                // this bad pixel contributes enough to "OR" in the bad bits
                *outAccessor.mask |= *imCol.mask;
            }
            
            imCol.nextCol();
            kCol.next_col();
        }
        imRow.nextRow();
        kRow.next_row();
    }
}

/**
 * \brief Convolve a MaskedImage with a Kernel, setting pixels of an existing image
 *
 * convolvedImage must be the same size as maskedImage.
 * It has a border in which the output pixels are just a copy of the input pixels
 * and the output mask bit edgeBit has been set. This border will have size:
 * * kernel.getCtrCol/Row() along the left/bottom edge
 * * kernel.getCols/Rows() - 1 - kernel.getCtrCol/Row() along the right/top edge
 *
 * \throw lsst::mwi::exceptions::InvalidParameter if convolvedImage is not the same size as maskedImage.
 * \throw lsst::mwi::exceptions::InvalidParameter if maskedImage is smaller (in colums or rows) than kernel.
 *
 * \ingroup fw
 */
template <typename ImageT, typename MaskT, typename KernelT>
void lsst::fw::kernel::convolve(
    lsst::fw::MaskedImage<ImageT, MaskT> &convolvedImage,       ///< convolved image
    lsst::fw::MaskedImage<ImageT, MaskT> const &maskedImage,    ///< image to convolve
    lsst::fw::Kernel<KernelT> const &kernel,    ///< convolution kernel
    KernelT threshold,  ///< if kernel pixel > threshold then corresponding maskedImage mask pixel is OR'd in
    int edgeBit         ///< mask bit to indicate pixel includes edge-extended data;
                        ///< if negative then no bit is set
) {
    typedef lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageAccessorType;
    typedef typename lsst::fw::Image<KernelT>::pixel_accessor kernelAccessorType;

    const unsigned int imCols = maskedImage.getCols();
    const unsigned int imRows = maskedImage.getRows();
    const unsigned int kCols = kernel.getCols();
    const unsigned int kRows = kernel.getRows();
    const unsigned int kCtrCol = kernel.getCtrCol();
    const unsigned int kCtrRow = kernel.getCtrRow();
    if ((convolvedImage.getCols() != imCols) || (convolvedImage.getRows() != imRows)) {
        throw lsst::mwi::exceptions::InvalidParameter("convolvedImage not the same size as maskedImage");
    }
    if ((imCols< kCols) || (imRows < kRows)) {
        throw lsst::mwi::exceptions::InvalidParameter("maskedImage smaller than kernel in columns and/or rows");
    }

    const unsigned int cnvCols = imCols + 1 - kernel.getCols();
    const unsigned int cnvRows = imRows + 1 - kernel.getRows();

    // create input and output image accessors
    // and advance output accessor to lower left pixel that is set by convolution
    imageAccessorType imRow(maskedImage);
    imageAccessorType outRow(convolvedImage);
    outRow.advance(kernel.getCtrCol(), kernel.getCtrRow());
    
    if (kernel.isSpatiallyVarying()) {
        lsst::fw::Image<KernelT> kernelImage(kernel.getCols(), kernel.getRows());
        kernelAccessorType kernelAccessor = kernelImage.origin();
        mwiu::Trace("lsst.fw.kernel.convolve", 1, "kernel is spatially varying");
        for (int row = 0; row < static_cast<int>(cnvRows); ++row) {
            double rowPos = lsst::fw::image::indexToPosition(row);
            imageAccessorType imCol = imRow;
            imageAccessorType outCol = outRow;
            for (int col = 0; col < static_cast<int>(cnvCols); ++col) {
                KernelT kSum;
                kernel.computeImage(
                    kernelImage, kSum, lsst::fw::image::indexToPosition(col), rowPos, false);
                KernelT adjThreshold = threshold * kSum;
                lsst::fw::kernel::apply(outCol, imCol, kernelAccessor, kCols, kRows, adjThreshold);
                *outCol.image /= static_cast<ImageT>(kSum);
                *outCol.variance /= static_cast<ImageT>(kSum * kSum);
                outCol.nextCol();
                imCol.nextCol();
            }
            outRow.nextRow();
            imRow.nextRow();
        }
    } else {
        mwiu::Trace("lsst.fw.kernel.convolve", 1, "kernel is spatially invariant");
        KernelT kSum;
        lsst::fw::Image<KernelT> kernelImage = kernel.computeNewImage(kSum, 0.0, 0.0, true);
        kernelAccessorType kernelAccessor = kernelImage.origin();
        for (unsigned int row = 0; row < cnvRows; ++row) {
            imageAccessorType imCol = imRow;
            imageAccessorType outCol = outRow;
            for (unsigned int col = 0; col < cnvCols; ++col) {
                lsst::fw::kernel::apply(outCol, imCol, kernelAccessor, kCols, kRows, threshold);
                outCol.nextCol();
                imCol.nextCol();
            }
            outRow.nextRow();
            imRow.nextRow();
        }
    }

    // set edge pixels
    MaskT edgeOrVal = edgeBit < 0 ? 0 : 1 << edgeBit;
    
    vw::BBox2i bottomEdge(0, 0, imCols, kCtrRow);
    lsst::fw::kernel::_copyRegion(convolvedImage, maskedImage, bottomEdge, edgeOrVal);
    
    vw::int32 numRows = kCols - (1 + kCtrRow);
    vw::BBox2i topEdge(0, imRows - numRows, imCols, numRows);
    lsst::fw::kernel::_copyRegion(convolvedImage, maskedImage, topEdge, edgeOrVal);

    vw::BBox2i leftEdge(0, kCtrRow, kCtrCol, imRows + 1 - kRows);
    lsst::fw::kernel::_copyRegion(convolvedImage, maskedImage, leftEdge, edgeOrVal);
    
    vw::int32 numCols = kCols - (1 + kernel.getCtrCol());
    vw::BBox2i rightEdge(imCols - numCols, kCtrRow, numCols, imRows + 1 - kRows);
    lsst::fw::kernel::_copyRegion(convolvedImage, maskedImage, rightEdge, edgeOrVal);
}

/**
 * \brief Convolve a MaskedImage with a Kernel, returning a new image.
 *
 * \return the convolved MaskedImage.
 *
 * The returned MaskedImage is the same size as maskedImage.
 * It has a border in which the output pixels are just a copy of the input pixels
 * and the output mask bit edgeBit has been set. This border will have size:
 * * kernel.getCtrCol/Row() along the left/bottom edge
 * * kernel.getCols/Rows() - 1 - kernel.getCtrCol/Row() along the right/top edge
 *
 * \throw lsst::mwi::exceptions::InvalidParameter if maskedImage is smaller (in colums or rows) than kernel.
 *
 * \ingroup fw
 */
template <typename ImageT, typename MaskT, typename KernelT>
lsst::fw::MaskedImage<ImageT, MaskT> lsst::fw::kernel::convolve(
    lsst::fw::MaskedImage<ImageT, MaskT> const &maskedImage,    ///< image to convolve
    lsst::fw::Kernel<KernelT> const &kernel,    ///< convolution kernel
    KernelT threshold,  ///< if kernel pixel > threshold then corresponding maskedImage mask pixel is OR'd in
    int edgeBit         ///< mask bit to indicate pixel includes edge-extended data;
                        ///< if negative then no bit is set
) {
    lsst::fw::MaskedImage<ImageT, MaskT> convolvedImage(maskedImage.getCols(), maskedImage.getRows());
    lsst::fw::kernel::convolve(convolvedImage, maskedImage, kernel, threshold, edgeBit);
    return convolvedImage;
}

/**
 * \brief Print the pixel values of a kernel to std::cout
 *
 * Rows increase upward and columns to the right; thus the lower left pixel is (0,0).
 *
 * \ingroup fw
 */
template <typename PixelT>
void lsst::fw::kernel::printKernel(
    lsst::fw::Kernel<PixelT> const &kernel,    ///< the kernel
    double x,   ///< x at which to evaluate kernel
    double y,   ///< y at which to evaluate kernel
    bool doNormalize,   ///< if true, normalize kernel
    string pixelFmt     ///< format for pixel values
) {
    typedef typename lsst::fw::Image<PixelT>::pixel_accessor imageAccessorType;
    PixelT kSum;
    lsst::fw::Image<PixelT> kImage = kernel.computeNewImage(kSum, x, y, doNormalize);
    imageAccessorType imRow = kImage.origin();
    imRow.advance(0, kImage.getRows()-1);
    for (unsigned int row=0; row < kImage.getRows(); ++row) {
        imageAccessorType imCol = imRow;
        for (unsigned int col = 0; col < kImage.getCols(); ++col) {
            std::cout << boost::format(pixelFmt) % (*imCol);
            imCol.next_col();
        }
        std::cout << std::endl;
        imRow.prev_row();
    }
    if (doNormalize && abs(static_cast<double>(kSum) - 1.0) > 1.0e-5) {
        std::cout << boost::format("Warning! Sum of all pixels = %9.5f != 1.0\n") % kSum;
    }
    std::cout << std::endl;
}

/**
 * \brief Private function to copy a rectangular region from one MaskedImage to another
 *
 * I hope eventually to replace this by calls to MaskedImage.getSubImage
 * and MaskedImage.replaceSubImage, but that is currently too messy
 * because getSubImage requires a shared pointer to the source image.
 *
 * \throw invalid_argument if the region extends off of either image.
 */
template <typename ImageT, typename MaskT>
inline void lsst::fw::kernel::_copyRegion(
    typename lsst::fw::MaskedImage<ImageT, MaskT> &destImage,           ///< destination MaskedImage
    typename lsst::fw::MaskedImage<ImageT, MaskT> const &sourceImage,   ///< source MaskedImage
    vw::BBox2i const &region,   ///< region to copy
    MaskT orMask    ///< data to "or" into the mask pixels
) {
    typedef lsst::fw::MaskedPixelAccessor<ImageT, MaskT> imageAccessorType;

    vw::math::Vector<vw::int32> const startColRow = region.min();
    vw::math::Vector<vw::int32> const numColRow = region.size();
    mwiu::Trace("lsst.fw.kernel._copyRegion", 1, str(boost::format(
        "_copyRegion: dest size %d, %d; src size %d, %d; region start=%d, %d; region size=%d, %d; orMask=%d")
        % destImage.getCols() % destImage.getRows() % sourceImage.getCols() % sourceImage.getRows()
        % startColRow[0] % startColRow[1]% numColRow[0] % numColRow[1] % orMask
    ));

    vw::math::Vector<vw::int32> const endColRow = region.max();
    if ((static_cast<unsigned int>(endColRow[0]) > min(destImage.getCols(), sourceImage.getCols()))
        || ((static_cast<unsigned int>(endColRow[1]) > min(destImage.getRows(), sourceImage.getRows())))) {
        throw lsst::mwi::exceptions::InvalidParameter("Region out of range");
    }
    imageAccessorType inRow(sourceImage);
    imageAccessorType outRow(destImage);
    inRow.advance(startColRow[0], startColRow[1]);
    outRow.advance(startColRow[0], startColRow[1]);
    for (int row = 0; row < numColRow[1]; ++row) {
        imageAccessorType inCol = inRow;
        imageAccessorType outCol = outRow;
        for (int col = 0; col < numColRow[0]; ++col) {
            *(outCol.image) = *(inCol.image);
            *(outCol.variance) = *(inCol.variance);
            *(outCol.mask) = *(inCol.mask) | orMask;
            inCol.nextCol();
            outCol.nextCol();
        }
        inRow.nextRow();
        outRow.nextRow();
    }
}

//
// Explicit instantiations
//
//template void lsst::fw::kernel::printKernel<float>(
//    lsst::fw::Kernel<PixelT> const &kernel,
//    double x,
//    double y,
//    bool doNormalize,
//    string pixelFmt);
//template void lsst::fw::kernel::printKernel<double>(
//    lsst::fw::Kernel<PixelT> const &kernel,
//    double x,
//    double y,
//    bool doNormalize,
//    string pixelFmt);
