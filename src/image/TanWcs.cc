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

#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>

#include "boost/format.hpp"

#include "wcslib/wcs.h"
#include "wcslib/wcsfix.h"
#include "wcslib/wcshdr.h"

#include "lsst/daf/base.h"
#include "lsst/daf/base/Citizen.h"
#include "lsst/afw/formatters/Utils.h"
#include "lsst/afw/formatters/TanWcsFormatter.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/WcsRecordFactory.h"

namespace lsst { namespace afw { namespace image {

const int lsstToFitsPixels = +1;
const int fitsToLsstPixels = -1;

TanWcs::TanWcs() :
    Wcs(),
    _hasDistortion(false),
    _sipA(), _sipB(), _sipAp(), _sipBp()
{}

geom::Angle TanWcs::pixelScale() const {
    // HACK -- assume "CD" elements are set (and are in degrees)
    double* cd = _wcsInfo->m_cd;
    assert(cd);
    return sqrt(fabs(cd[0]*cd[3] - cd[1]*cd[2])) * geom::degrees;
}

TanWcs::TanWcs(CONST_PTR(daf::base::PropertySet) const& fitsMetadata) :
    Wcs(fitsMetadata),
    _hasDistortion(false),
    _sipA(), _sipB(), _sipAp(), _sipBp() {

    //Internal params for wcslib. These should be set via policy - but for the moment...
    _relax = 1;
    _wcsfixCtrl = 2;
    _wcshdrCtrl = 2;

    //Check that the header isn't empty
    if(fitsMetadata->nameCount() == 0) {
        std::string msg = "Fits metadata contains no cards";
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, msg);
    }

    //Check for tangent plane projection
    std::string ctype1 = fitsMetadata->getAsString("CTYPE1");
    std::string ctype2 = fitsMetadata->getAsString("CTYPE2");

    if((ctype1.substr(5, 3) != "TAN") || (ctype2.substr(5, 3) != "TAN") ) {
        std::string msg = "One or more axes isn't in TAN projection (ctype1 = \"" + ctype1 + "\", ctype2 = \"" + ctype2 + "\")";
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, msg);
    }

    //Check for distorton terms. With two ctypes, there are 4 alternatives, only
    //two of which are valid.. Both have distortion terms or both don't.
    int nSip = (((ctype1.substr(8, 4) == "-SIP") ? 1 : 0) +
		((ctype2.substr(8, 4) == "-SIP") ? 1 : 0));

    switch (nSip) {
        case 0:
            _hasDistortion = false;
            break;
        case 1:
            {//Invalid case. Throw an exception
                std::string msg = "Distortion key found for only one CTYPE";
                throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, msg);
            }
            break;  //Not necessary, but looks naked without it.
        case 2:
            _hasDistortion = true;

            //Hide the distortion from wcslib
            //Make a copy that we can hack up
            PTR(daf::base::PropertySet) const& hackMetadata = fitsMetadata->deepCopy();
            hackMetadata->set<std::string>("CTYPE1", ctype1.substr(0,8));
            hackMetadata->set<std::string>("CTYPE2", ctype2.substr(0,8));

            //Save SIP information
            decodeSipHeader(*hackMetadata, "A", _sipA);
            decodeSipHeader(*hackMetadata, "B", _sipB);
            decodeSipHeader(*hackMetadata, "AP", _sipAp);
            decodeSipHeader(*hackMetadata, "BP", _sipBp);

            // this gets called in the Wcs (base class) constructor
            // We just changed fitsMetadata, so we have to re-init wcslib
            initWcsLibFromFits(hackMetadata);

            break;
    }

    //Check that the existence of forward sip matrices <=> existence of reverse matrices
    if (_hasDistortion) {
        if (_sipA.rows() <= 1 || _sipB.rows() <= 1) {
                std::string msg = "Existence of forward distorton matrices suggested, but not found";
                throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, msg);
        }

        if (_sipAp.rows() <= 1 || _sipBp.rows() <= 1) {
                std::string msg = "Forward distorton matrices present, but no reverse matrices";
                throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, msg);
        }
    }

}

void TanWcs::decodeSipHeader(
    daf::base::PropertySet const & fitsMetadata,
    std::string const& which,
    Eigen::MatrixXd & m
) {
    std::string header = which + "_ORDER";
    if (!fitsMetadata.exists(header)) return;
    int order = fitsMetadata.getAsInt(header);
    m.resize(order + 1, order + 1);
    boost::format format("%1%_%2%_%3%");
    for (int i = 0; i <= order; ++i) {
        for (int j = 0; j <= order; ++j) {
            header = (format % which % i % j).str();
            if (fitsMetadata.exists(header)) {
                m(i,j) = fitsMetadata.getAsDouble(header);
            }
            else {
                m(i, j) = 0.0;
            }
        }
    }
}



TanWcs::TanWcs(
    geom::Point2D const & crval,
    geom::Point2D const & crpix,
    Eigen::Matrix2d const & cd,
    double equinox,
    std::string const & raDecSys,
    std::string const & cunits1,
    std::string const & cunits2
) :
    Wcs(crval, crpix, cd, "RA---TAN", "DEC--TAN", equinox, raDecSys, cunits1, cunits2),
    _hasDistortion(false),
    _sipA(), _sipB(), _sipAp(), _sipBp()
{}

TanWcs::TanWcs(
    geom::Point2D const &crval,
    geom::Point2D const & crpix,
    Eigen::Matrix2d const & cd,
    Eigen::MatrixXd const & sipA,
    Eigen::MatrixXd const & sipB,
    Eigen::MatrixXd const & sipAp,
    Eigen::MatrixXd const & sipBp,
    double equinox, std::string const & raDecSys,
    std::string const & cunits1, std::string const & cunits2
) :
    Wcs(crval, crpix, cd, "RA---TAN", "DEC--TAN", equinox, raDecSys, cunits1, cunits2),
    _hasDistortion(true),
    //Sip's set by a dedicated method that does error checking
    _sipA(), _sipB(), _sipAp(), _sipBp()
{
    //Input checking is done constructor of base class, so don't need to do any here.

    //Set the distortion terms
    setDistortionMatrices(sipA, sipB, sipAp, sipBp);
}

TanWcs::TanWcs(TanWcs const & rhs) :
    Wcs(rhs),
    _hasDistortion(rhs._hasDistortion),
    _sipA(rhs._sipA),
    _sipB(rhs._sipB),
    _sipAp(rhs._sipAp),
    _sipBp(rhs._sipBp) {

}

bool TanWcs::_isSubset(Wcs const & rhs) const {
    if (!Wcs::_isSubset(rhs)) {
        return false;
    }
    // We only care about the derived-class part if we have a distortion; this could mean
    // a TanWcs with no distortion may be equal to a plain Wcs, but that doesn't happen
    // in practice because have different wcslib data structures.
    if (this->hasDistortion()) {
        TanWcs const * other = dynamic_cast<TanWcs const *>(&rhs);
        return other && other->_hasDistortion &&
            _sipA == other->_sipA &&
            _sipB == other->_sipB &&
            _sipAp == other->_sipAp &&
            _sipBp == other->_sipBp;
    }
    return true;
}

PTR(Wcs) TanWcs::clone(void) const {
    return PTR(Wcs)(new TanWcs(*this));
}

//
// Accessors
//
geom::Point2D TanWcs::skyToPixelImpl(
    geom::Angle sky1, // RA
    geom::Angle sky2  // Dec
) const {
    if(_wcsInfo == NULL) {
        throw(LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Wcs structure not initialised"));
    }

    double skyTmp[2];
    double imgcrd[2];
    double phi, theta;
    double pixTmp[2];

    //Estimate undistorted pixel coordinates
    int stat[1];
    int status = 0;

    skyTmp[_wcsInfo->lng] = sky1.asDegrees();
    skyTmp[_wcsInfo->lat] = sky2.asDegrees();

    status = wcss2p(_wcsInfo, 1, 2, skyTmp, &phi, &theta, imgcrd, pixTmp, stat);
    if (status > 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException,
            (boost::format("Error: wcslib returned a status code of %d at sky %s, %s deg: %s") %
            status % sky1.asDegrees() % sky2.asDegrees() % wcs_errmsg[status]).str());
    }


    //Correct for distortion. We follow the notation of Shupe et al. here, including
    //capitalisation
    if( _hasDistortion) {
        geom::Point2D pix = geom::Point2D(pixTmp[0], pixTmp[1]);
        geom::Point2D dpix = distortPixel(pix);
        pixTmp[0] = dpix[0];
        pixTmp[1] = dpix[1];
    }

    // wcslib assumes 1-indexed coords
    double offset = PixelZeroPos + fitsToLsstPixels;
    return geom::Point2D(pixTmp[0]+offset, pixTmp[1]+offset);

}

geom::Point2D TanWcs::undistortPixel(geom::Point2D const & pix) const {
    if (!_hasDistortion) {
        return geom::Point2D(pix);
    }
    //If the following assertions aren't true then something has gone seriously wrong.
    assert(_sipB.rows() > 0 );
    assert(_sipA.rows() == _sipA.cols());
    assert(_sipB.rows() == _sipB.cols());

    double u = pix[0] - _wcsInfo->crpix[0];  //Relative pixel coords
    double v = pix[1] - _wcsInfo->crpix[1];

    double f = 0;
    for(int i=0; i< _sipA.rows(); ++i) {
        for(int j=0; j< _sipA.cols(); ++j) {
            if (i+j>1 && i+j < _sipA.rows() ) {
                f += _sipA(i,j)* pow(u, i) * pow(v, j);
            }
        }
    }

    double g = 0;
    for(int i=0; i< _sipB.rows(); ++i) {
        for(int j=0; j< _sipB.cols(); ++j) {
            if (i+j>1 && i+j < _sipB.rows() ) {
                g += _sipB(i,j)* pow(u, i) * pow(v, j);
            }
        }
    }

    return geom::Point2D(pix[0] + f, pix[1] + g);
}

geom::Point2D TanWcs::distortPixel(geom::Point2D const & pix) const {
    if (!_hasDistortion) {
        return geom::Point2D(pix);
    }
    //If the following assertions aren't true then something has gone seriously wrong.
    assert(_sipBp.rows() > 0 );
    assert(_sipAp.rows() == _sipAp.cols());
    assert(_sipBp.rows() == _sipBp.cols());

    double U = pix[0] - _wcsInfo->crpix[0];  //Relative, undistorted pixel coords
    double V = pix[1] - _wcsInfo->crpix[1];

    double F = 0;
    for(int i=0; i< _sipAp.rows(); ++i) {
        for(int j=0; j< _sipAp.cols(); ++j) {
            F += _sipAp(i,j)* pow(U, i) * pow(V, j);
        }
    }

    double G = 0;
    for(int i=0; i< _sipBp.rows(); ++i) {
        for(int j=0; j< _sipBp.cols(); ++j) {
            G += _sipBp(i,j)* pow(U, i) * pow(V, j);
        }
    }
    return geom::Point2D(U + F + _wcsInfo->crpix[0],
                            V + G + _wcsInfo->crpix[1]);
}

/************************************************************************************************************/
/*
 * Worker routine for pixelToSky
 */
void
TanWcs::pixelToSkyImpl(double pixel1, double pixel2, geom::Angle sky[2]) const
{
    if(_wcsInfo == NULL) {
        throw(LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "Wcs structure not initialised"));
    }

    // wcslib assumes 1-indexed coordinates
    double pixTmp[2] = { pixel1 - PixelZeroPos + lsstToFitsPixels,
                         pixel2 - PixelZeroPos + lsstToFitsPixels};
    double imgcrd[2];
    double phi, theta;

    //Correct pixel positions for distortion if necessary
    if( _hasDistortion) {
        geom::Point2D pix = geom::Point2D(pixTmp[0], pixTmp[1]);
        geom::Point2D dpix = undistortPixel(pix);
        pixTmp[0] = dpix[0];
        pixTmp[1] = dpix[1];
    }

    int status = 0;
	double skyTmp[2];
    if (wcsp2s(_wcsInfo, 1, 2, pixTmp, imgcrd, &phi, &theta, skyTmp, &status) > 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
            (boost::format("Error: wcslib returned a status code of %d at pixel %s, %s: %s") %
            status % pixel1 % pixel2 % wcs_errmsg[status]).str());
    }
	sky[0] = skyTmp[0] * geom::degrees;
	sky[1] = skyTmp[1] * geom::degrees;
}

/************************************************************************************************************/

PTR(daf::base::PropertyList) TanWcs::getFitsMetadata() const {
    return formatters::TanWcsFormatter::generatePropertySet(*this);
}

//
// Mutators
//

void TanWcs::setDistortionMatrices(
    Eigen::MatrixXd const & sipA,
    Eigen::MatrixXd const & sipB,
    Eigen::MatrixXd const & sipAp,
    Eigen::MatrixXd const & sipBp
) {

    if (sipA.rows() != sipA.cols() ){
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Error: Matrix sipA must be square");
    }

    if (sipB.rows() != sipB.cols() ){
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Error: Matrix sipB must be square");
    }

    if (sipAp.rows() != sipAp.cols() ){
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Error: Matrix sipAp must be square");
    }

    if (sipBp.rows() != sipBp.cols() ){
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Error: Matrix sipBp must be square");
    }

    //Set the SIP terms
    _hasDistortion = true;
    _sipA = sipA;
    _sipB = sipB;
    _sipAp = sipAp;
    _sipBp = sipBp;
}

// -------------- Record Persistence ------------------------------------------------------------------------

/*
 *  We use the Wcs base class persistence to write one table, and then add another containing
 *  the SIP coefficients only if hasDistortion() is true.
 *
 *  The second table's schema depends on the SIP orders, so it will not necessarily be the same
 *  for all TanWcs objects.
 */
class TanWcsRecordOutputGenerator : public table::RecordOutputGenerator {
public:

    TanWcsRecordOutputGenerator(TanWcs const & wcs) :
        table::RecordOutputGenerator(table::Schema(), 1),
        _wcs(&wcs),
        _keyA(_schema.addField< table::Array<double> >(
                  "A", "x forward transform coefficients (column-major)", _wcs->_sipA.size()
              )),
        _keyB(_schema.addField< table::Array<double> >(
                  "B", "y forward transform coefficients (column-major)", _wcs->_sipB.size()
              )),
        _keyAp(_schema.addField< table::Array<double> >(
                   "Ap", "x reverse transform coefficients (column-major)", _wcs->_sipAp.size()
               )),
        _keyBp(_schema.addField< table::Array<double> >(
                   "Bp", "y reverse transform coefficients (column-major)", _wcs->_sipBp.size()
               ))
        {}

    virtual void fill(table::BaseRecord & record) {
        Eigen::Map<Eigen::MatrixXd> mapA(record[_keyA].getData(), _wcs->_sipA.rows(), _wcs->_sipA.cols());
        mapA = _wcs->_sipA;
        Eigen::Map<Eigen::MatrixXd> mapB(record[_keyB].getData(), _wcs->_sipB.rows(), _wcs->_sipB.cols());
        mapB = _wcs->_sipB;
        Eigen::Map<Eigen::MatrixXd> mapAp(record[_keyAp].getData(), _wcs->_sipAp.rows(), _wcs->_sipAp.cols());
        mapAp = _wcs->_sipAp;
        Eigen::Map<Eigen::MatrixXd> mapBp(record[_keyBp].getData(), _wcs->_sipBp.rows(), _wcs->_sipBp.cols());
        mapBp = _wcs->_sipBp;
    }

private:
    TanWcs const * _wcs;
    table::Key< table::Array<double> > _keyA;
    table::Key< table::Array<double> > _keyB;
    table::Key< table::Array<double> > _keyAp;
    table::Key< table::Array<double> > _keyBp;
};

afw::table::RecordOutputGeneratorSet TanWcs::writeToRecords() const {
    afw::table::RecordOutputGeneratorSet result = Wcs::writeToRecords();
    result.name = "TAN";
    if (hasDistortion()) {
        result.generators.push_back(
            boost::make_shared<TanWcsRecordOutputGenerator>(*this)
        );
    }
    return result;
}

class TanWcsRecordFactory : public WcsRecordFactory {
public:

    explicit TanWcsRecordFactory(std::string const & name) :
        WcsRecordFactory(name) {}

    virtual PTR(Wcs) operator()(table::RecordInputGeneratorSet const & inputs) const {
        CONST_PTR(table::BaseRecord) sipRecord;
        if (inputs.generators.size() > 1u) {
            sipRecord = inputs.generators.back()->next();
        }
        PTR(Wcs) result(new TanWcs(*inputs.generators.front()->next(), sipRecord));
        return result;
    }
};

TanWcs::TanWcs(
    afw::table::BaseRecord const & mainRecord,
    CONST_PTR(afw::table::BaseRecord) sipRecord
) : Wcs(mainRecord), _hasDistortion(sipRecord)
{
    if (_hasDistortion) {
        typedef afw::table::Array<double> T;
        afw::table::SchemaItem<T> sA = sipRecord->getSchema().find<T>("A");
        afw::table::SchemaItem<T> sB = sipRecord->getSchema().find<T>("B");
        afw::table::SchemaItem<T> sAp = sipRecord->getSchema().find<T>("Ap");
        afw::table::SchemaItem<T> sBp = sipRecord->getSchema().find<T>("Bp");
        // Adding 0.5 and truncating the result here guarantees we'll get the right answer
        // for small ints even when round-off error is involved.
        int nA = int(std::sqrt(sA.field.getSize() + 0.5));
        int nB = int(std::sqrt(sB.field.getSize() + 0.5));
        int nAp = int(std::sqrt(sAp.field.getSize() + 0.5));
        int nBp = int(std::sqrt(sBp.field.getSize() + 0.5));
        if (nA * nA != sA.field.getSize()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Forward X SIP matrix is not square."
            );
        }
        if (nB * nB != sB.field.getSize()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Forward Y SIP matrix is not square."
            );
        }
        if (nAp * nAp != sAp.field.getSize()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Reverse X SIP matrix is not square."
            );
        }
        if (nBp * nBp != sBp.field.getSize()) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Reverse Y SIP matrix is not square."
            );
        }
        Eigen::Map<Eigen::MatrixXd const> mapA((*sipRecord)[sA.key].getData(), nA, nA);
        _sipA = mapA;
        Eigen::Map<Eigen::MatrixXd const> mapB((*sipRecord)[sB.key].getData(), nB, nB);
        _sipB = mapB;
        Eigen::Map<Eigen::MatrixXd const> mapAp((*sipRecord)[sAp.key].getData(), nAp, nAp);
        _sipAp = mapAp;
        Eigen::Map<Eigen::MatrixXd const> mapBp((*sipRecord)[sBp.key].getData(), nBp, nBp);
        _sipBp = mapBp;
    }
}


namespace {

TanWcsRecordFactory registration("TAN");

} // anonymous

}}} // namespace lsst::afw::image
