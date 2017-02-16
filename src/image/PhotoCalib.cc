/*
 * LSST Data Management System
 * Copyright 2017 LSST Corporation.
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

#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/table/Source.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/pex/exceptions.h"

namespace lsst {
namespace afw {
namespace image {


PhotoCalib::PhotoCalib()
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

PhotoCalib::PhotoCalib(double fluxMag0, double fluxMag0Sigma)
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

PhotoCalib::PhotoCalib(std::shared_ptr<math::BoundedField> zeroPoint, double fluxMag0Sigma)
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

double PhotoCalib::countsToMaggies(double counts,
                       const afw::geom::Point<double, 2> &point,
                       bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair<double, double> PhotoCalib::countsToMaggies(double counts, double countsSigma,
                                          const afw::geom::Point<double, 2> &point,
                                          bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair<double, double> PhotoCalib::countsToMaggies(const afw::table::SourceRecord &sourceRecord,
                                          const std::string &fluxField,
                                          bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair< ndarray::Array<double,1>, ndarray::Array<double,1> >
    PhotoCalib::countsToMaggies(const afw::table::SourceCatalog &sourceCatalog,
                    const std::string &fluxField,
                    bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

void PhotoCalib::countsToMaggies(const afw::table::SourceCatalog &sourceCatalog,
                     const std::string &fluxField,
                     ndarray::Array<double,1> maggies,
                     ndarray::Array<double,1> maggiesSigma,
                     bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

void PhotoCalib::countsToMaggies(afw::table::SourceCatalog &sourceCatalog,
                     const std::string &fluxField,
                     const std::string &outField,
                     bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}


double PhotoCalib::countsToMagnitude(double counts,
                         const afw::geom::Point<double, 2> &point,
                         bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair<double, double> PhotoCalib::countsToMagnitude(double counts, double countsSigma,
                                            const afw::geom::Point<double, 2> &point,
                                            bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair<double, double> PhotoCalib::countsToMagnitude(const afw::table::SourceRecord &sourceRecord,
                                            const std::string &fluxField,
                                            bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::pair< ndarray::Array<double,1>, ndarray::Array<double,1> >
    PhotoCalib::countsToMagnitude(const afw::table::SourceCatalog &sourceCatalog,
                      const std::string &fluxField,
                      bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

void PhotoCalib::countsToMagnitude(const afw::table::SourceCatalog &sourceCatalog,
                       const std::string &fluxField,
                       ndarray::Array<double,1> mag,
                       ndarray::Array<double,1> magSigma,
                       bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

void PhotoCalib::countsToMagnitude(afw::table::SourceCatalog &sourceCatalog,
                       const std::string &fluxField,
                       const std::string &outField,
                       bool throwOnNegativeFlux) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

double PhotoCalib::magnitudeToCounts(double magnitude) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

double PhotoCalib::getFluxMag0() const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::shared_ptr<math::BoundedField> PhotoCalib::computeScaledZeroPoint() const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::shared_ptr<math::BoundedField> PhotoCalib::computeScalingTo(std::shared_ptr<PhotoCalib> other) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

bool PhotoCalib::operator==(PhotoCalib const& rhs) const
{
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

}}}
