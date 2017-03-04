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

/**
 * @brief      The photometric calibration of an exposure.
 *
 * A PhotoCalib is a BoundedField (a function with a specified domain) that converts between calibrated
 * counts-on-chip (ADU) to flux and magnitude. It is defined in terms of "maggies", which are a linear
 * unit defined in SDSS: http://www.sdss.org/dr12/algorithms/magnitudes/#nmgy
 *
 * PhotoCalib is immutable.
 *
 * The spatially varying flux/magnitude zero point is defined such that,
 * at a position (x,y) in the domain of the boundedField zeroPoint
 * and for a given measured source counts:
 *     zeroPoint(x,y) * counts = flux (in maggies)
 * while the errors (constant on the domain) are defined as:
 *     sqrt(countsSigma^2 + zeroPointSigma^2) = fluxSigma (in maggies)
 */

#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/table/Source.h"

namespace lsst { namespace afw { namespace image {

class PhotoCalib {
public:
    PhotoCalib(PhotoCalib const &) = delete;
    PhotoCalib(PhotoCalib &&) = delete;
    PhotoCalib & operator=(PhotoCalib const &) = delete;
    PhotoCalib & operator=(PhotoCalib &&) = delete;

    /**
     * @brief      Create a empty, zeroed calibration.
     */
    PhotoCalib(): _fluxMag0(0.0), _fluxMag0Sigma(0.0) {}

    /**
     * @brief      Create a non-spatially-varying calibration.
     *
     * @param[in]  fluxMag0       The constant flux/magnitude zero point (counts at magnitude 0).
     * @param[in]  fluxMag0Sigma  The error on the zero point.
     */
    explicit PhotoCalib(double fluxMag0, double fluxMag0Sigma=0):
        _fluxMag0(fluxMag0), _fluxMag0Sigma(fluxMag0Sigma) {}


    /**
     * @brief      Create a spatially-varying calibration.
     *
     * @param[in]  zeroPoint      The spatially varying photometric zero point.
     * @param[in]  fluxMag0Sigma  The error on the zero point.
     */
    explicit PhotoCalib(std::shared_ptr<afw::math::BoundedField> zeroPoint, double fluxMag0Sigma=0);

    //@{
    /**
     * @brief      Convert counts in ADU to maggies.
     *
     * If passed point, use the exact calculation at that point, otherwise, use the mean scaling factor.
     *
     * @param[in]  counts  The source flux in ADU.
     * @param[in]  point   The point that flux is measured at (must be within the domain of the
     *                     BoundedField of this PhotoCalib).
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The flux in maggies.
     */
    double countsToMaggies(double counts,
                           const afw::geom::Point<double, 2> &point,
                           bool throwOnNegativeFlux=false) const;
    double countsToMaggies(double counts,
                           bool throwOnNegativeFlux=false) const;
    //@}

    //@{}
    /**
     * @brief      Convert counts and error in counts (ADU) to maggies and maggies error.
     *
     * If passed point, use the exact calculation at that point, otherwise, use the mean scaling factor.
     *
     * @param[in]  counts       The source flux in ADU.
     * @param[in]  countsSigma  The counts error (sigma).
     * @param[in]  point        The point that flux is measured at (must be within the domain of the
     *                          BoundedField of this PhotoCalib).
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The flux in maggies and error (sigma).
     */
    std::pair<double, double> countsToMaggies(double counts, double countsSigma,
                                              const afw::geom::Point<double, 2> &point,
                                              bool throwOnNegativeFlux=false) const;
    std::pair<double, double> countsToMaggies(double counts, double countsSigma,
                                              bool throwOnNegativeFlux=false) const;
    //@}

    /**
     * @brief      Convert sourceRecord[fluxField] (ADU) at location
     *             (sourceRecord.get('x'), sourceRecord.get('y')) (pixels) to maggies and maggie error.
     *
     * @param[in]  sourceRecord  The source record to get flux and position from.
     * @param[in]  fluxField     The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                           For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The flux in maggies and error (sigma) for this source.
     */
    std::pair<double, double> countsToMaggies(const afw::table::SourceRecord &sourceRecord,
                                              const std::string &fluxField,
                                              bool throwOnNegativeFlux=false) const;

    //@{
    /**
     * @brief      Convert sourceCatalog[fluxField] (ADU) at locations
     *             (sourceCatalog.get('x'), sourceCatalog.get('y')) (pixels) to maggies.
     *
     * @param[in]  sourceCatalog  The source catalog to get flux and position from.
     * @param[in]  fluxField      The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                            For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[out] maggies, maggiesSigma  Pass pre-allocated arrays to put the output in, instead of
     *                                    returning it.
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The flux in maggies and error (sigma) for this source.
     */
    std::pair< ndarray::Array<double,1>, ndarray::Array<double,1> >
        countsToMaggies(const afw::table::SourceCatalog &sourceCatalog,
                        const std::string &fluxField,
                        bool throwOnNegativeFlux=false) const;
    void countsToMaggies(const afw::table::SourceCatalog &sourceCatalog,
                         const std::string &fluxField,
                         ndarray::Array<double,1> maggies,
                         ndarray::Array<double,1> maggiesSigma,
                         bool throwOnNegativeFlux=false) const;
    //@}

    /**
     * @brief      Convert sourceCatalog[fluxField_flux] (ADU) at locations
     *             (sourceCatalog.get('x'), sourceCatalog.get('x')) (pixels) to maggies
     *             and write the results back to sourceCatalog[outField_mag].
     *
     * @param[in]  sourceCatalog  The source catalog to get flux and position from.
     * @param[in]  fluxField      The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                            For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[in]  outField       The field to write the maggies and maggie errors to.
     *                            Keys of the form "*_flux" and "*_fluxSigma" must exist in the schema.
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     */
    void countsToMaggies(afw::table::SourceCatalog &sourceCatalog,
                         const std::string &fluxField,
                         const std::string &outField,
                         bool throwOnNegativeFlux=false) const;

    //@{
    /**
     * @brief      Convert counts in ADU to AB magnitude.
     *
     * If passed point, use the exact calculation at that point, otherwise, use the mean scaling factor.
     *
     * @param[in]  counts  The source flux in ADU.
     * @param[in]  point   The point that flux is measured at (must be within the domain of the
     *                     BoundedField of this PhotoCalib).
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The AB magnitude.
     */
    double countsToMagnitude(double counts,
                             const afw::geom::Point<double, 2> &point,
                             bool throwOnNegativeFlux=false) const;
    double countsToMagnitude(double counts,
                             bool throwOnNegativeFlux=false) const;
    //@}

    //@{
    /**
     * @brief      Convert counts and error in counts (ADU) to AB magnitude and magnitude error.
     *
     * If passed point, use the exact calculation at that point, otherwise, use the mean scaling factor.
     *
     * @param[in]  counts       The source flux in ADU.
     * @param[in]  countsSigma  The counts error (standard deviation).
     * @param[in]  point        The point that flux is measured at (must be within the domain of the
     *                          BoundedField of this PhotoCalib).
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The AB magnitude and error (sigma).
     */
    std::pair<double, double> countsToMagnitude(double counts, double countsSigma,
                                                const afw::geom::Point<double, 2> &point,
                                                bool throwOnNegativeFlux=false) const;
    std::pair<double, double> countsToMagnitude(double counts, double countsSigma,
                                                bool throwOnNegativeFlux=false) const;
    //@}

    /**
     * @brief      Convert sourceRecord[fluxField] (ADU) at location
     *             (sourceRecord.get('x'), sourceRecord.get('y')) (pixels) to AB magnitude.
     *
     * @param[in]  sourceRecord  The source record to get flux and position from.
     * @param[in]  fluxField     The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                           For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The magnitude and magnitude error for this source.
     */
    std::pair<double, double> countsToMagnitude(const afw::table::SourceRecord &sourceRecord,
                                                const std::string &fluxField,
                                                bool throwOnNegativeFlux=false) const;

    //@{
    /**
     * @brief      Convert sourceCatalog[fluxField] (ADU) at locations
     *             (sourceCatalog.get('x'), sourceCatalog.get('y')) (pixels) to AB magnitudes.
     *
     * @param[in]  sourceCatalog  The source catalog to get flux and position from.
     * @param[in]  fluxField      The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                            For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[out] mag, magSigma  Pass pre-allocated arrays to put the output in, instead of
     *                            returning it.
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     *
     * @return     The magnitudes and magnitude errors for the sources.
     */
    std::pair< ndarray::Array<double,1>, ndarray::Array<double,1> >
        countsToMagnitude(const afw::table::SourceCatalog &sourceCatalog,
                          const std::string &fluxField,
                          bool throwOnNegativeFlux=false) const;
    void countsToMagnitude(const afw::table::SourceCatalog &sourceCatalog,
                           const std::string &fluxField,
                           ndarray::Array<double,1> mag,
                           ndarray::Array<double,1> magSigma,
                           bool throwOnNegativeFlux=false) const;
    //@}

    /**
     * @brief      Convert sourceCatalog[fluxField_flux] (ADU) at locations
     *             (sourceCatalog.get('x'), sourceCatalog.get('x')) (pixels) to AB magnitudes
     *             and write the results back to sourceCatalog[outField_mag].
     *
     * @param[in]  sourceCatalog  The source catalog to get flux and position from.
     * @param[in]  fluxField      The flux field: Keys of the form "*_flux" and "*_fluxSigma" must exist.
     *                            For example: fluxField = "PsfFlux" -> "PsfFlux_flux", "PsfFlux_fluxSigma"
     * @param[in]  outField       The field to write the magnitudes and magnitude errors to.
     *                            Keys of the form "*_mag" and "*_magSigma" must exist in the schema.
     * @param[in]  throwOnNegativeFlux  Raise an exception when passed negative counts.
     */
    void countsToMagnitude(afw::table::SourceCatalog &sourceCatalog,
                           const std::string &fluxField,
                           const std::string &outField,
                           bool throwOnNegativeFlux=false) const;

    /**
     * @brief      Convert AB magnitude to counts (ADU), using the mean flux/magnitude scaling factor.
     *
     * @param[in]  magnitude  The AB magnitude to convert.
     *
     * @return     Source counts in ADU.
     */
    double magnitudeToCounts(double magnitude) const;

    /**
     * @brief      Get the mean flux/magnitude zero point.
     *
     * This value is defined, for counts at (x,y), such that:
     *   getFluxMag0() * counts * computeScaledZeroPoint()(x,y) = countsToMaggies(counts, (x,y))
     *
     * @see PhotoCalib::computeScaledZeroPoint(), getFluxMag0Sigma()
     *
     * @return     The flux magnitude zero point.
     */
    double getFluxMag0() const {return _fluxMag0;}

    /**
     * @brief      Get the mean flux/magnitude zero point error.
     *
     * This value is defined such that for some countsSigma:
     *     sqrt(countsSigma^2 + zeroPointSigma^2) = fluxSigma (in maggies)
     *
     * @see PhotoCalib::computeScaledZeroPoint(), getFluxMag0()
     *
     * @return     The flux magnitude zero point error.
     */
    double getFluxMag0Sigma() const {return _fluxMag0Sigma;}

    /**
     * @brief      Calculates the spatially-variable zero point, normalized by the mean in the valid domain.
     *
     * This value is defined, for counts at (x,y), such that:
     *   getFluxMag0() * counts * computeScaledZeroPoint()(x,y) = countsToMaggies(counts, (x,y))
     *
     * @see PhotoCalib::getFluxMag0()
     *
     * @return     The normalized spatially-variable zero point.
     */
    std::shared_ptr<afw::math::BoundedField> computeScaledZeroPoint() const;

    /**
     * @brief      Calculates the scaling between this PhotoCalib and another PhotoCalib.
     *
     * With:
     *   c = counts at position (x,y)
     *   this = this PhotoCalib
     *   other = other PhotoCalib
     *   return = BoundedField returned by this method
     * the return value from this method is defined as:
     *   this.countsToMaggies(c, (x,y)) * return(x, y) = other.countsToMaggies(c, (x,y))
     *
     * @param[in]  other  The PhotoCalib to scale to.
     *
     * @return     The BoundedField as defined above.
     */
    std::shared_ptr<afw::math::BoundedField> computeScalingTo(std::shared_ptr<PhotoCalib> other) const;

    /**
     * @brief      Compare two PhotoCalibs for equality.
     *
     * @param      rhs   The PhotoCalib to compare to.
     *
     * @return     True if both PhotoCalibs have identical representations.
     */
    bool operator==(PhotoCalib const& rhs) const;

    /**
     * @brief      Compare two PhotoCalibs for non-equality.
     *
     * @param      rhs   The PhotoCalib to compare to.
     *
     * @return     True if both PhotoCalibs do not have identical representations.
     */
    bool operator!=(PhotoCalib const& rhs) const { return !(*this == rhs); }

private:
    std::shared_ptr<afw::math::BoundedField> _zeroPoint;

    // The "mean" zero point, defined as the geometric mean of _zeroPoint evaluated over _zeroPoint's bbox.
    // Computed on instantiation as a convinience.
    // Also, the actual zeroPoint for a spatially-constant calibration.
    double _fluxMag0;

    double _fluxMag0Sigma;
};

}}}
