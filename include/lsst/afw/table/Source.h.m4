changecom(`###')dnl
// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014, 2011 LSST Corporation.
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

// THIS FILE IS AUTOMATICALLY GENERATED from Source.h.m4, AND WILL BE OVERWRITTEN IF EDITED MANUALLY.

define(`m4def', defn(`define'))dnl
m4def(`DECLARE_SLOT_GETTERS',
`/// @brief Get the value of the $1$2 slot measurement.
    $2SlotDefinition::MeasValue get$1$2() const;

    /// @brief Get the uncertainty on the $1$2 slot measurement.
    $2SlotDefinition::ErrValue get$1$2Err() const;

    /// @brief Return true if the measurement in the $1$2 slot failed.
    bool get$1$2Flag() const;
')dnl
m4def(`DEFINE_SLOT_GETTERS',
`inline $2SlotDefinition::MeasValue SourceRecord::get$1$2() const {
    return this->get(getTable()->get$1$2Slot().getMeasKey());
}

inline $2SlotDefinition::ErrValue SourceRecord::get$1$2Err() const {
    return this->get(getTable()->get$1$2Slot().getErrKey());
}

inline bool SourceRecord::get$1$2Flag() const {
    return this->get(getTable()->get$1$2Slot().getFlagKey());
}
')dnl
m4def(`DECLARE_SLOT_DEFINERS',
`
    $2SlotDefinition const & get$1$2Slot() const { return _slots.def$1$2; }

    /**
     *  @brief Set the measurement used for the $1$2 slot.
     *
     *  The definitions for slots are actually managed by the Schema object, and its associated
     *  AliasMap, so this simply sets the "slot_$1$2" alias
     *  to point to the given field name prefix.  See $2SlotDefinition for more information.
     */
    void define$1$2(std::string const & name) {
        getSchema().getAliasMap()->set(get$1$2Slot().getAlias(), name);
    }

    /**
     *  @brief Return the name of the field used for the $1$2 slot.
     *
     *  @throw NotFoundError if the slot is not defined.
     *
     *  @deprecated in favor of
     *  @code
     *  getSchema().getAliasMap()->get("slot_$1$2")
     *  @endcode
     */
    std::string get$1$2Definition() const {
        return getSchema().getAliasMap()->get(get$1$2Slot().getAlias());
    }

    /**
     *  @brief Return true if the $1$2 slot corresponds to a valid field.
     *
     *  @deprecated in favor of get$1$2Slot().isValid().
     */
    bool has$1$2Slot() const {
        return get$1$2Slot().isValid();
    }

    /**
     *  @brief Return the key used for the $1$2 slot measurement value.
     *
     *  @deprecated in favor of get$1$2Slot().getMeasKey().
     */
    $2SlotDefinition::MeasKey get$1$2Key() const {
        return get$1$2Slot().getMeasKey();
    }

    /**
     *  @brief Return the key used for the $1$2 slot uncertainty.
     *
     *  @deprecated in favor of get$1$2Slot().getErrKey().
     */
    $2SlotDefinition::ErrKey get$1$2ErrKey() const {
        return get$1$2Slot().getErrKey();
    }

    /**
     *  @brief Return the key used for the $1$2 slot failure flag.
     *
     *  @deprecated in favor of get$1$2Slot().getFlagKey().
     */
    Key<Flag> get$1$2FlagKey() const {
        return get$1$2Slot().getFlagKey();
    }
')dnl
define(`m4def', defn(`define'))dnl
m4def(`DEFINE_FLUX_COLUMN_GETTERS',
`/// @brief Get the value of the $1Flux slot measurement.
    ndarray::Array<double,1> get$1Flux() const {
        return this->operator[](this->getTable()->get$1FluxSlot().getMeasKey());
    }
    /// @brief Get the uncertainty on the $1Flux slot measurement.
    ndarray::Array<double,1> get$1FluxErr() const {
        return this->operator[](this->getTable()->get$1FluxSlot().getErrKey());
    }
')dnl
#ifndef AFW_TABLE_Source_h_INCLUDED
#define AFW_TABLE_Source_h_INCLUDED

#include "boost/array.hpp"
#include "boost/type_traits/is_convertible.hpp"

#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/table/Simple.h"
#include "lsst/afw/table/aggregates.h"
#include "lsst/afw/table/IdFactory.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/slots.h"
#include "lsst/afw/table/io/FitsWriter.h"

namespace lsst { namespace afw {

namespace image {
class Wcs;
} // namespace image

namespace table {

/**
 *  @brief Bitflags to be passed to SourceCatalog::readFits and SourceCatalog::writeFits
 *
 *  Note that these flags may also be passed when reading/writing SourceCatalogs via the Butler,
 *  by passing a "flags" key/value pair as part of the data ID.
 */
enum SourceFitsFlags {
    SOURCE_IO_NO_FOOTPRINTS = 0x1,       ///< Do not read/write footprints at all
    SOURCE_IO_NO_HEAVY_FOOTPRINTS = 0x2  ///< Read/write heavy footprints as non-heavy footprints
};

typedef lsst::afw::detection::Footprint Footprint;

class SourceRecord;
class SourceTable;

template <typename RecordT> class SourceColumnViewT;

/**
 *  @brief Record class that contains measurements made on a single exposure.
 *
 *  Sources provide four additions to SimpleRecord / SimpleRecord:
 *   - Specific fields that must always be present, with specialized getters.
 *     The schema for a SourceTable should always be constructed by starting with the result of
 *     SourceTable::makeMinimalSchema.
 *   - A shared_ptr to a Footprint for each record.
 *   - A system of aliases (called slots) in which a SourceTable instance stores keys for particular
 *     measurements (a centroid, a shape, and a number of different fluxes) and SourceRecord uses
 *     this keys to provide custom getters and setters.  These are not separate fields, but rather
 *     aliases that can point to custom fields.  See the SlotDefinition hierarchy for more information.
 */
class SourceRecord : public SimpleRecord {
public:

    typedef SourceTable Table;
    typedef SourceColumnViewT<SourceRecord> ColumnView;
    typedef SortedCatalogT<SourceRecord> Catalog;
    typedef SortedCatalogT<SourceRecord const> ConstCatalog;

    PTR(Footprint) getFootprint() const { return _footprint; }

    void setFootprint(PTR(Footprint) const & footprint) { _footprint = footprint; }

    CONST_PTR(SourceTable) getTable() const {
        return std::static_pointer_cast<SourceTable const>(BaseRecord::getTable());
    }

    //@{
    /// @brief Convenience accessors for the keys in the minimal source schema.
    RecordId getParent() const;
    void setParent(RecordId id);
    //@}

    DECLARE_SLOT_GETTERS(`Psf', `Flux')
    DECLARE_SLOT_GETTERS(`Model', `Flux')
    DECLARE_SLOT_GETTERS(`Ap', `Flux')
    DECLARE_SLOT_GETTERS(`Inst', `Flux')
    DECLARE_SLOT_GETTERS(`Calib', `Flux')
    DECLARE_SLOT_GETTERS(`', `Centroid')
    DECLARE_SLOT_GETTERS(`', `Shape')

    /// @brief Return the centroid slot x coordinate.
    double getX() const;

    /// @brief Return the centroid slot y coordinate.
    double getY() const;

    /// @brief Return the shape slot Ixx value.
    double getIxx() const;

    /// @brief Return the shape slot Iyy value.
    double getIyy() const;

    /// @brief Return the shape slot Ixy value.
    double getIxy() const;

    /// @brief Update the coord field using the given Wcs and the field in the centroid slot.
    void updateCoord(image::Wcs const & wcs);

    /// @brief Update the coord field using the given Wcs and the image center from the given key.
    void updateCoord(image::Wcs const & wcs, PointKey<double> const & key);

protected:

    SourceRecord(PTR(SourceTable) const & table);

    virtual void _assign(BaseRecord const & other);

private:
    PTR(Footprint) _footprint;
};

/**
 *  @brief Table class that contains measurements made on a single exposure.
 *
 *  @copydetails SourceRecord
 */
class SourceTable : public SimpleTable {
public:

    typedef SourceRecord Record;
    typedef SourceColumnViewT<SourceRecord> ColumnView;
    typedef SortedCatalogT<Record> Catalog;
    typedef SortedCatalogT<Record const> ConstCatalog;

    /**
     *  @brief Construct a new table.
     *
     *  @param[in] schema            Schema that defines the fields, offsets, and record size for the table.
     *  @param[in] idFactory         Factory class to generate record IDs when they are not explicitly given.
     *                               If null, record IDs will default to zero.
     *
     *  Note that not passing an IdFactory at all will call the other override of make(), which will
     *  set the ID factory to IdFactory::makeSimple().
     */
    static PTR(SourceTable) make(Schema const & schema, PTR(IdFactory) const & idFactory);

    /**
     *  @brief Construct a new table.
     *
     *  @param[in] schema            Schema that defines the fields, offsets, and record size for the table.
     *
     *  This overload sets the ID factory to IdFactory::makeSimple().
     */
    static PTR(SourceTable) make(Schema const & schema) { return make(schema, IdFactory::makeSimple()); }

    /**
     *  @brief Return a minimal schema for Source tables and records.
     *
     *  The returned schema can and generally should be modified further,
     *  but many operations on sources will assume that at least the fields
     *  provided by this routine are present.
     *
     *  Keys for the standard fields added by this routine can be obtained
     *  from other static member functions of the SourceTable class.
     */
    static Schema makeMinimalSchema() {
        Schema r = getMinimalSchema().schema;
        r.disconnectAliases();
        return r;
    }

    /**
     *  @brief Return true if the given schema is a valid SourceTable schema.
     *
     *  This will always be true if the given schema was originally constructed
     *  using makeMinimalSchema(), and will rarely be true otherwise.
     */
    static bool checkSchema(Schema const & other) {
        return other.contains(getMinimalSchema().schema);
    }

    /// @brief Key for the parent ID.
    static Key<RecordId> getParentKey() { return getMinimalSchema().parent; }

    /// @copydoc BaseTable::clone
    PTR(SourceTable) clone() const { return std::static_pointer_cast<SourceTable>(_clone()); }

    /// @copydoc BaseTable::makeRecord
    PTR(SourceRecord) makeRecord() { return std::static_pointer_cast<SourceRecord>(_makeRecord()); }

    /// @copydoc BaseTable::copyRecord
    PTR(SourceRecord) copyRecord(BaseRecord const & other) {
        return std::static_pointer_cast<SourceRecord>(BaseTable::copyRecord(other));
    }

    /// @copydoc BaseTable::copyRecord
    PTR(SourceRecord) copyRecord(BaseRecord const & other, SchemaMapper const & mapper) {
        return std::static_pointer_cast<SourceRecord>(BaseTable::copyRecord(other, mapper));
    }

    DECLARE_SLOT_DEFINERS(`Psf', `Flux')
    DECLARE_SLOT_DEFINERS(`Model', `Flux')
    DECLARE_SLOT_DEFINERS(`Ap', `Flux')
    DECLARE_SLOT_DEFINERS(`Inst', `Flux')
    DECLARE_SLOT_DEFINERS(`Calib', `Flux')
    DECLARE_SLOT_DEFINERS(`', `Centroid')
    DECLARE_SLOT_DEFINERS(`', `Shape')

protected:

    SourceTable(Schema const & schema, PTR(IdFactory) const & idFactory);

    SourceTable(SourceTable const & other);

    virtual void handleAliasChange(std::string const & alias);

private:

    // Struct that holds the minimal schema and the special keys we've added to it.
    struct MinimalSchema {
        Schema schema;
        Key<RecordId> parent;

        MinimalSchema();
    };

    // Return the singleton minimal schema.
    static MinimalSchema & getMinimalSchema();

    friend class io::FitsWriter;
    friend class SourceRecord;

     // Return a writer object that knows how to save in FITS format.  See also FitsWriter.
    virtual PTR(io::FitsWriter) makeFitsWriter(fits::Fits * fitsfile, int flags) const;

    SlotSuite _slots;
};

template <typename RecordT>
class SourceColumnViewT : public ColumnViewT<RecordT> {
public:

    typedef RecordT Record;
    typedef typename RecordT::Table Table;

    // See the documentation for BaseColumnView for an explanation of why these
    // accessors *appear* to violate const-correctness.

    DEFINE_FLUX_COLUMN_GETTERS(`Psf')
    DEFINE_FLUX_COLUMN_GETTERS(`Ap')
    DEFINE_FLUX_COLUMN_GETTERS(`Model')
    DEFINE_FLUX_COLUMN_GETTERS(`Inst')
    DEFINE_FLUX_COLUMN_GETTERS(`Calib')

    ndarray::Array<double,1> const getX() const {
        return this->operator[](this->getTable()->getCentroidKey().getX());
    }
    ndarray::Array<double,1> const getY() const {
        return this->operator[](this->getTable()->getCentroidKey().getY());
    }

    ndarray::Array<double,1> const getIxx() const {
        return this->operator[](this->getTable()->getShapeKey().getIxx());
    }
    ndarray::Array<double,1> const getIyy() const {
        return this->operator[](this->getTable()->getShapeKey().getIyy());
    }
    ndarray::Array<double,1> const getIxy() const {
        return this->operator[](this->getTable()->getShapeKey().getIxy());
    }

    /// @brief @copydoc BaseColumnView::make
    template <typename InputIterator>
    static SourceColumnViewT make(PTR(Table) const & table, InputIterator first, InputIterator last) {
        return SourceColumnViewT(BaseColumnView::make(table, first, last));
    }

protected:
    explicit SourceColumnViewT(BaseColumnView const & base) : ColumnViewT<RecordT>(base) {}
};

typedef SourceColumnViewT<SourceRecord> SourceColumnView;

#ifndef SWIG

DEFINE_SLOT_GETTERS(`Psf', `Flux')
DEFINE_SLOT_GETTERS(`Model', `Flux')
DEFINE_SLOT_GETTERS(`Ap', `Flux')
DEFINE_SLOT_GETTERS(`Inst', `Flux')
DEFINE_SLOT_GETTERS(`Calib', `Flux')
DEFINE_SLOT_GETTERS(`', `Centroid')
DEFINE_SLOT_GETTERS(`', `Shape')

inline RecordId SourceRecord::getParent() const { return get(SourceTable::getParentKey()); }
inline void SourceRecord::setParent(RecordId id) { set(SourceTable::getParentKey(), id); }
inline double SourceRecord::getX() const {
    return get(getTable()->getCentroidKey().getX());
}
inline double SourceRecord::getY() const {
    return get(getTable()->getCentroidKey().getY());
}
inline double SourceRecord::getIxx() const {
    return get(getTable()->getShapeKey().getIxx());
}
inline double SourceRecord::getIyy() const {
    return get(getTable()->getShapeKey().getIyy());
}
inline double SourceRecord::getIxy() const {
    return get(getTable()->getShapeKey().getIxy());
}

#endif // !SWIG

}}} // namespace lsst::afw::table

#endif // !AFW_TABLE_Source_h_INCLUDED
