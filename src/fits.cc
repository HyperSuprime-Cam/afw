// -*- lsst-c++ -*-

#include <cstdio>
#include <complex>
#include <sstream>

#include "fitsio.h"
extern "C" {
#include "fitsio2.h"
}

#include "boost/regex.hpp"
#include "boost/filesystem.hpp"
#include "boost/preprocessor/seq/for_each.hpp"
#include "boost/cstdint.hpp"
#include "boost/format.hpp"
#include "boost/scoped_array.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/geom/Angle.h"

namespace lsst { namespace afw { namespace fits {

// ----------------------------------------------------------------------------------------------------------
// ---- Miscellaneous utilities -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

namespace {

// Strip leading and trailing single quotes and whitespace from a string.
std::string strip(std::string const & s) {
    if (s.empty()) return s;
    std::size_t i1 = s.find_first_not_of(" '");
    std::size_t i2 = s.find_last_not_of(" '");
    return s.substr(i1, (i1 == std::string::npos) ? 0 : 1 + i2 - i1);
}

// ---- FITS binary table format codes for various C++ types. -----------------------------------------------

char getFormatCode(bool*) { return 'X'; }
char getFormatCode(boost::uint8_t*) { return 'B'; }
char getFormatCode(boost::int16_t*) { return 'I'; }
char getFormatCode(boost::uint16_t*) { return 'U'; }
char getFormatCode(boost::int32_t*) { return 'J'; }
char getFormatCode(boost::uint32_t*) { return 'V'; }
char getFormatCode(boost::int64_t*) { return 'K'; }
char getFormatCode(float*) { return 'E'; }
char getFormatCode(double*) { return 'D'; }
char getFormatCode(std::complex<float>*) { return 'C'; }
char getFormatCode(std::complex<double>*) { return 'M'; }
char getFormatCode(lsst::afw::geom::Angle*) { return 'D'; }

// ---- Create a TFORM value for the given type and size ----------------------------------------------------

template <typename T>
std::string makeColumnFormat(int size = 1) {
    if (size > 0) {
        return (boost::format("%d%c") % size % getFormatCode((T*)0)).str();
    } else if (size < 0) {
        // variable length, max size given as -size
        return (boost::format("1P%c(%d)") % getFormatCode((T*)0) % (-size)).str();
    } else {
        // variable length, max size unknown
        return (boost::format("1P%c") % getFormatCode((T*)0)).str();
    }
}

// ---- Traits class to get cfitsio type constants from templates -------------------------------------------

template <typename T> struct FitsType;

template <> struct FitsType<bool> { static int const CONSTANT = TLOGICAL; };
template <> struct FitsType<unsigned char> { static int const CONSTANT = TBYTE; };
template <> struct FitsType<short> { static int const CONSTANT = TSHORT; };
template <> struct FitsType<unsigned short> { static int const CONSTANT = TUSHORT; };
template <> struct FitsType<int> { static int const CONSTANT = TINT; };
template <> struct FitsType<unsigned int> { static int const CONSTANT = TUINT; };
template <> struct FitsType<long> { static int const CONSTANT = TLONG; };
template <> struct FitsType<unsigned long> { static int const CONSTANT = TULONG; };
template <> struct FitsType<LONGLONG> { static int const CONSTANT = TLONGLONG; };
template <> struct FitsType<float> { static int const CONSTANT = TFLOAT; };
template <> struct FitsType<double> { static int const CONSTANT = TDOUBLE; };
template <> struct FitsType<lsst::afw::geom::Angle> { static int const CONSTANT = TDOUBLE; };
template <> struct FitsType< std::complex<float> > { static int const CONSTANT = TCOMPLEX; };
template <> struct FitsType< std::complex<double> > { static int const CONSTANT = TDBLCOMPLEX; };

// We use TBIT when writing booleans to table cells, but TLOGICAL in headers.
template <typename T> struct FitsTableType : public FitsType<T> {};
template <> struct FitsTableType<bool> { static int const CONSTANT = TBIT; };

template <typename T> struct FitsBitPix;

template <> struct FitsBitPix<unsigned char> { static int const CONSTANT = BYTE_IMG; };
template <> struct FitsBitPix<short> { static int const CONSTANT = SHORT_IMG; };
template <> struct FitsBitPix<unsigned short> { static int const CONSTANT = USHORT_IMG; };
template <> struct FitsBitPix<int> { static int const CONSTANT = LONG_IMG; }; // not a typo!
template <> struct FitsBitPix<unsigned int> { static int const CONSTANT = ULONG_IMG; };
template <> struct FitsBitPix<long> { static int const CONSTANT = LONG_IMG; };
template <> struct FitsBitPix<unsigned long> { static int const CONSTANT = ULONG_IMG; };
template <> struct FitsBitPix<LONGLONG> { static int const CONSTANT = LONGLONG_IMG; };
template <> struct FitsBitPix<float> { static int const CONSTANT = FLOAT_IMG; };
template <> struct FitsBitPix<double> { static int const CONSTANT = DOUBLE_IMG; };

} // anonymous

// ----------------------------------------------------------------------------------------------------------
// ---- Implementations for stuff in fits.h -----------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

std::string makeErrorMessage(std::string const & fileName, int status, std::string const & msg) {
    std::ostringstream os;
    os << "cfitsio error";
    if (fileName != "") {
        os << " (" << fileName << ")";
    }
    if (status != 0) {
        char fitsErrMsg[FLEN_ERRMSG];
        fits_get_errstatus(status, fitsErrMsg);
        os << ": " << fitsErrMsg << " (" << status << ")";
    }
    if (msg != "") {
        os << " : " << msg;
    }
    return os.str();
}

std::string makeErrorMessage(void * fptr, int status, std::string const & msg) {
    std::string fileName = "";
    fitsfile * fd = reinterpret_cast<fitsfile*>(fptr);
    if (fd != 0 && fd->Fptr != 0 && fd->Fptr->filename != 0) {
        fileName = fd->Fptr->filename;
    }
    return makeErrorMessage(fileName, status, msg);
}

void MemFileManager::reset() {
    if (_managed) std::free(_ptr);
    _ptr = 0;
    _len = 0;
    _managed = true;
}

void MemFileManager::reset(std::size_t len) {
    reset();
    _ptr = std::malloc(len);
    _len = len;
    _managed = true;
}

// ----------------------------------------------------------------------------------------------------------
// ---- Implementations for Fits class ----------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

int Fits::getHdu() {
    int n = 0;
    fits_get_hdu_num(reinterpret_cast<fitsfile*>(fptr), &n);
    return n;
}

void Fits::setHdu(int hdu) {
    if (hdu != 0) {
        fits_movabs_hdu(reinterpret_cast<fitsfile*>(fptr), hdu, 0, &status);
    }
    if (hdu == 0 && getHdu() == 1 && getImageDim() == 0) {
        // want a silent failure here
        int tmpStatus = status;
        fits_movrel_hdu(reinterpret_cast<fitsfile*>(fptr), 1, 0, &tmpStatus);
    }
    if (behavior & AUTO_CHECK) LSST_FITS_CHECK_STATUS(*this, boost::format("Moving to HDU %d") % hdu);
}

int Fits::countHdus() {
    int n = 0;
    fits_get_num_hdus(reinterpret_cast<fitsfile*>(fptr), &n, &status);
    if (behavior & AUTO_CHECK) LSST_FITS_CHECK_STATUS(*this, "Getting number of HDUs in file.");
    return n;
}

// ---- Writing and updating header keys --------------------------------------------------------------------

namespace {

// Impl functions in the anonymous namespace do special handling for strings, bools, and IEEE fp values.

template <typename T>
void updateKeyImpl(Fits & fits, char const * key, T const & value, char const * comment) {
    fits_update_key(
        reinterpret_cast<fitsfile*>(fits.fptr),
        FitsType<T>::CONSTANT,
        const_cast<char*>(key),
        const_cast<T *>(&value),
        const_cast<char*>(comment),
        &fits.status
    );    
}

void updateKeyImpl(Fits & fits, char const * key, std::string const & value, char const * comment) {
    fits_update_key_longstr(
        reinterpret_cast<fitsfile*>(fits.fptr),
        const_cast<char*>(key),
        const_cast<char*>(value.c_str()),
        const_cast<char*>(comment),
        &fits.status
    );    
}

void updateKeyImpl(Fits & fits, char const * key, bool const & value, char const * comment) {
    int v = value;
    fits_update_key(
        reinterpret_cast<fitsfile*>(fits.fptr),
        TLOGICAL,
        const_cast<char*>(key),
        &v,
        const_cast<char*>(comment),
        &fits.status
    );
}

template <typename T>
void writeKeyImpl(Fits & fits, char const * key, T const & value, char const * comment) {
    fits_write_key(
        reinterpret_cast<fitsfile*>(fits.fptr),
        FitsType<T>::CONSTANT,
        const_cast<char*>(key),
        const_cast<T *>(&value),
        const_cast<char*>(comment),
        &fits.status
    );    
}

void writeKeyImpl(Fits & fits, char const * key, std::string const & value, char const * comment) {
    if (strncmp(key, "COMMENT", 7) == 0) {
        fits_write_comment(
            reinterpret_cast<fitsfile*>(fits.fptr),
            const_cast<char*>(value.c_str()),
            &fits.status
        );
    } else if (strncmp(key, "HISTORY", 7) == 0) {
        fits_write_history(
            reinterpret_cast<fitsfile*>(fits.fptr),
            const_cast<char*>(value.c_str()),
            &fits.status
        );
    } else {
        fits_write_key_longstr(
            reinterpret_cast<fitsfile*>(fits.fptr),
            const_cast<char*>(key),
            const_cast<char*>(value.c_str()),
            const_cast<char*>(comment),
            &fits.status
        );
    }
}

void writeKeyImpl(Fits & fits, char const * key, bool const & value, char const * comment) {
    int v = value;
    fits_write_key(
        reinterpret_cast<fitsfile*>(fits.fptr),
        TLOGICAL,
        const_cast<char*>(key),
        &v,
        const_cast<char*>(comment),
        &fits.status
    );
}

} // anonymous

template <typename T>
void Fits::updateKey(std::string const & key, T const & value, std::string const & comment) {
    updateKeyImpl(*this, key.c_str(), value, comment.c_str());
    if (behavior & AUTO_CHECK) 
        LSST_FITS_CHECK_STATUS(*this, boost::format("Updating key '%s': '%s'") % key % value);
}

template <typename T>
void Fits::writeKey(std::string const & key, T const & value, std::string const & comment) {
    writeKeyImpl(*this, key.c_str(), value, comment.c_str());
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Writing key '%s': '%s'") % key % value);
}

template <typename T>
void Fits::updateKey(std::string const & key, T const & value) {
    updateKeyImpl(*this, key.c_str(), value, 0);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Updating key '%s': '%s'") % key % value);
}

template <typename T>
void Fits::writeKey(std::string const & key, T const & value) {
    writeKeyImpl(*this, key.c_str(), value, 0);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Writing key '%s': '%s'") % key % value);
}

template <typename T>
void Fits::updateColumnKey(std::string const & prefix, int n, T const & value, std::string const & comment) {
    updateKey((boost::format("%s%d") % prefix % (n + 1)).str(), value, comment);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Updating key '%s%d': '%s'") % prefix % (n+1) % value);
}

template <typename T>
void Fits::writeColumnKey(std::string const & prefix, int n, T const & value, std::string const & comment) {
    writeKey((boost::format("%s%d") % prefix % (n + 1)).str(), value, comment);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Writing key '%s%d': '%s'") % prefix % (n+1) % value);
}

template <typename T>
void Fits::updateColumnKey(std::string const & prefix, int n, T const & value) {
    updateKey((boost::format("%s%d") % prefix % (n + 1)).str(), value);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Updating key '%s%d': '%s'") % prefix % (n+1) % value);
}

template <typename T>
void Fits::writeColumnKey(std::string const & prefix, int n, T const & value) {
    writeKey((boost::format("%s%d") % prefix % (n + 1)).str(), value);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Writing key '%s%d': '%s'") % prefix % (n+1) % value);
}

// ---- Reading header keys ---------------------------------------------------------------------------------

namespace {

template <typename T>
void readKeyImpl(Fits & fits, char const * key, T & value) {
    fits_read_key(
        reinterpret_cast<fitsfile*>(fits.fptr), 
        FitsType<T>::CONSTANT,
        const_cast<char*>(key),
        &value,
        0,
        &fits.status
    );
}

void readKeyImpl(Fits & fits, char const * key, std::string & value) {
    char * buf = 0;
    fits_read_key_longstr(
        reinterpret_cast<fitsfile*>(fits.fptr),
        const_cast<char*>(key),
        &buf,
        0,
        &fits.status
    );
    if (buf) {
        value = strip(buf);
        free(buf);
    }
}

} // anonymous

template <typename T>
void Fits::readKey(std::string const & key, T & value) {
    readKeyImpl(*this, key.c_str(), value);
    if (behavior & AUTO_CHECK) LSST_FITS_CHECK_STATUS(*this, boost::format("Reading key '%s'") % key);
}

void Fits::forEachKey(HeaderIterationFunctor & functor) {
    char key[81];                       // allow for terminating NUL
    char value[81];
    char comment[81];
    int nKeys = 0;
    fits_get_hdrspace(reinterpret_cast<fitsfile*>(fptr), &nKeys, 0, &status);
    std::string keyStr;
    std::string valueStr;
    std::string commentStr;
    int i = 1;
    while (i <= nKeys) {
        fits_read_keyn(reinterpret_cast<fitsfile*>(fptr), i, key, value, comment, &status);
        keyStr = key;
        valueStr = value;
        commentStr = comment;
        ++i;
        while (valueStr.size() > 2 && valueStr[valueStr.size() - 2] == '&' && i <= nKeys) {
            // we're using key to hold the entire record here; the actual key is safe in keyStr
            fits_read_record(reinterpret_cast<fitsfile*>(fptr), i, key, &status);
            if (!strncmp(key, "CONTINUE", 8) == 0) {
                // require both trailing '&' and CONTINUE to invoke long-string handling
                break;
            }
            std::string card = key;
            valueStr.erase(valueStr.size() - 2);
            std::size_t firstQuote = card.find('\'');
            if (firstQuote == std::string::npos) {
                throw LSST_EXCEPT(
                    FitsError,
                    makeErrorMessage(
                        fptr, status,
                        boost::format("Invalid CONTINUE at header key %d: \"%s\".") % i % card
                    )
                );
            }
            std::size_t lastQuote = card.find('\'', firstQuote + 1);
            if (lastQuote == std::string::npos) {
                throw LSST_EXCEPT(
                    FitsError,
                    makeErrorMessage(
                        fptr, status,
                        boost::format("Invalid CONTINUE at header key %d: \"%s\".") % i % card
                    )
                );
            }
            valueStr += card.substr(firstQuote + 1, lastQuote - firstQuote);
            std::size_t slash = card.find('/', lastQuote + 1);
            if (slash != std::string::npos) {
                commentStr += strip(card.substr(slash + 1));
            }
            ++i;
        }
        if (behavior & AUTO_CHECK) LSST_FITS_CHECK_STATUS(*this, boost::format("Reading key '%s'") % keyStr);
        functor(keyStr, valueStr, commentStr);
    }
}

// ---- Reading and writing PropertySet/PropertyList --------------------------------------------------------

namespace {

bool isKeyIgnored(std::string const & key) {
    return key == "SIMPLE" || key == "BITPIX" || key == "EXTEND" ||
        key == "GCOUNT" || key == "PCOUNT" || key == "XTENSION" || key == "TFIELDS" ||
        key == "BSCALE" || key == "BZERO" || key.compare(0, 5, "NAXIS") == 0;
}

class MetadataIterationFunctor : public HeaderIterationFunctor {
public:

    virtual void operator()(
        std::string const & key,
        std::string const & value,
        std::string const & comment
    );

    template <typename T>
    void add(std::string const & key, T value, std::string const & comment) {
        if (list) {
            list->add(key, value, comment);
        } else {
            set->add(key, value);
        }
    }

    bool strip;
    daf::base::PropertySet * set;
    daf::base::PropertyList * list;
};

void MetadataIterationFunctor::operator()(
    std::string const & key,
    std::string const & value,
    std::string const & comment
) {
    static boost::regex const boolRegex("[tTfF]");
    static boost::regex const intRegex("[+-]?[0-9]+");
    static boost::regex const doubleRegex("[+-]?([0-9]*\\.[0-9]+|[0-9]+\\.[0-9]*)([eE][+-]?[0-9]+)?");
    static boost::regex const fitsStringRegex("'(.*?) *'");
    boost::smatch matchStrings;

    if (strip && isKeyIgnored(key)) {
        return;
    }

    std::istringstream converter(value);
    if (boost::regex_match(value, boolRegex)) {
        // convert the string to an bool
        add(key, bool(value == "T" || value == "t"), comment);
    } else if (boost::regex_match(value, intRegex)) {
        // convert the string to an int
        boost::int64_t val;
        converter >> val;
        if (val < (1LL << 31) && val > -(1LL << 31)) {
            add(key, static_cast<int>(val), comment);
        } else {
            add(key, val, comment);
        }
    } else if (boost::regex_match(value, doubleRegex)) {
        // convert the string to a double
        double val;
        converter >> val;
        add(key, val, comment);
    } else if (boost::regex_match(value, matchStrings, fitsStringRegex)) {
        // strip off the enclosing single quotes and return the string
        add(key, matchStrings[1].str(), comment);
    } else if (value.empty()) {
        // do nothing for empty values
    } else if (key == "HISTORY" ||
               (key == "COMMENT" &&
                comment != "  FITS (Flexible Image Transport System) format is defined in 'Astronomy" &&
                comment != "  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")) {
        if (list) {
            list->add(key, comment);
        } else {
            set->add(key, comment);
        }
    } else {
        throw LSST_EXCEPT(
            afw::fits::FitsError,
            (boost::format("Could not parse header value for key '%s': '%s'") % key % value).str()
        );
    }
}

void writeKeyFromProperty(
    Fits & fits, daf::base::PropertySet const & metadata, std::string const & key, char const * comment=0
) {
    std::type_info const & valueType = metadata.typeOf(key); 
    if (valueType == typeid(bool)) {
        if (metadata.isArray(key)) {
            std::vector<bool> tmp = metadata.getArray<bool>(key);
            // work around unfortunate specialness of std::vector<bool>
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), static_cast<bool>(tmp[i]), comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<bool>(key), comment);
        }
    } else if (valueType == typeid(int)) {
        if (metadata.isArray(key)) {
            std::vector<int> tmp = metadata.getArray<int>(key);
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), tmp[i], comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<int>(key), comment);
        }
    } else if (valueType == typeid(long)) {
        if (metadata.isArray(key)) {
            std::vector<long> tmp = metadata.getArray<long>(key);
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), tmp[i], comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<long>(key), comment);
        }
    } else if (valueType == typeid(boost::int64_t)) {
        if (metadata.isArray(key)) {
            std::vector<boost::int64_t> tmp = metadata.getArray<boost::int64_t>(key);
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), tmp[i], comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<boost::int64_t>(key), comment);
        }
    } else if (valueType == typeid(double)) {
        if (metadata.isArray(key)) {
            std::vector<double> tmp = metadata.getArray<double>(key);
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), tmp[i], comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<double>(key), comment);
        }
    } else if (valueType == typeid(std::string)) {
        if (metadata.isArray(key)) {
            std::vector<std::string> tmp = metadata.getArray<std::string>(key);
            for (std::size_t i = 0; i != tmp.size(); ++i) {
                writeKeyImpl(fits, key.c_str(), tmp[i], comment);
            }
        } else {
            writeKeyImpl(fits, key.c_str(), metadata.get<std::string>(key), comment);
        }
    } else {
        // FIXME: inherited this error handling from fitsIo.cc; need a better option.
        pex::logging::Log::getDefaultLog().log(
            pex::logging::Log::WARN,
            makeErrorMessage(
                fits.fptr, fits.status,
                boost::format("In %s, unknown type '%s' for key '%s'.")
                % BOOST_CURRENT_FUNCTION % valueType.name() % key
            )
        );
    }
    if (fits.behavior & Fits::AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(fits, boost::format("Writing key '%s'") % key);
}

} // anonymous

void Fits::readMetadata(daf::base::PropertySet & metadata, bool strip) {
    MetadataIterationFunctor f;
    f.strip = strip;
    f.set = &metadata;
    f.list = dynamic_cast<daf::base::PropertyList*>(&metadata);
    forEachKey(f);
}

void Fits::writeMetadata(daf::base::PropertySet const & metadata) {
    typedef std::vector<std::string> NameList;
    daf::base::PropertyList const * pl = dynamic_cast<daf::base::PropertyList const*>(&metadata);
    NameList paramNames;
    if (pl) {
        paramNames = pl->getOrderedNames();
    } else {
        paramNames = metadata.paramNames(false);
    }
    for (NameList::const_iterator i = paramNames.begin(); i != paramNames.end(); ++i) {
        if (!isKeyIgnored(*i)) {
            if (pl) {
                writeKeyFromProperty(*this, metadata, *i, pl->getComment(*i).c_str());
            } else {
                writeKeyFromProperty(*this, metadata, *i);
            }
        }
    }
}

// ---- Manipulating tables ---------------------------------------------------------------------------------

void Fits::createTable() {
    char * ttype = 0;
    char * tform = 0;
    fits_create_tbl(reinterpret_cast<fitsfile*>(fptr), BINARY_TBL, 0, 0, &ttype, &tform, 0, 0, &status);
    if (behavior & AUTO_CHECK) LSST_FITS_CHECK_STATUS(*this, "Creating binary table");
}

template <typename T>
int Fits::addColumn(std::string const & ttype, int size) {
    int nCols = 0;
    fits_get_num_cols(
        reinterpret_cast<fitsfile*>(fptr),
        &nCols,
        &status
    );
    std::string tform = makeColumnFormat<T>(size);
    fits_insert_col(
        reinterpret_cast<fitsfile*>(fptr),
        nCols + 1,
        const_cast<char*>(ttype.c_str()),
        const_cast<char*>(tform.c_str()),
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Adding column '%s' with size %d") % ttype % size);
    return nCols;
}

template <typename T>
int Fits::addColumn(std::string const & ttype, int size, std::string const & comment) {
    int nCols = addColumn<T>(ttype, size);
    updateColumnKey("TTYPE", nCols, ttype, comment);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Adding column '%s' with size %d") % ttype % size);
    return nCols;
}

std::size_t Fits::addRows(std::size_t nRows) {
    long first = 0;
    fits_get_num_rows(
        reinterpret_cast<fitsfile*>(fptr),
        &first,
        &status
    );
    fits_insert_rows(
        reinterpret_cast<fitsfile*>(fptr),
        first,
        nRows,
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Adding %d rows to binary table") % nRows);
    return first;
}

std::size_t Fits::countRows() {
    long r = 0;
    fits_get_num_rows(
        reinterpret_cast<fitsfile*>(fptr),
        &r,
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, "Checking how many rows are in table");
    return r;
}

template <typename T>
void Fits::writeTableArray(std::size_t row, int col, int nElements, T const * value) {
    fits_write_col(
        reinterpret_cast<fitsfile*>(fptr), 
        FitsTableType<T>::CONSTANT, 
        col + 1, row + 1, 
        1, nElements,
        const_cast<T*>(value),
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Writing value at table cell (%d, %d)") % row % col);
}

template <typename T>
void Fits::readTableArray(std::size_t row, int col, int nElements, T * value) {
    int anynul = false;
    fits_read_col(
        reinterpret_cast<fitsfile*>(fptr), 
        FitsTableType<T>::CONSTANT, 
        col + 1, row + 1, 
        1, nElements,
        0,
        value,
        &anynul,
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Reading value at table cell (%d, %d)") % row % col);
}

long Fits::getTableArraySize(int col) {
    int typecode = 0;
    long result = 0;
    long width = 0;
    fits_get_coltype(
        reinterpret_cast<fitsfile*>(fptr),
        col + 1,
        &typecode,
        &result,
        &width,
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Looking up array size for column %d") % col);
    return result;
}

long Fits::getTableArraySize(std::size_t row, int col) {
    long result = 0;
    long offset = 0;
    fits_read_descript(
        reinterpret_cast<fitsfile*>(fptr),
        col + 1,
        row + 1,
        &result,
        &offset,
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Looking up array size for cell (%d, %d)") % row % col);
    return result;
}

// ---- Manipulating images ---------------------------------------------------------------------------------

void Fits::createEmpty() {
    long naxes = 0;
    fits_create_img(reinterpret_cast<fitsfile*>(fptr), 8, 0, &naxes, &status);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, "Creating empty image HDU");    
}

template <typename T>
void Fits::createImageImpl(int naxis, long * naxes) {
    fits_create_img(reinterpret_cast<fitsfile*>(fptr), FitsBitPix<T>::CONSTANT, naxis, naxes, &status);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, "Creating new image HDU");
}

template <typename T>
void Fits::writeImageImpl(T const * data, int nElements) {
    fits_write_img(
        reinterpret_cast<fitsfile*>(fptr),
        FitsType<T>::CONSTANT,
        1, nElements,
        const_cast<T*>(data),
        &status
    );
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, "Writing image");
}

int Fits::getImageDim() {
    int nAxis = 0;
    fits_get_img_dim(reinterpret_cast<fitsfile*>(fptr), &nAxis, &status);
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, "Getting NAXIS");
    return nAxis;
}

// ---- Manipulating files ----------------------------------------------------------------------------------

Fits::Fits(std::string const & filename, std::string const & mode, int behavior_)
    : fptr(0), status(0), behavior(behavior_)
{
    if (mode == "r" || mode == "rb") {
        fits_open_file(
            reinterpret_cast<fitsfile**>(&fptr),
            const_cast<char*>(filename.c_str()), 
            READONLY,
            &status
        );
    } else if (mode == "w" || mode == "wb") {
        boost::filesystem::remove(filename); // cfitsio doesn't like over-writing files
        fits_create_file(
            reinterpret_cast<fitsfile**>(&fptr), 
            const_cast<char*>(filename.c_str()),
            &status
        );
    } else if (mode == "a" || mode == "ab") {
        fits_open_file(
            reinterpret_cast<fitsfile**>(&fptr),
            const_cast<char*>(filename.c_str()), 
            READWRITE,
            &status
        );
        int nHdu = 0;
        fits_get_num_hdus(reinterpret_cast<fitsfile*>(fptr), &nHdu, &status);
        fits_movabs_hdu(reinterpret_cast<fitsfile*>(fptr), nHdu, NULL, &status);
        if ((behavior & AUTO_CHECK) && (behavior & AUTO_CLOSE) && (status) && (fptr)) {
            // We're about to throw an exception, and the destructor won't get called
            // because we're in the constructor, so cleanup here first.
            int tmpStatus = 0;
            fits_close_file(reinterpret_cast<fitsfile*>(fptr), &tmpStatus);
        }
    } else {
        throw LSST_EXCEPT(
            FitsError,
            (boost::format("Invalid mode '%s' given when opening file '%s'") % mode % filename).str()
        );
    }
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(*this, boost::format("Opening file '%s' with mode '%s'") % filename % mode);
}

Fits::Fits(MemFileManager & manager, std::string const & mode, int behavior_)
    : fptr(0), status(0), behavior(behavior_)
{
    typedef void * (*Reallocator)(void *, std::size_t);
    // It's a shame this logic is essentially a duplicate of above, but the innards are different enough
    // we can't really reuse it.
    if (mode == "r" || mode == "rb") {
        fits_open_memfile(
            reinterpret_cast<fitsfile**>(&fptr),
            "unused",
            READONLY,
            &manager._ptr, &manager._len,
            0, 0, // no reallocator or deltasize necessary for READONLY
            &status
        );
    } else if (mode == "w" || mode == "wb") {
        if (!manager._ptr) manager.reset(2880);
        Reallocator reallocator = 0;
        if (manager._managed) reallocator = &std::realloc;
        fits_create_memfile(
            reinterpret_cast<fitsfile**>(&fptr),
            &manager._ptr, &manager._len,
            0, reallocator, // use default deltasize
            &status
        );
    } else if (mode == "a" || mode == "ab") {
        Reallocator reallocator = 0;
        if (manager._managed) reallocator = &std::realloc;
        fits_open_memfile(
            reinterpret_cast<fitsfile**>(&fptr),
            "unused",
            READWRITE,
            &manager._ptr, &manager._len,
            0, reallocator,
            &status
        );
        int nHdu = 0;
        fits_get_num_hdus(reinterpret_cast<fitsfile*>(fptr), &nHdu, &status);
        fits_movabs_hdu(reinterpret_cast<fitsfile*>(fptr), nHdu, NULL, &status);
        if ((behavior & AUTO_CHECK) && (behavior & AUTO_CLOSE) && (status) && (fptr)) {
            // We're about to throw an exception, and the destructor won't get called
            // because we're in the constructor, so cleanup here first.
            int tmpStatus = 0;
            fits_close_file(reinterpret_cast<fitsfile*>(fptr), &tmpStatus);
        }
    } else {
        throw LSST_EXCEPT(
            FitsError,
            (boost::format("Invalid mode '%s' given when opening memory file at '%s'")
             % mode % manager._ptr).str()
        );
    }
    if (behavior & AUTO_CHECK)
        LSST_FITS_CHECK_STATUS(
            *this, boost::format("Opening memory file at '%s' with mode '%s'") % manager._ptr % mode
        );
}

void Fits::closeFile() {
    fits_close_file(reinterpret_cast<fitsfile*>(fptr), &status);
}

#define INSTANTIATE_KEY_OPS(r, data, T)                                \
    template void Fits::updateKey(std::string const &, T const &, std::string const &); \
    template void Fits::writeKey(std::string const &, T const &, std::string const &); \
    template void Fits::updateKey(std::string const &, T const &);      \
    template void Fits::writeKey(std::string const &, T const &);       \
    template void Fits::updateColumnKey(std::string const &, int, T const &, std::string const &); \
    template void Fits::writeColumnKey(std::string const &, int, T const &, std::string const &); \
    template void Fits::updateColumnKey(std::string const &, int, T const &); \
    template void Fits::writeColumnKey(std::string const &, int, T const &); \
    template void Fits::readKey(std::string const &, T &);

#define INSTANTIATE_IMAGE_OPS(r, data, T)                        \
    template void Fits::createImageImpl<T>(int, long *);         \
    template void Fits::writeImageImpl(T const *, int);

#define INSTANTIATE_TABLE_OPS(r, data, T)                               \
    template int Fits::addColumn<T>(std::string const & ttype, int size); \
    template int Fits::addColumn<T>(std::string const & ttype, int size, std::string const & comment); \
    template void Fits::writeTableArray(std::size_t row, int col, int nElements, T const * value); \
    template void Fits::readTableArray(std::size_t row, int col, int nElements, T * value);

// ----------------------------------------------------------------------------------------------------------
// ---- Explicit instantiation ------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

#define KEY_TYPES                                                       \
    (bool)(unsigned char)(short)(unsigned short)(int)(unsigned int)(long)(unsigned long)(LONGLONG) \
    (float)(double)(std::complex<float>)(std::complex<double>)(std::string)

#define COLUMN_TYPES                            \
    (bool)(boost::uint8_t)(boost::int16_t)(boost::uint16_t)(boost::int32_t)(boost::uint32_t) \
    (boost::int64_t)(float)(double)(lsst::afw::geom::Angle)(std::complex<float>)(std::complex<double>)

#define IMAGE_TYPES                                                       \
    (unsigned char)(short)(unsigned short)(int)(unsigned int)(long)(unsigned long)(LONGLONG) \
    (float)(double)

BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_KEY_OPS, _, KEY_TYPES)
BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_TABLE_OPS, _, COLUMN_TYPES)
BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_IMAGE_OPS, _, IMAGE_TYPES)

}}} // namespace lsst::afw::fits
