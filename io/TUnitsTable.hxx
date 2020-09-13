#ifndef TUnitsTable_h
#define TUnitsTable_h

#include <string>
#include <map>

namespace unit {
    class TUnitsTable;
};

class TVector3;
class TLorentzVector;

namespace unit {
    /// Take a value with a unit type and an uncertainty on the value and
    /// format it into a string.  This handles values that are free, or fixed
    /// based on the magnituded of the uncertainty (as defined by
    /// unit::IsFixed() and unit::IsFree()).  This
    /// knows about the following unit types: length, time, direction, angle,
    /// momentum, and charge.  If the uncertainty is negative, then only the
    /// value is printed. This provides a simplified interface to the
    /// unit::TUnitsTable value conversion routines.
    std::string AsString(double val, double sig, std::string type = "");

    /// Take a value with a unit type and format it into a string.  This knows
    /// about the following unit types: length, time, direction, angle,
    /// momentum, and charge.  This provides a simplified interface to the
    /// unit::TUnitsTable value conversion routines.
    std::string AsString(double val, std::string type);

    /// Take a TVector3 with a unit type and format it into a string.  This
    /// knows about the following unit types: length, time, direction, angle,
    /// momentum, and charge.  If the type is empty, then the vector is
    /// converted as a pure number.  This provides a simplified interface to
    /// the unit::TUnitsTable value conversion routines.
    std::string AsString(const TVector3& val, std::string type = "");

    /// Take a TLorentzVector with a unit type and format it into a string.
    /// This knows about the following unit types: length, time, direction,
    /// angle, momentum, and charge.  If the type is empty, then the vector is
    /// converted as a pure number.  This provides a simplified interface to
    /// the unit::TUnitsTable value conversion routines.
    std::string AsString(const TLorentzVector& val, std::string type = "");

    /// Check to see if a variance value means that the parameter is free.
    bool IsFree(double variance);

    /// Check to see if a variance value means that the parameter is fixed.
    bool IsFixed(double variance);

};

/// This class provides a method for converting a string like "1.5 cm" into a
/// double with the appropriate unit.  To do so it defines a set of units,
/// using the same base units as in HEPUnits.hxx (e.g. mm, ns, MeV).  Only a
/// fairly limited set of units is defined.  Care must be taken to keep this
/// in sync with HEPUnits.hxx.
class unit::TUnitsTable {
public:
    /// Constructor.  Creates list of units.
    ~TUnitsTable();

    ///  Get a reference to the singleton instance of runtime parameters.  If
    ///  this is first attempt at reference then singleton is instantiated and
    ///  parameters are read from text files.
    static unit::TUnitsTable& Get(void) {
        if (!fTUnitsTable)
            fTUnitsTable = new TUnitsTable();
        return *fTUnitsTable;
    }

    /// Converts a string like "1.5 cm" into a number represented by a string.
    /// For example, since the basic units are millimeters, "1.5 cm" is
    /// converted to "15"
    std::string ConvertWithUnit(std::string line);

    /// Convert a energy value (i.e. in MeV) into a string for printing.  The
    /// output will be formated as "xxx uu" where xxx is a number between -999
    /// and 999, and uu is a unit.
    std::string ConvertEnergy(double enr);

    /// Convert a length value (i.e. in mm) into a string for printing.  The
    /// output will be formated as "xxx uu" where xxx is a number between -999
    /// and 999, and uu is a unit.
    std::string ConvertLength(double len);

    /// Convert a charge value (i.e. in eplus) into a string for printing.  The
    /// output will be formated as "xxx uu" where xxx is a number between -999
    /// and 999, and uu is a unit.
    std::string ConvertCharge(double q);

    /// Convert a time value (i.e. in ns) into a string for printing.  The
    /// output will be formated as "xxx uu" where xxx is a number between -999
    /// and 999, and uu is a unit.
    std::string ConvertTime(double tim);

    /// Prints all the defined units.
    void PrintListOfUnits();

private:
    TUnitsTable();

    /// The table mapping known unit names to multipliers.
    std::map<std::string, double> units;

    /// The static pointer to the singleton instance.
    static TUnitsTable* fTUnitsTable;
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
