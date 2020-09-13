#include "TUnitsTable.hxx"

#include "CubeUnits.hxx"
#include "CubeLog.hxx"

#include <TVector3.h>
#include <TLorentzVector.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

unit::TUnitsTable*
unit::TUnitsTable::fTUnitsTable = NULL;

unit::TUnitsTable::TUnitsTable() {
    // Define units.  Add to this list as needed.

    // Length
    units["km"] = unit::kilometer;
    units["m"] = unit::meter;
    units["cm"] = unit::centimeter ;
    units["mm"] = unit::millimeter ;
    units["nm"] = unit::nanometer ;

    // Area
    units["m2"] = unit::meter*unit::meter ;
    units["cm2"] = unit::centimeter*unit::centimeter ;
    units["mm2"] = unit::millimeter*unit::millimeter ;
    units["barn"] = unit::barn ;
    units["mbarn"] = unit::millibarn ;
    units["mubarn"] = unit::microbarn  ;
    units["nbarn"] = unit::nanobarn ;
    units["pbarn"] = unit::picobarn ;

    // Volume
    units["m3"] = unit::meter*unit::meter*unit::meter ;
    units["cm3"] = unit::cm3  ;
    units["mm3"] = unit::mm3  ;

    // Degree
    units["rad"] = unit::radian ;
    units["mrad"] = unit::milliradian  ;
    units["sr"] = unit::steradian ;
    units["deg"] = unit::degree ;

    // Time
    units["s"] = unit::second ;
    units["ms"] = unit::millisecond ;
    units["us"] = unit::microsecond ;
    units["ns"] = unit::nanosecond ;

    // Frequency
    units["Hz"] = unit::hertz ;
    units["kHz"] = unit::kilohertz ;
    units["MHz"] = unit::megahertz  ;

    // Electric Charge
    units["electrons"] = unit::eplus  ;
    units["e+"] = unit::eplus  ;
    units["C"] = unit::coulomb ;
    units["muC"] = unit::microcoulomb ;
    units["fC"] = unit::femtocoulomb ;

    // Energy
    units["eV"] = unit::eV ;
    units["keV"] = unit::keV  ;
    units["MeV"] = unit::MeV ;
    units["GeV"] = unit::GeV  ;
    units["TeV"] = unit::TeV  ;
    units["J"] = unit::joule ;

    // Energy/Length
    units["GeV/cm"] = units["GeV"]/units["cm"]  ;
    units["MeV/cm"] = units["MeV"]/units["cm"]  ;
    units["keV/cm"] = units["keV"]/units["cm"]  ;
    units["eV/cm"] = units["eV"]/units["cm"] ;

    // Mass
    units["mg"] = unit::milligram  ;
    units["g"] = unit::gram ;
    units["kg"] = unit::kilogram  ;

    // Volumic Mass
    units["g/cm3"] = units["g"]/units["cm3"] ;
    units["mg/cm3"] = units["mg"]/units["cm3"]  ;
    units["kg/m3"] = units["kg"]/units["m3"] ;

    units["g/cm2"] = units["g"]/units["cm2"] ;
    units["mg/cm2"] = units["mg"]/units["cm2"] ;
    units["kg/m2"] = units["kg"]/units["m2"]  ;

    // Power
    units["W"] = unit::watt  ;

    // Pressure
    units["Pa"] = unit::pascal ;
    units["bar"] = unit::bar  ;
    units["atm"] = unit::atmosphere  ;

    // Electric current
    units["A"] = unit::ampere ;
    units["mA"] = unit::milliampere ;
    units["muA"] = unit::microampere  ;
    units["nA"] = unit::nanoampere ;

    // Electric potential
    units["V"] = unit::volt ;
    units["kV"] = unit::kilovolt  ;
    units["MV"] = unit::megavolt ;

    // Magnetic flux density
    units["T"] = unit::tesla ;
    units["kG"] = unit::kilogauss  ;
    units["G"] = unit::gauss ;

    // Speed
    units["cm/us"] = units["cm"]/units["us"]  ;
    units["cm/ns"] = units["cm"]/units["ns"]  ;
    units["mm/us"] = units["mm"]/units["us"]  ;
    units["mm/ns"] = units["mm"]/units["ns"]  ;

    // Length/Energy
    units["mm/MeV"] = units["mm"]/units["MeV"];
    units["mm/keV"] = units["mm"]/units["keV"];
    units["cm/MeV"] = units["cm"]/units["MeV"];
    units["cm/keV"] = units["cm"]/units["keV"];

    // Dummy units for diffusion coefficient
    units["cm2/s"] = units["cm2"]/(units["s"]);

}

unit::TUnitsTable::~TUnitsTable() {}

std::string unit::TUnitsTable::ConvertWithUnit(std::string input) {

    double value;
    std::string unit;

    std::istringstream line(input);
    if (!(line >> value >> unit)) {
        PrintListOfUnits();
        CUBE_ERROR << "Badly formatted input string." << std::endl;
        std::exit(1);
    }

    // Check if requested unit is in map.
    if (units.find(unit) == units.end()) {
        PrintListOfUnits();
        CUBE_ERROR << "Tequested unit '"
                   << unit << "' not found."
                   << std::endl;
        std::exit(1);
    }

    value = value * units[unit];
    std::ostringstream result;
    result << value;

    return result.str();
}

std::string unit::TUnitsTable::ConvertEnergy(double enr) {
    std::stringstream output;

    output << std::setprecision(4);
    double val = enr/unit::PeV;
    if (std::abs(val) >= 1) {
        output << val << " PeV";
        return output.str();
    }

    val = enr/unit::TeV;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " TeV";
        return output.str();
    }

    val = enr/unit::GeV;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " GeV";
        return output.str();
    }

    val = enr/unit::MeV;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " MeV";
        return output.str();
    }

    val = enr/unit::keV;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " keV";
        return output.str();
    }

    val = enr/unit::eV;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " eV";
        return output.str();
    }

    val = enr;
    if (val < 0.001*unit::eV) val = 0.0;
    output << val/unit::MeV << " eV";
    return output.str();
}

std::string unit::TUnitsTable::ConvertLength(double len) {
    std::stringstream output;

    output << std::setprecision(4);
    double val = len/unit::km;
    if (std::abs(val) >= 1) {
        output << val << " km";
        return output.str();
    }

    val = len/unit::m;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " m";
        return output.str();
    }

    val = len/unit::cm;
    if (std::abs(val) <= 1200 && std::abs(val) >= 100) {
        output << val << " cm";
        return output.str();
    }

    val = len/unit::mm;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " mm";
        return output.str();
    }

    val = len/unit::micrometer;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " um";
        return output.str();
    }

    val = len/unit::nm;
    if (std::abs(val) <= 1200 && std::abs(val) >= 1) {
        output << val << " nm";
        return output.str();
    }

    val = len;
    if (val < 0.001*unit::mm) val = 0.0;
    output << val/unit::mm << " mm";
    return output.str();
}

std::string unit::TUnitsTable::ConvertCharge(double q) {
    std::stringstream output;

    output << std::setprecision(4);
    double val = q/unit::coulomb;
    if (std::abs(val) >= 1) {
        output << val << " C";
        return output.str();
    }

    val = q/unit::millicoulomb;
    if (std::abs(val) >= 1) {
        output << val << " mC";
        return output.str();
    }

    val = q/unit::microcoulomb;
    if (std::abs(val) >= 1) {
        output << val << " uC";
        return output.str();
    }

    val = q/unit::picocoulomb;
    if (std::abs(val) >= 1) {
        output << val << " pC";
        return output.str();
    }

    val = q/unit::femtocoulomb;
    if (std::abs(val) >= 1) {
        output << val << " fC";
        return output.str();
    }

    val = q/unit::attocoulomb;
    output << val << " aC";
    return output.str();
}

std::string unit::TUnitsTable::ConvertTime(double tim) {
    std::stringstream output;
    double val;

    val = tim/(365*24*3600*unit::second);
    if (std::abs(val) >= 100) {
        output << val << " yr";
        return output.str();
    }

    output << std::setprecision(4);
    val = tim/(7*24*3600*unit::second);
    if (std::abs(val) >= 100) {
        output << val << " wk";
        return output.str();
    }

    val = tim/(24*3600*unit::second);
    if (std::abs(val) >= 3) {
        output << val << " day";
        return output.str();
    }

    val = tim/(3600*unit::second);
    if (std::abs(val) >= 15) {
        output << val << " hr";
        return output.str();
    }

    val = tim/(60*unit::second);
    if (std::abs(val) > 15) {
        output << val << " min";
        return output.str();
    }

    val = tim/unit::second;
    if (std::abs(val) >= 1) {
        output << val << " s";
        return output.str();
    }

    val = tim/unit::millisecond;
    if (std::abs(val) >= 1) {
        output << val << " ms";
        return output.str();
    }

    val = tim/unit::microsecond;
    if (std::abs(val) >= 1) {
        output << val << " us";
        return output.str();
    }

    val = tim/unit::nanosecond;
    if (std::abs(val) >= 1) {
        output << val << " ns";
        return output.str();
    }

    val = tim;
    if (val < 0.001*unit::ns) val = 0.0;
    output << val/unit::ns << " ns";

    return output.str();
}

void unit::TUnitsTable::PrintListOfUnits() {
    std::cout << std::endl;
    std::cout << "***** List of available units *****" << std::endl;
    for (std::map<std::string, double>::iterator unit = units.begin();
         unit != units.end(); unit++) {
        std::cout << "    Unit Name: " << (*unit).first
                  << "    Multiplier: " << (*unit).second
                  <<  std::endl;
    }
    std::cout << std::endl;
}

std::string unit::AsString(double val, double sig, std::string type) {
    unit::TUnitsTable& u = unit::TUnitsTable::Get();
    std::ostringstream measure;
    if (type == "length") {
        measure << u.ConvertLength(val);
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << u.ConvertLength(sig);
        }
    }
    else if (type == "time") {
        measure << u.ConvertTime(val);
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << u.ConvertTime(sig);
        }
    }
    else if (type == "charge") {
        measure << u.ConvertCharge(val);
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << u.ConvertCharge(sig);
        }
    }
    else if (type == "momentum") {
        measure << u.ConvertEnergy(val) << "/c";
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << u.ConvertEnergy(sig) << "/c";
        }
    }
    else if (type == "energy") {
        measure << u.ConvertEnergy(val);
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << u.ConvertEnergy(sig);
        }
    }
    else if (type == "direction") {
        measure << std::fixed;
        measure << std::setprecision(3);
        measure << val;
        if (sig >= 0) {
            if (unit::IsFree(sig)) measure << " (free)";
            else if (unit::IsFixed(sig)) measure << " (fixed)";
            else measure << "+-" << sig;
        }
    }
    else if (type == "angle") {
        measure << std::fixed;
        measure << std::setprecision(1);
        measure << val/unit::degree;
        if (sig<0) measure << " deg";
        else if (unit::IsFree(sig)) measure << " deg (free)";
        else if (unit::IsFixed(sig)) measure << " deg (fixed)";
        else measure << "+-" << sig/unit::degree << " deg";
    }
    else if (type == "electrons" || type == "pe") {
        measure << std::fixed;
        measure << std::setprecision(1);
        measure << val;
        if (sig<0) measure << " e";
        else if (unit::IsFree(sig)) measure << " e (free)";
        else if (unit::IsFixed(sig)) measure << " e (fixed)";
        else measure << "+-" << sig << " e";
    }
    else if (type != "") {
        measure << std::fixed;
        measure << std::setprecision(3);
        measure << val;
        if (sig<0) measure << " <" << type << ">";
        else if (unit::IsFree(sig)) {
            measure << " <" << type << "> (free)";
        }
        else if (unit::IsFixed(sig)) {
            measure << " <" << type << "> (fixed)";
        }
        else measure << "+-" << sig << " <" << type << ">";
    }
    else {
        measure << std::fixed;
        measure << std::setprecision(4);
        measure << val;
        if (sig<0) measure << "";
        else if (unit::IsFree(sig)) {
            measure << " (free)";
        }
        else if (unit::IsFixed(sig)) {
            measure << " (fixed)";
        }
        else measure << "+-" << sig;
    }

    return measure.str();
}

std::string unit::AsString(double val, std::string type) {
    return AsString(val,-1.0,type);
}

std::string unit::AsString(const TVector3& val, std::string type) {
    std::ostringstream measure;
    if (type != "") {
        measure << "(" << unit::AsString(val.X(),type)
                << ", " << unit::AsString(val.Y(),type)
                << ", " << unit::AsString(val.Z(),type) << ")";
    }
    else {
        measure << "(" << unit::AsString(val.X(),type)
                << ", " << unit::AsString(val.Y(),type)
                << ", " << unit::AsString(val.Z(),type)
                << ")";
    }
    return measure.str();
}

std::string unit::AsString(const TLorentzVector& val, std::string type) {
    std::ostringstream measure;
    if (type =="length") {
        measure << "(" << unit::AsString(val.X(),"length")
                << ", " << unit::AsString(val.Y(),"length")
                << ", " << unit::AsString(val.Z(),"length")
                << ", " << unit::AsString(val.T(),"time")
                << ")";
    }
    else if (type =="momentum") {
        measure << "(" << unit::AsString(val.X(),"momentum")
                << ", " << unit::AsString(val.Y(),"momentum")
                << ", " << unit::AsString(val.Z(),"momentum")
                << ", " << unit::AsString(val.T(),"energy")
                << ")";
    }
    else {
        measure << "(" << unit::AsString(val.X(),type)
                << ", " << unit::AsString(val.Y(),type)
                << ", " << unit::AsString(val.Z(),type)
                << ", " << unit::AsString(val.T(),type)
                << ")";
    }
    return measure.str();
}

/// Check to see if a variance value means that the parameter is free.
bool unit::IsFree(double variance) {
    return (std::sqrt(std::numeric_limits<float>::max())<variance);
}

/// Check to see if a variance value means that the parameter is fixed.
bool unit::IsFixed(double variance) {
    return (variance<1.0/std::sqrt(std::numeric_limits<float>::max()));
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
