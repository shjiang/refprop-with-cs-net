using System;
//using System.Collections.Generic;
//using System.Linq;
using System.Text;

namespace RefpropCSNET
{
    //public class ThermoProperties
    //{
    //    //public TemperatureUnits Temperature;
    //    //public PressureUnits Pressure;
    //    double density;
    //    double internalenergy;
    //    double enthalpy;
    //    double entropy;
    //    double quality;
    //    double isochoric_specificheat;
    //    double isobaric_specificheat;
    //    double soundspeed;
    //    double density_liquidphase;				//molar density [mol/l] of the liquid phase
    //    double density_vaporphase;				//molar density [mol/l] of the vapor phase
    //    double[] xliq;							//composition of liquid phase [array of mol frac]
    //    double[] xvap;							//composition of vapor phase [array of mol frac]

    //    public double specific_volume
    //    {
    //        get { return 1 / density; }
    //    }
    //    public void showvalues()
    //    {
    //        //system.console.writeline("temperature = " + temperature.tostring());
    //        //system.console.writeline("pressure = " + pressure.tostring());
    //        //system.console.writeline("density = " + density.tostring());
    //        //system.console.writeline("enthalpy = " + enthalpy.tostring());
    //        //system.console.writeline("entropy = " + entropy.tostring());
    //        //system.console.writeline("internal energy = " + internalenergy.tostring());
    //        //system.console.writeline("quality = " + quality.tostring());

    //    }
    //}


    public enum ReferenceState
    {
        DEF,
        NBP,        //h,s = 0 at pure component normal boiling point(s)
        ASH = 3,        //h,s = 0 for sat liquid at -40 C (ASHRAE convention)
        IIR = 4         //:  h = 200, s = 1.0 for sat liq at 0 C (IIR convention)
    }
    public enum RefrigerantCategory
    {
        PureFluid,
        PredefinedMixture,
        NewMixture,
        PseudoPureFluid
    }

    public enum SaturationPoint
    {
        Bubble_Point = 1,	   //liquid composition (bubble point)
        Dew_Point = 2,	       //vapor composition (dew point)
        Freezing_Point = 3,    //freezing point
        Sublimation_Point = 4  //it is the melting point in case of solid-liquid phase change 

    };

    public enum SubstancePhase
    {
        Solid,
        Subcooled,
        Saturated_Liquid,
        Saturated_Vapor,
        Superheated = 998,
        Superheated_T_Over_Critical,
        Superheated_SuperCritical
    };

    //public struct SaturatedProperties
    //{
    //    double Temperature;
    //    double Pressure;
    //    double Density_SaturatedLiquid;							//rhol--molar density [mol/L] of saturated liquid
    //    double Density_SaturatedVapor;						//rhov--molar density [mol/L] of saturated vapor
    //};

    /*
      T = temperature [K]
      P = pressure [kPa]
      D = density [mol/L]  L=dm^3
      E = internal energy [J/mol]
      H = enthalpy [J/mol]
      S = entropy [J/mol-K]
      Q = vapor quality [moles vapor/total moles]
          or [kg vapor/total kg] depending on the
          value of the input flag kq

    */
    public enum ThermoProperties
    {
        Temperature,
        Pressure,		// 00000010
        Density,		// 00000100
        SpecificVolume,
        IsochoricHeatCapacity,
        IsobaricHeatCapacity,
        InternalEnergy,
        Enthalpy,       // 00010000
        Entropy,       // 00100000
        SpeedofSound,
        Viscosity,
        ThermalConductivity,
        SurfaceTension,
        Quality        // 01000000
    };
    public enum UnitTypes
    {
        T,    //         Temperature                            K
        P,    //         Pressure                               Pa
        D,    //         Density || specific volume         mol/m^3 || kg/m^3 (|| m^3/mol || m^3/kg)
        H,    //         Enthalpy || specific energy        J/mol || J/kg
        S,    //         Entropy || heat capacity           J/mol-K || J/kg-K
        W,    //         Speed of sound                         m/s
        U,    //         Viscosity                              Pa-s
        K,    //         Thermal conductivity                   W/m-K
        JT,   //         Joule Thompson                         K/Pa
        L,    //         Length                                 m
        A,    //         Area                                   m,2
        V,    //         Volume                                 m^3
        M,    //         Mass                                   kg
        F,    //         F||ce                                  N
        E,    //         Energy                                 J
        Q,    //         Power                                  W
        N,     //         Surface tension                        N/m
        X
    };

    /*		
		temperature                     K
		pressure, fugacity              kPa
		density                         mol/L
		composition                     mole fraction
		quality                         mole basis (moles vapor/total moles)
		enthalpy, internal energy       J/mol
		Gibbs, Helmholtz free energy    J/mol
		entropy, heat capacity          J/(mol.K)
		speed of sound                  m/s
		Joule-Thompson coefficient      K/kPa
		d(p)/d(rho)                     kPa.L/mol
		d2p)/d(rho)2                    kPa.(L/mol)^2
		viscosity                       microPa.s (10^-6 Pa.s)
		thermal conductivity            W/(m.K)
		dipole moment                   debye
		surface tension                 N/m
		*/
    public static class UnitConversion
    {
        public static double CtoK = 273.15;               //Exact conversion
        public static double FtoR = 459.67;            //Exact conversion
        public static double RtoK = 5.0 / 9.0;               //Exact conversion
        public static double HtoS = 3600.0;               //Exact conversion
        public static double ATMtoMPa = 0.101325;            //Exact conversion
        public static double BARtoMPA = 0.1;           //Exact conversion
        public static double KGFtoN = 98.0665 / 10.0;          //Exact conversion
        public static double INtoM = 0.0254;         //Exact conversion
        public static double FTtoM = 12.0 * INtoM;        //Exact conversion
        public static double LBMtoKG = 0.45359237;       //Exact conversion
        public static double CALtoJ = 4.184;      //Exact conversion (tc)
        //public static double CALtoJ = 4.1868         ;      //Exact conversion (IT) (Use this one only with pure water)
        public static double MMHGtoMPA = ATMtoMPa / 760.0;    //Exact conversion
        public static double INH2OtoMPA = 0.000249082;

        public static double BTUtoKJ = CALtoJ * LBMtoKG * RtoK;
        public static double LBFtoN = LBMtoKG * KGFtoN;
        public static double IN3toM3 = INtoM * INtoM * INtoM;
        public static double FT3toM3 = FTtoM * FTtoM * FTtoM;
        public static double GALLONtoM3 = IN3toM3 * 231.0;
        public static double PSIAtoMPA = LBMtoKG / INtoM / INtoM * KGFtoN / 1000000.0;
        public static double PSIAtoKPA = LBMtoKG / INtoM / INtoM * KGFtoN / 1000.0;
        public static double FTLBFtoJ = FTtoM * LBFtoN;
        public static double HPtoW = 550.0 * FTLBFtoJ;
        public static double BTUtoW = BTUtoKJ * 1000.0;
        public static double LBFTtoNM = LBFtoN / FTtoM;
        public static string[] T = new string[] { "K", "K", "C", "K", "K", "K", "F", "C" };
        public static string[] P = new string[] { "kPa", "MPa", "MPa", "MPa", "kPa", "MPa", "psia", "bar" };
        public static string[] D = new string[] { "mol/dm^3", "kg/m^3", "kg/m^3", "mol/dm^3", "kg/m^3", "g/cm^3", "lbm/ft^3", "g/cm^3", "dm^3/mol", "m^3/kg", "m^3/kg", "dm^3/mol", "m^3/kg", "cm^3/g", "ft^3/lbm", "cm^3/g" };
        //public static string[] V = new string[] {  };
        public static string[] H = new string[] { "J/mol", "kJ/kg", "kJ/kg", "J/mol", "kJ/kg", "J/g", "Btu/lbm", "J/g" };
        public static string[] S = new string[] { "J/mol-K", "kJ/kg-K", "kJ/kg-K", "J/mol-K", "kJ/kg-K", "J/g-K", "Btu/lbm-R", "J/g-K" };
        public static string[] W = new string[] { "m/s", "m/s", "m/s", "m/s", "m/s", "cm/s", "ft/s", "cm/s" };
        public static string[] U = new string[] { "uPa-s", "uPa-s", "uPa-s", "uPa-s", "uPa-s", "uPa-s", "lbm/ft-s", "centipoise" };
        public static string[] K = new string[] { "W/m-K", "mW/m-K", "mW/m-K", "mW/m-K", "W/m-K", "mW/m-K", "Btu/h-ft-F", "mW/m-K" };
        public static string[] N = new string[] { "N/m", "mN/m", "mN/m", "mN/m", "mN/m", "dyn/cm", "lbf/ft", "mN/m" };
        public static string[] X = new string[] { "-", "-", "-", "-", "-", "-", "-", "-" };

    }

    public enum UnitSystems
    {
        Refprop = 0,
        SI = 1,
        SIwithC = 2,
        MolarSI = 3,
        mks = 4,
        cgs = 5,
        E = 6,
        ME = 7
        //User = 8
    };
    public static class TemperatureUnitCatalog
    {
        public static string Kelvin = "K";
        public static string Celsius = "C";
        public static string Fahrenheit = "F";
        public static string Reaumur = "R";
    }
    public static class PressureUnitCatalog
    {
        public static string Pa = "Pa";
        public static string kPa = "kPa";
        public static string MPa = "MPa";
        public static string GPa = "GPa";
        public static string Bar = "Bar";
        public static string KBar = "KBar";
        public static string ATM = "ATM";
        public static string KGFperM2 = "KGF/m^2";
        public static string PSI = "psi";
        public static string PSF = "psf";
        public static string mmHG = "mmHG";
        public static string cmHG = "cmHG";
        public static string inHG = "inHG";
        public static string PSIG = "psig";
        public static string inH2O = "inH2O";
    }
    public static class DensityUnitCatalog
    {
        public static string MOLperDM3 = "MOL/DM^3";
        public static string MOLperL = "MOL/L";
        public static string KMOLperM3 = "KMOL/M^3";
        public static string MOLperCM3 = "MOL/CM^3";
        public static string MOLperCC = "MOL/CC";
        public static string MOLperM3 = "MOL/M^3";
        public static string KGperM3 = "KG/M^3";
        public static string KGperDM3 = "KG/DM^3";
        public static string KGperL = "KG/L";
        public static string GperDM3 = "G/DM^3";
        public static string GperL = "G/L";
        public static string GperCC = "G/CC";
        public static string GperCM3 = "G/CM^3";
        public static string GperML = "G/ML";
        public static string LBMperFT3 = "LBM/FT^3";
        public static string LBperFT3 = "LB/FT^3";
        public static string LBMOLperFT3 = "LBMOL/FT^3";
        public static string SLUGperFT3 = "SLUG/FT^3";
        public static string LBperGAL = "LB/GAL";
    }
    public static class SpecificVolumeUnitCatalog
    {
        public static string DM3perMOL = "DM^3/MOL";
        public static string LperMOL = " L/MOL";
        public static string M3perKMOL = "M^3/KMOL";
        public static string CM3perMOL = "CM^3/MOL";
        public static string CCperMOL = "CC/MOL";
        public static string M3perMOL = "M^3/MOL";
        public static string M3perKG = "M^3/KG";
        public static string DM3perKG = "DM^3/KG ";
        public static string LperKG = "L/KG";
        public static string DM3perG = "DM^3/G";
        public static string LperG = "L/G";
        public static string CCperG = "CC/G";
        public static string CM3perG = "CM^3/G";
        public static string MLperG = "ML/G";
        public static string FT3perLBM = "FT^3/LBM";
        public static string FT3perLB = " FT^3/LB";
        public static string FT3perLBMOL = "FT^3/LBMOL";
        public static string FT3perSLUG = "FT^3/SLUG ";
        public static string GALperLB = " GAL/LB";
    }
    public static class SpecificEnergyUnitCatalog
    {
        public static string JperMOL = " J/MOL ";
        public static string KJperKMOL = " KJ/KMOL ";
        public static string KJperMOL = " KJ/MOL ";
        public static string MJperMOL = " MJ/MOL ";
        public static string KJperKG = " KJ/KG ";
        public static string JperG = " J/G ";
        public static string JperKG = " J/KG ";
        public static string M2perS2 = " M^2/S^2 ";
        public static string FT2perS2 = " FT^2/S^2 ";
        public static string CALperMOL = " CAL/MOL ";
        public static string KCALperKMOL = " KCAL/KMOL ";
        public static string KCALperKG = " KCAL/KG ";
        public static string CALperG = " CAL/G ";
        public static string BTUperLBM = " BTU/LBM ";
        public static string BTUperLBMOL = " BTU/LBMOL ";
        public static string BTUperLB = " BTU/LB ";
    }
    public static class EnthalpyUnitCatalog
    {
        public static string JperMOL = " J/MOL ";
        public static string KJperKMOL = " KJ/KMOL ";
        public static string KJperMOL = " KJ/MOL ";
        public static string MJperMOL = " MJ/MOL ";
        public static string KJperKG = " KJ/KG ";
        public static string JperG = " J/G ";
        public static string JperKG = " J/KG ";
        public static string M2perS2 = " M^2/S^2 ";
        public static string FT2perS2 = " FT^2/S^2 ";
        public static string CALperMOL = " CAL/MOL ";
        public static string KCALperKMOL = " KCAL/KMOL ";
        public static string KCALperKG = " KCAL/KG ";
        public static string CALperG = " CAL/G ";
        public static string BTUperLBM = " BTU/LBM ";
        public static string BTUperLBMOL = " BTU/LBMOL ";
        public static string BTUperLB = " BTU/LB ";
    }
    public static class EntropyUnitCatalog
    {
        public static string JperMOLK = " J/MOL-K ";
        public static string KJperKMOLK = " KJ/KMOL-K ";
        public static string KJperMOLK = " KJ/MOL-K ";
        public static string JperGK = " J/G-K ";
        public static string KJperKGK = " KJ/KG-K ";
        public static string JperKGK = " J/KG-K ";
        public static string BTUperLBR = " BTU/LB-R ";
        public static string BTUperLBMR = " BTU/LBM-R ";
        public static string BTUperLBMOLR = " BTU/LBMOL-R ";
        public static string CALperGK = " CAL/G-K ";
        public static string CALperMOLK = " CAL/MOL-K ";
        public static string CALperGC = " CAL/G-C ";
        public static string KCALperKGK = " KCAL/KG-K ";
        public static string KCALperKGC = " KCAL/KG-C ";
        public static string CALperMOLC = " CAL/MOL-C ";
        public static string FTLBFperLBMOLR = " FT-LBF/LBMOL-R";
        public static string CPperR = " CP/R";
    }
    public static class HeatCApacityUnitCatalog
    {
        public static string JperMOLK = " J/MOL-K ";
        public static string KJperKMOLK = " KJ/KMOL-K ";
        public static string KJperMOLK = " KJ/MOL-K ";
        public static string JperGK = " J/G-K ";
        public static string KJperKGK = " KJ/KG-K ";
        public static string JperKGK = " J/KG-K ";
        public static string BTUperLBR = " BTU/LB-R ";
        public static string BTUperLBMR = " BTU/LBM-R ";
        public static string BTUperLBMOLR = " BTU/LBMOL-R ";
        public static string CALperGK = " CAL/G-K ";
        public static string CALperMOLK = " CAL/MOL-K ";
        public static string CALperGC = " CAL/G-C ";
        public static string KCALperKGK = " KCAL/KG-K ";
        public static string KCALperKGC = " KCAL/KG-C ";
        public static string CALperMOLC = " CAL/MOL-C ";
        public static string FTLBFperLBMOLR = " FT-LBF/LBMOL-R";
        public static string CPperR = " CP/R";
    }
    public static class SpeedofSoundUnitCatalog
    {
        public static string W = " W ";
        public static string MperS = " M/S ";
        public static string M2perS2 = " M^2/S^2 ";
        public static string CMperS = " CM/S ";
        public static string KMperH = " KM/H ";
        public static string FTperS = " FT/S ";
        public static string INperS = " IN/S ";
        public static string MILEperH = " MILE/H ";
        public static string MPH = " MPH ";
        public static string KNOT = " KNOT ";
        public static string MACH = " MACH ";
    }
    public static class ViscosityUnitCatalog
    {
        public static string PAS = " PA-S ";
        public static string MPAS = " MPA-S ";
        public static string KGperMS = " KG/M-S ";
        public static string UPAS = " UPA-S ";
        public static string POISE = " POISE ";
        public static string GperCMS = " G/CM-S ";
        public static string MPOISE = " MPOISE ";
        public static string CENTIPOISE = " CENTIPOISE ";
        public static string UPOISE = " UPOISE ";
        public static string MILLIPOISE = " MILLIPOISE ";
        public static string LBperFTS = " LB/FT-S ";
        public static string MICROPOISE = " MICROPOISE ";
        public static string LBperFTH = " LB/FT-H ";
        public static string LBMperFTS = " LBM/FT-S ";
        public static string LBMperFTH = " LBM/FT-H ";
        public static string LBFSperFT2 = " LBF-S/FT^2 ";
    }
    public static class ThermalConductivityUnitCatalog
    {
        public static string MWperMK = " MW/M-K ";
        public static string WperMK = " W/M-K ";
        public static string GCMperS3K = " G-CM/S^3-K ";
        public static string KGMperS3K = " KG-M/S^3-K ";
        public static string CALperSCMK = " CAL/S-CM-K ";
        public static string KCALperHRMK = " KCAL/HR-M-K ";
        public static string LBMFTperS3F = " LBM-FT/S^3-F ";
        public static string LBFTperS3F = " LB-FT/S^3-F ";
        public static string LBFperSF = " LBF/S-F ";
        public static string BTUperHFTF = " BTU/H-FT-F ";
    }
    public static class JouleThomsonUnitCatalog
    {
        public static string KperMPA = " K/MPA ";
        public static string KperKPA = " K/KPA ";
        public static string CperMPA = " C/MPA ";
        public static string KperPA = " K/PA ";
        public static string CperKPA = " C/KPA ";
        public static string CperATM = " C/ATM ";
        public static string CperPA = " C/PA ";
        public static string CperBAR = " C/BAR ";
        public static string KperPSIA = " K/PSIA ";
        public static string KperPSI = " K/PSI ";
        public static string FperPSIA = " F/PSIA ";
        public static string FperPSI = " F/PSI ";
        public static string RperPSIA = " R/PSIA ";
    }
    public static class LengthUnitCatalog
    {
        public static string L = " L ";
        public static string METER = " METER ";
        public static string DM = " DM ";
        public static string M = " M ";
        public static string CM = " CM ";
        public static string IN = " IN ";
        public static string MM = " MM ";
        public static string FT = " FT ";
        public static string KM = " KM ";
        public static string YD = " YD ";
        public static string INCH = " INCH ";
        public static string MI = " MI ";
        public static string FOOT = " FOOT ";
        public static string YARD = " YARD ";
        public static string MILE = " MILE ";
        public static string LIGHTYEAR = " LIGHT YEAR ";
        public static string ANGSTROM = " ANGSTROM ";
        public static string FATHOM = " FATHOM ";
        public static string MIL = " MIL ";
        public static string ROD = " ROD ";
        public static string PARSEC = " PARSEC ";
    }
    public static class AreaUnitCatalog
    {
        public static string METER2 = " METER^2 ";
        public static string CM2 = " CM^2 ";
        public static string M2 = " M^2 ";
        public static string MM2 = " MM^2 ";
        public static string IN2 = " IN^2 ";
        public static string KM2 = " KM^2 ";
        public static string FT2 = " FT^2 ";
        public static string INCH2 = " INCH^2 ";
        public static string YD2 = " YD^2 ";
        public static string FOOT2 = " FOOT^2 ";
        public static string MI2 = " MI^2 ";
        public static string YARD2 = " YARD^2 ";
        public static string MILE2 = " MILE^2 ";
        public static string ACRE = " ACRE ";
        public static string BARN = " BARN ";
        public static string HECTARE = " HECTARE ";
    }
    public static class VolumeUnitCatalog
    {
        public static string METER3 = " METER^3 ";
        public static string CM3 = " CM^3 ";
        public static string M3 = " M^3 ";
        public static string LITER = " LITER ";
        public static string L = " L ";
        public static string INCH3 = " INCH^3 ";
        public static string DM3 = " DM^3 ";
        public static string FOOT3 = " FOOT^3 ";
        public static string IN3 = " IN^3 ";
        public static string YARD3 = " YARD^3 ";
        public static string FT3 = " FT^3 ";
        public static string GALLON = " GALLON ";
        public static string YD3 = " YD^3 ";
        public static string QUART = " QUART ";
        public static string GAL = " GAL ";
        public static string PINT = " PINT ";
        public static string QT = " QT ";
        public static string CUP = " CUP ";
        public static string PT = " PT ";
        public static string OUNCE = " OUNCE ";
        public static string TBSP = " TBSP ";
        public static string TABLESPOON = " TABLESPOON ";
        public static string TSP = " TSP ";
        public static string TEASPOON = " TEASPOON ";
        public static string CORD = " CORD ";
        public static string BARREL = " BARREL ";
        public static string BOARDFOOT = " BOARD FOOT ";
        public static string BUSHEL = " BUSHEL ";
    }
    public static class MassUnitCatalog
    {
        public static string KG = " KG ";
        public static string G = " G ";
        public static string MG = " MG ";
        public static string LB = " LB ";
        public static string LBM = " LBM ";
        public static string GRAIN = " GRAIN ";
        public static string SLUG = " SLUG ";
        public static string TON = " TON ";
        public static string TONNE = " TONNE ";

    }
    public static class ForceUnitCatalog
    {
        public static string NEWTON = " NEWTON ";
        public static string MN = " MN ";
        public static string N = " N ";
        public static string KGF = " KGF ";
        public static string DYNE = " DYNE ";
        public static string LBF = " LBF ";
        public static string POUNDAL = " POUNDAL ";
        public static string OZF = " OZF ";
    }
    public static class EnergyUnitCatalog
    {
        public static string JOULE = " JOULE ";
        public static string KJ = " KJ ";
        public static string J = " J ";
        public static string MJ = " MJ ";
        public static string KWH = " KW-H ";
        public static string CAL = " CAL ";
        public static string KCAL = " KCAL ";
        public static string ERG = " ERG ";
        public static string BTU = " BTU ";
        public static string FTLBF = " FT-LBF ";
    }
    public static class PowerUnitCatalog
    {
        public static string WATT = " WATT ";
        public static string KWATT = " KWATT ";
        public static string W = " W ";
        public static string BTUperS = " BTU/S ";
        public static string KW = " KW ";
        public static string BTUperMIN = " BTU/MIN ";
        public static string BTUperH = " BTU/H ";
        public static string CALperS = " CAL/S ";
        public static string KCALperS = " KCAL/S ";
        public static string CALperMIN = " CAL/MIN ";
        public static string KCALperMIN = " KCAL/MIN ";
        public static string FTLBFperS = " FT-LBF/S ";
        public static string FTLBFperMIN = " FT-LBF/MIN ";
        public static string FTLBFperH = " FT-LBF/H ";
        public static string HP = " HP ";
    }
    public static class SurfaceTensionUnitCatalog
    {
        public static string NperM = " N/M ";
        public static string MNperM = " MN/M ";
        public static string DYNperCM = " DYN/CM ";
        public static string DYNEperCM = " DYNE/CM ";
        public static string LBFperFT = " LBF/FT ";
    }

    public static class RefpropErrorHandler
    {
        public static void ErrorHandler(object source, string err, int errcode)
        {
            string msg = source.ToString() + "Error:[" + err + "],Code:[" + errcode.ToString() + "]";
            throw new Exception(msg);
            //Console.WriteLine(String.Format("Source:{0},Error:{1},Code:{2}", source.ToString(), err, errcode));
            //Console.WriteLine("== Please Press Any Key to Exit ==");
            //Console.ReadKey();
            //Environment.Exit(errcode);
        }
        public static void ErrorHandler(object source, string err)
        {
            string msg = source.ToString() + "Error:[" + err + "]";
            throw new Exception(msg);
            //Console.WriteLine(String.Format("Source:{0},Error:{1}", source.ToString(), err.ToString()));
            //Console.WriteLine("== Please Press Any Key to Exit ==");
            //Console.ReadKey();
            //Environment.Exit(0);
        }
        public static void ErrorHandler(string err)
        {
            string msg = "Error:[" + err + "]";
            throw new Exception(msg);
            //Console.WriteLine(String.Format("Error:{0}", err.ToString()));
            //Console.WriteLine("== Please Press Any Key to Exit ==");
            //Console.ReadKey();
            //Environment.Exit(0);
        }
    }

    public static class Global
    {
        //public static bool isMassQ { get; set; }
        public const double NearZero = 1.0E-6;
        //public const double AsoluteZero = -273.15;
        public static UnitSystems UnitSystem;
    }
}
