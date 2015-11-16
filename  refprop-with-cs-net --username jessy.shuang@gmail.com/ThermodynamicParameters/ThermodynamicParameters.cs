using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace RefpropCSNET
{
    public class ThermodynamicParameter
    {
        //private string currentUnits;
        //private double currentValue;
        private UnitSystems unitSystem;
        private string unit;
        //private Refrigerant currentRefrigerant = new Refrigerant();
        public double MolecularWeight { get; set; }
        //private ThermoProperties ThermoProperty { set; get; }
        public UnitTypes UnitType { private set; get; } //不可更改
        /// <summary>
        /// 
        /// unit只能从unit system获得，而不能由外部直接修改
        /// </summary>
        public string Unit
        {
            protected set
            {
                unit = value;
            }
            get { return unit; }
        }
        public double Value { set; get; }

        /// <summary>
        /// 根据UnitSystem 和UnitType确定Unit.
        /// </summary>
        public UnitSystems UnitSystem
        {
            get
            {
                return this.unitSystem;
            }
            set
            {
                if (this.unitSystem != value)
                {
                    try
                    {
                        //利用映射
                        string[] temp = (string[])(typeof(UnitConversion).GetField(UnitType.ToString()).GetValue(null));
                        this.Value = UnitConvert(this.Value, UnitType.ToString(), this.Unit, temp[(int)value], MolecularWeight);
                        this.Unit = temp[(int)value];
                    }
                    catch
                    {
                        RefpropErrorHandler.ErrorHandler(this, "Unknown unit", 002);
                    }

                    this.unitSystem = value;
                }
            }
        }
        public ThermodynamicParameter()
        {
        }
        public ThermodynamicParameter(double value)  //????
            : this(0, UnitTypes.X, value, UnitSystems.SI)
        {
        }
        public ThermodynamicParameter(UnitTypes type, UnitSystems unitsys)
            : this(0, type, 0, unitsys)
        {
        }
        public ThermodynamicParameter(UnitTypes type, double value, UnitSystems unitsys)
            : this(0, type, 0, unitsys)
        {
        }

        public ThermodynamicParameter(double wm, UnitTypes utp, double value, UnitSystems unitsys)
        {
            this.MolecularWeight = wm;
            this.UnitType = utp;
            this.Value = value;
            try
            {
                //利用映射
                string[] temp = (string[])(typeof(UnitConversion).GetField(UnitType.ToString()).GetValue(null));
                this.Unit = temp[(int)unitsys];
            }
            catch
            {
                RefpropErrorHandler.ErrorHandler(this, "Unknown unit", 1000);
            }

            this.UnitSystem = unitsys;
        }
        public ThermodynamicParameter(double wm, UnitTypes utp, double value, UnitSystems unitsys, List<ThermodynamicParameter> List)
        {
            this.MolecularWeight = wm;
            this.UnitType = utp;
            this.Value = value;
            try
            {
                //利用映射
                string[] temp = (string[])(typeof(UnitConversion).GetField(UnitType.ToString()).GetValue(null));
                this.Unit = temp[(int)unitsys];
            }
            catch
            {
                RefpropErrorHandler.ErrorHandler(this, "Unknown unit", 1000);
            }

            this.UnitSystem = unitsys;
            List.Add(this);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="unit"></param>
        /// <returns></returns>
        public bool CheckUnit(string unit)
        {
            unit = unit.Replace(" ", "");
            //利用映射
            string[] temp = (string[])(typeof(UnitConversion).GetField(UnitType.ToString()).GetValue(null));
            for (int i = 0; i < temp.Length; i++)
            {
                if (temp[i].Equals(unit, StringComparison.OrdinalIgnoreCase))
                {
                    //this.CurrentUnitSystem = (UnitSystem)i;
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="InputValue">InputValue is the value to be converted from OldUnits to NewUnits</param>
        /// <param name="UnitType">UnitType is one of the following letters (one character only in most cases):</param>
        /// //'UnitType     Unit name                          SI units
        ///'  T         Temperature                            K
        ///'  P         Pressure                               Pa
        ///'  D         Density || specific volume         mol/m^3 || kg/m^3 (|| m^3/mol || m^3/kg)
        ///'  H         Enthalpy || specific energy        J/mol || J/kg
        ///'  S         Entropy || heat capacity           J/mol-K || J/kg-K
        ///'  W         Speed of sound                         m/s
        ///'  U         Viscosity                              Pa-s
        ///'  K         Thermal conductivity                   W/m-K
        ///'  JT        Joule Thompson                         K/Pa
        ///'  L         Length                                 m
        ///'  A         Area                                   m,2
        ///'  V         Volume                                 m^3
        ///'  M         Mass                                   kg
        ///'  F         F||ce                                  N
        ///'  E         Energy                                 J
        ///'  Q         Power                                  W
        ///'  N         Surface tension                        N/m
        ///' Gage pressures can be used by adding "_g" to the unit, e.g., "MPa_g"      
        /// <param name="OldUnits"></param>
        /// <param name="NewUnits"></param>
        /// <param name="wm"></param>
        /// <returns></returns>
        public double UnitConvert(double InputValue, string UnitType, string OldUnits, string NewUnits, double wm)
        {
            double Value, MolWt, Rgas;
            string Tpe, Unit1, Unit2;
            Int32 Drct, Gage = 0, Vacm = 0;

            Value = InputValue;
            Tpe = UnitType.ToUpper().Trim();
            Unit1 = OldUnits.ToUpper().Trim();
            Unit2 = NewUnits.ToUpper().Trim();

            Rgas = 8.314472;
            //Call WMOLdll(x(1), wm)
            //if CompFlag = 1  Call WMOLdll(xliq(1), wm)
            //if CompFlag = 2  Call WMOLdll(xvap(1), wm)
            MolWt = wm;

            for (Drct = 1; Drct > -2; Drct -= 2)  // To -1 Step -2
            {
                //'-----------------------------------------------------------------------
                //'   Temperature Conversion
                //'-----------------------------------------------------------------------
                if (Tpe == "T")
                {
                    if (Unit1 == "K") { }
                    else if (Unit1 == "C") Value = Value + Drct * UnitConversion.CtoK;
                    else if (Unit1 == "R") Value = Value * Math.Pow(UnitConversion.RtoK, Drct);
                    else if (Unit1 == "F")
                    {
                        if (Drct == 1)
                            //Value = RtoK * (Value + Ft||)    'Does not give exactly zero at 32 F
                            Value = (Value - 32) * UnitConversion.RtoK + UnitConversion.CtoK;
                        else
                            //Value = Value / RtoK - Ft||      'Does not give exactly 32 at 273.15 K
                            Value = (Value - UnitConversion.CtoK) / UnitConversion.RtoK + 32.0;
                    }
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Pressure Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "P")
                {
                    Gage = Unit1.IndexOf("GAGE");
                    Vacm = Unit1.IndexOf("VACM");
                    if (Gage == -1) Gage = Unit1.IndexOf("_G");
                    if (Vacm == -1) Vacm = Unit1.IndexOf("_V");
                    if (Gage != -1 && Drct == -1) Value = Value - UnitConversion.ATMtoMPa;
                    if (Vacm != -1 && Drct == -1) Value = UnitConversion.ATMtoMPa - Value;
                    if (Gage != -1) Unit1 = Unit1.Substring(0, Gage).Trim();
                    if (Vacm != -1) Unit1 = Unit1.Substring(0, Vacm).Trim();
                    if (Unit1 == "PA") Value = Value / Math.Pow(1000000, Drct);
                    else if (Unit1 == "KPA") Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "MPA") { }
                    else if (Unit1 == "GPA") Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "BAR") Value = Value * Math.Pow(UnitConversion.BARtoMPA, Drct);
                    else if (Unit1 == "KBAR") Value = Value * Math.Pow((UnitConversion.BARtoMPA * 1000), Drct);
                    else if (Unit1 == "ATM") Value = Value * Math.Pow(UnitConversion.ATMtoMPa, Drct);
                    else if (Unit1 == "KGF/CM^2" || Unit1 == "KG/CM^2" || Unit1 == "ATA" || Unit1 == "AT" || Unit1 == "ATMA")
                        Value = Value * Math.Pow((UnitConversion.KGFtoN / 100), Drct);
                    else if (Unit1 == "PSI" || Unit1 == "PSIA")
                        Value = Value * Math.Pow(UnitConversion.PSIAtoMPA, Drct);
                    else if (Unit1 == "PSF")
                        Value = Value * Math.Pow((UnitConversion.PSIAtoMPA / 144), Drct);
                    else if (Unit1 == "MMHG" || Unit1 == "TORR")
                        Value = Value * Math.Pow(UnitConversion.MMHGtoMPA, Drct);
                    else if (Unit1 == "CMHG")
                        Value = Value * Math.Pow((UnitConversion.MMHGtoMPA * 10), Drct);
                    else if (Unit1 == "INHG")
                        Value = Value * Math.Pow((UnitConversion.MMHGtoMPA * UnitConversion.INtoM * 1000), Drct);
                    else if (Unit1 == "INH2O")
                        Value = Value * Math.Pow(UnitConversion.INH2OtoMPA, Drct);
                    else if (Unit1 == "PSIG")
                    {
                        if (Drct == 1)
                            Value = UnitConversion.PSIAtoMPA * Value + UnitConversion.ATMtoMPa;
                        else
                            Value = (Value - UnitConversion.ATMtoMPa) / UnitConversion.PSIAtoMPA;
                    }
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                    if (Gage != -1 && Drct == 1) Value = Value + UnitConversion.ATMtoMPa;
                    if (Vacm != -1 && Drct == 1) Value = UnitConversion.ATMtoMPa - Value;
                }
                //'-----------------------------------------------------------------------
                //'   Density Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "D")
                {
                    if (Value == 0) Value = 1E-50;
                    if (Unit1 == "MOL/DM^3" || Unit1 == "MOL/L" || Unit1 == "KMOL/M^3") { }
                    else if (Unit1 == "MOL/CM^3" || Unit1 == "MOL/CC")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "MOL/M^3")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "KG/M^3")
                        Value = Value / Math.Pow(MolWt, Drct);
                    else if (Unit1 == "KG/DM^3" || Unit1 == "KG/L")
                        Value = Value * Math.Pow((1000 / MolWt), Drct);
                    else if (Unit1 == "G/DM^3" || Unit1 == "G/L")
                        Value = Value * Math.Pow((1 / MolWt), Drct);
                    else if (Unit1 == "G/CC" || Unit1 == "G/CM^3" || Unit1 == "G/ML")
                        Value = Value * Math.Pow((1000 / MolWt), Drct);
                    else if (Unit1 == "G/DM^3")
                        Value = Value * Math.Pow((1 / MolWt), Drct);
                    else if (Unit1 == "LBM/FT^3" || Unit1 == "LB/FT^3")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.FT3toM3 / MolWt), Drct);
                    else if (Unit1 == "LBMOL/FT^3")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.FT3toM3), Drct);
                    else if (Unit1 == "SLUG/FT^3")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.FT3toM3 / MolWt * UnitConversion.KGFtoN / UnitConversion.FTtoM), Drct);
                    else if (Unit1 == "LB/GAL")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.GALLONtoM3 / MolWt), Drct);

                //'-----------------------------------------------------------------------
                    //'   Specific Volume Conversion
                    //'-----------------------------------------------------------------------
                    else if (Unit1 == "DM^3/MOL" || Unit1 == "L/MOL" || Unit1 == "M^3/KMOL")
                        Value = 1 / Value;
                    else if (Unit1 == "CM^3/MOL" || Unit1 == "CC/MOL" || Unit1 == "ML/MOL")
                        Value = 1000 / Value;
                    else if (Unit1 == "M^3/MOL")
                        Value = 1 / Value / 1000;
                    else if (Unit1 == "M^3/KG")
                        Value = 1 / Value / MolWt;
                    else if (Unit1 == "DM^3/KG" || Unit1 == "L/KG")
                        Value = 1000 / Value / MolWt;
                    else if (Unit1 == "CC/G" || Unit1 == "CM^3/G" || Unit1 == "ML/G")
                        Value = 1000 / Value / MolWt;
                    else if (Unit1 == "DM^/G")
                        Value = 1 / Value / MolWt;
                    else if (Unit1 == "FT^3/LBM" || Unit1 == "FT^3/LB")
                        Value = 1 / Value * (UnitConversion.LBMtoKG / UnitConversion.FT3toM3 / MolWt);
                    else if (Unit1 == "FT^3/LBMOL")
                        Value = 1 / Value * (UnitConversion.LBMtoKG / UnitConversion.FT3toM3);
                    else if (Unit1 == "FT^3/SLUG")
                        Value = 1 / Value * (UnitConversion.LBMtoKG / UnitConversion.FT3toM3 / MolWt * UnitConversion.KGFtoN / UnitConversion.FTtoM);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                    if (Math.Abs(Value) < 1E-30) Value = 0;
                }
                //'-----------------------------------------------------------------------
                //'   Specific Energy && Enthalpy Conversions
                //'-----------------------------------------------------------------------
                else if (Tpe == "H")
                {
                    if (Unit1 == "J/MOL" || Unit1 == "KJ/KMOL") { }
                    else if (Unit1 == "KJ/MOL")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "MJ/MOL")
                        Value = Value * Math.Pow(1000000, Drct);
                    else if (Unit1 == "KJ/KG" || Unit1 == "J/G")
                        Value = Math.Pow(MolWt, Drct) * Value;
                    else if (Unit1 == "J/KG")
                        Value = Math.Pow((MolWt / 1000), Drct) * Value;
                    else if (Unit1 == "M^2/S^2")
                        Value = Math.Pow((MolWt / 1000), Drct) * Value;
                    else if (Unit1 == "FT^/S^2")
                        Value = Math.Pow((MolWt / 1000 * Math.Pow(UnitConversion.FTtoM, 2)), Drct) * Value;
                    else if (Unit1 == "CAL/MOL" || Unit1 == "KCAL/KMOL")
                        Value = Math.Pow(UnitConversion.CALtoJ, Drct) * Value;
                    else if (Unit1 == "CAL/G" || Unit1 == "KCAL/KG")
                        Value = Math.Pow((UnitConversion.CALtoJ * MolWt), Drct) * Value;
                    else if (Unit1 == "BTU/LBM" || Unit1 == "BTU/LB")
                        Value = Math.Pow((UnitConversion.BTUtoKJ / UnitConversion.LBMtoKG * MolWt), Drct) * Value;
                    else if (Unit1 == "BTU/LBMOL")
                        Value = Math.Pow((UnitConversion.BTUtoKJ / UnitConversion.LBMtoKG), Drct) * Value;
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Entropy && Heat Capacity Conversions
                //'-----------------------------------------------------------------------
                else if (Tpe == "S")
                {
                    if (Unit1 == "J/MOL-K" || Unit1 == "KJ/KMOL-K") { }
                    else if (Unit1 == "KJ/MOL-K")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "KJ/KG-K" || Unit1 == "J/G-K")
                        Value = Math.Pow(MolWt, Drct) * Value;
                    else if (Unit1 == "J/KG-K")
                        Value = Math.Pow((MolWt / 1000), Drct) * Value;
                    else if (Unit1 == "BTU/LBM-R" || Unit1 == "BTU/LB-R")
                        Value = Math.Pow((UnitConversion.BTUtoKJ / UnitConversion.LBMtoKG / UnitConversion.RtoK * MolWt), Drct) * Value;
                    else if (Unit1 == "BTU/LBMOL-R")
                        Value = Math.Pow((UnitConversion.BTUtoKJ / UnitConversion.LBMtoKG / UnitConversion.RtoK), Drct) * Value;
                    else if (Unit1 == "CAL/G-K" || Unit1 == "CAL/G-C" || Unit1 == "KCAL/KG-K" || Unit1 == "KCAL/KG-C")
                        Value = Math.Pow((UnitConversion.CALtoJ * MolWt), Drct) * Value;
                    else if (Unit1 == "CAL/MOL-K" || Unit1 == "CAL/MOL-C")
                        Value = Math.Pow(UnitConversion.CALtoJ, Drct) * Value;
                    else if (Unit1 == "FT-LBF/LBMOL-R")
                        Value = Math.Pow((UnitConversion.FTLBFtoJ / UnitConversion.LBMtoKG / UnitConversion.RtoK / 1000), Drct) * Value;
                    else if (Unit1 == "CP/R")
                        Value = Math.Pow(Rgas, Drct) * Value * 1000;
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Speed of Sound Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "W")
                {
                    if (Unit1 == "M/S") { }
                    else if (Unit1 == "M^2/S^2")
                        Value = Math.Sqrt(Value);
                    else if (Unit1 == "CM/S")
                        Value = Math.Pow(Value / 100, Drct);
                    else if (Unit1 == "KM/H")
                        Value = Math.Pow(Value * (1000 / UnitConversion.HtoS), Drct);
                    else if (Unit1 == "FT/S")
                        Value = Value * Math.Pow(UnitConversion.FTtoM, Drct);
                    else if (Unit1 == "IN/S")
                        Value = Value * Math.Pow(UnitConversion.INtoM, Drct);
                    else if (Unit1 == "MILE/H" || Unit1 == "MPH")
                        Value = Value * Math.Pow((UnitConversion.INtoM * 63360 / UnitConversion.HtoS), Drct);
                    else if (Unit1 == "KNOT")
                        Value = Value * Math.Pow(0.5144444444, Drct);
                    else if (Unit1 == "MACH")
                        Value = Value * Math.Pow(Math.Sqrt(1.4 * 298.15 * 8314.51 / 28.95853816), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Viscosity Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "U")
                {
                    if (Unit1 == "PA-S" || Unit1 == "KG/M-S") { }
                    else if (Unit1 == "MPA-S")   //    'Note:  This is milliPa-s, not MPa-s
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "UPA-S")
                        Value = Value / Math.Pow(1000000, Drct);
                    else if (Unit1 == "G/CM-S" || Unit1 == "POISE")
                        Value = Value / Math.Pow(10, Drct);
                    else if (Unit1 == "CENTIPOISE")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "MILLIPOISE" || Unit1 == "MPOISE")
                        Value = Value / Math.Pow(10000, Drct);
                    else if (Unit1 == "MICROPOISE" || Unit1 == "UPOISE")
                        Value = Value / Math.Pow(10000000, Drct);
                    else if (Unit1 == "LBM/FT-S" || Unit1 == "LB/FT-S")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.FTtoM), Drct);
                    else if (Unit1 == "LBF-S/FT,2")
                        Value = Value * Math.Pow((UnitConversion.LBFtoN / Math.Pow(UnitConversion.FTtoM, 2)), Drct);
                    else if (Unit1 == "LBM/FT-H" || Unit1 == "LB/FT-H")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / UnitConversion.FTtoM / UnitConversion.HtoS), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Thermal Conductivity Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "K")
                {
                    if (Unit1 == "MW/M-K") { }
                    else if (Unit1 == "W/M-K")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "G-CM/S^3-K")
                        Value = Value / Math.Pow(100, Drct);
                    else if (Unit1 == "KG-M/S^3-K")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "CAL/S-CM-K")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ * 100000), Drct);
                    else if (Unit1 == "KCAL/HR-M-K")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ * 100000 * 1000 / 100 / 3600), Drct);
                    else if (Unit1 == "LBM-FT/S^3-F" || Unit1 == "LB-FT/S^3-F")
                        Value = Value * Math.Pow((1000 * UnitConversion.LBMtoKG * UnitConversion.FTtoM / UnitConversion.RtoK), Drct);
                    else if (Unit1 == "LBF/S-F")
                        Value = Value * Math.Pow((1000 * UnitConversion.LBFtoN / UnitConversion.RtoK), Drct);
                    else if (Unit1 == "BTU/H-FT-F")
                        Value = Value * Math.Pow((1000 * UnitConversion.BTUtoW / UnitConversion.HtoS / UnitConversion.FTtoM / UnitConversion.RtoK), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Joule-Thomson Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "JT")
                {
                    if (Unit1 == "K/MPA" || Unit1 == "C/MPA") { }
                    else if (Unit1 == "K/KPA" || Unit1 == "C/KPA")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "K/PA" || Unit1 == "C/PA")
                        Value = Value * Math.Pow(1000000, Drct);
                    else if (Unit1 == "C/ATM")
                        Value = Value / Math.Pow(UnitConversion.ATMtoMPa, Drct);
                    else if (Unit1 == "C/BAR")
                        Value = Value / Math.Pow(UnitConversion.BARtoMPA, Drct);
                    else if (Unit1 == "K/PSI" || Unit1 == "K/PSIA")
                        Value = Value / Math.Pow(UnitConversion.PSIAtoMPA, Drct);
                    else if (Unit1 == "F/PSI" || Unit1 == "F/PSIA" || Unit1 == "R/PSIA")
                        Value = Value / Math.Pow((UnitConversion.PSIAtoMPA / UnitConversion.RtoK), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Length Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "L")
                {
                    if (Unit1 == "METER" || Unit1 == "M") { }
                    else if (Unit1 == "DM")
                        Value = Value / Math.Pow(10, Drct);
                    else if (Unit1 == "CM")
                        Value = Value / Math.Pow(100, Drct);
                    else if (Unit1 == "MM")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "KM")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "INCH" || Unit1 == "IN")
                        Value = Value * Math.Pow(UnitConversion.INtoM, Drct);
                    else if (Unit1 == "FOOT" || Unit1 == "FT")
                        Value = Value * Math.Pow(UnitConversion.FTtoM, Drct);
                    else if (Unit1 == "YARD" || Unit1 == "YD")
                        Value = Value * Math.Pow((UnitConversion.INtoM * 36), Drct);
                    else if (Unit1 == "MILE" || Unit1 == "MI")
                        Value = Value * Math.Pow((UnitConversion.INtoM * 63360), Drct);
                    else if (Unit1 == "LIGHT YEAR")
                        Value = Value * Math.Pow(9.46055E+15, Drct);
                    else if (Unit1 == "ANGSTROM")
                        Value = Value / Math.Pow(10000000000, Drct);
                    else if (Unit1 == "FATHOM")
                        Value = Value * Math.Pow((UnitConversion.FTtoM * 6), Drct);
                    else if (Unit1 == "MIL")
                        Value = Value * Math.Pow((UnitConversion.INtoM / 1000), Drct);
                    else if (Unit1 == "ROD")
                        Value = Value * Math.Pow((UnitConversion.INtoM * 16.5 * 12), Drct);
                    else if (Unit1 == "PARSEC")
                        Value = Value * Math.Pow((30837400000000 * 1000), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }

              //'-----------------------------------------------------------------------
                //'   Area Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "A")
                {
                    if (Unit1 == "METER^2" || Unit1 == "M^2") { }
                    else if (Unit1 == "CM^2")
                        Value = Value / Math.Pow(10000, Drct);
                    else if (Unit1 == "MM^2")
                        Value = Value / Math.Pow(1000000, Drct);
                    else if (Unit1 == "KM^2")
                        Value = Value * Math.Pow(1000000, Drct);
                    else if (Unit1 == "INCH^2" || Unit1 == "IN^2")
                        Value = Value * Math.Pow(Math.Pow(UnitConversion.INtoM, 2), Drct);
                    else if (Unit1 == "FOOT^2" || Unit1 == "FT^2")
                        Value = Value * Math.Pow(Math.Pow(UnitConversion.FTtoM, 2), Drct);
                    else if (Unit1 == "YARD^2" || Unit1 == "YD^2")
                        Value = Value * Math.Pow(Math.Pow((UnitConversion.INtoM * 36), 2), Drct);
                    else if (Unit1 == "MILE^2" || Unit1 == "MI^2")
                        Value = Value * Math.Pow(Math.Pow((UnitConversion.INtoM * 63360), 2), Drct);
                    else if (Unit1 == "ACRE")
                        Value = Value * Math.Pow(Math.Pow((UnitConversion.INtoM * 36), 2) * 4840, Drct);
                    else if (Unit1 == "BARN")
                        Value = Value * Math.Pow(1E-28, Drct);
                    else if (Unit1 == "HECTARE")
                        Value = Value * Math.Pow(10000, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Volume Conversion (Note: not specific volume)
                //'-----------------------------------------------------------------------
                else if (Tpe == "V")
                {
                    if (Unit1 == "METER^3" || Unit1 == "M^3") { }
                    else if (Unit1 == "CM^3")
                        Value = Value / Math.Pow(1000000, Drct);
                    else if (Unit1 == "LITER" || Unit1 == "L" || Unit1 == "DM^3")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "INCH^3" || Unit1 == "IN^3")
                        Value = Value * Math.Pow(UnitConversion.IN3toM3, Drct);
                    else if (Unit1 == "FOOT^3" || Unit1 == "FT^3")
                        Value = Value * Math.Pow((UnitConversion.IN3toM3 * Math.Pow(12, 3)), Drct);
                    else if (Unit1 == "YARD^3" || Unit1 == "YD^3")
                        Value = Value * Math.Pow((UnitConversion.IN3toM3 * Math.Pow(36, 3)), Drct);
                    else if (Unit1 == "GALLON" || Unit1 == "GAL")
                        Value = Value * Math.Pow(UnitConversion.GALLONtoM3, Drct);
                    else if (Unit1 == "QUART" || Unit1 == "QT")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 4), Drct);
                    else if (Unit1 == "PINT" || Unit1 == "PT")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 8), Drct);
                    else if (Unit1 == "CUP")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 16), Drct);
                    else if (Unit1 == "OUNCE")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 128), Drct);
                    else if (Unit1 == "TABLESPOON" || Unit1 == "TBSP")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 256), Drct);
                    else if (Unit1 == "TEASPOON" || Unit1 == "TSP")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 / 768), Drct);
                    else if (Unit1 == "C||D")
                        Value = Value * Math.Pow((UnitConversion.FT3toM3 * 128), Drct);
                    else if (Unit1 == "BARREL")
                        Value = Value * Math.Pow((UnitConversion.GALLONtoM3 * 42), Drct);
                    else if (Unit1 == "BOARD FOOT")
                        Value = Value * Math.Pow((UnitConversion.IN3toM3 * 144), Drct);
                    else if (Unit1 == "BUSHEL")
                        Value = Value * Math.Pow(0.03523907016688, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }

              //'-----------------------------------------------------------------------
                //'   Mass Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "M")
                {
                    if (Unit1 == "KG") { }
                    else if (Unit1 == "G")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "MG")    //        'milligram
                        Value = Value / Math.Pow(1000000, Drct);
                    else if (Unit1 == "LBM" || Unit1 == "LB")
                        Value = Value * Math.Pow(UnitConversion.LBMtoKG, Drct);
                    else if (Unit1 == "GRAIN")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG / 7000), Drct);
                    else if (Unit1 == "SLUG")
                        Value = Value * Math.Pow((UnitConversion.KGFtoN * UnitConversion.LBMtoKG / UnitConversion.FTtoM), Drct);
                    else if (Unit1 == "TON")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG * 2000), Drct);
                    else if (Unit1 == "TONNE")
                        Value = Value * Math.Pow(1000, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }

              //'-----------------------------------------------------------------------
                //'   F||ce Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "F")
                {
                    if (Unit1 == "NEWTON" || Unit1 == "N") { }
                    else if (Unit1 == "MN")    //'milliNewtons
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "KGF")
                        Value = Value * Math.Pow(UnitConversion.KGFtoN, Drct);
                    else if (Unit1 == "DYNE")
                        Value = Value / Math.Pow(100000, Drct);
                    else if (Unit1 == "LBF")
                        Value = Value * Math.Pow(UnitConversion.LBFtoN, Drct);
                    else if (Unit1 == "POUNDAL")
                        Value = Value * Math.Pow((UnitConversion.LBMtoKG * UnitConversion.FTtoM), Drct);
                    else if (Unit1 == "OZF")
                        Value = Value * Math.Pow((UnitConversion.LBFtoN / 16), Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Energy Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "E")
                {
                    if (Unit1 == "JOULE" || Unit1 == "J") { }
                    else if (Unit1 == "KJ")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "MJ")
                        Value = Value * Math.Pow(1000000, Drct);
                    else if (Unit1 == "KW-H")
                        Value = Value * Math.Pow((UnitConversion.HtoS * 1000), Drct);
                    else if (Unit1 == "CAL")
                        Value = Math.Pow(UnitConversion.CALtoJ, Drct) * Value;
                    else if (Unit1 == "KCAL")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ * 1000), Drct);
                    else if (Unit1 == "ERG")
                        Value = Value / Math.Pow(10000000, Drct);
                    else if (Unit1 == "BTU")
                        Value = Value * Math.Pow((UnitConversion.BTUtoKJ * 1000), Drct);
                    else if (Unit1 == "FT-LBF")
                        Value = Value * Math.Pow(UnitConversion.FTLBFtoJ, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                //'-----------------------------------------------------------------------
                //'   Power Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "Q")
                {
                    if (Unit1 == "WATT" || Unit1 == "W") { }
                    else if (Unit1 == "KWATT" || Unit1 == "KW")
                        Value = Value * Math.Pow(1000, Drct);
                    else if (Unit1 == "BTU/S")
                        Value = Value * Math.Pow(UnitConversion.BTUtoW, Drct);
                    else if (Unit1 == "BTU/MIN")
                        Value = Value * Math.Pow((UnitConversion.BTUtoW / 60), Drct);
                    else if (Unit1 == "BTU/H")
                        Value = Value * Math.Pow((UnitConversion.BTUtoW / UnitConversion.HtoS), Drct);
                    else if (Unit1 == "CAL/S")
                        Value = Value * Math.Pow(UnitConversion.CALtoJ, Drct);
                    else if (Unit1 == "KCAL/S")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ * 1000), Drct);
                    else if (Unit1 == "CAL/MIN")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ / 60), Drct);
                    else if (Unit1 == "KCAL/MIN")
                        Value = Value * Math.Pow((UnitConversion.CALtoJ / 60 * 1000), Drct);
                    else if (Unit1 == "FT-LBF/S")
                        Value = Value * Math.Pow(UnitConversion.FTLBFtoJ, Drct);
                    else if (Unit1 == "FT-LBF/MIN")
                        Value = Value * Math.Pow((UnitConversion.FTLBFtoJ / 60), Drct);
                    else if (Unit1 == "FT-LBF/H")
                        Value = Value * Math.Pow((UnitConversion.FTLBFtoJ / UnitConversion.HtoS), Drct);
                    else if (Unit1 == "HP")
                        Value = Value * Math.Pow(UnitConversion.HPtoW, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");

                }
                //'-----------------------------------------------------------------------
                //'   Surface Tension Conversion
                //'-----------------------------------------------------------------------
                else if (Tpe == "N")
                {
                    if (Unit1 == "N/M") { }
                    else if (Unit1 == "MN/M")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "DYNE/CM" || Unit1 == "DYN/CM")
                        Value = Value / Math.Pow(1000, Drct);
                    else if (Unit1 == "LBF/FT")
                        Value = Value * Math.Pow(UnitConversion.LBFTtoNM, Drct);
                    else
                        RefpropErrorHandler.ErrorHandler(this, "Undefined input unit");
                }
                Unit1 = Unit2;
                OldUnits = NewUnits;  //新加
            }
            return Value;
        }
        public override string ToString()
        {
            return //this.ThermoProperty.ToString() + " : " + 
               this.Value.ToString("F4") + " "
                + "[" + this.Unit + "]";
        }
        public void OutPut()
        {
            Console.WriteLine(this.ToString());
        }

        public static implicit operator ThermodynamicParameter(double value)
        {
            return new ThermodynamicParameter(value);
        }

        //public static implicit operator ThermodynamicParameter(string d)
        //{
        //    string[] s = d.Split(new char[] { ' ' });
        //    string unit = s[1];
        //    double value = s[0];
        //    return new ThermodynamicParameter(d);
        //}

        public static implicit operator double(ThermodynamicParameter tp)
        {
            return tp.Value;
        }
    }

    public class Temperature : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Temperature(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Temperature (UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.T, 0.0, unitSystem,List) { }
        public Temperature(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.T, value, unitSystem) { }

        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Temperature(double value)
        {
            return new Temperature(value,UnitSystems.SI);
        }
    }

    public class Pressure : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Pressure(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Pressure(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.P, 0.0, unitSystem,List) { }
        public Pressure(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.P, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Pressure(double value)
        {
            return new Pressure(value, UnitSystems.SI);
        }
    }
    public class Density : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Density(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Density(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.D, 0.0, unitSystem,List) { }
        public Density(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.D, value, unitSystem) { }
        public Density(double _wm, double value, UnitSystems unitSystem)
            : base(_wm, UnitTypes.D, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Density(double value)
        {
            return new Density(value, UnitSystems.SI);
        }
    }

    public class Enthalpy : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Enthalpy(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Enthalpy(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.H, 0.0, unitSystem,List) { }
        public Enthalpy(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.H, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Enthalpy(double value)
        {
            return new Enthalpy(value, UnitSystems.SI);
        }
    }

    public class Entropy : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Entropy(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Entropy(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.S, 0.0, unitSystem,List) { }
        public Entropy(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.S, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Entropy(double value)
        {
            return new Entropy(value, UnitSystems.SI);
        }
    }

    public class SpecificHeat : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public SpecificHeat(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public SpecificHeat(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.S, 0.0, unitSystem,List) { }
        public SpecificHeat(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.S, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator SpecificHeat(double value)
        {
            return new SpecificHeat(value, UnitSystems.SI);
        }
    }

    public class Viscosity : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Viscosity(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Viscosity(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.U, 0.0, unitSystem,List) { }
        public Viscosity(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.U, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Viscosity(double value)
        {
            return new Viscosity(value, UnitSystems.SI);
        }
    }

    public class ThermalConductivity : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public ThermalConductivity(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public ThermalConductivity(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.K, 0.0, unitSystem,List) { }
        public ThermalConductivity(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.K, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator ThermalConductivity(double value)
        {
            return new ThermalConductivity(value, UnitSystems.SI);
        }
    }

    public class SurfaceTension : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public SurfaceTension(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public SurfaceTension(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.N, 0.0, unitSystem,List) { }
        public SurfaceTension(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.N, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator SurfaceTension(double value)
        {
            return new SurfaceTension(value, UnitSystems.SI);
        }
    }
    public class PrantdlNumber : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public PrantdlNumber(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public PrantdlNumber(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.X, 0.0, unitSystem,List) { }
        public PrantdlNumber(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.X, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator PrantdlNumber(double value)
        {
            return new PrantdlNumber(value, UnitSystems.SI);
        }
    }
   public class Quality : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public Quality(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public Quality(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.X, 0.0, unitSystem,List) { }
        public Quality(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.X, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator Quality(double value)
        {
            return new Quality(value, UnitSystems.SI);
        }
    }


   public class CompressibilityFactor : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public CompressibilityFactor(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public CompressibilityFactor(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.X, 0.0, unitSystem,List) { }
        public CompressibilityFactor(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.X, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator CompressibilityFactor(double value)
        {
            return new CompressibilityFactor(value, UnitSystems.SI);
        }
    }
   public class InternalEnergy : ThermodynamicParameter
    {
        protected new UnitTypes UnitType { private set; get; }
        public InternalEnergy(UnitSystems unitSystem)
            : this(0.0, unitSystem) { }
        public InternalEnergy(UnitSystems unitSystem, List<ThermodynamicParameter> List)
            : base(0.0, UnitTypes.H, 0.0, unitSystem,List) { }
        public InternalEnergy(double value, UnitSystems unitSystem)
            : base(0.0, UnitTypes.H, value, unitSystem) { }
        /// <summary>
        /// value with SI units
        /// </summary>
        /// <param name="value">value</param>
        /// <returns></returns>
        public static implicit operator InternalEnergy(double value)
        {
            return new InternalEnergy(value, UnitSystems.SI);
        }
    }

}