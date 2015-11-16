using System;
//using System.Collections.Generic;
//using System.Linq;
using System.Text;

namespace RefpropCSNET
{
    public sealed partial class Refrigerant
    {       
        public ThermodynamicParameter Func_T(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.T;
        }
        public  Temperature Func_T3(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.T as Temperature;
        }

        public  ThermodynamicParameter Func_P(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.P;
        }

        public  ThermodynamicParameter Func_H(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.H;
        }

        public  ThermodynamicParameter Func_S(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.S;
        }
        
        public  ThermodynamicParameter Func_D(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.D;
        }
        public  ThermodynamicParameter Func_CP(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.CP;
        }
        public  ThermodynamicParameter Func_CV(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.CV;
        }
        public  ThermodynamicParameter Func_E(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.E;
        }
        
        public  ThermodynamicParameter Func_X(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.X;
        }

        public  ThermodynamicParameter Func_U(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.U;
        }

        public ThermodynamicParameter Func_K(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = fluidName;
            this.FindState(InpCode, units, Prop1, Prop2);
            return this.K;
        }

        public  void WriteErrMsg(string ErrMsg)
        {
            Console.WriteLine(ErrMsg);
            Console.WriteLine("== Please Press Any Key to Exit ==");
            Console.ReadKey();
            Environment.Exit(0);
        }

        public double Func_T(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_T(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_P(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_P(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_H(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_H(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_S(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_S(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_D(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_D(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_CP(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CP(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_CV(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CV(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_X(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_X(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_IE(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_E(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_U(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_U(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_K(string fluidName, string InpCode, double Prop1, double Prop2){ return Func_K(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }


        public double Func_T(string InpCode, double Prop1, double Prop2) { return Func_T(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_P(string InpCode, double Prop1, double Prop2) { return Func_P(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_H(string InpCode, double Prop1, double Prop2) { return Func_H(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_S(string InpCode, double Prop1, double Prop2) { return Func_S(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_D(string InpCode, double Prop1, double Prop2) { return Func_D(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_CP(string InpCode, double Prop1, double Prop2) { return Func_CP(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_CV(string InpCode, double Prop1, double Prop2) { return Func_CV(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_X(string InpCode, double Prop1, double Prop2) { return Func_X(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_IE(string InpCode, double Prop1, double Prop2) { return Func_E(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_U(string InpCode, double Prop1, double Prop2) { return Func_U(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        public double Func_K(string InpCode, double Prop1, double Prop2) { return Func_K(Name, InpCode, CurrentUnitSystem, Prop1, Prop2); }


        public double Func_T(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_T(Name, InpCode, units, Prop1, Prop2); }
        public double Func_P(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_P(Name, InpCode, units, Prop1, Prop2); }
        public double Func_H(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_H(Name, InpCode, units, Prop1, Prop2); }
        public double Func_S(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_S(Name, InpCode, units, Prop1, Prop2); }
        public double Func_D(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_D(Name, InpCode, units, Prop1, Prop2); }
        public double Func_CP(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_CP(Name, InpCode, units, Prop1, Prop2); }
        public double Func_CV(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_CV(Name, InpCode, units, Prop1, Prop2); }
        public double Func_X(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_X(Name, InpCode, units, Prop1, Prop2); }
        public double Func_IE(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_E(Name, InpCode, units, Prop1, Prop2); }
        public double Func_U(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_U(Name, InpCode, units, Prop1, Prop2); }
        public double Func_K(string InpCode, UnitSystems units, double Prop1, double Prop2) { return Func_K(Name, InpCode, units, Prop1, Prop2); }

        //public double Temperature(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_T(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Pressure(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_P(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Enthalpy(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_H(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Entropy(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_S(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Density(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_D(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double IsobaricHeatCapacity(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CP(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double IsochoricHeatCapacity(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CV(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Quality(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_X(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double InternalEnergy(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_E(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public double Viscosity(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_U(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Temperature2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_T(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Pressure2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_P(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Enthalpy2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_H(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Entropy2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_S(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Density2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_D(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter IsobaricHeatCapacity2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CP(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter IsochoricHeatCapacity2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_CV(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Quality2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_X(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter InternalEnergy2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_E(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }
        //public ThermodynamicParameter Viscosity2(string fluidName, string InpCode, double Prop1, double Prop2) { return Func_U(fluidName, InpCode, CurrentUnitSystem, Prop1, Prop2); }

    }
}
