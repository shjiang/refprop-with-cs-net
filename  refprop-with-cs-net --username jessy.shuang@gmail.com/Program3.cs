using System;
//using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace RefpropCSNET
{
    partial class Program3
    {
        static void Main(string[] args)
        {
            //Refrigerant REFPROP = Refrigerant.GetInstance("R22=.85,ISOBUTAN=.15", UnitSystems.SI);
            //Refrigerant REFPROP = Refrigerant.GetInstance("R410A", UnitSystems.SI);
            // REFPROP.Reference = ReferenceState.DEF;//ASH
            Refrigerant REFPROP = Refrigerant.GetInstance("R717");//R22=.85,ISOBUTAN=.15
            REFPROP.Reference = ReferenceState.DEF;//ASH
            Console.WriteLine(REFPROP.Func_P("TQ", 197.6833, 1.0));
            REFPROP.Name = "R32";
            Console.WriteLine(REFPROP.Func_H("PS", 0.70905, 2.34416));
            REFPROP.Name = "R22";
            Console.WriteLine(REFPROP.Func_H("PS", 0.70905, 2.34416));
            REFPROP.Name = "R32";
            Console.WriteLine(REFPROP.Func_H("PS", 0.70905, 2.34416));
            Pause();
            REFPROP.DisplayInformation();
            Console.WriteLine(REFPROP.CriticalTemperature.ToString());
            Console.WriteLine(REFPROP.CriticalPressure.ToString());
            Console.WriteLine(REFPROP.CriticalDensity.ToString());
            Pause();
            //Console.WriteLine(REFPROP.Func_P("TQ", 10 + 273.15, 1.0));
            //Pause();
            for (int i = 210; i < 370; i += 10)
            {
                Console.WriteLine(i.ToString() + ":" + REFPROP.Func_P("TQ", i, 0.5).ToString());
                //Console.WriteLine(REFPROP.Func_CP("PT", 3 * REFPROP.P, i * Math.Pow(3, (0.2 / 1.2))));
            }
            Pause();

            var myPressure = Pressure2("CO2", "TQ", UnitSystems.ME, -8, 0); //or Pressure("CO2", "TQ", UnitSystem.SIwithC, -8, 0)
            Console.WriteLine(myPressure);
            Console.WriteLine(Temperature2("CO2", "PS", UnitSystems.ME, 39.5, 2.093).ToString());
            Console.WriteLine(Entropy2("CO2", "PT", UnitSystems.ME, 15.3, -18).ToString());
            Console.WriteLine(Entropy2("CO2", "PH", UnitSystems.ME, 39.5, 511.62).ToString());
            Console.WriteLine(Quality2("CO2", "PH", UnitSystems.ME, 39.5, 296.42).ToString());
            Console.WriteLine(Enthalpy2("CO2", "PT", UnitSystems.ME, 15.3, -18).ToString());
            Console.WriteLine(Density2("CO2", "PT", UnitSystems.ME, 15.3, -15).ToString());
            Console.WriteLine(Viscosity2("CO2", "PT", UnitSystems.ME, 15.3, -15).ToString());

            Console.ReadKey();


            REFPROP.Name = "R410A";
            REFPROP.Reference = ReferenceState.DEF;//ASH
            for (int i = 210; i < 345; i += 5)
            {
                Console.WriteLine(REFPROP.Func_P("TQ", i, 1.0));
            }
            Pause("R410A");
            REFPROP.Name = "R404A";
            for (int i = 210; i < 345; i += 5)
            {
                Console.WriteLine(REFPROP.Func_P("TQ", i, 1.0));
            }
            Pause("R404A");
            REFPROP.Name = "R410A";

            for (int i = 210; i < 345; i += 5)
            {
                Console.WriteLine(REFPROP.Func_P("TQ", i, 1.0));
            }
            Pause("R410A");


            Timing tObj = new Timing();
            tObj.StartTime();
            //CalPm();
            //REFPROP.DisableCalThermalQuantities();
            //for (int i = 0; i < 100000; i++) put(REFPROP.Func_D("R 4 10a ", "TP", UnitSystems.SI, 288.15, 0.1010));
            tObj.StopTime();
            Console.WriteLine("Running time (.NET): {0} Seconds", tObj.Result().TotalSeconds);
            Pause(0);

            Temperature T0 = 273.0;
            T0.OutPut();
            double d0 = T0;
            Console.WriteLine(d0);
            Pause(0);
            #region test code
            Refrigerant RTest = Refrigerant.GetInstance("nitrogen", ReferenceState.DEF, UnitSystems.cgs);
            RTest.FindState("TP", 100, 1.000);
            RTest.DisplayThermoDynamicState();
            Pause(1);

            Temperature T1 = new Temperature(500, UnitSystems.SI);
            T1.OutPut();
            T1.UnitSystem = UnitSystems.SIwithC;
            T1.OutPut();
            //T.UnitType = UnitTypes.P; //不应该让子类能修改
            T1.UnitSystem = UnitSystems.E;
            T1.OutPut();

            T1 = REFPROP.Func_T3("CO2", "PQ", UnitSystems.ME, 39.5, 0.6);
            REFPROP.DisplayThermoDynamicState();
            Pause(11);
            double satT, satP;
            for (int i = 0; i < 11; i++)
            {
                satP = 0.101325 + i * 0.4;
                satT = REFPROP.Func_T("R410a", "PQ", UnitSystems.SIwithC, satP, 0);
                Console.WriteLine("satT:{0},satP:{1},u:{2}", satT, satP, REFPROP.Func_U("R410a", "PT", UnitSystems.SIwithC, satP, satT - 3));
            }
            //REFPROP.D.CheckUnit("kg/cm^3");
            Pause(2);

            //Refrigerant r0 = new Refrigerant();
            //r0.UserUnitSystem = UnitSystem.SI;
            //method to evaluate the calculation time
            //delete the corresponding codes if you didn't need it.
            Pressure P1 = new Pressure(UnitSystems.E);
            P1.Value = 3626;
            P1.OutPut();
            Console.WriteLine(P1.UnitConvert(1.0, "P", "Psig", "pa", 0.0));
            P1.UnitSystem = UnitSystems.SI;
            P1.OutPut();
            Pause(3);

            var myPressure1 = REFPROP.Func_P("CO2", "TQ", UnitSystems.ME, -8, 0); //or Pressure("CO2", "TQ", UnitSystem.SIwithC, -8, 0)
            Console.WriteLine(myPressure1);
            Console.WriteLine(REFPROP.Func_T("CO2", "PS", UnitSystems.ME, 39.5, 2.093));
            Console.WriteLine(REFPROP.Func_S("CO2", "PT", UnitSystems.ME, 15.3, -18));
            Console.WriteLine(REFPROP.Func_S("CO2", "PH", UnitSystems.ME, 39.5, 511.62));
            Console.WriteLine(REFPROP.Func_X("CO2", "PH", UnitSystems.ME, 39.5, 296.42));
            Console.WriteLine(REFPROP.Func_H("CO2", "PT", UnitSystems.ME, 15.3, -18));
            Console.WriteLine(REFPROP.Func_D("CO2", "PT", UnitSystems.ME, 15.3, -15));
            Pause(4);

            //the refrigerant name can be either in upper case or lower case, with or without spaces, with or without suffix
            //names of predefined mixture: R410A,R502
            //names of pure refrigerant: CO2,nitrogen
            //names of new defined mixture:R125=0.50,R32=0.50. Here, R125 and R32 are compositions of the mixutre, 
            //and the the following number is the fraction of the composition, with Mass fraction as default.
            //if we want to use mole fraction, we can change it with "REFPROP.IsMassQ = false;". If you like another way, please let me know.

            //First way to use the code
            Console.WriteLine("\n===First way to use the code===");
            double h = REFPROP.Func_H("R 4 10a.M IX ", "TP", UnitSystems.Refprop, 288.15, 101.0); //The Enthalpy function with "TP"
            Console.WriteLine(h.ToString());
            h = REFPROP.Func_H("R 4 10a.M IX ", "PT", UnitSystems.Refprop, 101.0, 288.15); //with "PT"
            Console.WriteLine(h.ToString());
            Pause(5);


            //The second way to use the code
            //in this way, we can output the units
            Console.WriteLine("\n===The second way to use the code===");
            ThermodynamicParameter d = REFPROP.Func_D("R 4 10a ", "TP", UnitSystems.SI, 288.15, 0.1010);
            d.OutPut();  //both value and unit are output 
            Console.WriteLine(REFPROP.Func_D("R 4 10a ", "TP", UnitSystems.SI, 288.15, 0.1010));  //with the same output
            Console.WriteLine("value:" + d.Value);  //get its value
            Console.WriteLine("uint:" + d.Unit);  //get its unit
            Pause(6);


            //Another way to use the code
            //It is more efficient to get more thermodynamic parameters at one point.
            //For example, if we determine a state of the refrigerant with Pressure and Temperature,
            //all the other parameters can be obtained at the same time.
            Console.WriteLine("\n===Another way to use the code===");
            Refrigerant r = Refrigerant.GetInstance("R125=.50,R32=.50", UnitSystems.Refprop);  //first define a refrigerant
            //also, it is not necessary create "r" here. we can use "RERPROP" instead of "r" 
            //REFPROP.Name="R125=50,R32=50"; REFPROP.UserUnitSystem=UnitSystem.Refprop;
            r.FindState("TP", 288.15, 0.1010); //determine its state with two parameters
            r.H.OutPut();   //output enthalpy
            r.P.OutPut();   //output pressure
            r.T.OutPut();   //output temperaute
            r.X.OutPut();   //output quality
            r.S.OutPut();   //output entropy
            r.CP.OutPut();   //output CP
            r.CV.OutPut();   //output CV
            r.E.OutPut();   //output INTERNAL ENERGY
            Console.WriteLine("==Refrigerant Name is changed to R22==");
            r.Name = "R22"; //Change its Name, i.e., to another refrigerant            
            Console.WriteLine(r.MolecularWeight); //show the changes
            r.FindState("TP", 288.15, 0.1010);
            r.DisplayThermoDynamicState(UnitSystems.SI);
            Console.WriteLine("\n===give H another unit system===");
            r.H.UnitSystem = UnitSystems.E;  //give it another unit system
            r.H.OutPut(); //output with another unit
            Console.WriteLine("\n===display all the parameters with specified unit system===");
            r.DisplayThermoDynamicState(UnitSystems.SIwithC);  //display all the parameters with specified unit system
            Console.WriteLine("\n===give all another unit system===");
            r.CurrentUnitSystem = UnitSystems.SIwithC;  //give all another unit system
            r.H.OutPut();   //output enthalpy
            r.P.OutPut();   //output pressure
            r.T.OutPut();   //output temperaute
            r.X.OutPut();   //output quality
            r.S.OutPut();   //output entropy
            r.CP.OutPut();   //output CP
            r.CV.OutPut();   //output CV
            r.E.OutPut();   //output INTERNAL ENERGY
            Pause(7);

            //REFPROP.IsMassQ = false;
            ////Console.WriteLine(Quality2("R125=0.3024,R32=0.6976", "TD", UnitSystem.SI, 300, 200));
            //Console.WriteLine(Quality2("R125=0.3024,R32=0.6976", "TD", UnitSystem.SI, 300, 200));
            //REFPROP.DisplayInformation();
            //REFPROP.IsMassQ = true;
            //Console.WriteLine(Quality2("R125=50,R32=50", "TD", UnitSystem.SI, 300, 200));
            //REFPROP.DisplayInformation();
            //Console.ReadKey();    

            //The test examples are from the file "REFPROP.XLS"
            //Pure fluid calculations
            Console.WriteLine("\n=== ure fluid calculations ===");
            REFPROP.Name = "nitrogen";
            REFPROP.CurrentUnitSystem = UnitSystems.MolarSI;
            REFPROP.FindState("TP", UnitSystems.MolarSI, 100.0, 2.0);
            REFPROP.DisplayThermoDynamicState(UnitSystems.MolarSI);
            REFPROP.FindState("PS", UnitSystems.MolarSI, 1, 100);
            REFPROP.DisplayThermoDynamicState(UnitSystems.MolarSI);
            REFPROP.FindState("PH", UnitSystems.MolarSI, 1, -4000);
            REFPROP.DisplayThermoDynamicState(UnitSystems.MolarSI);
            Console.WriteLine(REFPROP.Func_P("nitrogen", "TH", UnitSystems.SI, 100, -40));
            Console.WriteLine(REFPROP.Func_P("nitrogen", "TH<", UnitSystems.SI, 100, -40));
            Console.WriteLine(REFPROP.Func_X("nitrogen", "TH", UnitSystems.SI, 100, -40));
            Console.WriteLine(REFPROP.Func_X("nitrogen", "TH<", UnitSystems.SI, 100, -40));
            Pause();

            //Two-Phase Calculations
            Console.WriteLine("\n=== Two-Phase Calculations ===");
            REFPROP.FindState("TD", UnitSystems.MolarSI, 100, 15);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("PQ", UnitSystems.MolarSI, 0.3, 0.0);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("TQ", UnitSystems.MolarSI, 100, 0.0);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("TQ", UnitSystems.MolarSI, 100, 1.0);
            REFPROP.DisplayThermoDynamicState();
            Pause();


            //!Calculate properties given only the saturation temperature or saturation pressure are not implemented.
            //If you want these functions, please show me the "function form" you needed.

            //Mixture Calculations (maximum number of components in a mixture is 20)
            Console.WriteLine("\n=== Mixture Calculations ===");
            Console.WriteLine(REFPROP.Func_D("nitrogen=0.7557,argon=0.0127,oxygen=0.2316", "TP", UnitSystems.SI, 100, 0.1));
            Console.WriteLine(REFPROP.Func_D("nitrogen=0.7557,argon=0.0127,oxygen=0.2316", "TP", UnitSystems.SI, 100, 1));
            Console.WriteLine("Name:{0}", REFPROP.Name);
            Console.WriteLine("MolecularWeight:{0}", REFPROP.MolecularWeight);
            Pause();


            //Predefined Mixture Calculations
            Console.WriteLine("\n=== Predefined Mixture Calculations===");
            REFPROP.Name = "R410A.mix";
            Console.WriteLine("Name:{0}", REFPROP.Name);
            Console.WriteLine("MolecularWeight:{0}", REFPROP.MolecularWeight);
            REFPROP.FindState("TP", UnitSystems.SI, 300, 10);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("TP", UnitSystems.E, 25, 1000);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("TD", UnitSystems.SI, 300, 200);
            REFPROP.DisplayThermoDynamicState();
            REFPROP.FindState("TQ", UnitSystems.SI, 300, 1);
            REFPROP.DisplayThermoDynamicState();
            Pause();

            #endregion
        }
        static void Pause()
        {
            Pause("");
        }
        static void Pause(string i)
        {
            Console.WriteLine(i + "==> " + "Press SPACE to QUIT, other keys to continue <==");
            if (Console.ReadKey().KeyChar.Equals(' ')) System.Environment.Exit(0);
        }
        static void Pause(int i)
        {
            Pause(i.ToString());
        }
        static void put(Object o)
        { Console.WriteLine(o.ToString()); }


        //计算在指定过冷度的情况下，保证制冷剂不汽化而允许的最大压降。
        static void CalMaxPDrop(Refrigerant r, double deltaT)
        {
            r.Name = "R410A";
            r.CurrentUnitSystem = UnitSystems.SIwithC;
            double satP, Tsub, deltaP;
            for (int i = 0; i < 60; i++)
            {
                satP = r.Func_P("TQ", i, 0);
                Tsub = i - deltaT;
                deltaP = satP - r.Func_P("TQ", Tsub, 0);
                Console.WriteLine("{0}", deltaP);
            }
        }


        readonly static Refrigerant REFPROP = Refrigerant.GetInstance();

        public static double Temperature(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.T.Value;
        }
        public static ThermodynamicParameter Temperature2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            REFPROP.Name = fluidName;
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.T;
        }
        public static double Pressure(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.P.Value;
        }
        public static ThermodynamicParameter Pressure2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.P;
        }
        public static double Enthalpy(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.H.Value;
        }
        public static ThermodynamicParameter Enthalpy2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.H;
        }
        public static double Entropy(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.S.Value;
        }
        public static ThermodynamicParameter Entropy2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.S;
        }

        public static double Density(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.D.Value;
        }
        public static ThermodynamicParameter Density2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.D;
        }
        public static double CP(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.CP.Value;
        }
        public static ThermodynamicParameter CP2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.CP;
        }
        public static double CV(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.CV.Value;
        }
        public static ThermodynamicParameter CV2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.CV;
        }

        public static double InternalEnergy(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.E.Value;
        }
        public static ThermodynamicParameter InternalEnergy2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.E;
        }

        public static double Quality(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.X.Value;
        }
        public static ThermodynamicParameter Quality2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.X;
        }

        public static ThermodynamicParameter Viscosity2(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.U;
        }

        public static double Viscosity(string fluidName, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            //Refrigerant r = new Refrigerant(fluidName);
            REFPROP.Name = fluidName;
            REFPROP.FindState(InpCode, units, Prop1, Prop2);
            return REFPROP.U.Value;
        }

        public static void WriteErrMsg(string ErrMsg)
        {
            Console.WriteLine(ErrMsg);
            Console.WriteLine("== Please Press Any Key to Exit ==");
            Console.ReadKey();
            Environment.Exit(0);
        }




    }

    public class Timing
    {
        TimeSpan startingTime;
        TimeSpan duration;
        public Timing()
        {
            startingTime = new TimeSpan(0);
            duration = new TimeSpan(0);
        }
        public void StopTime()
        {
            duration =
                //Process.GetCurrentProcess().Threads[0].UserProcessorTime.Subtract(startingTime);
            Process.GetCurrentProcess().UserProcessorTime.Subtract(startingTime);
        }
        public void StartTime()
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            //startingTime = Process.GetCurrentProcess().Threads[0].UserProcessorTime;
            startingTime = Process.GetCurrentProcess().UserProcessorTime;
        }
        public TimeSpan Result()
        {
            return duration;
        }
    }
}
