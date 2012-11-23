using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace sc.net
{
    class Refrigerant : Refprop
    {
        //public RefrigerantCategory Category { get; set; }
        public string Name { get; set; }
        
        /// <summary>
        /// Molecular Weight: g/mol
        /// </summary>
        public double MolecularWeight { get { return base.wm; } }

        public double CriticalTemperature { get { return base.tc; } }
        
        /// <summary>
        /// Critical Pressure: kPa
        /// </summary>
        public double CriticalPressure { get { return base.pc; } }
       
        /// <summary>
        /// Critical Density: g/L
        /// </summary>
        public double CriticalDensity { get { return base.Dc * MolecularWeight; } }

        public new double[]  MoleFractions { get { return base.MoleFractions; } }

        public new double[]  MassFractions { get { return base.MassFractions; } }

        /// <summary>
        /// Number of Components 
        /// </summary>
        public Int32 NumberOfComponents { get { return base.nc; } }

        Comonent[] components;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="refrigerantCategory">purefluid(e.g. R22),predefined mixture(e.g. R410A)</param>
        /// <param name="refrigerantName">e.g. R22</param>

        /// <summary>
        /// temperature of the current refrigerant : ℃
        /// </summary>
        public double Temperature { get { return base.temperature; } }
        /// <summary>
        ///  pressure of the current refrigerant : kPa
        /// </summary>
        public double Pressure { get { return base.pressure; } }
        /// <summary>
        /// overall molar density of the current refrigerant : mol/L * g/mol
        /// </summary>
        public double Density { get { return base.density * MolecularWeight; } }
        /// <summary>
        /// ??????
        /// </summary>
        public double VaporQuality { get { return base.quality; } }
        /// <summary>
        /// Density of the current refrigerant : J/mol
        /// </summary>
        public double InternalEnergy { get { return base.internalenergy; } }
        /// <summary>
        /// Enthalpy of the current refrigerant : kJ/kg
        /// </summary>
        public double Enthalpy { get { return base.enthalpy / MolecularWeight; } }
        /// <summary>
        /// Entropy of the current refrigerant : kJ/(kg-K)
        /// </summary>
        public double Entropy { get { return base.entropy / MolecularWeight; } }
        /// <summary>
        /// isochoric (constant V) heat capacity [J/mol-K]
        /// </summary>
        public double Cv { get { return base.cv; } }
        /// <summary>
        /// isobaric (constant p) heat capacity [J/mol-K]
        /// </summary>
        public double Cp { get { return base.cp; } }
        /// <summary>
        /// isenthalpic Joule-Thompson coefficient [K/kPa]
        /// </summary>
        public double HJT { get { return base.hjt; } }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="refrigerantCategory">PureFluid, PredefinedMixture, NewMixture, PseudoPureFluid</param>
        /// <param name="refrigerantName">e.g. :  R22,or R410A,or "R22=0.5,R134a=0.5"(0.5 is the mole composition</param>
        public Refrigerant(RefrigerantCategory refrigerantCategory, string refrigerantName,ReferenceState reference)
        {
            //this.Category = refrigerantCategory;

            this.Name = refrigerantName;
            if (refrigerantCategory == RefrigerantCategory.PureFluid)
            {
                base.LoadFluid(Name + ".fld", 1, new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },reference.ToString());
            }
            else if (refrigerantCategory == RefrigerantCategory.PredefinedMixture)
            {
                base.LoadPredefinedFluid(Name + ".mix",reference.ToString());
            }
            else if (refrigerantCategory == RefrigerantCategory.NewMixture)
            {
                string[] a = refrigerantName.Split(',');
                double[] b = new double[20];
                string c = "";
                for (int i = 0; i < a.Length; i++)
                {
                    c += a[i].Split('=')[0];
                    b[i] = Convert.ToDouble((a[i].Split('='))[1]);
                    c += ".fld|";
                }
                base.LoadFluid(c, a.Length, b,reference.ToString());
            }
            else
            {

            }
            if (NumberOfComponents < 1) throw new Exception("Number of Components is less than 1");
            components = new Comonent[NumberOfComponents];
            for (int i = 0; i < NumberOfComponents; i++)
            {
                double wm = 0.0, ttp = 0.0, tnbp = 0.0, tc = 0.0, pc = 0.0, Dc = 0.0, Zc = 0.0, acf = 0.0, dip = 0.0, Rgas = 0.0;
                base.GetComponentInfo(i + 1, ref wm, ref ttp, ref tnbp, ref tc, ref pc, ref Dc, ref Zc, ref acf, ref dip, ref Rgas);
                components[i] = new Comonent();
                components[i].MolecularWeight = wm;
                components[i].TriplePointTemperature = ttp;
                components[i].NormalBoilingPointTemperature = tnbp;
                components[i].CriticalTemperature = tc;
                components[i].CriticalPressure = pc;
                components[i].CriticalDensity = Dc;
                components[i].Compressibility = Zc;
                components[i].AccentricFactor = acf;
                components[i].DipoleMoment = dip;
                components[i].GasConstant = Rgas;
                components[i].Mass_Fraction = MassFractions[i];
                components[i].MoleFraction = MoleFractions[i];
                if (NumberOfComponents == 1) components[i].Name = refrigerantName;
                else if(refrigerantCategory == RefrigerantCategory.PredefinedMixture)
                    components[i].Name = base.hfiles.ToString().Split('|')[i];
                else if(refrigerantCategory == RefrigerantCategory.NewMixture)
                    components[i].Name = refrigerantName.Split(',', '=')[i * 2];
            }
        }

        /// <summary>
        /// Calculate the saturated state given the staturated temperature
        /// </summary>
        /// <param name="temperature">staturated temperature</param>
        /// <param name="PhaseFlag">saturated liquid ? saturated vapor</param>
        public void FindSaturatedStateWithTemperature(double temperature, SaturationPoint PhaseFlag)
        {
            base.GetSaturatedPressure(temperature, (Int32)PhaseFlag);
        }

        /// <summary>
        /// Calculate thermodynamic state given two independent variables(temperature and pressure) plus composition
        /// </summary>
        /// <param name="T">temperture of the refrigerant: k</param>
        /// <param name="P">pressure of the refrigerant: Kpa</param>
        public void FindStateWithTP(double T, double P)
        {
            base.TPFLSH(T, P);
        }
        /// <summary>
        /// Calculate thermodynamic state given two independent variables(temperature and density) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="D"></param>
        public void FindStateWithTD(double T, double D)
        {
            base.TDFLSH(T, D);
        }
        /// <summary>
        ///Calculate thermodynamic state given two independent variables(temperature and enthalpy) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="H"></param>
        /// <param name="root"></param>
        public void FindStateWithTH(double T, double H, int root)
        {
            base.THFLSH(T, H, root);
        }
        /// <summary>
        /// Calculate thermodynamic state given two independent variables(temperature and energy) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="E"></param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        public void FindStateWithTE(double T, double E, int root)
        {
            base.TEFLSH(T, E, root);
        }

        public void FindStatueWithPD(double P, double D)
        {
            base.PDFLSH(P, D);
        }

        //      subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine PHFLSH (p,h,z,t,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        public void FindStatueWithPH(double P, double H)
        {
            base.PHFLSH(P, H);
        }
        //      subroutine PSFLSH (p,s,z,t,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        public void FindStatueWithPS(double P, double S)
        {
            base.PSFLSH(P, S);
        }
        //      subroutine PEFLSH (p,e,z,t,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
        public void FindStatueWithPE(double P, double E)
        {
            base.PEFLSH(P,E);
        }
        //      subroutine HSFLSH (h,s,z,t,p,D,Dl,Dv,x,y,q,e,cv,cp,w,ierr,herr)
        public void FindStatueWithHS(double H, double S)
        {
            base.HSFLSH(H, S);
        }

        //      subroutine ESFLSH (e,s,z,t,p,D,Dl,Dv,x,y,q,h,cv,cp,w,ierr,herr)
        public void FindStatueWithES(double E, double S)
        {
            base.ESFLSH(E, S);
        }
        //      subroutine DHFLSH (D,h,z,t,p,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        public void FindStatueWithDH(double D, double H)
        {
            base.DHFLSH(D, H);
        }
        //      subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        public void FindStatueWithDS(double D, double S)
        {
            base.DSFLSH(D, S);
        }
        //      subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)

        public void FindStatueWithDE(double D, double E)
        {
            base.DEFLSH(D, E);
        }
        //      subroutine TQFLSH (t,q,z,kq,p,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        public void FindStatueWithTQ(double T, double Q,int kq)
        {
            base.TQFLSH(T,Q,kq);
        }
        //      subroutine PQFLSH (p,q,z,kq,t,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        public void FindStatueWithPQ(double P, double Q, int kq)
        {
            base.PQFLSH(P, Q, kq);
        }

        public void DisplayThermoDynamicState()
        {
            Console.WriteLine("===================Thermodynamic State========================");
            Console.WriteLine("Temperature            :" + this.Temperature.ToString());
            Console.WriteLine("Pressure               :" + this.Pressure.ToString());
            Console.WriteLine("Enthalpy               :" + this.Enthalpy.ToString());
            Console.WriteLine("Entropy                :" + this.Entropy.ToString());
            Console.WriteLine("Density                :" + this.Density.ToString());
            Console.WriteLine("VaporQuality           :" + this.VaporQuality.ToString());
            Console.WriteLine("Internal Energy        :" + this.InternalEnergy.ToString());
            Console.WriteLine("Isochoric SpecificHeat :" + this.Cv.ToString());
            Console.WriteLine("Isobaric SpecificHeat  :" + this.Cp.ToString());
            Console.WriteLine("==============================================================");
        }

        public void Display()
        {
            Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            Console.WriteLine("Refrigerant Name       :"+this.Name);
            Console.WriteLine("Critical Density       :"+this.CriticalDensity);
            Console.WriteLine("Critical Temperature   :"+this.CriticalTemperature);
            Console.WriteLine("Critical Pressure      :"+this.CriticalPressure);
            Console.WriteLine("Molecular Weight       :"+this.MolecularWeight);
            Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            Console.WriteLine("There are {0} components in this refrigerent", this.NumberOfComponents);
            Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            for (int i = 0; i < NumberOfComponents; i++)
            {
                Console.WriteLine("Component {0}#", i);
                
                var properities = System.ComponentModel.TypeDescriptor.GetProperties(components[i]);
                foreach (System.ComponentModel.PropertyDescriptor propertyDescriptror in properities)
                {
                    Console.WriteLine(string.Format("{0}:{1}", propertyDescriptror.Name, propertyDescriptror.GetValue(components[i])));
                }
                Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            }
        }
    }
}
