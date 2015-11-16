/* Refprop计算中使用的单位
temperature     	 	        K
pressure/fugacity   	        kPa
density          		        mol/L
composition     		        mole fraction
quality        	                mole basis (moles vapor/total moles)
enthalpy/ internal energy 	    J/mol
Gibbs/ Helmholtz free energy  	J/mol
Entropy/heat capacity        	J/(mol.K)
speed of sound             	    m/s
Joule-Thompson coefficient   	K/kPa
d(p)/d(rho)        		        kPa.L/mol
d2(p)/d(rho)2   		        kPa.(L/mol)^2
viscosity   			        microPa.s (10^-6 Pa.s)
thermal conductivity 		    W/(m.K)
dipole moment   		        debye
surface tension  		        N/m
 * 2013年1月6日
 * 1.为了便于计算Refrigerants.cs中的属性仍然保留Refprop默认的单位。Refrigerants.cs中的UserUnitSystem暂时没有实质作用。
 * 2.使用DisplayThermoDynamicState(UnitSystem unitSystem)来输出设定单位制的数值。
 * 2014.08
 * 1.ThermodynamicParameters of T P H S D CP CV E are added. And unit conversion function is implemented.
 * 2014.09
 * 1.加入了焓熵参考点的设定功能。在更改新的参考点的时候，会自动重新载入制冷剂。
 * 2.在更改制冷剂名称的时候，会自动重新载入制冷剂。
 * 2014.10
 * 1.将check name 和check category分开
 * 2.修正错误: R40和R41为纯制冷剂，而不是R40和R50
 * 3.加入R600A的判断
 * 4.修正错误：在载入new mixture的时候，对每个组分重新check name
 *2015.5
 * 1.加入R600/R744/R717/NH3的判断
 * 2.激活components[]
 * 3.问题：
 * 1)如何定义components修饰符，使得该变量的成员只能在Refrigerant类内被修改和访问。Refrigerant类的继承类只能访问components的成员，不能修改。
 * 2)似乎在某种情况下，通过修改fluid name后，计算新的fluid物性时，存在计算错误的问题。
 *
 *2015.11
 * 1.增加public void FindState(string _name, string InpCode, UnitSystems units, double Prop1, double Prop2)。考虑实用性。
 * 2. Func_* 多态化
 */


using System;
using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.Reflection;

namespace RefpropCSNET
{
    public sealed partial class Refrigerant : Refprop
    {
        private static volatile Refrigerant _instance;
        private static object _lockHelper = new object();

        public const double NearZero = 1.0E-6;
        private UnitSystems currentUnitSystem;
        public UnitSystems CurrentUnitSystem
        {
            get { return currentUnitSystem; }
            set
            {
                if (this.currentUnitSystem != value)
                {
                    currentUnitSystem = value;
                    ConvertUnitSystem(value);
                }
            }
        }  //
        //public RefrigerantCategory Category { get; set; }
        private string name;
        private string currentRef;
        public string Name
        {
            get { return name; }
            set
            {
                if (CheckName(value) != this.name)
                {
                    this.name = value;
                    EnsureCurrentFluid();   //如果制冷剂更改，自动载入新的制冷剂。
                }
            }
        }
        private bool isMassQ = true;    //default 
        public bool IsMassQ { get { return isMassQ; } set { this.isMassQ = value; } }
        /// <summary>
        /// Molecular Weight: g/mol
        /// </summary>
        public double MolecularWeight { get { return base.wm; } }
        /// <summary>
        /// Triple Point Temperature: K
        /// </summary>
        //public double TriplePointTemperature { get { return base.ttp - 273.15; } }
        /// <summary>
        /// Normal Boiling Point Temperature:K
        /// </summary>
        /// 
        //public double NormalBoilingPointTemperature { get { return base.tnbp - 273.15; } private set { } }
        /// <summary>
        /// Critical Temperature
        /// </summary>
        public Temperature CriticalTemperature
        {
            get
            {
                Temperature T = new Temperature(base.tc, UnitSystems.Refprop);
                T.UnitSystem = currentUnitSystem;
                return T;
            }
        }
        /// <summary>
        /// Critical Pressure
        /// </summary>
        public Pressure CriticalPressure
        {
            get
            {
                Pressure P = new Pressure(base.pc, UnitSystems.Refprop);
                P.UnitSystem = currentUnitSystem;
                return P;
            }
        }
        /// <summary>
        /// Critical Density
        /// </summary>
        public Density CriticalDensity
        {
            get
            {
                Density D = new Density(base.wm, base.Dc, UnitSystems.Refprop);
                D.UnitSystem = currentUnitSystem;
                return D;
            }
        }
        public new double[] MoleFractions { get { return base.MoleFractions; } }
        public new double[] MassFractions { get { return base.MassFractions; } }

        /// <summary>
        /// Critical Point Compressibility: pc/(Rgas*Tc*Dc)
        /// </summary>
        //public double CriticalPointCompressibility { get { return base.Zc; } private set { } }
        ///// <summary>
        ///// Accentric Factor
        ///// </summary>
        //public double AccentricFactor { get { return base.acf; } private set { } }
        ///// <summary>
        ///// Dipole Moment
        ///// </summary>
        //public double DipoleMoment { get { return base.dip; } private set { } }
        ///// <summary>
        ///// Gas Constatnt :J/mol-K
        ///// </summary>
        //public double GasConstatnt_R { get { return base.Rgas; } private set { } }

        /// <summary>
        /// Number of Components 
        /// </summary>
        public Int32 NumberOfComponents { get { return base.nc; } }

        public Comonent[] components;  //如何定义components修饰符，使得该变量的成员只能在Refrigerant类内被修改和访问。Refrigerant类的继承类只能访问，不能修改。
        /// <summary>
        /// Temperature
        /// </summary>
        public ThermodynamicParameter T { get; set; }
        /// <summary>
        /// Pressure
        /// </summary>
        public ThermodynamicParameter P { get; set; }
        /// <summary>
        /// Enthalpy
        /// </summary>
        public ThermodynamicParameter H { get; set; }
        /// <summary>
        /// Entropy
        /// </summary>
        public ThermodynamicParameter S { get; set; }
        /// <summary>
        /// Density
        /// </summary>
        public ThermodynamicParameter D { get; set; }
        /// <summary>
        /// isobaric (constant p) heat capacity
        /// </summary>
        public ThermodynamicParameter CP { get; set; }
        /// <summary>
        /// isochoric (constant V) heat capacity
        /// </summary>
        public ThermodynamicParameter CV { get; set; }
        /// <summary>
        /// Internal energy
        /// </summary>
        public ThermodynamicParameter E { get; set; }
        /// <summary>
        /// Quality
        /// </summary>
        public ThermodynamicParameter X { get; set; }
        /// <summary>
        /// Viscosity [microPa.s (10^-6 Pa.s)]  10^-6 kg*m/s^2/m^2.s=10^-6 kg/s.m
        /// </summary>
        public ThermodynamicParameter U { get; set; }
        /// <summary>
        /// Thermal conductivity
        /// </summary>
        public ThermodynamicParameter K { get; set; }
        List<ThermodynamicParameter> ThermodynamicParameters = new List<ThermodynamicParameter>();
        /// <summary>
        /// temperature of the current refrigerant : K
        /// </summary>
        //protected double Temperature { get { return base.temperature; } }//- 273.15
        //public TemperatureUnits Temperature=new TemperatureUnits(;

        /// <summary>
        ///  pressure of the current refrigerant : kPa
        /// </summary>
        //protected double Pressure { get { return base.pressure; } }

        /// <summary>
        /// overall molar density of the current refrigerant : mol/L  // * MolecularWeight >g/L
        /// </summary>
        public double MoleDensity { get { return base.density; } }

        /// <summary>
        /// overall mass density of the current refrigerant : g/L
        /// </summary>
        public double MassDensity { get { return base.density * MolecularWeight; } }

        /// <summary>
        /// moles vapor/total moles
        /// </summary>
        public double MoleQuality { get { return base.molequality; } }
        /// <summary>
        /// mass vapor / total mass
        /// </summary>
        public double MassQuality { get { return base.massquality; } }

        //public double Quality { get { return (IsMassQ ? MassQuality : MoleQuality); } }   //right??

        /// <summary>
        /// Internal Energy of the current refrigerant : J/mol
        /// </summary>
        //protected double InternalEnergy { get { return base.internalenergy; } }

        /// <summary>
        /// Enthalpy of the current refrigerant : J/mol //  /MolecularWeight>kJ/kg
        /// </summary>
        //protected double Enthalpy { get { return base.enthalpy; } } // / MolecularWeight

        /// <summary>
        /// Entropy of the current refrigerant :J/mol-K   /// kJ/(kg-K)
        /// </summary>
        //protected double Entropy { get { return base.entropy; } } // / MolecularWeight

        /// <summary>
        /// isochoric (constant V) heat capacity [J/mol-K]
        /// </summary>
        //protected double Cv { get { return base.cv; } }

        /// <summary>
        /// isobaric (constant p) heat capacity [J/mol-K]
        /// </summary>
        //protected double Cp { get { return base.cp; } }

        /// <summary>
        /// isenthalpic Joule-Thompson coefficient [K/kPa]
        /// </summary>
        //public double HJT { get { return base.hjt; } }

        /// <summary>
        /// Viscosity [microPa.s (10^-6 Pa.s)]  10^-6 kg*m/s^2/m^2.s=10^-6 kg/s.m
        /// </summary>
        //public double Viscosity { get { return base.eta; } }

        public double IsothermalCompressibility { get { return base.xkappa; } }
        public double IsentropicExpansionCofficient { get { return base.xisenk; } }
        public double IsothermalExpansionCoefficient { get { return base.xkt; } }
        public double AdiabaticCompressibility { get { return base.betas; } }
        public double AdiabaticBulkModulus { get { return base.bs; } }
        public double IsothermalBulkModulus { get { return base.xkkt; } }
        public double IsothermalThrottlingCoefficient { get { return base.thrott; } }
        public double InternalPressure { get { return base.thrott; } }
        public double SpecificHeatInput { get { return base.spht; } }

        /// <summary>
        /// 动力粘度Kinematic Viscosity [10^-6 m^2/s]
        /// </summary>
        public double KinematicViscosity { get { return base.eta / base.density / base.wm; } }
        /// <summary>
        /// thermal conductivity  [W/(m.K)]
        /// </summary> 
        //public double ThermalConductivity { get { return base.tcx; } }

        /// <summary>
        /// Thermal Diffusivity [m^2/s]
        /// </summary>
        public double ThermalDiffusivity { get { return 1.0E-3 * base.tcx / base.density / base.cp; } }
        /// <summary>
        /// volume expansivity [1/K]
        /// </summary>
        public double VolumeExpansivity { get { return base.beta; } }
        /// <summary>
        /// 普朗特数Prandtl Number
        /// </summary>
        public double PrandtlNumber
        {
            get
            {
                if (base.tcx < NearZero) return -10000;  //ThermalConductivity
                if (base.molequality > 0 && base.molequality < 1) return -20000;
                return base.cp / base.wm * base.eta / base.tcx / 1000;  //eta:运动粘度
            }
        }
        public RefrigerantCategory Category { get; set; }

        public ReferenceState? Reference
        {
            get { return (ReferenceState)Enum.Parse(typeof(ReferenceState), base.hrf, true); }
            set
            {
                if (!value.ToString().Equals(base.hrf, StringComparison.OrdinalIgnoreCase))
                {
                    base.hrf = value.ToString();
                    this.Load(Name, value);
                }
            }
        }

        //private Refrigerant(): this("") { }
        /// <summary>
        /// Constructor: The refrigerant category can be determined by its name
        /// </summary>
        /// <param name="refrigerantName"></param>
        /// <param name="reference"></param>
        /// <param name="currentUnits"></param>
        private Refrigerant(string refrigerantName, ReferenceState? reference, UnitSystems currentUnits)
        {
            Load(refrigerantName, reference);
            initThermDynamicParameters();
            this.CurrentUnitSystem = currentUnits;
        }
        //private Refrigerant(string refrigerantName)
        //    : this(refrigerantName, ReferenceState.DEF, UnitSystems.Refprop) { }
        //private Refrigerant(string refrigerantName, UnitSystems currentUnits)
        //    : this(refrigerantName, ReferenceState.DEF, currentUnits) { }

        public static Refrigerant GetInstance(string refrigerantName, ReferenceState? reference = ReferenceState.DEF, UnitSystems currentUnits = UnitSystems.SI)
        {
            if (_instance == null)
            {
                lock (_lockHelper)
                {
                    if (_instance == null)
                    {
                        _instance = new Refrigerant(refrigerantName, reference, currentUnits);
                    }
                }
            }
            else
            {
                _instance.Name = refrigerantName;
                if (reference != null) _instance.Reference = reference;
                _instance.CurrentUnitSystem = currentUnits;
            }
            return _instance;
        }
        public static Refrigerant GetInstance(string refrigerantName, UnitSystems currentUnits)
        {
            return GetInstance(refrigerantName, ReferenceState.DEF, currentUnits);
        }
        public static Refrigerant GetInstance()
        {
            return GetInstance("");
        }

        private void initThermDynamicParameters()
        {
            //这里有简便的方法么？？？？
            //this.T = new ThermodynamicParameter(this.wm, UnitTypes.T, 0.0, this.UserUnitSystem);
            this.T = new Temperature(this.CurrentUnitSystem, ThermodynamicParameters);
            this.P = new Pressure(this.CurrentUnitSystem, ThermodynamicParameters);
            this.H = new Enthalpy(this.CurrentUnitSystem, ThermodynamicParameters);
            this.S = new Entropy(this.CurrentUnitSystem, ThermodynamicParameters);
            this.D = new Density(this.CurrentUnitSystem, ThermodynamicParameters);
            this.CP = new SpecificHeat(this.CurrentUnitSystem, ThermodynamicParameters);
            this.CV = new SpecificHeat(this.CurrentUnitSystem, ThermodynamicParameters);
            this.E = new InternalEnergy(this.CurrentUnitSystem, ThermodynamicParameters);
            this.X = new Quality(this.CurrentUnitSystem, ThermodynamicParameters);
            this.U = new Viscosity(this.CurrentUnitSystem, ThermodynamicParameters);
            this.K = new ThermalConductivity(this.CurrentUnitSystem, ThermodynamicParameters);
            //this.Z = new CompressibilityFactor(this.UserUnitSystem);
            //this.Prt = new PrantdlNumber(this.UserUnitSystem);
        }
        public string CheckName(string FluidName)
        {
            //if (String.IsNullOrEmpty(FluidName)) RefpropErrorHandler.ErrorHandler("Invalid fluid name");
            FluidName = FluidName.Replace(" ", "");  //remove all the spaces
            FluidName = FluidName.Replace("-", "");  //Remove "-" if R-410A etc. is given as input
            FluidName = FluidName.ToUpper();   //to Uppercase 
            if (FluidName == "AIR") FluidName = "nitrogen=7812,argon=0092,oxygen=2096";
            else if (FluidName == "BUTENE") FluidName = "1BUTENE";
            else if (FluidName == "NH3" || FluidName == "R717") FluidName = "AMMONIA";
            else if (FluidName == "CARBONDIOXIDE" || FluidName == "R744") FluidName = "CO2";
            else if (FluidName == "CARBONMONOXIDE") FluidName = "CO";
            else if (FluidName == "CARBONYLSULFIDE") FluidName = "COS";
            else if (FluidName == "R290") FluidName = "PROPANE";   //丙烷R290
            else if (FluidName == "CIS-BUTENE") FluidName = "C2BUTENE";
            else if (FluidName == "CYCLOHEXANE") FluidName = "CYCLOHEX";
            else if (FluidName == "CYCLOPENTANE") FluidName = "CYCLOPEN";
            else if (FluidName == "CYCLOPROPANE" || FluidName == "RC270") FluidName = "CYCLOPRO";
            else if (FluidName == "DEUTERIUM") FluidName = "D2";
            else if (FluidName == "DIMETHYLCARBONATE") FluidName = "DMC";
            else if (FluidName == "DIMETHYLETHER") FluidName = "DME";
            else if (FluidName == "DIETHYLETHER") FluidName = "DEE";
            else if (FluidName == "DODECANE") FluidName = "C12";
            else if (FluidName == "ETHYLBENZENE") FluidName = "EBENZENE";
            else if (FluidName == "HEAVYWATER") FluidName = "D2O";
            else if (FluidName == "HYDROGENCHLORIDE") FluidName = "HCL";
            else if (FluidName == "HYDROGENSULFIDE") FluidName = "H2S";
            else if (FluidName == "R600") FluidName = "BUTANE"; //丁烷R-600
            else if (FluidName == "R600A") FluidName = "ISOBUTAN"; //异丁烷R-600a 
            else if (FluidName == "IBUTANE") FluidName = "ISOBUTAN";
            else if (FluidName == "ISOBUTANE") FluidName = "ISOBUTAN";
            else if (FluidName == "ISOBUTENE") FluidName = "IBUTENE";
            else if (FluidName == "ISOHEXANE") FluidName = "IHEXANE";
            else if (FluidName == "ISOPENTANE") FluidName = "IPENTANE";
            else if (FluidName == "ISOOCTANE") FluidName = "IOCTANE";
            else if (FluidName == "METHYLCYCLOHEXANE") FluidName = "C1CC6";
            else if (FluidName == "METHYLLINOLENATE") FluidName = "MLINOLEN";
            else if (FluidName == "METHYLLINOLEATE") FluidName = "MLINOLEA";
            else if (FluidName == "METHYLOLEATE") FluidName = "MOLEATE";
            else if (FluidName == "METHYLPALMITATE") FluidName = "MPALMITA";
            else if (FluidName == "METHYLSTEARATE") FluidName = "MSTEARAT";
            else if (FluidName == "NEOPENTANE") FluidName = "NEOPENTN";
            else if (FluidName == "NITROGENTRelse ifLUORIDE") FluidName = "NF3";
            else if (FluidName == "NITROUSOXIDE") FluidName = "N2O";
            else if (FluidName == "ORTHOHYDROGEN") FluidName = "ORTHOHYD";
            else if (FluidName == "PARAHYDROGEN") FluidName = "PARAHYD";
            else if (FluidName == "PERFLUOROBUTANE") FluidName = "C4F10";
            else if (FluidName == "PERFLUOROPENTANE") FluidName = "C5F12";
            else if (FluidName == "PROPYLCYCLOHEXANE") FluidName = "C3CC6";
            else if (FluidName == "PROPYLENE") FluidName = "PROPYLEN";
            else if (FluidName == "SULFUR HEXAFLUORIDE") FluidName = "SF6";
            else if (FluidName == "TRANS-BUTENE") FluidName = "T2BUTENE";
            else if (FluidName == "TRIFLUOROIODOMETHANE") FluidName = "CF3I";
            else if (FluidName == "SULFURDIOXIDE") FluidName = "SO2";
            else if (FluidName == "SULFURHEXAFLUORIDE") FluidName = "SF6";
            else if (FluidName == "UNDECANE") FluidName = "C11";
            else if (FluidName == "STEAM") FluidName = "WATER";
            else if (FluidName == "H2O") FluidName = "WATER";
            else if (FluidName == "NATURALGAS") FluidName = "NGSAMPLE.MIX";
            //FluidName = FluidName.Replace(".FLD", "");
            //FluidName = FluidName.Replace(".MIX", "");   
            return FluidName;
        }

        private void CheckCategory(string FluidName)
        {
            //除了R40和R41外，以R4和R5开头的制冷剂均为混合制冷剂
            if (FluidName.EndsWith(".fld", true, null) || FluidName.EndsWith(".ppf", true, null)) this.Category = RefrigerantCategory.PureFluid; //以扩展名确定
            else if (FluidName.EndsWith(".mix", true, null)) this.Category = RefrigerantCategory.PredefinedMixture;//以扩展名确定
            else if (FluidName.Contains("=") && FluidName.Contains(",")) this.Category = RefrigerantCategory.NewMixture;
            else if (FluidName.StartsWith("R4") || FluidName.StartsWith("R5"))
            {
                if (FluidName != "R40" && FluidName != "R41") this.Category = RefrigerantCategory.PredefinedMixture;
                else this.Category = RefrigerantCategory.PureFluid;
            }
            else this.Category = RefrigerantCategory.PureFluid;  //默认为单一制冷剂            
        }

        private void Load(string refrigerantName, ReferenceState? reference)
        {
            refrigerantName = CheckName(refrigerantName);
            if (String.IsNullOrEmpty(refrigerantName)) { currentRef = refrigerantName; return; }
            CheckCategory(refrigerantName);
            if (Category == RefrigerantCategory.PureFluid)
            {
                if (refrigerantName.EndsWith(".fld", true, null)
                    || refrigerantName.EndsWith(".ppf", true, null))
                    base.LoadFluid(refrigerantName, 1, new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, reference.ToString(), isMassQ);
                else
                    base.LoadFluid(refrigerantName + ".fld", 1, new double[] { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, reference.ToString(), isMassQ);
            }
            else if (Category == RefrigerantCategory.PredefinedMixture)
            {
                if (refrigerantName.EndsWith(".mix", true, null))
                    base.LoadPredefinedFluid(refrigerantName, reference.ToString());
                else base.LoadPredefinedFluid(refrigerantName + ".mix", reference.ToString());
            }
            else if (Category == RefrigerantCategory.NewMixture)
            {
                string[] a = refrigerantName.Split(',');
                double[] b = new double[20];
                string c = "";
                double sum = 0.0;

                System.Globalization.NumberFormatInfo provider = new System.Globalization.NumberFormatInfo();
                provider.NumberDecimalSeparator = ".";
                provider.NumberGroupSeparator = ",";

                for (int i = 0; i < a.Length; i++)
                {
                    string s = " ";
                    c += CheckName(a[i].Split('=')[0]);  //
                    s = a[i].Split('=')[1];
                    try
                    {
                        b[i] = Convert.ToDouble(s, provider);
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Unable to convert '{0}' to a Double.", s);
                    }
                    catch (OverflowException)
                    {
                        Console.WriteLine("'{0}' is outside the range of a Double.", s);
                    }
                    c += ".fld|";
                    sum += b[i];
                }
                for (int i = 0; i < a.Length; i++) b[i] = b[i] / sum;

                base.LoadFluid(c, a.Length, b, reference.ToString(), isMassQ);
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
                components[i].MassFraction = MassFractions[i];
                components[i].MoleFraction = MoleFractions[i];
                if (NumberOfComponents == 1) components[i].Name = refrigerantName;
                else if (Category == RefrigerantCategory.PredefinedMixture)
                    components[i].Name = base.hfiles.ToString().Split('|')[i];
                else if (Category == RefrigerantCategory.NewMixture)
                    components[i].Name = refrigerantName.Split(',', '=')[i * 2];
                string temps = components[i].Name;
                components[i].Name = temps.Substring(temps.LastIndexOf('/') + 1).Replace(".FLD", "");
            }
            currentRef = refrigerantName;
            this.name = refrigerantName;
            base.hrf = reference.ToString();   //好奇怪的做法啊
        }

        private void SetNumOfComp(int n)
        {
            base.SetNC(n);
        }


        /// <summary>
        /// use UserUnitSystem as defaut if no UnitSystem is specified
        /// </summary>
        /// <param name="InpCode"></param>
        /// <param name="Prop1"></param>
        /// <param name="Prop2"></param>
        public void FindState(string InpCode, double Prop1, double Prop2)
        {
            this.FindState(InpCode, CurrentUnitSystem, Prop1, Prop2);
        }

        /// <summary>
        /// Find the state with two or three parameters. For expample, P and T.
        /// </summary>
        /// <param name="InpCode">Function Code: PT, PH, etc.</param>
        /// <param name="units">UnitSystem: SI, SIwithC, etc.</param>
        /// <param name="Prop1"></param>
        /// <param name="Prop2"></param>
        public void FindState(string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.FindState(this.Name, InpCode, units, Prop1, Prop2);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_name">refrigerant Name</param>
        /// <param name="InpCode">Function Code: PT, PH, etc.</param>
        /// <param name="units">UnitSystem: SI, SIwithC, etc.</param>
        /// <param name="Prop1">value 1</param>
        /// <param name="Prop2">value 2</param>
        public void FindState(string _name, string InpCode, UnitSystems units, double Prop1, double Prop2)
        {
            this.Name = _name;
            EnsureCurrentFluid();
            this.CurrentUnitSystem = units;
            InpCode = InpCode.Replace(" ", "");
            InpCode = InpCode.ToUpper();
            if (InpCode.Contains("P"))
            {
                this.P.UnitSystem = units;
                if (InpCode.IndexOf('P') == 0) this.P.Value = Prop1;
                else this.P.Value = Prop2;
            }
            if (InpCode.Contains("T"))
            {
                this.T.UnitSystem = units;
                if (InpCode.IndexOf('T') == 0) this.T.Value = Prop1;
                else this.T.Value = Prop2;
            }
            if (InpCode.Contains("H"))
            {
                this.H.UnitSystem = units;
                if (InpCode.IndexOf('H') == 0) this.H.Value = Prop1;
                else this.H.Value = Prop2;
            }
            if (InpCode.Contains("D"))
            {
                this.D.UnitSystem = units;
                if (InpCode.IndexOf('D') == 0) this.D.Value = Prop1;
                else this.D.Value = Prop2;
            }
            if (InpCode.Contains("S"))
            {
                this.S.UnitSystem = units;
                if (InpCode.IndexOf('S') == 0) this.S.Value = Prop1;
                else this.S.Value = Prop2;
            }
            if (InpCode.Contains("E"))
            {
                this.E.UnitSystem = units;
                if (InpCode.IndexOf('E') == 0) this.E.Value = Prop1;
                else this.E.Value = Prop2;
            }
            if (InpCode.Contains("Q"))
            {
                this.X.UnitSystem = units;
                if (InpCode.IndexOf('Q') == 0) this.X.Value = Prop1;
                else this.X.Value = Prop2;
            }
            ConvertUnitSystem(UnitSystems.Refprop); //convert to the Refprop unit system before calculation
            if (InpCode == "PH" || InpCode == "HP") base.PHFLSH(P.Value, H.Value);
            else if (InpCode == "PE" || InpCode == "EP") base.PEFLSH(P.Value, E.Value);
            else if (InpCode == "PQ" || InpCode == "QP") base.PQFLSH(P.Value, X.Value, isMassQ ? 2 : 1);  //2:mass fraction 1:mole fraction
            else if (InpCode == "PS" || InpCode == "SP") base.PSFLSH(P.Value, S.Value);
            else if (InpCode == "TD" || InpCode == "DT") base.TDFLSH(T.Value, D.Value);
            else if (InpCode.StartsWith("TE") || InpCode.StartsWith("ET"))
            {
                if (InpCode.EndsWith("<")) base.TEFLSH(T.Value, E.Value, 1); //1 = return lower density root;
                else base.TEFLSH(T.Value, E.Value, 2);  //2 = return higher density root
            }
            else if (InpCode == "TD" || InpCode == "TD") base.TDFLSH(T.Value, D.Value);
            else if (InpCode == "PT" || InpCode == "TP") base.TPFLSH(T.Value, P.Value);
            else if (InpCode == "TQ" || InpCode == "QT") base.TQFLSH(T.Value, X.Value, isMassQ ? 2 : 1);//2:mass fraction 1:mole fraction
            else if (InpCode.StartsWith("TH") || InpCode.StartsWith("HT"))
            {
                if (InpCode.EndsWith("<")) base.THFLSH(T.Value, H.Value, 1);  // 1 = return lower density root;
                else base.THFLSH(T.Value, H.Value, 2);   //2 = return higher density root.
            }
            else if (InpCode == "DE" || InpCode == "ED") base.DEFLSH(D.Value, E.Value);
            else if (InpCode == "DS" || InpCode == "SD") base.DSFLSH(D.Value, S.Value);
            else if (InpCode == "ES" || InpCode == "SE") base.ESFLSH(E.Value, S.Value);
            else RefpropErrorHandler.ErrorHandler(this, "No such function");
            T.Value = base.temperature;
            P.Value = base.pressure;
            H.Value = base.enthalpy;
            S.Value = base.entropy;
            D.Value = base.density;
            CP.Value = base.cp;
            CV.Value = base.cv;
            E.Value = base.internalenergy;
            //X.Value = this.Quality;
            U.Value = base.eta;
            K.Value = base.tcx;
            //Prt.Value = this.PrandtlNumber;
            //Z.Value = base.z;
            //
            if (!InpCode.Contains("Q"))
            {
                if (isMassQ && this.molequality >= 0 && this.molequality <= 1.0) { base.QMASS(); X.Value = base.massquality; }
                else X.Value = base.molequality;
            }
            ConvertUnitSystem(CurrentUnitSystem); //convert to the user defined unit system
        }
        private void ConvertUnitSystem(UnitSystems unitsystem)
        {
            for (int i = 0; i < ThermodynamicParameters.Count; i++)
            {
                ThermodynamicParameters[i].MolecularWeight = this.MolecularWeight;
                ThermodynamicParameters[i].UnitSystem = unitsystem;
            }
        }
        private void ConvertUnit(ThermodynamicParameter tp, UnitSystems unit)
        {
            tp.UnitSystem = unit;
        }


        #region Functions without unit conversion, i.e.,using the Refprop units
        /// <summary>
        /// Calculate the saturated state given the staturated temperature
        /// </summary>
        /// <param name="temperature">staturated temperature</param>
        /// <param name="PhaseFlag">saturated liquid ? saturated vapor</param>
        private void FindSaturatedStateWithTemperature(double temperature, SaturationPoint PhaseFlag)
        {
            EnsureCurrentFluid();
            base.GetSaturatedPressure(temperature, (Int32)PhaseFlag);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pressure"></param>
        /// <param name="PhaseFlag"></param>
        private void FindSaturatedStateWithPressure(double pressure, SaturationPoint PhaseFlag)
        {
            EnsureCurrentFluid();
            base.GetSaturatedTemperature(pressure, (Int32)PhaseFlag);
        }

        /// <summary>
        /// a determination of the thermodynamic state given two independent variables(temperature and pressure) plus composition
        /// </summary>
        /// <param name="T">temperture of the refrigerant: k</param>
        /// <param name="P">pressure of the refrigerant: Kpa</param>
        private void FindStateWithTP(double T, double P)
        {
            EnsureCurrentFluid();
            base.TPFLSH(T, P);
        }
        /// <summary>
        /// a determination of the thermodynamic state given two independent variables(temperature and density) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="D"></param>
        private void FindStateWithTD(double T, double D)
        {
            EnsureCurrentFluid();
            base.TDFLSH(T, D);
        }
        /// <summary>
        /// a determination of the thermodynamic state given two independent variables(temperature and enthalpy) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="H"></param>
        /// <param name="root"></param>
        private void FindStateWithTH(double T, double H, int root)
        {
            EnsureCurrentFluid();
            base.THFLSH(T, H, root);
        }
        /// <summary>
        /// a determination of the thermodynamic state given two independent variables(temperature and energy) plus composition
        /// </summary>
        /// <param name="T"></param>
        /// <param name="E"></param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        private void FindStateWithTE(double T, double E, int root)
        {
            EnsureCurrentFluid();
            base.TEFLSH(T, E, root);
        }

        private void FindStateWithPD(double P, double D)
        {
            EnsureCurrentFluid();
            base.PDFLSH(P, D);
        }

        //      subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine PHFLSH (p,h,z,t,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        private void FindStateWithPH(double P, double H)
        {
            EnsureCurrentFluid();
            base.PHFLSH(P, H);
        }
        //      subroutine PSFLSH (p,s,z,t,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        private void FindStateWithPS(double P, double S)
        {
            EnsureCurrentFluid();
            base.PSFLSH(P, S);
        }
        //      subroutine PEFLSH (p,e,z,t,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
        private void FindStateWithPE(double P, double E)
        {
            EnsureCurrentFluid();
            base.PEFLSH(P, E);
        }
        //      subroutine HSFLSH (h,s,z,t,p,D,Dl,Dv,x,y,q,e,cv,cp,w,ierr,herr)
        private void FindStateWithHS(double H, double S)
        {
            EnsureCurrentFluid();
            base.HSFLSH(H, S);
        }

        //      subroutine ESFLSH (e,s,z,t,p,D,Dl,Dv,x,y,q,h,cv,cp,w,ierr,herr)
        private void FindStateWithES(double E, double S)
        {
            EnsureCurrentFluid();
            base.ESFLSH(E, S);
        }
        //      subroutine DHFLSH (D,h,z,t,p,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        private void FindStateWithDH(double D, double H)
        {
            EnsureCurrentFluid();
            base.DHFLSH(D, H);
        }
        //      subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        private void FindStateWithDS(double D, double S)
        {
            EnsureCurrentFluid();
            base.DSFLSH(D, S);
        }
        //      subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)

        private void FindStateWithDE(double D, double E)
        {
            EnsureCurrentFluid();
            base.DEFLSH(D, E);
        }
        //      subroutine TQFLSH (t,q,z,kq,p,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        private void FindStateWithTQ(double T, double Q, int kq)
        {
            EnsureCurrentFluid();
            base.TQFLSH(T, Q, kq);
        }
        //      subroutine PQFLSH (p,q,z,kq,t,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        //kq=1 mole base  kq=2 mass base
        private void FindStateWithPQ(double P, double Q, int kq)
        {
            EnsureCurrentFluid();
            base.PQFLSH(P, Q, kq);
        }


        /// <summary>
        /// 已知单相状态的T P计算焓值
        /// </summary>
        /// <param name="T">温度[K]</param>
        /// <param name="P">压力[kPa]</param>
        /// <param name="phase">相态</param>
        /// <returns></returns>
        private double EnthalpyCal(double T, double P, SubstancePhase phase)
        {
            EnsureCurrentFluid();
            Int32 phaseflag;
            if ((phase == SubstancePhase.Saturated_Liquid) || (phase == SubstancePhase.Subcooled)) phaseflag = 1;
            else if ((phase == SubstancePhase.Saturated_Vapor) || (phase == SubstancePhase.Superheated)) phaseflag = 2;
            else
            {
                RefpropErrorHandler.ErrorHandler(this, "错误的相态");
                return 0;
            }
            base.S_TPRHO(T, P, this.MoleFractions, phaseflag);
            base.EnthalpyCal(T, this.MoleDensity, this.MoleFractions);
            return base.enthalpy;
        }
        /// <summary>
        ///  已知单相状态的T P计算熵值
        /// </summary>
        /// <param name="T">温度[K]</param>
        /// <param name="P">压力[kPa]</param>
        /// <param name="phase">相态</param>
        /// <returns></returns>
        private double EntropyCal(double T, double P, SubstancePhase phase)
        {
            EnsureCurrentFluid();
            Int32 phaseflag = 0;
            if ((phase == SubstancePhase.Saturated_Liquid) || (phase == SubstancePhase.Subcooled)) phaseflag = 1;
            else if ((phase == SubstancePhase.Saturated_Vapor) || (phase == SubstancePhase.Superheated)) phaseflag = 2;
            else
            {
                RefpropErrorHandler.ErrorHandler(this, "错误的相态");
                return 0;
            }
            base.S_TPRHO(T, P, this.MassFractions, phaseflag);
            base.EntropyCal(T, this.MoleDensity, this.MoleFractions);
            return base.entropy;
        }


        private double PressureCal(double T, double D)
        {
            EnsureCurrentFluid();
            base.PressureCal(T, D, this.MoleFractions);
            return base.pressure;
        }

        /// <summary>
        /// 已知单相状态的P H计算状态点
        /// </summary>
        /// <param name="P"></param>
        /// <param name="H"></param>
        /// <param name="phase"></param>
        private void FindSPhaseStateWithPH(double P, double H, SubstancePhase phase)
        {
            EnsureCurrentFluid();
            Int32 phaseflag = 0;
            if ((phase == SubstancePhase.Saturated_Liquid) || (phase == SubstancePhase.Subcooled)) phaseflag = 1;
            else if ((phase == SubstancePhase.Saturated_Vapor) || (phase == SubstancePhase.Superheated)) phaseflag = 2;
            else
            {
                RefpropErrorHandler.ErrorHandler(this, "错误的相态");
            }
            base.S_PHFL1(P, H, phaseflag);
            //base.Thermal2Cal(this.Temperature, this.MoleDensity, this.MoleFractions);
            //base.TRNPRPCal(this.temperature, this.density, this.MoleFractions);
        }
        /// <summary>
        /// 已知单相状态的T P计算状态点
        /// </summary>
        /// <param name="T">温度[K]</param>
        /// <param name="P">压力[kPa]</param>
        /// <param name="phase">相态</param>
        private void FindSPhaseStateWithTP(double T, double P, SubstancePhase phase)
        {
            EnsureCurrentFluid();
            Int32 phaseflag = 0;
            if ((phase == SubstancePhase.Saturated_Liquid) || (phase == SubstancePhase.Subcooled)) phaseflag = 1;
            else if ((phase == SubstancePhase.Saturated_Vapor) || (phase == SubstancePhase.Superheated)) phaseflag = 2;
            else
            {
                RefpropErrorHandler.ErrorHandler(this, "错误的状态");
            }
            base.S_TPRHO(T, P, this.MoleFractions, phaseflag);
            //base.Thermal2Cal(this.Temperature, this.MoleDensity, this.MoleFractions);
            //base.TRNPRPCal(this.temperature, this.density, this.MoleFractions);
        }

        /// <summary>
        /// 调用base.Thermal3Cal计算一些热动力学参数
        /// 此函数需要单独调用
        /// </summary>
        /// <param name="T"></param>
        /// <param name="D"></param>
        private void CalMiscellaneousThermodynamic(double T, double D)
        {
            EnsureCurrentFluid();
            base.Thermal3Cal(T, D, this.MoleFractions);
        }
        #endregion

        /// <summary>
        /// Displays all the thermodynamicParameters with SI unit system
        /// </summary>
        public void DisplayThermoDynamicState()
        {
            DisplayThermoDynamicState(CurrentUnitSystem);
        }
        /// <summary>
        /// Displays all the thermodynamicParameters with specified unit system
        /// </summary>
        /// <param name="unitSystem">Unit system</param>
        public void DisplayThermoDynamicState(UnitSystems unitSystem)
        {
            Console.WriteLine("<--");
            for (int i = 0; i < ThermodynamicParameters.Count; i++)
            {
                ThermodynamicParameters[i].UnitSystem = unitSystem;
                Console.WriteLine(ThermodynamicParameters[i].GetType().Name + " : " + ThermodynamicParameters[i].ToString());
            }
            Console.WriteLine("-->");
        }

        /// <summary>
        /// Display information of the current refrigerant
        /// </summary>
        public void DisplayInformation()
        {
            Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            Console.WriteLine("Refrigerant Name       :" + this.Name);
            Console.WriteLine("CriticalDensity        :" + this.CriticalDensity.ToString());
            Console.WriteLine("CriticalTemperature    :" + this.CriticalTemperature.ToString());
            Console.WriteLine("CriticalPressure       :" + this.CriticalPressure.ToString());
            Console.WriteLine("MolecularWeight        :" + this.MolecularWeight.ToString());
            //Console.WriteLine("CriticalDensity        :" + this.CriticalDensity.ToString());
            Console.WriteLine("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            Console.WriteLine("There are {0} components in this refrigerent", this.NumberOfComponents);
            for (int i = 0; i < NumberOfComponents; i++)
            {
                Console.WriteLine("Component #{0}", i + 1);
                Console.WriteLine(components[i].Name);
                Console.WriteLine(components[i].AccentricFactor.ToString());
                Console.WriteLine(components[i].Compressibility);
                Console.WriteLine(components[i].CriticalDensity);
                Console.WriteLine(components[i].CriticalPressure);
                Console.WriteLine(components[i].CriticalTemperature);
                Console.WriteLine(components[i].DipoleMoment);
                Console.WriteLine(components[i].GasConstant);
                Console.WriteLine(components[i].MassFraction);
                Console.WriteLine(components[i].MolecularWeight);
                //var properities = System.ComponentModel.TypeDescriptor.GetProperties(components[i]);
            }
        }

        /// <summary> 
        /// </summary>
        public void EnsureCurrentFluid()
        {
            if ((!string.Equals(currentRef, this.Name, StringComparison.OrdinalIgnoreCase)))
            {
                this.Load(this.Name, this.Reference);
                for (int i = 0; i < ThermodynamicParameters.Count; i++)
                {
                    ThermodynamicParameters[i].MolecularWeight = this.MolecularWeight;
                }
            }
        }

        public void EnableCalThermalQuantities() { base.calThermalQuantities = true; }
        public void DisableCalThermalQuantities() { base.calThermalQuantities = false; }
    }


    public class Comonent
    {
        public string Name { get; set; }
        public double MoleFraction { set; get; }
        public double MassFraction { get; set; }
        public double MolecularWeight { get; set; }
        public double TriplePointTemperature { get; set; }
        public double NormalBoilingPointTemperature { get; set; }
        public double CriticalTemperature { get; set; }
        public double CriticalPressure { get; set; }
        public double CriticalDensity { get; set; }
        public double Compressibility { get; set; }
        public double AccentricFactor { get; set; }
        public double DipoleMoment { get; set; }
        public double GasConstant { get; set; }
        public Comonent()
        {

        }
    }
}
