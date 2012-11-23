#define setup_use_string

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace sc.net
{
    /// <summary>
    /// Data Passed between C# and Fortran is by Reference
    /// A Example:
    /// SUBROUTINE TESTPROC( INTPARM, REALPARM,CHARPARM, CHARARRAYPARM )
    /// cDEC$ ATTRIBUTES DLLEXPORT :: TESTPROC
    /// INTEGER INTPARM
    /// REAL(8) REALPARM
    /// CHARACTER*255 CHARPRAM
    /// CHARACTER*255 CHARARRAYPARM(20)
    /// END SUBROUTINE
    /// C# CODE LIKE THIS:
    /// [DllImport(DLL_path, EntryPoint = "TESTPROC", SetLastError = true)]
    ///    public static extern void TESTPROC(ref Int32 intparm,ref Double realparm,
    ///         [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder charparm, 
    ///         [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder chararrayparm, Int32 l1, Int32 l2);
    /// when called:
    /// Int32 intparm;
    /// Double realparm;
    /// StringBuilder charparm = new StringBuilder(255);
    /// StringBuilder chararrayparm = new StringBuilder(5100);
    /// TESTPROC(intparm, realparm, charparm, chararrayparm,charparm.length, chararrayparm.length);
    /// Noticed that CHARPARM & CHARARRAYPARM are passed from Fortran to C#
    /// When a string or string Array is passed from c# to Fortran:
    /// C#                         FORTRAN
    /// String L                   CHARACTER*255  La
    /// String LL                  CHARACTER*255  LLa(20)
    /// LL="abc|def|ghi|"   →     LLa(1)=abc       LLa(2)=def     LLa(3)=ghi
    /// </summary>
    class Refprop
    {
        private const Int32 MaxComosition = 20;
        private const Int32 RefpropFluidPathLength = (RefpropCharLength + 1) * MaxComosition; //+1 for '|' between file names
        private const Int32 RefpropCharLength = 255;
        private const Int32 FilePathLength = 255;
        private const Int32 LengthOfReference = 3;
        private const Int32 ErrorMessageLength = 255;

        public Int32 ierr;
        public StringBuilder herr = new StringBuilder(ErrorMessageLength);

        UnitsBasis CurrentUnitsBasis;

        protected Int32 nc;   // number of components
        protected double[] MoleFractions = new double[MaxComosition]; // this is array of components  x[nc]
        protected double[] MassFractions = new double[MaxComosition];   

        private char[] hfld = new char[RefpropFluidPathLength]; // the full address of the fluid file
        private string hfmix = "HMX.BNC";
        private string hrf = "DEF";

        //private bool IsPureFluid;
        //private bool IsMixedFluid;
        //private bool Reload;

        //private string FluidName;
        //private double[] x;//Fractions;

        protected double wm, tc, pc, Dc;	//  INFO variables  cached 
        //protected double wm, ttp, tnbp, tc, pc, Dc, Zc, acf, dip, Rgas;	//  INFO variables  cached 

        protected double temperature, pressure, density, quality, internalenergy, enthalpy, entropy, 
            cv, cp, speedofsound, rhol, rhov,densityofvaporphase, densityofliquidphase,hjt;

        //saturation variables
        protected double[] XLIQ = new double[MaxComosition];
        protected double[] XVAP = new double[MaxComosition];
        protected double[] XOVERALL = new double[MaxComosition];  //z

        public StringBuilder hfiles = new StringBuilder(5100);

        //public double MolecularWeight { get { return this.wm; } private set { } }
        //public double TriplePointTemperature { get { return this.ttp; } private set { } }
        //public double NormalBoilingPointTemperature { get { return this.tnbp; } private set { } }
        //public double CriticalTemperature { get { return this.tc; } private set { } }
        //public double CriticalPressure { get { return this.pc; } private set { } }
        //public double CriticalDensity { get { return this.Dc; } private set { } }
        //public double CriticalPointCompressibility { get { return this.Zc; } private set { } }
        //public double AccentricFactor { get { return this.acf; } private set { } }
        //public double DipoleMoment { get { return this.dip; } private set { } }
        //public double GasConstatnt_R { get { return this.Rgas; } private set { } }

        //c     rhol--molar density [mol/L] of saturated liquid
        //c     rhov--molar density [mol/L] of saturated vapor
        //c     xliq--liquid phase composition [array of mol frac]
        //c     xvap--vapor phase composition [array of mol frac]

        /*---------------------------------------------------------------------*/

        //for flash and other routines

        public Refprop()
        {

        }


        #region import Functions of the REFPROP PROGRAM

        //const string FluidsDirectory = "../../fluids/";
        //const string MixturesDirectory = "../../mixtures/";
        const string refpropDLL_path = "../../Refprop.dll";

        // subroutine SETUP (nc,hfiles,hfmix,hrf,ierr,herr)
        //c
        //c  define models and initialize arrays
        //c
        //c  A call to this routine is required.
        //c
        //c  inputs:
        //c       nc--number of components (1 for pure fluid) [integer]
        //c   hfiles--array of file names specifying fluid/mixture components
        //c           [character*255 variable] for each of the nc components;
        //c           e.g., :fluids:r134a.fld (Mac) or fluids\r134a.fld (DOS) or
        //c           [full_path]/fluids/r134a.fld (UNIX)
        //c    hfmix--mixture coefficients [character*255]
        //c           file name containing coefficients for mixture model,
        //c           if applicable
        //c           e.g.,  fluids\hmx.bnc
        //c      hrf--reference state for thermodynamic calculations [character*3]
        //c           'DEF':  default reference state as specified in fluid file
        //c                   is applied to each pure component
        //c           'NBP':  h,s = 0 at pure component normal boiling point(s)
        //c           'ASH':  h,s = 0 for sat liquid at -40 C (ASHRAE convention)
        //c           'IIR':  h = 200, s = 1.0 for sat liq at 0 C (IIR convention)
        //c           other choices are possible, but these require a separate
        //c           call to SETREF
        //c  outputs:
        //c     ierr--error flag:  0 = successful
        //c                      101 = error in opening file
        //c                      102 = error in file or premature end of file
        //c                     -103 = unknown model encountered in file
        //c                      104 = error in setup of model
        //c                      105 = specified model not found
        //c                      111 = error in opening mixture file
        //c                      112 = mixture file of wrong type
        //c                      114 = nc<>nc from setmod
        //c     herr--error string (character*255 variable if ierr<>0)
        //c     [fluid parameters, etc. returned via various common blocks]
        [DllImport(refpropDLL_path, EntryPoint = "SETUPdll", SetLastError = true)]
#if (setup_use_string)
        //public static extern void SETUPdll(ref Int32 NumberOfComponents, string HFILES, string HFMIX, string HRF, [In, Out] ref Int32 ierr, [In, Out] char[] HERR, Int32 l1, Int32 l2, Int32 l3, Int32 l4);
        public static extern void SETUPdll(ref Int32 NumberOfComponents, string HFILES, string HFMIX, string HRF, [In, Out] ref Int32 ierr,
            [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 l1, Int32 l2, Int32 l3, Int32 l4);
#else
        public static extern void SETUPdll(ref Int32 NumberOfComponents, char[] HFILES, char[] HFMIX, char[] HRF, [In, Out] ref Int32 ierr, [In, Out] char[] HERR, Int32 l1, Int32 l2, Int32 l3, Int32 l4);
#endif
        [DllImport(refpropDLL_path, EntryPoint = "SETMIXdll", SetLastError = true)]
        public static extern void SETMIXdll(string hmxnme, string hfmix, string hrf, ref Int32 ncc,
             [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder hfiles, double[] x, ref Int32 ierr,
             [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 l1, Int32 l2, Int32 l3, Int32 l4, Int32 l5);

        [DllImport(refpropDLL_path, EntryPoint = "INFOdll", SetLastError = true)]
        public static extern void INFOdll(ref Int32 icomp, ref double wmm, ref  double ttp, ref double tnbp, ref double tc, ref double pc, ref double dc, ref double zc, ref double acf, ref double dip, ref double rgas);

#region SATURATION-STATE SUBROUTINES
        [DllImport(refpropDLL_path, EntryPoint = "SATTdll", SetLastError = true)]
        public static extern void SATTdll(ref double TK, double[] X, ref Int32 KPH, [In, Out] ref double PkPa,
            [In, Out] ref double RHOF, [In, Out] ref double RHOG, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "SATPdll", SetLastError = true)]
        public static extern void SATPdll(ref double PkPa, double[] X, ref Int32 KPH, [In, Out] ref double TK,
            [In, Out] ref double RHOF, [In, Out] ref double RHOG, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        // subroutine SATD (rho,x,kph,kr,t,p,rhol,rhov,xliq,xvap,ierr,herr)
        //c
        //c  iterate for temperature and pressure given a density along the
        //c  saturation boundary and the composition
        //c
        //c  inputs:
        //c      rho--molar density [mol/L]
        //c        x--composition [array of mol frac]
        //c      kph--flag specifying desired root for multi-valued inputs
        //c           has meaning only for water at temperatures close to its triple point
        //c          -1 = return middle root (between 0 and 4 C)
        //c           1 = return highest temperature root (above 4 C)
        //c           3 = return lowest temperature root (along freezing line)
        //c  outputs:
        //c        t--temperature [K]
        //c        p--pressure [kPa]
        //c     rhol--molar density [mol/L] of saturated liquid
        //c     rhov--molar density [mol/L] of saturated vapor
        //c     xliq--liquid phase composition [array of mol frac]
        //c     xvap--vapor phase composition [array of mol frac]
        //c       kr--phase flag: 1 = input state is liquid
        //c                       2 = input state is vapor in equilibrium with liq
        //c                       3 = input state is liquid in equilibrium with solid
        //c                       4 = input state is vapor in equilibrium with solid
        //c     ierr--error flag:   0 = successful
        //c                         2 = D > Dmax
        //c                         8 = x out of range
        //c                        10 = D and x out of range
        //c                       160 = CRITP did not converge
        //c                       161 = SATD did not converge
        //c     herr--error string (character*255 variable if ierr<>0)
        //c
        //c  N.B. kr = 3,4 presently working only for pure components
        //c
        //c  either (rhol,xliq) or (rhov,xvap) will correspond to the input state
        //c  with the other pair corresponding to the other phase in equilibrium
        //c  with the input state
        [DllImport(refpropDLL_path, EntryPoint = "SATDdll", SetLastError = true)]
        public static extern void SATDdll(ref double rhom, double[] X, ref Int32 KPH, [In, Out] ref Int32 kr, [In, Out] ref double tk,
            [In, Out] ref double pk, [In, Out] ref double RHOL, [In, Out] ref double RHOV, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);  //255
        //subroutine SATH (h,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
        //c
        //c  iterate for temperature, pressure, and density given enthalpy along
        //c  the saturation boundary and the composition
        //c
        //c  inputs:
        //c        h--molar enthalpy [J/mol]
        //c        x--composition [array of mol frac]
        //c      kph--flag specifying desired root
        //c           0 = return all roots along the liquid-vapor line
        //c           1 = return only liquid VLE root
        //c           2 = return only vapor VLE roots
        //c           3 = return liquid SLE root (melting line)
        //c           4 = return vapor SVE root (sublimation line)
        //c  outputs:
        //c    nroot--number of roots.  Set to one for kph=1,3,4 if ierr=0
        //c       k1--phase of first root (1-liquid, 2-vapor, 3-melt, 4-subl)
        //c       t1--temperature of first root [K]
        //c       p1--pressure of first root [kPa]
        //c       d1--molar density of first root [mol/L]
        //c       k2--phase of second root (1-liquid, 2-vapor, 3-melt, 4-subl)
        //c       t2--temperature of second root [K]
        //c       p2--pressure of second root [kPa]
        //c       d2--molar density of second root [mol/L]
        //c     ierr--error flag:   0 = successful
        //c                         2 = h < hmin
        //c                         4 = h > hmax
        //c                         8 = h > htrp (for subl input)
        //c                       160 = CRITP did not converge
        //c                       161 = SATH did not converge for one root
        //c                       162 = SATH did not converge for both roots
        //c     herr--error string (character*255 variable if ierr<>0)
        //c
        //c  The second root is always set as the root in the vapor at temperatures
        //c  below the maximum enthalpy on the vapor saturation line.  If kph is
        //c  set to 2, and only one root is found in the vapor (this occurs when h<hcrit)
        //c  the state point will be placed in k2,t2,p2,d2.  If kph=0 and this situation
        //c  occurred, the first root (k1,t1,p1,d1) would be in the liquid (k1=1, k2=2).
        //c
        //c  N.B. kph = 3,4 presently working only for pure components
        [DllImport(refpropDLL_path, EntryPoint = "SATHdll", SetLastError = true)]
        public static extern void SATHdll(ref double hm, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);  //255

        //  subroutine SATE (e,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,ierr,herr)
        //c
        //c  iterate for temperature, pressure, and density given energy along
        //c  the saturation boundary and the composition
        //c
        //c  inputs:
        //c        e--molar energy [J/mol]
        //c        x--composition [array of mol frac]
        //c      kph--flag specifying desired root
        //c           0 = return all roots along the liquid-vapor line
        //c           1 = return only liquid VLE root
        //c           2 = return only vapor VLE roots
        //c           3 = return liquid SLE root (melting line)
        //c           4 = return vapor SVE root (sublimation line)
        //c  outputs:
        //c    see SATH for description of outputs
        [DllImport(refpropDLL_path, EntryPoint = "SATEdll", SetLastError = true)]
        public static extern void SATEdll(ref double Em, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);  //255
        // subroutine SATS (s,x,kph,nroot,k1,t1,p1,d1,k2,t2,p2,d2,
        //     &                 k3,t3,p3,d3,ierr,herr)
        //c
        //c  iterate for temperature, pressure, and density given an entropy along
        //c  the saturation boundary and the composition
        //c
        //c  inputs:
        //c        s--molar entropy [J/mol-K]
        //c        x--composition [array of mol frac]
        //c      kph--flag specifying desired root
        //c           0 = return all roots along the liquid-vapor line
        //c           1 = return only liquid VLE root
        //c           2 = return only vapor VLE roots
        //c           3 = return liquid SLE root (melting line)
        //c           4 = return vapor SVE root (sublimation line)
        //c  outputs:
        //c    nroot--number of roots.  Set to one for kph=1,3,4 if ierr=0
        //c       k1--phase of first root (1-liquid, 2-vapor, 3-melt, 4-subl)
        //c       t1--temperature of first root [K]
        //c       p1--pressure of first root [kPa]
        //c       dl--molar density of first root [mol/L]
        //c       k2--phase of second root (1-liquid, 2-vapor, 3-melt, 4-subl)
        //c       t2--temperature of second root [K]
        //c       p2--pressure of second root [kPa]
        //c       d2--molar density of second root [mol/L]
        //c       k3--phase of third root (1-liquid, 2-vapor, 3-melt, 4-subl)
        //c       t3--temperature of third root [K]
        //c       p3--pressure of third root [kPa]
        //c       d3--molar density of third root [mol/L]
        //c     ierr--error flag:   0 = successful
        //c                         2 = s < smin
        //c                         4 = s > smax
        //c                         8 = s > strp (for subl input)
        //c                       160 = CRITP did not converge
        //c                       161 = SATS did not converge for one root
        //c                       162 = SATS did not converge for two roots
        //c                       163 = SATS did not converge for all roots
        //c     herr--error string (character*255 variable if ierr<>0)
        //c
        //c  The second root is always set as the root in the vapor at temperatures
        //c  below the maximum entropy on the vapor saturation line.  If kph is
        //c  set to 2, and only one root is found in the vapor (this occurs when s<scrit)
        //c  the state point will be placed in k2,t2,p2,d2.  If kph=0 and this situation
        //c  occurred, the first root (k1,t1,p1,d1) would be in the liquid (k1=1, k2=2).
        //c
        //c  The third root is the root with the lowest temperature.  For fluids
        //c  with multiple roots:  When only one root is found in the vapor phase
        //c  (this happens only at very low temperatures past the region where three
        //c  roots are located), the value of the root is still placed in
        //c  k3,t3,p3,d3.  For fluids that never have more than one root (when there
        //c  is no maximum entropy along the saturated vapor line), the value of the
        //c  root is always placed in k1,t1,p1,d1.
        //c
        //c  N.B. kph = 3,4 presently working only for pure components
        [DllImport(refpropDLL_path, EntryPoint = "SATSdll", SetLastError = true)]
        public static extern void SATSdll(ref double Sm, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 K3, [In, Out] ref double TK3, [In, Out] ref double PK3, [In, Out] ref double DM3,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);  //255
#endregion

        //subroutine CSATK (icomp,t,kph,p,rho,csat,ierr,herr)
        //c
        //c  compute the heat capacity along the saturation line as a function of
        //c  temperature for a given component
        //c
        //c  csat can be calculated two different ways:
        //c     Csat = Cp - T(DvDT)(DPDTsat)
        //c     Csat = Cp - beta/rho*hvap/(vliq - vvap)
        //c     where beta is the volume expansivity
        //c
        //c  inputs:
        //c    icomp--component number in mixture (1..nc); 1 for pure fluid
        //c        t--temperature [K]
        //c      kph--phase flag: 1 = liquid calculation
        //c                       2 = vapor calculation
        //c  outputs:
        //c        p--saturation pressure [kPa]
        //c      rho--saturation molar density [mol/L]
        //c     csat--saturation heat capacity [J/mol-K]

        //[DllImport(refpropDLL_path, EntryPoint = "CSATKdll", SetLastError = true)]

        //subroutine DPTSATK (icomp,t,kph,p,rho,csat,dpt,ierr,herr)
        //c
        //c  compute the heat capacity and dP/dT along the saturation line as a
        //c  function of temperature for a given component.  See also subroutine CSATK.
        //c
        //c  inputs:
        //c    icomp--component number in mixture (1..nc); 1 for pure fluid
        //c        t--temperature [K]
        //c      kph--phase flag: 1 = liquid calculation
        //c                       2 = vapor calculation
        //c  outputs:
        //c        p--saturation pressure [kPa]
        //c      rho--saturation molar density [mol/L]
        //c     csat--saturation heat capacity [J/mol-K] (same as that called from CSATK)
        //c      dpt--dP/dT along the saturation line [kPa/K]
        //c           (this is not dP/dT "at" the saturation line for the single phase
        //c            state, but the change in saturated vapor pressure as the
        //c            saturation temperature changes.)

        //      subroutine CV2PK (icomp,t,rho,cv2p,csat,ierr,herr)
        //c
        //c  compute the isochoric heat capacity in the two phase (liquid+vapor)
        //c  region
        //c
        //c  inputs:
        //c    icomp--component number in mixture (1..nc); 1 for pure fluid
        //c        t--temperature [K]
        //c      rho--density [mol/l] if known
        //c           If rho=0, then a saturated liquid state is assumed.
        //c
        //c  outputs:
        //c     cv2p--isochoric two-phase heat capacity [J/mol-K]
        //c     csat--saturation heat capacity [J/mol-K]
        //c           (Although there is already a csat routine in REFPROP,
        //c            it is also returned here.  However, the calculation
        //c            speed is slower than csat.)

        //Two similar routines are provided for calculating surface tension.
        //SURTEN is more efficient if the liquid and vapor density and composition
        //are known (e.g., from a previous call to SATT).  If these are not known,
        //then SURFT may be used.

        //      subroutine SURFT (t,rhol,xl,sigma,ierr,herr)
        //c
        //c  compute surface tension
        //c
        //c  inputs:
        //c        t--temperature [K]
        //c       xl--composition of liquid phase [array of mol frac]
        //c  outputs:
        //c     rhol--molar density of liquid phase [mol/L]
        //c           if rho > 0 use as input value
        //c                  < 0 call SATT to find density
        //c    sigma--surface tension [N/m]
        //c     ierr--error flag:   0 = successful
        //c                         1 = T < Tmin
        //c                         8 = x out of range
        //c                         9 = T and x out of range
        //c                       120 = CRITP did not converge
        //c                       121 = T > Tcrit
        //c                       122 = TPRHO-liquid did not converge in SATT
        //c                       123 = TPRHO-vapor did not converge in SATT
        //c                       124 = SATT pure fluid iteration did not converge
        //c                       128 = SATT mixture iteration did not converge
        //c     herr--error string if ierr<>0 (character*255)
        [DllImport(refpropDLL_path, EntryPoint = "SURFT", SetLastError = true)]
        public static extern void SURFdll(ref double Tk, ref double rholm, double[] xl, ref double sigma,
            [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        //         subroutine SURTEN (t,rhol,rhov,xl,xv,sigma,ierr,herr)
        //c
        //c  compute surface tension
        //c
        //c  inputs:
        //c        t--temperature [K]
        //c     rhol--molar density of liquid phase [mol/L]
        //c     rhov--molar density of vapor phase [mol/L]
        //c           if either rhol or rhov < 0 call SATT to find densities
        //c       xl--composition of liquid phase [array of mol frac]
        //c       xv--composition of liquid phase [array of mol frac]
        //c           (xv is optional input if rhol < 0 or rhov < 0)
        //c  outputs:
        //c    sigma--surface tension [N/m]
        //c     ierr--error flag:   0 = successful
        //c                         1 = T < Tmin
        //c                         8 = x out of range
        //c                         9 = T and x out of range
        //c                       120 = CRITP did not converge
        //c                       121 = T > Tcrit
        //c                       122 = TPRHO-liquid did not converge in SATT
        //c                       123 = TPRHO-vapor did not converge in SATT
        //c                       124 = SATT pure fluid iteration did not converge
        //c                       128 = SATT mixture iteration did not converge
        //c     herr--error string if ierr<>0 (character*255)

        [DllImport(refpropDLL_path, EntryPoint = "SURTEN", SetLastError = true)]
        public static extern void SURFTENdll(ref double Tk, ref double rholm, ref double rhovm, double[] xl, double[] xv, ref double sigma,
            [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        //c  compute entropy as a function of temperature, density and composition
        //c  inputs:
        //c        t--temperature [K]
        //c      rho--molar density [mol/L]
        //c        x--composition [array of mol frac]
        //c  output:
        //c        s--entropy [J/mol-K]
        //         subroutine ENTRO (t,rho,x,s)
        //c
        //c  compute entropy as a function of temperature, density and composition
        //c  using core functions (temperature derivative of Helmholtz free energy
        //c  and ideal gas integrals)
        //c
        //c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
        //c  equations A5, A19 - A26
        //c
        //c  inputs:
        //c        t--temperature [K]
        //c      rho--molar density [mol/L]
        //c        x--composition [array of mol frac]
        //c  output:
        //c        s--entropy [J/mol-K]
        [DllImport(refpropDLL_path, EntryPoint = "ENTROdll", SetLastError = true)]
        public static extern void ENTROdll(ref double tk, ref double rho, double[] X, ref double s);
        //subroutine ENTHAL (t,rho,x,h)
        //c
        //c  compute enthalpy as a function of temperature, density, and
        //c  composition using core functions (Helmholtz free energy and ideal
        //c  gas integrals)
        //c
        //c  based on derivations in Younglove & McLinden, JPCRD 23 #5, 1994,
        //c  equations A7, A18, A19
        //c
        //c  inputs:
        //c        t--temperature [K]
        //c      rho--molar density [mol/L]
        //c        x--composition [array of mol frac]
        //c  output:
        //c        h--enthalpy [J/mol]
        [DllImport(refpropDLL_path, EntryPoint = "ENTHALdll", SetLastError = true)]
        public static extern void ENTHALdll(ref double tk, ref double rho, double[] X, ref double h);

        #region  GENERAL FLASH SUBROUTINES

        //For cases where the phase is not known, the following routines are available.

        //      subroutine TPFLSH (t,p,z,D,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine TDFLSH (t,D,z,p,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine THFLSH (t,h,z,kr,p,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        //      subroutine TSFLSH (t,s,z,kr,p,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        //      subroutine TEFLSH (t,e,z,kr,p,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
        //      subroutine PDFLSH (p,D,z,t,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine PHFLSH (p,h,z,t,D,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        //      subroutine PSFLSH (p,s,z,t,D,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        //      subroutine PEFLSH (p,e,z,t,D,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
        //      subroutine HSFLSH (h,s,z,t,p,D,Dl,Dv,x,y,q,e,cv,cp,w,ierr,herr)
        //      subroutine ESFLSH (e,s,z,t,p,D,Dl,Dv,x,y,q,h,cv,cp,w,ierr,herr)
        //      subroutine DHFLSH (D,h,z,t,p,Dl,Dv,x,y,q,e,s,cv,cp,w,ierr,herr)
        //      subroutine DSFLSH (D,s,z,t,p,Dl,Dv,x,y,q,e,h,cv,cp,w,ierr,herr)
        //      subroutine DEFLSH (D,e,z,t,p,Dl,Dv,x,y,q,h,s,cv,cp,w,ierr,herr)
        //      subroutine TQFLSH (t,q,z,kq,p,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        //      subroutine PQFLSH (p,q,z,kq,t,D,Dl,Dv,x,y,e,h,s,cv,cp,w,ierr,herr)
        //c
        //c  flash calculation given two independent variables and bulk composition
        //c
        //c  These routines accept both single-phase and two-phase states as the
        //c  input; if the phase is known, the specialized routines are faster
        //c
        //c  inputs--two of the following as indicated by the first two letters of
        //c          the subroutine name:
        //c        t--temperature [K]
        //c        p--pressure [kPa]
        //c        e--internal energy [J/mol]
        //c        h--enthalpy [J/mol]
        //c        s--entropy [[J/mol-K]
        //c        q--vapor quality [basis specified by kq]
        //c           q = 0 indicates saturated liquid
        //c           q = 1 indicates saturated vapor
        //c           q < 0 or q > 1 are not allowed and will result in warning
        //c
        //c  additional input--required for all routines
        //c        z--overall (bulk) composition [array of mol frac]
        //c
        //c  additional input--only for TQFLSH and PQFLSH
        //c       kq--flag specifying units for input quality
        //c           kq = 1 quality on MOLAR basis [moles vapor/total moles]
        //c           kq = 2 quality on MASS basis [mass vapor/total mass]
        //c
        //c  outputs--one, two, or all of the following, depending on the inputs:
        //c        t--temperature [K]
        //c        p--pressure [kPa]
        //c        D--overall (bulk) molar density [mol/L]
        //c
        //c  additional outputs--common to all routines
        //c       Dl--molar density [mol/L] of the liquid phase
        //c       Dv--molar density [mol/L] of the vapor phase
        //c           if only one phase is present, Dl = Dv = D
        //c        x--composition of liquid phase [array of mol frac]
        //c        y--composition of vapor phase [array of mol frac]
        //c           if only one phase is present, x = y = z
        //c
        //c  additional output--common to all routines except TQFLSH and PQFLSH
        //c        q--vapor quality on a MOLAR basis [moles vapor/total moles]
        //c           q < 0 indicates subcooled (compressed) liquid
        //c           q = 0 indicates saturated liquid
        //c           q = 1 indicates saturated vapor
        //c           q > 1 indicates superheated vapor
        //c           q = 998 superheated vapor, but quality not defined (t > Tc)
        //c           q = 999 indicates supercritical state (t > Tc) and (p > Pc)
        //c
        //c  additional outputs--common to all routines, except that input
        //c                      quantities are not repeated
        //c        e--overall (bulk) internal energy [J/mol]
        //c        h--overall (bulk) enthalpy [J/mol]
        //c        s--overall (bulk) entropy [J/mol-K]
        //c       cv--isochoric (constant V) heat capacity [J/mol-K]
        //c       cp--isobaric (constant p) heat capacity [J/mol-K]
        //c        w--speed of sound [m/s]
        //c           Cp, w are not defined for 2-phase states
        //c           in such cases, a flag = -9.99998d6 is returned
        //c     ierr--error flag:  0 = successful
        //c                        1 = Tin < Tmin
        //c                        4 = Pin < 0
        //c                        5 = T and P out of range
        //c                        8 = x out of range (component and/or sum < 0
        //c                            or > 1)
        //c                        9 = x and T out of range
        //c                       12 = x out of range and P < 0
        //c                       13 = x and T and P out of range
        //c     herr--error string (character*255 variable if ierr<>0)

        //subroutine TPFLSH (t,p,z,D,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        [DllImport(refpropDLL_path, EntryPoint = "TPFLSHdll", SetLastError = true)]
        public static extern void TPFLSHdll(ref double tk, ref double pk, double[] Z, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TDFLSHdll", SetLastError = true)]
        public static extern void TDFLSHdll(ref double tk, ref double D, double[] Z, ref double P, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "THFLSHdll", SetLastError = true)]
        public static extern void THFLSHdll(ref double tk, ref double H, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double e, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TSFLSHdll", SetLastError = true)]
        public static extern void TSFLSHdll(ref double tk, ref double S, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double e, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TEFLSHdll", SetLastError = true)]
        public static extern void TEFLSHdll(ref double tk, ref double E, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PDFLSHdll", SetLastError = true)]
        public static extern void PDFLSHdll(ref double Pk, ref double D, double[] Z, ref double Tk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PHFLSHdll", SetLastError = true)]
        public static extern void PHFLSHdll(ref double Pk, ref double H, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PSFLSHdll", SetLastError = true)]
        public static extern void PSFLSHdll(ref double Pk, ref double S, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PEFLSHdll", SetLastError = true)]
        public static extern void PEFLSHdll(ref double Pk, ref double E, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "HSFLSHdll", SetLastError = true)]
        public static extern void HSFLSHdll(ref double H, ref double S, double[] Z, ref double Tk, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "ESFLSHdll", SetLastError = true)]
        public static extern void ESFLSHdll(ref double E, ref double S, double[] Z, ref double Tk, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DHFLSHdll", SetLastError = true)]
        public static extern void DHFLSHdll(ref double D, ref double H, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double S, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DSFLSHdll", SetLastError = true)]
        public static extern void DSFLSHdll(ref double D, ref double S, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DEFLSHdll", SetLastError = true)]
        public static extern void DEFLSHdll(ref double D, ref double E, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TQFLSHdll", SetLastError = true)]
        public static extern void TQFLSHdll(ref double tk, ref double Q, double[] Z, ref Int32 kq, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PQFLSHdll", SetLastError = true)]
        public static extern void PQFLSHdll(ref double Pk, ref double Q, double[] Z, ref Int32 kq, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        #endregion

        //STRANGE!!! WOML IS A FUNCTION IN FORTRAN, EXPRESSED AS "function WMOL (x)"
        // BUT HERE IT IS NOT IMPORTED AS "public static extern double WMOLdll(double[] X);"
        // BUT AS THE FOLLOWING:
        [DllImport(refpropDLL_path, EntryPoint = "WMOLdll", SetLastError = true)]
        public static extern void WMOLdll(double[] X, ref double wm);
        //converts composition on a mole fraction basis to mass fraction
        [DllImport(refpropDLL_path, EntryPoint = "XMASSdll", SetLastError = true)]
        public static extern void XMASSdll(double[] xmol, double[] xkg,ref double wmix);
        //converts composition on a mass fraction basis to mole fraction
        [DllImport(refpropDLL_path, EntryPoint = "XMOLEdll", SetLastError = true)]
        public static extern void XMOLEdll(double[] xkg, double[] xmol, ref double wmix);
        
        [DllImport(refpropDLL_path, EntryPoint = "CRITPdll", SetLastError = true)]
        public static extern void CRITPdll(double[]xmol,ref double crit,ref double pcrit,ref double Dcrit,
            ref Int32 ierr,[Out, MarshalAs(UnmanagedType.LPStr)]StringBuilder HERR, Int32 LengthHERR);
        //[DllImport(refpropDLL_path, EntryPoint = "PSFL1dll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "TPRHOdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "PRESSdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "GIBBSdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "CVCPdll", SetLastError = true)]
//c  inputs:
//c        t--temperature [K]
//c      rho--molar density [mol/L]
//c        x--composition [array of mol frac]
//c  outputs:
//c        p--pressure [kPa]
//c        e--internal energy [J/mol]
//c        h--enthalpy [J/mol]
//c        s--entropy [J/mol-K]
//c       Cv--isochoric heat capacity [J/mol-K]
//c       Cp--isobaric heat capacity [J/mol-K]
//c        w--speed of sound [m/s]
//c      hjt--isenthalpic Joule-Thompson coefficient [K/kPa]
        [DllImport(refpropDLL_path, EntryPoint = "THERMdll", SetLastError = true)]
        public static extern void THERMdll(ref double tk, ref double rho, double[] x, ref double pk, ref double e, 
            ref double h, ref double s, ref double cv, ref double cp, ref double w, ref double hjt);
        //[DllImport(refpropDLL_path, EntryPoint = "THERM2dll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "TRNPRPdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DHDTdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DPDDdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DPDD2dll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DPDTdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DDDPdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "DDDTdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "FGCTYdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "LIMITXdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "SETKTVdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "GETKTVdll", SetLastError = true)]
        //[DllImport(refpropDLL_path, EntryPoint = "GETFIJdll", SetLastError = true)]
        #endregion


        //public void LoadFluid(string FileName, Int32 NumberOfComponents, double[] Fraction)
        /// <summary>
        /// Load a fluid with name, number of components,and the mass composition array
        /// </summary>
        /// <param name="FileName"></param>
        /// <param name="NumberOfComponents"></param>
        /// <param name="mole_x">mole composition of the refrigerant</param>
        protected void LoadFluid(string FileName, Int32 NumberOfComponents, double[] mole_x,string hrf)
        {
            try
            {
                this.MoleFractions = mole_x;
                this.nc = NumberOfComponents;
                this.hrf = hrf;
                //XMOLEdll(mass_x, MoleFractions, ref this.wm);
                //SETUPdll(ref nc, "R32.fld|R125.fld|", hfmix, "def", ref ierr, herr, 5100, 255, 3, 255);
                SETUPdll(ref NumberOfComponents, FileName, hfmix, hrf, ref ierr, herr, 5100, 255, 3, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
                CRITPdll(MoleFractions, ref this.tc, ref this.pc, ref this.Dc, ref ierr, herr, 255);
                //WMOLdll(MoleFractions, ref this.wm);
                XMASSdll(MoleFractions, MassFractions, ref this.wm);
                //INFOdll(ref nc, ref wm, ref ttp, ref tnbp, ref tc, ref pc, ref Dc, ref Zc, ref acf, ref dip, ref Rgas);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="FileName">name of the mixture, e.g.: "R410A.mix"</param>
        /// <param name="hfiles">array of file names specifying mixture components that were used to call setup.</param>
        public void LoadPredefinedFluid(string FileName,string hrf)
        {
            this.hrf = hrf;
            try
            {
                SETMIXdll(FileName, hfmix, hrf, ref nc, hfiles, MoleFractions, ref ierr, herr, 255, 255, 3, hfiles.Length, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
                //INFOdll(ref nc, ref wm, ref ttp, ref tnbp, ref tc, ref pc, ref Dc, ref Zc, ref acf, ref dip, ref Rgas);
                //WMOLdll(MoleFractions, ref wm);
                CRITPdll(MoleFractions, ref this.tc, ref this.pc, ref this.Dc, ref ierr, herr, 255);
                XMASSdll(MoleFractions, MassFractions, ref this.wm);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="temperature_value"></param>
        /// <param name="PhaseFlag"></param>
        protected void GetSaturatedPressure(double temperature_value, Int32 PhaseFlag)
        {
            this.temperature = temperature_value;
            try
            {
                //SATPdll(ref temperature, x, ref PhaseFlag, ref SaturatedPressure, ref rhol, ref rhov, XLIQ, XVAP, ref ierr, herr, 255);
                SATTdll(ref temperature_value, MoleFractions, ref PhaseFlag, ref pressure, ref rhol, ref rhov, XLIQ, XVAP, ref ierr, herr, 255);
                if (PhaseFlag == (Int32)SaturationPoint.Bubble_Point)
                {
                    THERMdll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                    //this.density = rhol;
                }
                else if (PhaseFlag == (Int32)SaturationPoint.Dew_Point)
                {
                    THERMdll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                    this.density = rhov;
                }
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="pressure_value"></param>
        /// <param name="PhaseFlag"></param>
        protected void GetSaturatedTemperature(double pressure_value, Int32 PhaseFlag)
        {
            this.pressure = pressure_value;
            try
            {
                SATPdll(ref pressure_value, MoleFractions, ref PhaseFlag, ref temperature, ref rhol, ref rhov, XLIQ, XVAP, ref ierr, herr, 255);
                if (PhaseFlag == (Int32)SaturationPoint.Bubble_Point)
                    THERMdll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                else if (PhaseFlag == (Int32)SaturationPoint.Dew_Point)
                    THERMdll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="i"></param>
        /// <param name="wm0"></param>
        /// <param name="ttp0"></param>
        /// <param name="tnbp0"></param>
        /// <param name="tc0"></param>
        /// <param name="pc0"></param>
        /// <param name="Dc0"></param>
        /// <param name="Zc0"></param>
        /// <param name="acf0"></param>
        /// <param name="dip0"></param>
        /// <param name="Rgas0"></param>
        protected void GetComponentInfo(int i, ref double wm0, ref double ttp0, ref double tnbp0, ref double tc0, 
            ref double pc0, ref double Dc0, ref double Zc0, ref double acf0, ref double dip0, ref double Rgas0)
        {
            try
            {
                INFOdll(ref i, ref wm0, ref ttp0, ref tnbp0, ref tc0, ref pc0, ref Dc0, ref Zc0, ref acf0, ref dip0, ref Rgas0);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }


        /// <summary>
        /// Calculate Thermdynamic values with Pressure and Vapor Quality
        /// </summary>
        /// <param name="press">pressure(kPa)</param>
        /// <param name="quality">Vapor Quality(basis specified by kq)</param>
        /// <param name="kq">flag specifying units for input quality
        /// kq = 1 quality on MOLAR basis [moles vapor/total moles]
        /// kq = 2 quality on MASS basis [mass vapor/total mass]</param>
        protected void PQFLSH(double press, double quality, Int32 kq)
        {
            this.pressure = press;
            this.quality = quality;
            try
            {
                TQFLSHdll(ref pressure, ref quality, MoleFractions, ref kq, ref temperature, ref density, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref internalenergy, ref entropy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        /// <summary>
        ///  Calculate Thermdynamic values with Temperature and Vapor Quality
        /// </summary>
        /// <param name="temp">temperature(K)</param>
        /// <param name="quality">Vapor Quality(basis specified by kq)</param>
        /// <param name="kq">flag specifying units for input quality
        /// kq = 1 quality on MOLAR basis [moles vapor/total moles]
        /// kq = 2 quality on MASS basis [mass vapor/total mass]</param>
        protected void TQFLSH(double temp, double quality, Int32 kq)
        {
            this.temperature = temp;
            this.quality = quality;
            try
            {
                TQFLSHdll(ref temperature, ref quality, MoleFractions, ref kq, ref pressure, ref density,ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref internalenergy, ref entropy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        /// <summary>
        /// Calculate Thermdynamic values with Density and Internal Energy
        /// </summary>
        /// <param name="densi">overall (bulk) molar density (mol/L)</param>
        /// <param name="energy">Internal Energy(J/mol)</param>
        protected void DEFLSH(double densi, double energy)
        {
            this.density = densi;
            this.internalenergy = energy;
            try
            {
                DEFLSHdll(ref density, ref internalenergy, MoleFractions, ref temperature, ref pressure, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref entropy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }

        /// <summary>
        /// Calculate Thermdynamic values with Density and Entropy
        /// </summary>
        /// <param name="densi">overall (bulk) molar density [mol/L]</param>
        /// <param name="entro">Entropy(J/mol-K)</param>
        protected void DSFLSH(double densi, double entro)
        {
            this.density = densi;
            this.entropy = entro;
            try
            {
                DSFLSHdll(ref density, ref entropy, MoleFractions, ref temperature, ref pressure, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy,ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        
        /// <summary>
        /// Calculate Thermdynamic values with Density and Enthalhy
        /// </summary>
        /// <param name="densi">overall (bulk) molar density [mol/L]</param>
        /// <param name="entha">Enthapy(J/mol-K)</param>
        protected void DHFLSH(double densi, double entha)
        {
            this.density = densi;
            this.enthalpy = entha;
            try
            {
                DHFLSHdll(ref density,  ref enthalpy,MoleFractions, ref temperature, ref pressure, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }

        /// <summary>
        /// Calculate Thermdynamic values with InternalEnergy and Entropy
        /// </summary>
        /// <param name="energy">Internal Energy(J/mol)</param>
        /// <param name="entro">Entropy(J/mol-K)</param>
        protected void ESFLSH(double energy, double entro)
        {
            this.internalenergy = energy;
            this.entropy = entro;
            try
            {
                ESFLSHdll(ref internalenergy, ref entropy, MoleFractions, ref temperature, ref pressure, ref density, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        /// <summary>
        /// Calculate Thermdynamic values with Enthalpy and Entropy
        /// </summary>
        /// <param name="entha">Enthalpy(J/mol)</param>
        /// <param name="entro">Entropy(J/mol-K)</param>
        protected void HSFLSH(double entha, double entro)
        {
            this.enthalpy = entha;
            this.entropy = entro;
            try
            {
                HSFLSHdll(ref enthalpy, ref entropy, MoleFractions, ref temperature, ref pressure, ref density, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);

                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
        }
        /// <summary>
        /// Calculate Thermdynamic values with Pressure and Internal Energy
        /// </summary>
        /// <param name="press">Pressure(kPa)</param>
        /// <param name="energ">Internal Energy(J/mol)</param>
        protected void PEFLSH(double press, double energ)
        {
            this.pressure = press;
            this.internalenergy = energ;
            try
            {
                PEFLSHdll(ref pressure, ref internalenergy, MoleFractions, ref temperature, ref density, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
         /// <summary>
         /// Calculate Thermdynamic values with Pressure and Entropy
         /// </summary>
         /// <param name="press">Pressure(kPa)</param>
         /// <param name="entro">Entropy(J/mol-k)</param>
        protected void PSFLSH(double press, double entro)
        {
            this.pressure = press;
            this.entropy = entro;
            try
            {
                PSFLSHdll(ref pressure, ref entropy, MoleFractions, ref temperature, ref density, ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
        /// <summary>
        /// Calculate Thermdynamic values with Pressure and Enthalpy
        /// </summary>
        /// <param name="press">Pressure(kPa)</param>
        /// <param name="enth">Enthalpy(J/mol)</param>
        protected void PHFLSH(double press, double enth)
        {
            this.pressure = press;
            this.enthalpy = enth;
            try
            {
                PHFLSHdll(ref pressure, ref enthalpy, MoleFractions, ref temperature, ref density,ref densityofliquidphase, ref densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);        
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
        /// <summary>
        /// Calculate Thermdynamic values with Pressure and Density
        /// </summary>
        /// <param name="press">Pressure(kPa)</param>
        /// <param name="den">overall (bulk) molar density [mol/L]</param>
        protected void PDFLSH(double press, double den)
        {
            this.pressure = press;
            this.density = den;
            try
            {
                PDFLSHdll(ref pressure, ref density, MoleFractions, ref temperature, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
        
        /// <summary>
        /// Calculate Thermdynamic values with Pressure and Temperature
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="pres">Pressure(kPa)</param>
        protected void TPFLSH(double temp, double pres)
        {
            this.temperature = temp;
            this.pressure = pres;
            try
            {
                TPFLSHdll(ref temperature, ref pressure, MoleFractions, ref density, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }
            
        }
        //TDFLSH (t,D,z,p,Dl,Dv,x,y,q,e,h,s,cv,cp,w,ierr,herr)
        /// <summary>
        /// Calculate Thermdynamic values with Temperature and Density 
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="dens">overall (bulk) molar density [mol/L](</param>
        protected void TDFLSH(double temp, double dens)
        {
            this.temperature = temp;
            this.density = dens;
            try
            {
                TDFLSHdll(ref temperature, ref density, MoleFractions, ref pressure, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
        /// <summary>
        /// Calculate Thermdynamic values with Temperature and Enthalpy
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="h">Enthalpy(J/mol)</param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        protected void THFLSH(double temp, double h,int root)
        {
            this.temperature = temp;
            this.enthalpy = h;
            try
            {
                THFLSHdll(ref temperature, ref enthalpy, MoleFractions, ref root, ref pressure, ref density, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }

        /// <summary>
        /// Calculate Thermdynamic values with Temperature and Entropy
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="ss">Entropy(J/mol-K)</param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        protected void TSFLSH(double temp, double ss, int root)
        {
            this.temperature = temp;
            this.entropy = ss;
            try
            {
                THFLSHdll(ref temperature, ref entropy, MoleFractions, ref root, ref pressure, ref density, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref internalenergy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }

        /// <summary>
        /// Calculate Thermdynamic values with Temperature and Internal Energy
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="ee">Internal Energy(J/mol)</param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        protected void TEFLSH(double temp, double ee, int root)
        {
            this.temperature = temp;
            this.internalenergy = ee;
            try
            {
                THFLSHdll(ref temperature, ref internalenergy, MoleFractions, ref root, ref pressure, ref density, ref densityofliquidphase, ref this.densityofvaporphase,
                            XLIQ, XVAP, ref quality, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr, 255);
                if (ierr > 0) throw new Exception(ierr.ToString());
            }
            catch (Exception e)
            {
                ErrorDisplay.DisplayError(this, e.ToString());
            }

        }
        /// <summary>
        /// Calculate 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="wm"></param>
        protected void WMOL(double[]x,ref double wm)
        {
            WMOLdll(MoleFractions,ref wm);
        }

        public void EnsureCurrentFluid()
        {
        }

    }
}
