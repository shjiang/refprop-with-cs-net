
/*Refrigerants.cs
 * 
 * ---2013年1月5日-----
 * (1)偶尔出现“the value of the file specifier in an open statement is invalid  unit=12”错误,未解决
 * (2)引入REFPROP中的计算单相状态的函数，以提高计算速度。
 * ---2013年1月6日----
 * (1)“the value of the file specifier in an open statement is invalid  unit=12”错误，通过设置文件路径解决.
 * ---2014年5月7日----
 * (1)发现Thermdll Therm2dll计算出来的压力与原来的压力不符。原因未知。去掉所有关于它们的计算。
 * 计算的例子：
 * Refrigerant RTest = new Refrigerant(RefrigerantCategory.PredefinedMixture, "R410A", ReferenceState.DEF, UnitSystem.Refprop);
 * RTest.FindStateWithPQ(1600, 0.7, 2);
 * ---2014年9月---
 * (1)发现在英文系统下用StringBuilder，在载入混合制冷剂的时候会出问题，原因未知。改为byte[]数组后问题解决
 * (2)找到Thermdll Therm2dll的错误原因在于不能在两相区进行计算。因此改为判断是否处于两相区进行计算。
 * (3)9.1的dll文件计算速度快很多，但是在调试模式下(64位系统)关闭窗口的时候会出错。
 * ---2014年10月---
 * (1)修正了PQFlsh和TQFlsh中计算组分的错误。
 */

//#define Win64
using System;
//using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.Diagnostics;
//using System.Collections;
namespace RefpropCSNET
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
    public abstract class Refprop
    {
        const string FluidsDirectory = @"../fluids/";
        const string MixturesDirectory = @"../mixtures/";
#if(Win64)                
            const string refpropDLL_path = @"../REFPRP64.dll";
#else
        const string refpropDLL_path = @"../Refprop.dll";
#endif
        private const Int32 MaxComosition = 20;
        private const Int32 RefpropFluidPathLength = (RefpropCharLength + 1) * MaxComosition; //+1 for '|' between file names
        private const Int32 RefpropCharLength = 255;
        private const Int32 FilePathLength = 255;
        private const Int32 LengthOfReference = 3;
        private const Int32 ErrorMessageLength = 255;

        public Int32 ierr;

        private byte[] herr_c = new byte[255];
        private byte[] hfiles_c = new byte[10200];
        protected string herr { get { return System.Text.Encoding.ASCII.GetString(herr_c); } set { ;} }
        protected string hfiles { get { return System.Text.Encoding.ASCII.GetString(hfiles_c); } set { ;} }

        protected Int32 nc;   // number of components
        protected double[] MoleFractions = new double[MaxComosition]; // this is array of components  x[nc]
        protected double[] MassFractions = new double[MaxComosition];

        private char[] hfld = new char[RefpropFluidPathLength]; // the full address of the fluid file
        private string hfmix = FluidsDirectory + "HMX.BNC";
        protected string hrf = "DEF";

        protected double wm, tc, pc, Dc;	//  INFO variables  cached 
        //protected double wm, ttp, tnbp, tc, pc, Dc, Zc, acf, dip, Rgas;	//  INFO variables  cached 

        protected double temperature, pressure, density, molequality, massquality, internalenergy, enthalpy, entropy,
            cv, cp, speedofsound, rhov, rhol, z, hjt, eta, tcx, wliq, wvap, wmix,
            A, G, xkappa, beta, dPdD, d2PdD2, dPdT, dDdT, dDdP, d2PT2, d2PdTD, spare1, spare2, spare3, spare4,
            xisenk, xkt, betas, bs, xkkt, thrott, pi, spht;   //thermal3
        //protected int NRoots;
        protected double[] xlkg = new double[MaxComosition];
        protected double[] xvkg = new double[MaxComosition];
        //saturation variables
        protected double[] XLIQ = new double[MaxComosition];
        protected double[] XVAP = new double[MaxComosition];
        protected double[] XOVERALL = new double[MaxComosition];  //z

        protected bool calThermalQuantities = true;
        //c     rhol--molar density [mol/L] of saturated liquid
        //c     rhov--molar density [mol/L] of saturated vapor
        //c     xliq--liquid phase composition [array of mol frac]
        //c     xvap--vapor phase composition [array of mol frac]

        public Refprop()
        {

        }

        #region import Functions of the REFPROP PROGRAM
        [DllImport(refpropDLL_path, EntryPoint = "SETUPdll", CallingConvention = CallingConvention.StdCall, SetLastError = true)]
        public static extern void SETUPdll(ref Int32 NumberOfComponents, string HFILES, string HFMIX, string HRF, [In, Out] ref Int32 ierr,
            [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 l1, Int32 l2, Int32 l3, Int32 l4);
        [DllImport(refpropDLL_path, EntryPoint = "SETMIXdll", CallingConvention = CallingConvention.StdCall, SetLastError = true, CharSet = CharSet.Ansi)]
        public static extern void SETMIXdll(string hmxnme, string hfmix, string hrf, ref Int32 ncc,
             [Out, MarshalAs(UnmanagedType.LPArray)]byte[] hfiles, double[] x, ref Int32 ierr,
             [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 l1, Int32 l2, Int32 l3, Int32 l4, Int32 l5);

        [DllImport(refpropDLL_path, EntryPoint = "SATSPLNdll", SetLastError = true)]
        public static extern void SATSPLNdll(double[] x, [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 l4);

        [DllImport(refpropDLL_path, EntryPoint = "INFOdll", SetLastError = true)]
        public static extern void INFOdll(ref Int32 icomp, ref double wmm, ref  double ttp, ref double tnbp, ref double tc, ref double pc, ref double dc, ref double zc, ref double acf, ref double dip, ref double rgas);
        [DllImport(refpropDLL_path, EntryPoint = "SETNCdll", SetLastError = true)]
        public static extern void SETNCdll(ref Int32 icomp);
        [DllImport(refpropDLL_path, EntryPoint = "PUREFLDdll", SetLastError = true)]
        public static extern void PUREFLDdll(ref Int32 icomp);




        //converts quality and composition on a mass basis to a molar basis
        //subroutine QMOLE (qkg,xlkg,xvkg,qmol,xl,xv,wliq,wvap,ierr,herr)
        [DllImport(refpropDLL_path, EntryPoint = "QMOLEdll", SetLastError = true)]
        public static extern void QMOLEdll(ref double qkg, double[] xlkg, double[] xvkg, ref double qmol, double[] xl, double[] xv, ref double wliq, ref double wvap,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);

        //converts quality and composition on a mole basis to a mass basis
        [DllImport(refpropDLL_path, EntryPoint = "QMASSdll", SetLastError = true)]
        public static extern void QMASSdll(ref double qmol, double[] xl, double[] xv, ref double qkg, double[] xlkg, double[] xvkg, ref double wliq, ref double wvap,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);

        protected void QMASS()
        {
            QMASSdll(ref molequality, XLIQ, XVAP, ref massquality, xlkg, xvkg, ref wliq, ref wvap, ref ierr, herr_c, 255);
        }
        //V8.0 无此函数
        protected void SetNC(int n)
        {
            SETNCdll(ref n);
        }

        #region SATURATION-STATE SUBROUTINES
        [DllImport(refpropDLL_path, EntryPoint = "SATTdll", SetLastError = true)]
        public static extern void SATTdll(ref double TK, double[] X, ref Int32 KPH, [In, Out] ref double PkPa,
            [In, Out] ref double RHOF, [In, Out] ref double RHOG, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "SATPdll", SetLastError = true)]
        public static extern void SATPdll(ref double PkPa, double[] X, ref Int32 KPH, [In, Out] ref double TK,
            [In, Out] ref double RHOF, [In, Out] ref double RHOG, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "SATDdll", SetLastError = true)]
        public static extern void SATDdll(ref double rhom, double[] X, ref Int32 KPH, [In, Out] ref Int32 kr, [In, Out] ref double tk,
            [In, Out] ref double pk, [In, Out] ref double RHOL, [In, Out] ref double RHOV, [Out] double[] XLIQ, [Out] double[] XVAP,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);  //255
        [DllImport(refpropDLL_path, EntryPoint = "SATHdll", SetLastError = true)]
        public static extern void SATHdll(ref double hm, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);  //255
        [DllImport(refpropDLL_path, EntryPoint = "SATEdll", SetLastError = true)]
        public static extern void SATEdll(ref double Em, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);  //255
        [DllImport(refpropDLL_path, EntryPoint = "SATSdll", SetLastError = true)]
        public static extern void SATSdll(ref double Sm, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
            [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
            [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
            [In, Out] ref Int32 K3, [In, Out] ref double TK3, [In, Out] ref double PK3, [In, Out] ref double DM3,
            [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);  //255
        #endregion

        #region THERMODYNAMIC PROPERTY SUBROUTINES AS F(T,rho,X)
        [DllImport(refpropDLL_path, EntryPoint = "SURFTdll", SetLastError = true)]
        public static extern void SURFdll(ref double Tk, ref double rholm, double[] xl, ref double sigma,
            [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "SURTENdll", SetLastError = true)]
        public static extern void SURFTENdll(ref double Tk, ref double rholm, ref double rhovm, double[] xl, double[] xv, ref double sigma,
            [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "ENTROdll", SetLastError = true)]
        public static extern void ENTROdll(ref double tk, ref double rho, double[] X, ref double s);
        [DllImport(refpropDLL_path, EntryPoint = "ENTHALdll", SetLastError = true)]
        public static extern void ENTHALdll(ref double tk, ref double rho, double[] X, ref double h);
        [DllImport(refpropDLL_path, EntryPoint = "CRITPdll", SetLastError = true)]
        public static extern void CRITPdll(double[] xmol, ref double crit, ref double pcrit, ref double Dcrit,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PRESSdll", SetLastError = true)]
        public static extern void PRESSdll(ref double t, ref double rho, double[] x, ref double p);
        //[DllImport(refpropDLL_path, EntryPoint = "GIBBSdll", SetLastError = true)]
        [DllImport(refpropDLL_path, EntryPoint = "CVCPdll", SetLastError = true)]
        public static extern void CVCPdll(ref double t, ref double rho, double[] x, ref double cv, ref double cp);
        [DllImport(refpropDLL_path, EntryPoint = "THERMdll", SetLastError = true)]
        public static extern void THERMdll(ref double tk, ref double rho, double[] x, ref double pk, ref double e,
            ref double h, ref double s, ref double cv, ref double cp, ref double w, ref double hjt);
        [DllImport(refpropDLL_path, EntryPoint = "THERM2dll", SetLastError = true)]
        public static extern void THERM2dll(ref double tk, ref double rho, double[] x,
            ref double pk, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref double Z, ref double hjt, ref double A, ref double G, ref double xkappa, ref double beta, ref double dPdD,
            ref double d2PdD2, ref double dPdT, ref double dDdT, ref double dDdP, ref double spare1, ref double spare2,
            ref double spare3, ref double spare4);
        [DllImport(refpropDLL_path, EntryPoint = "THERM3dll", SetLastError = true)]
        public static extern void THERM3dll(ref double tk, ref double rho, double[] x,
            ref double xkappa, ref double beta, ref double xisenk, ref double xkt, ref double betas, ref double bs,
            ref double xkkt, ref double thrott, ref double pi, ref double spht);
        [DllImport(refpropDLL_path, EntryPoint = "TRNPRPdll", SetLastError = true)]
        public static extern void TRNPRPdll(ref double tk, ref double rho, double[] x, ref double eta, ref double tex,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        protected void PressureCal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                PRESSdll(ref temperature, ref density, MoleFractions, ref pressure);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void EnthalpyCal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                ENTHALdll(ref temperature, ref density, MoleFractions, ref enthalpy);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void EntropyCal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                ENTROdll(ref temperature, ref density, MoleFractions, ref entropy);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void CVCPCal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                CVCPdll(ref temperature, ref density, MoleFractions, ref cv, ref cp);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void Thermal2Cal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                THERM2dll(ref temperature, ref density, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                        ref entropy, ref cv, ref cp, ref speedofsound, ref z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                        ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref spare1, ref spare2, ref spare3, ref spare4);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void Thermal3Cal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                THERM3dll(ref temperature, ref density, MoleFractions,
            ref xkappa, ref beta, ref xisenk, ref xkt, ref betas, ref bs, ref xkkt, ref thrott, ref pi, ref spht);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void TRNPRPCal(double t, double d, double[] x)
        {
            this.temperature = t;
            this.density = d;
            this.MoleFractions = x;
            try
            {
                TRNPRPdll(ref temperature, ref density, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);

            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        #endregion

        #region  GENERAL FLASH SUBROUTINES
        [DllImport(refpropDLL_path, EntryPoint = "TPFLSHdll", SetLastError = true)]
        public static extern void TPFLSHdll(ref double tk, ref double pk, double[] Z, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TDFLSHdll", SetLastError = true)]
        public static extern void TDFLSHdll(ref double tk, ref double D, double[] Z, ref double P, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "THFLSHdll", SetLastError = true)]
        public static extern void THFLSHdll(ref double tk, ref double H, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double e, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TSFLSHdll", SetLastError = true)]
        public static extern void TSFLSHdll(ref double tk, ref double S, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double e, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TEFLSHdll", SetLastError = true)]
        public static extern void TEFLSHdll(ref double tk, ref double E, double[] Z, ref Int32 kr, ref double P, ref double D, ref double Dl,
            ref double Dv, double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PDFLSHdll", SetLastError = true)]
        public static extern void PDFLSHdll(ref double Pk, ref double D, double[] Z, ref double Tk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PHFLSHdll", SetLastError = true)]
        public static extern void PHFLSHdll(ref double Pk, ref double H, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PSFLSHdll", SetLastError = true)]
        public static extern void PSFLSHdll(ref double Pk, ref double S, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PEFLSHdll", SetLastError = true)]
        public static extern void PEFLSHdll(ref double Pk, ref double E, double[] Z, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "HSFLSHdll", SetLastError = true)]
        public static extern void HSFLSHdll(ref double H, ref double S, double[] Z, ref double Tk, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "ESFLSHdll", SetLastError = true)]
        public static extern void ESFLSHdll(ref double E, ref double S, double[] Z, ref double Tk, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DHFLSHdll", SetLastError = true)]
        public static extern void DHFLSHdll(ref double D, ref double H, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double S, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DSFLSHdll", SetLastError = true)]
        public static extern void DSFLSHdll(ref double D, ref double S, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double E, ref double H, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DEFLSHdll", SetLastError = true)]
        public static extern void DEFLSHdll(ref double D, ref double E, double[] Z, ref double Tk, ref double Pk, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double q, ref double H, ref double S, ref double cv, ref double cp, ref double w,
           ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TQFLSHdll", SetLastError = true)]
        public static extern void TQFLSHdll(ref double tk, ref double Q, double[] Z, ref Int32 kq, ref double Pk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PQFLSHdll", SetLastError = true)]
        public static extern void PQFLSHdll(ref double Pk, ref double Q, double[] Z, ref Int32 kq, ref double Tk, ref double D, ref double Dl, ref double Dv,
            double[] x, double[] y, ref double e, ref double h, ref double s, ref double cv, ref double cp, ref double w,
            ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);

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
            //this.quality = quality;
            try
            {
                PQFLSHdll(ref pressure, ref quality, MoleFractions, ref kq, ref temperature, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (kq == 1)
                {
                    this.molequality = quality;
                    this.massquality = quality * rhov / density;
                }
                else if (kq == 2)
                {
                    //this.molequality = quality * density / rhov;
                    //this.massquality = quality;
                    XMASSdll(XLIQ, xlkg, ref wmix);
                    XMASSdll(XVAP, xvkg, ref wmix);
                    QMOLEdll(ref quality, xlkg, xvkg, ref molequality, XLIQ, XVAP, ref wliq, ref wvap, ref ierr, herr_c, 255);
                }
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                Debug.WriteLine("EEE:::" + this.molequality.ToString());
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
            //this.quality = quality;
            try
            {
                TQFLSHdll(ref temperature, ref quality, MoleFractions, ref kq, ref pressure, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (kq == 1)
                {
                    this.molequality = quality;
                    this.massquality = quality * rhov / density;
                }
                else if (kq == 2)
                {
                    //this.molequality = quality * density / rhov;
                    //this.massquality = quality;
                    XMASSdll(XLIQ, xlkg, ref wmix);
                    XMASSdll(XVAP, xvkg, ref wmix);
                    QMOLEdll(ref quality, xlkg, xvkg, ref molequality, XLIQ, XVAP, ref wliq, ref wvap, ref ierr, herr_c, 255);
                }
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                DEFLSHdll(ref density, ref internalenergy, MoleFractions, ref temperature, ref pressure, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                DSFLSHdll(ref density, ref entropy, MoleFractions, ref temperature, ref pressure, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                DHFLSHdll(ref density, ref enthalpy, MoleFractions, ref temperature, ref pressure, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                ESFLSHdll(ref internalenergy, ref entropy, MoleFractions, ref temperature, ref pressure, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                HSFLSHdll(ref enthalpy, ref entropy, MoleFractions, ref temperature, ref pressure, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                PEFLSHdll(ref pressure, ref internalenergy, MoleFractions, ref temperature, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                PSFLSHdll(ref pressure, ref entropy, MoleFractions, ref temperature, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                PHFLSHdll(ref pressure, ref enthalpy, MoleFractions, ref temperature, ref density, ref rhol, ref rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                PDFLSHdll(ref pressure, ref density, MoleFractions, ref temperature, ref rhol, ref this.rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                TPFLSHdll(ref temperature, ref pressure, MoleFractions,
                    ref density, ref rhol, ref this.rhov, XLIQ, XVAP, ref molequality,
                    ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                TDFLSHdll(ref temperature, ref density, MoleFractions, ref pressure, ref rhol, ref this.rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        /// <summary>
        /// Calculate Thermdynamic values with Temperature and Enthalpy
        /// </summary>
        /// <param name="temp">Temperature(K)</param>
        /// <param name="h">Enthalpy(J/mol)</param>
        /// <param name="root">Often in the liquid, two solutions exist, one of them in the two phase. 
        /// 1 = return lower density root;2 = return higher density root. Use 1 in most cases.</param>
        protected void THFLSH(double temp, double h, int root)
        {
            this.temperature = temp;
            this.enthalpy = h;
            try
            {
                THFLSHdll(ref temperature, ref enthalpy, MoleFractions, ref root, ref pressure, ref density, ref rhol, ref this.rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                THFLSHdll(ref temperature, ref entropy, MoleFractions, ref root, ref pressure, ref density, ref rhol, ref this.rhov,
                            XLIQ, XVAP, ref molequality, ref internalenergy, ref enthalpy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);

                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                THFLSHdll(ref temperature, ref internalenergy, MoleFractions, ref root, ref pressure, ref density, ref rhol, ref this.rhov,
                            XLIQ, XVAP, ref molequality, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound,
                            ref ierr, herr_c, 255);
                if (calThermalQuantities) CalThermalQuantities();
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }

        }
        protected void CalThermalQuantities()
        {
            if (molequality > 1 || molequality < 0)
                THERM2dll(ref temperature, ref density, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                        ref entropy, ref cv, ref cp, ref speedofsound, ref z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                        ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref d2PT2, ref d2PdTD, ref spare3, ref spare4);
            TRNPRPdll(ref temperature, ref density, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);
        }


        /// <summary>
        /// Calculate 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="wm"></param>

        #endregion


        #region SINGLE-PHASE FLASH SUBROUTINES
        //These routines accept only single-phase states as inputs.  They will be
        //faster than the corresponding general routines, but will fail if called
        //with an incorrect phase specification.  The phase-specific subroutines
        //also do not check limits, so may fail if called outside the range of the
        //equation of state.  The following single-phase routines are available.
        //      subroutine THFL1 (t,h,x,Dmin,Dmax,D,ierr,herr)
        //      subroutine TSFL1 (t,s,x,Dmin,Dmax,D,ierr,herr)
        //      subroutine TEFL1 (t,e,x,Dmin,Dmax,D,ierr,herr)
        //      subroutine PDFL1 (p,rho,x,t,ierr,herr)
        //      subroutine PHFL1 (p,h,x,kph,t,D,ierr,herr)
        //      subroutine PSFL1 (p,s,x,kph,t,D,ierr,herr)
        //      subroutine PEFL1 (p,e,x,kph,t,D,ierr,herr)
        //      subroutine HSFL1 (h,s,x,Dmin,Dmax,t,D,ierr,herr)
        //      subroutine DHFL1 (rho,h,x,t,ierr,herr)
        //      subroutine DSFL1 (rho,s,x,t,ierr,herr)
        //      subroutine DEFL1 (rho,e,x,t,ierr,herr)
        //c
        //c  inputs--two of the following as indicated by the first two letters of
        //c          the subroutine name:
        //c        t--temperature [K]
        //c        e--internal energy [J/mol]
        //c        h--enthalpy [J/mol]
        //c        s--entropy [[J/mol-K]
        //c      rho--molar density [mol/L]
        //c
        //c  additional inputs
        //c        x--overall (bulk) composition [array of mol frac]
        //c     Dmin--lower bound on density [mol/L]
        //c     Dmax--upper bound on density [mol/L]
        //c      kph--phase flag:  1 = liquid
        //c                        2 = vapor
        //c
        //c  outputs:
        //c        t--temperature [K] (present only for HSFL1)
        //c        D--molar density [mol/L]
        //c     ierr--error flag:  0 = successful
        //c     herr--error string (character*255 variable if ierr<>0)
        //The single-phase temperature-pressure flash is called many times by
        //other routines, and has been optimized for speed; it requires a specific
        //calling sequence.

        //      subroutine TPRHO (t,p,x,kph,kguess,rho,ierr,herr)
        //c
        //c  iterate for density as a function of temperature, pressure, and
        //c  composition for a specified phase
        //c
        //c***********************************************************************
        //c  WARNING:
        //c  Invalid densities will be returned for T & P outside range of validity,
        //c  i.e., pressure > melting pressure, pressure less than saturation
        //c  pressure for kph=1, etc.
        //c
        //c***********************************************************************
        //c  inputs:
        //c        t--temperature [K]
        //c        p--pressure [kPa]
        //c        x--composition [array of mol frac]
        //c      kph--phase flag:  1 = liquid
        //c                        2 = vapor
        //c                 N.B.:  0 = stable phase--NOT ALLOWED (use TPFLSH)
        //c                            (unless an initial guess is supplied for rho)
        //c                       -1 = force the search in the liquid phase
        //c                       -2 = force the search in the vapor phase
        //c   kguess--input flag:  1 = first guess for rho provided
        //c                        0 = no first guess provided
        //c      rho--first guess for molar density [mol/L], only if kguess = 1
        //c
        //c  outputs:
        //c      rho--molar density [mol/L]
        //c     ierr--error flag:  0 = successful
        //c                      200 = CRITP did not converge
        //c                      201 = illegal input (kph <= 0)
        //c                      202 = liquid-phase iteration did not converge
        //c                      203 = vapor-phase iteration did not converge
        //c     herr--error string (character*255 variable if ierr<>0)
        [DllImport(refpropDLL_path, EntryPoint = "THFL1dll", SetLastError = true)]
        public static extern void THFL1dll(ref double t, ref double h, double[] X, ref double Dmin, ref double Dmax,
            ref double wm, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TSFL1dll", SetLastError = true)]
        public static extern void TSFL1dll(ref double t, ref double s, double[] X, ref double Dmin, ref double Dmax,
            ref double wm, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "TEFL1dll", SetLastError = true)]
        public static extern void TEFL1dll(ref double t, ref double e, double[] X, ref double Dmin, ref double Dmax,
            ref double wm, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PDFL1dll", SetLastError = true)]
        public static extern void PDFL1dll(ref double p, ref double rho, double[] X,
            ref double t, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PHFL1dll", SetLastError = true)]
        public static extern void PHFL1dll(ref double p, ref double h, double[] X, ref Int32 kph,
            ref double t, ref double D, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PSFL1dll", SetLastError = true)]
        public static extern void PSFL1dll(ref double p, ref double s, double[] X, ref Int32 kph,
            ref double t, ref double D, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "PEFL1dll", SetLastError = true)]
        public static extern void PEFL1dll(ref double p, ref double e, double[] X, ref Int32 kph,
            ref double t, ref double D, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "HSFL1dll", SetLastError = true)]
        public static extern void HSFL1dll(ref double h, ref double s, double[] X, ref double Dmin, ref double Dmax,
            ref double t, ref double wm, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DHFL1dll", SetLastError = true)]
        public static extern void DHFL1dll(ref double rho, ref double h, double[] X,
            ref double t, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DSFL1dll", SetLastError = true)]
        public static extern void DSFL1dll(ref double rho, ref double s, double[] X,
            ref double t, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);
        [DllImport(refpropDLL_path, EntryPoint = "DEFL1dll", SetLastError = true)]
        public static extern void DEFL1dll(ref double rho, ref double e, double[] X,
            ref double t, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);

        [DllImport(refpropDLL_path, EntryPoint = "TPRHOdll", SetLastError = true)]
        public static extern void TPRHOdll(ref double t, ref double p, double[] X, ref Int32 kph, ref Int32 kguess,
            ref double rho, ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);

        /// <summary>
        /// iterate for density as a function of temperature, pressure, and composition for a specified phase
        /// </summary>
        /// <param name="t">温度:[K] (输入)</param>
        /// <param name="p">压力:[kPa](输入)</param>
        /// <param name="x">摩尔组分(输入)</param>
        /// <param name="kq">相态(输入)</param>
        protected void S_TPRHO(double t, double p, double[] x, Int32 kq)
        {
            this.temperature = t;
            this.pressure = p;
            this.MoleFractions = x;
            Int32 kguess = 0;
            TPRHOdll(ref temperature, ref pressure, MoleFractions, ref kq, ref kguess, ref density, ref ierr, herr_c, 255);
        }
        protected void S_THFL1(double t, double h)
        {
            this.temperature = t;
            this.enthalpy = h;
            double Dmin = 0.0;
            double Dmax = 10000.0;//????
            try
            {
                THFL1dll(ref temperature, ref enthalpy, MoleFractions, ref Dmin, ref Dmax, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_TSFL1(double t, double s)
        {
            this.temperature = t;
            this.entropy = s;
            double Dmin = 0.0;
            double Dmax = 10000.0;//????
            try
            {
                TSFL1dll(ref temperature, ref entropy, MoleFractions, ref Dmin, ref Dmax, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_TEFL1(double t, double ie)
        {
            this.temperature = t;
            this.internalenergy = ie;
            double Dmin = 0.0;
            double Dmax = 10000.0;//????
            try
            {
                TEFL1dll(ref temperature, ref internalenergy, MoleFractions, ref Dmin, ref Dmax, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_PDFL1(double p, double d)
        {
            this.pressure = p;
            this.density = d;
            try
            {
                PDFL1dll(ref pressure, ref density, MoleFractions, ref temperature, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }

        protected void S_PHFL1(double p, double h, Int32 kph)
        {
            this.pressure = p;
            this.enthalpy = h;
            try
            {
                PHFL1dll(ref pressure, ref enthalpy, MoleFractions, ref kph, ref temperature, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }

        protected void S_PSFL1(double p, double s, Int32 kph)
        {
            this.pressure = p;
            this.entropy = s;
            try
            {
                PSFL1dll(ref pressure, ref entropy, MoleFractions, ref kph, ref temperature, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_PEFL1(double p, double ie, Int32 kph)
        {
            this.pressure = p;
            this.internalenergy = ie;
            try
            {
                PEFL1dll(ref pressure, ref internalenergy, MoleFractions, ref kph, ref temperature, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_HSFL1(double h, double s, Int32 kph)
        {
            this.enthalpy = h;
            this.entropy = s;
            double Dmin = 0.0;
            double Dmax = 10000.0;//????
            try
            {
                HSFL1dll(ref enthalpy, ref entropy, MoleFractions, ref Dmin, ref Dmax, ref temperature, ref density, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_DHFL1(double d, double h)
        {
            this.enthalpy = h;
            this.density = d;
            try
            {
                DHFL1dll(ref density, ref enthalpy, MoleFractions, ref temperature, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_DSFL1(double d, double s)
        {
            this.entropy = s;
            this.density = d;
            try
            {
                DSFL1dll(ref density, ref entropy, MoleFractions, ref temperature, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        protected void S_DEFL1(double d, double ie)
        {
            this.internalenergy = ie;
            this.density = d;
            try
            {
                DEFL1dll(ref density, ref internalenergy, MoleFractions, ref temperature, ref ierr, herr_c, 255);
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        #endregion

        // WOML IS A FUNCTION IN FORTRAN, EXPRESSED AS "function WMOL (x)"
        // BUT HERE IT IS NOT IMPORTED AS "public static extern double WMOLdll(double[] X);"
        // BUT AS THE FOLLOWING:
        [DllImport(refpropDLL_path, EntryPoint = "WMOLdll", SetLastError = true)]
        public static extern void WMOLdll(double[] X, ref double wm);
        //converts composition on a mole fraction basis to mass fraction
        [DllImport(refpropDLL_path, EntryPoint = "XMASSdll", SetLastError = true)]
        public static extern void XMASSdll(double[] xmol, double[] xkg, ref double wmix);
        //converts composition on a mass fraction basis to mole fraction
        [DllImport(refpropDLL_path, EntryPoint = "XMOLEdll", SetLastError = true)]
        public static extern void XMOLEdll(double[] xkg, double[] xmol, ref double wmix);

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

        public double RefpropDLLVersionNumber()
        {
            int nc2 = -1;
            SETUPdll(ref nc2, "".PadRight(255), hfmix, hrf, ref ierr, herr_c, (Int32)10000, (Int32)255, (Int32)3, (Int32)255);
            double RefpropDLLVersionNumber = ierr;
            if (ierr == 900) RefpropDLLVersionNumber = 9;
            return RefpropDLLVersionNumber / 10000;
        }

        //public void LoadFluid(string FileName, Int32 NumberOfComponents, double[] Fraction)
        /// <summary>
        /// Load a fluid with name, number of components,and the mass composition array
        /// </summary>
        /// <param name="FileName"></param>
        /// <param name="NumberOfComponents"></param>
        /// <param name="mole_x">mole composition of the refrigerant</param>
        public void LoadFluid(string FileName, Int32 NumberOfComponents, double[] x, string hrf, bool isMass)
        {
            //Marshal the managed string to unmanaged memory.
            //IntPtr FileNamePointer = (IntPtr)Marshal.StringToHGlobalAnsi(FileName);
            Debug.WriteLine("FileName in LoadFluid(*,*,*,*):" + FileName);
            FileName = FluidsDirectory + FileName;
            FileName = FileName.PadRight(255);

            try
            {
                this.nc = NumberOfComponents;
                this.hrf = hrf;
                //XMOLEdll(mass_x, MoleFractions, ref this.wm);
                //SETUPdll(ref nc, "R32.fld|R125.fld|", hfmix, "def", ref ierr, herr_c, 5100, 255, 3, 255);
                //SETUPdll(ref NumberOfComponents, FileName, hfmix, hrf, ref ierr, herr_c, 
                //    FileName.Length,hfmix.Length, hrf.Length, herr.Length);
                SETUPdll(ref NumberOfComponents, FileName, hfmix, hrf, ref ierr, herr_c,
                    FileName.Length, hfmix.Length, hrf.Length, herr.Length);
                //SETUPdll(ref NumberOfComponents, FileName, hfmix, hrf, ref ierr, herr_c, 5100, 255, 3, 255);
                if (ierr > 0) throw new Exception(herr.ToString());
                if (isMass)
                {
                    this.MassFractions = x;
                    XMOLEdll(MassFractions, MoleFractions, ref wm);
                }
                else
                {
                    this.MoleFractions = x;
                    XMASSdll(MoleFractions, MassFractions, ref wm);
                }
                //WMOLdll(MoleFractions, ref wm);
                //SATSPLNdll(x, ref ierr, herr_c, 255);  //有问题
                CRITPdll(MoleFractions, ref this.tc, ref this.pc, ref this.Dc, ref ierr, herr_c, 255);
                //THERMdll(ref tc, ref Dc, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                //WMOLdll(MoleFractions, ref this.wm);                
                //INFOdll(ref nc, ref wm, ref ttp, ref tnbp, ref tc, ref pc, ref Dc, ref Zc, ref acf, ref dip, ref Rgas);
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="FileName">name of the mixture, e.g.: "R410A.mix"</param>
        /// <param name="hfiles">array of file names specifying mixture components that were used to call setup.</param>
        public void LoadPredefinedFluid(string FileName, string hrf)
        {
            FileName = MixturesDirectory + FileName;
            FileName = FileName.PadRight(255);
            //Debug.WriteLine("FileName in LoadPredefinedFluid(*,*):"+"\n" + FileName +"\n"+ hfmix);
            this.hrf = hrf;
            //IntPtr FileNamePointer = (IntPtr)Marshal.StringToHGlobalAnsi(FileName);
            try
            {
                SETMIXdll(FileName, hfmix, hrf, ref nc, hfiles_c, MoleFractions, ref ierr, herr_c,
                    FileName.Length, hfmix.Length, 3, hfiles.Length, herr.Length);
                //SATSPLNdll(MoleFractions, ref ierr, herr_c,255); //有问题
                if (ierr > 0) throw new Exception(herr.ToString());
                //INFOdll(ref nc, ref wm, ref ttp, ref tnbp, ref tc, ref pc, ref Dc, ref Zc, ref acf, ref dip, ref Rgas);
                //WMOLdll(MoleFractions, ref wm);
                CRITPdll(MoleFractions, ref this.tc, ref this.pc, ref this.Dc, ref ierr, herr_c, 255);
                XMASSdll(MoleFractions, MassFractions, ref this.wm);
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                SATTdll(ref temperature_value, MoleFractions, ref PhaseFlag, ref pressure, ref rhol, ref rhov, XLIQ, XVAP, ref ierr, herr_c, 255);
                //if (PhaseFlag == (Int32)SaturationPoint.Bubble_Point)
                //{
                //    //THERMdll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                //    THERM2dll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                //        ref entropy, ref cv, ref cp, ref speedofsound, ref Z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                //        ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref spare1, ref spare2, ref spare3, ref spare4);
                //    this.density = rhol;
                //    TRNPRPdll(ref temperature, ref rhol, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);
                //}
                //else if (PhaseFlag == (Int32)SaturationPoint.Dew_Point)
                //{
                //    //THERMdll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                //    THERM2dll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                //        ref entropy, ref cv, ref cp, ref speedofsound, ref Z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                //        ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref spare1, ref spare2, ref spare3, ref spare4);
                //    this.density = rhov;
                //    TRNPRPdll(ref temperature, ref rhov, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);
                //}
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
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
                SATPdll(ref pressure_value, MoleFractions, ref PhaseFlag, ref temperature, ref rhol, ref rhov, XLIQ, XVAP, ref ierr, herr_c, 255);
                //if (PhaseFlag == (Int32)SaturationPoint.Bubble_Point)
                //{
                //    //THERMdll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                //    THERM2dll(ref temperature, ref rhol, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                //                         ref entropy, ref cv, ref cp, ref speedofsound, ref Z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                //                         ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref spare1, ref spare2, ref spare3, ref spare4);
                //    this.density = rhol;
                //    TRNPRPdll(ref temperature, ref density, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);
                //}
                //else if (PhaseFlag == (Int32)SaturationPoint.Dew_Point)
                //{
                //    //THERMdll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy, ref entropy, ref cv, ref cp, ref speedofsound, ref hjt);
                //    THERM2dll(ref temperature, ref rhov, MoleFractions, ref pressure, ref internalenergy, ref enthalpy,
                //            ref entropy, ref cv, ref cp, ref speedofsound, ref Z, ref hjt, ref A, ref G, ref xkappa, ref beta, ref dPdD,
                //            ref d2PdD2, ref dPdT, ref dDdT, ref  dDdP, ref spare1, ref spare2, ref spare3, ref spare4);
                //    this.density = rhov;
                //    TRNPRPdll(ref temperature, ref density, MoleFractions, ref eta, ref tcx, ref ierr, herr_c, 255);
                //}
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr);
            }
        }
        //public static extern void SATHdll(ref double hm, double[] X, ref Int32 KPH, [In, Out] ref Int32 NROOT,
        //  [In, Out] ref Int32 K1, [In, Out] ref double TK1, [In, Out] ref double PK1, [In, Out] ref double DM1,
        //  [In, Out] ref Int32 K2, [In, Out] ref double TK2, [In, Out] ref double PK2, [In, Out] ref double DM2,
        //  [In, Out] ref Int32 ierr, [Out, MarshalAs(UnmanagedType.LPArray)]byte[] herr_c, Int32 LengthHERR);  //255

        //protected void SATH(double hh, Int32 kph)
        //{
        //    SATHdll(ref hh,this.MoleFractions,ref kph,ref NRoots,
        //}
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
                if (ierr > 0) throw new Exception(herr.ToString());
            }
            catch (Exception e)
            {
                RefpropErrorHandler.ErrorHandler(this, e.ToString(), ierr); ;
            }
        }


        protected void WMOL(double[] x, ref double wm)
        {
            WMOLdll(MoleFractions, ref wm);
        }

    }
}
