using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace sc.net
{
    public interface IFluid
    {
        double MolecularWeight { get; }
        double TriplePointTemperature { get; }
        double NormalBoilingPointTemperature { get; }
        double CriticalTemperature { get; }
        double CriticalPressure { get; }
        double CriticalDensity { get; }
        double CriticalPointCompressibility { get; }
        double AccentricFactor { get; }
        double DipoleMoment { get; }
        double GasConstatnt_R { get; }

        //public UnitsBasis CurrentUnitsBasis;

        //void GetSaturatedPressure(ref double SaturatedPressure, ref double temperature, ref Int32 PhaseFlag);
        //void GetSaturatedTemperature(double pressure);
    }
}
