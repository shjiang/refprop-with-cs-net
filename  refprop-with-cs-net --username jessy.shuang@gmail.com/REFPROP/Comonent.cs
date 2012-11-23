using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace sc.net
{
    class Comonent
    {
        public string Name { get; set; }
        public double MoleFraction { set; get; }
        public double Mass_Fraction { get; set; }
        public double MolecularWeight { get; set; }
        public double TriplePointTemperature { get; set; }
        public double NormalBoilingPointTemperature{get;set;}
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
