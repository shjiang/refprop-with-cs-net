using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

namespace sc.net
{
    static class Program
    {
        /// <summary>
        /// 应用程序的主入口点。
        /// </summary>
        [STAThread]
        static void Main()
        {
            //Application.EnableVisualStyles();
            //Application.SetCompatibleTextRenderingDefault(false);
            //mf_net_wa__refprop refcs=new mf_net_wa__refprop();
            
            Refrigerant myRefrigerant = new Refrigerant(RefrigerantCategory.PredefinedMixture, "R410A",ReferenceState.DEF);

            //Refrigerant myRefrigerant = new Refrigerant(RefrigerantCategory.PureFluid, "R22", ReferenceState.DEF);
            
            //Refrigerant myRefrigerant = new Refrigerant(RefrigerantCategory.NewMixture, "R32=0.69761,R125=0.30239",ReferenceState.DEF);
            
            Console.WriteLine(myRefrigerant.MolecularWeight);
            //myRefrigerant.FindSaturatedStateWithTemperature(273, SaturationPoint.Dew_Point);
            //myRefrigerant.DisplayThermoDynamicState();
            
            //myRefrigerant.FindStateWithTP(300.0, 1097.0);// 在两相区，TP计算出来的是饱和液相的状态
            //myRefrigerant.DisplayThermoDynamicState();

            myRefrigerant.FindStateWithTD(300, 40 / myRefrigerant.MolecularWeight);
            myRefrigerant.DisplayThermoDynamicState();
            
            //myRefrigerant.Display();

            System.Console.ReadKey();
        }
    }
}
