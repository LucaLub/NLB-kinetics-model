within NLBDynamics;
package Components "Lumped nodes and links of the thermal network"
  extends Modelica.Icons.Package;
  import SI = Modelica.Units.SI;
  import Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a;
  import Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b;

  model FuelPlasma "Fuel plasma node of one cavity (internal energy U_F, mass M_F)"
    parameter SI.Mass M_F0 "Reference fuel mass in the cavity";
    parameter SI.SpecificHeatCapacity cv "Effective heat capacity (e = cv T at Levels 0-2)";
    parameter SI.Temperature T0 "Start value of the volume-averaged temperature";
    parameter SI.Temperature T_inj "Temperature of the injected fuel";
    parameter SI.DecayConstant lambda_F "Fuel extraction constant";

    input SI.HeatFlowRate Q_fis "Fission heat deposited in the plasma";

    SI.Mass M_F "Fuel mass (fixed at Level 1, state from Level 2)";
    SI.Energy U_F(start=M_F0*cv*T0) "Internal energy";
    SI.Temperature T_F(start=T0) "Volume-averaged fuel temperature";
    SI.MassFlowRate mdot_out "Extracted fuel flow";
    SI.MassFlowRate mdot_in "Injected fuel flow";

    HeatPort_a port "Volume-averaged fuel temperature; radiation and convection leave here"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  equation
    // ---- Level 1 closures (redeclared / replaced at Level 2) ----
    M_F = M_F0;
    mdot_out = lambda_F*M_F;
    mdot_in = mdot_out;
    // ---- balances ----
    U_F = M_F*cv*T_F;
    der(U_F) = Q_fis + port.Q_flow - mdot_out*U_F/M_F + mdot_in*cv*T_inj;
    port.T = T_F;
  end FuelPlasma;

  model NeonBuffer "Well-mixed neon buffer of one cavity, flow-through"
    parameter SI.Mass M "Neon hold-up";
    parameter SI.SpecificHeatCapacity cp "Neon specific heat";
    parameter SI.Temperature T0 "Start value";

    input SI.MassFlowRate mdot "Neon flow (boundary)";
    input SI.Temperature T_in "Neon inlet temperature (boundary)";

    SI.Temperature T(start=T0) "Buffer (and outlet) temperature";
    HeatPort_a port annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  equation
    M*cp*der(T) = port.Q_flow - mdot*cp*(T - T_in);
    port.T = T;
  end NeonBuffer;

  model CooledSolid
    "Solid node + well-mixed coolant node; hA calibrated from the reference temperature"
    parameter Data.SolidData d "Node data";
    parameter SI.ThermalConductance hA(fixed=false, start=1e4, min=0)
      "Solid-to-coolant conductance (found at initialisation from d.T0)";

    input SI.MassFlowRate mdot_c "Coolant flow (boundary)";

    SI.Temperature T(start=d.T0) "Solid temperature";
    SI.Temperature Tc(start=0.5*(d.T0 + d.Tc_in)) "Coolant (outlet) temperature";
    SI.HeatFlowRate Q_c "Heat to the coolant";

    HeatPort_a port "Heat deposited in the solid" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  equation
    Q_c = hA*(T - Tc);
    d.m*d.cp*der(T) = port.Q_flow - Q_c;
    d.m_c*d.cp_c*der(Tc) = Q_c - mdot_c*d.cp_c*(Tc - d.Tc_in);
    port.T = T;
    assert(hA > 0, "CooledSolid: calibrated hA <= 0; the coolant cannot remove the reference heat "
      + "at the reference temperature. Lower d.Tc_in or raise d.mdot_c.");

  initial equation
    T = d.T0;
    der(T) = 0;
    der(Tc) = 0;
  end CooledSolid;

  model DuctCell "Propellant duct cell, equilibrium H2/H from Cantera tables"
    parameter SI.Volume V "Cell volume";
    parameter String fileName "Property table file";
    parameter SI.SpecificEnthalpy h_start "Start value of the cell enthalpy";

    input SI.Pressure p "Duct pressure (fixed at Level 1)";
    input SI.MassFlowRate mdot "Flow through the cell";
    input SI.SpecificEnthalpy h_in "Enthalpy of the inflow (upwind cell)";

    SI.SpecificEnthalpy h(start=h_start) "Cell (and outflow) enthalpy";
    SI.Temperature T "Cell temperature";
    SI.Density rho "Cell density";
    SI.Mass m "Hydrogen hold-up";

    HeatPort_a port "Radiation and direct nuclear heating" annotation (Placement(transformation(extent={{-10,90},{10,110}})));

    Modelica.Blocks.Tables.CombiTable2Ds tabT(tableOnFile=true, tableName="T_ph",
      fileName=fileName, smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
    Modelica.Blocks.Tables.CombiTable2Ds tabRho(tableOnFile=true, tableName="rho_ph",
      fileName=fileName, smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

  equation
    tabT.u1 = p;
    tabT.u2 = h;
    tabRho.u1 = p;
    tabRho.u2 = h;
    T = tabT.y;
    rho = tabRho.y;
    m = rho*V;
    // Level 1: constant mass flow through the cell (mass storage neglected)
    m*der(h) = mdot*(h_in - h) + port.Q_flow;
    port.T = T;
  end DuctCell;

  model RadiationLink
    "Net grey radiation fuel -> duct cell, Turner form: sigma eps A (eta (T_F/f_avg)^4 - T_k^4)"
    parameter SI.Area A "Radiating fuel surface facing the cell";
    parameter Real eps "Effective emissivity";
    parameter Real f_avg "Volume-average to effective radiating temperature ratio";
    parameter Real eta "Calibration factor (set by the cavity at initialisation)";

    SI.HeatFlowRate Q "Net radiation fuel -> cell";
    HeatPort_a port_a "Fuel side" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    HeatPort_b port_b "Propellant cell side" annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  equation
    Q = Modelica.Constants.sigma*eps*A*(eta*(port_a.T/f_avg)^4 - port_b.T^4);
    port_a.Q_flow = Q;
    port_b.Q_flow = -Q;
  end RadiationLink;

  model Conductance "Linear thermal conductance with a calibratable value"
    parameter SI.ThermalConductance G "Conductance";
    SI.HeatFlowRate Q "Heat a -> b";
    HeatPort_a port_a annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    HeatPort_b port_b annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  equation
    Q = G*(port_a.T - port_b.T);
    port_a.Q_flow = Q;
    port_b.Q_flow = -Q;
  end Conductance;
end Components;
