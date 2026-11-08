within NLBDynamics;
model ReactorFMU "Reactor with its main results as outputs: export THIS model as an FMU"
  extends Reactor;

  Modelica.Blocks.Interfaces.RealOutput y_P(unit="W") = kin.P "Fission power"
    annotation (Placement(transformation(extent={{100,70},{120,90}})));
  Modelica.Blocks.Interfaces.RealOutput y_T_F_central(unit="K") = central.fuel.T_F
    "Central fuel temperature" annotation (Placement(transformation(extent={{100,40},{120,60}})));
  Modelica.Blocks.Interfaces.RealOutput y_T_F_peripheral(unit="K") = peripheral.fuel.T_F
    "Peripheral fuel temperature" annotation (Placement(transformation(extent={{100,10},{120,30}})));
  Modelica.Blocks.Interfaces.RealOutput y_T_out_central(unit="K") = central.T_out
    "Central duct outlet temperature" annotation (Placement(transformation(extent={{100,-20},{120,0}})));
  Modelica.Blocks.Interfaces.RealOutput y_T_out_peripheral(unit="K") = peripheral.T_out
    "Peripheral duct outlet temperature" annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
  Modelica.Blocks.Interfaces.RealOutput y_T_BeO(unit="K") = beo.T
    "BeO temperature" annotation (Placement(transformation(extent={{100,-80},{120,-60}})));
  Modelica.Blocks.Interfaces.RealOutput y_rho_fb = reac.rho - rho_ext
    "Total feedback reactivity (Delta k/k)" annotation (Placement(transformation(extent={{100,-110},{120,-90}})));

  annotation (Documentation(info="<html>
<p>FMU inputs: rho_ext, mdotP_frac, T_P_in, mdotNe_frac, coolMod_frac, coolSec_frac.
Give them their nominal values (0, 1, 2250 K, 1, 1, 1) for the steady state.</p>
<p>OMEdit: select this model, File &gt; Export &gt; FMU (FMI 2.0, Model Exchange and
Co-Simulation). Script: <code>buildModelFMU(NLBDynamics.ReactorFMU, version=\"2.0\",
fmuType=\"me_cs\")</code>.</p>
<p>Check that the FMU's resources folder contains h2_props.txt. If it does not,
set the parameter fileName to the table's absolute path on the target machine.</p>
</html>"));
end ReactorFMU;
