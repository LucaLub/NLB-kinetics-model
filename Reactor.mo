within NLBDynamics;
model Reactor "NLB reactor, Level 1: two cavity groups, shared solids, point kinetics"
  import SI = Modelica.Units.SI;

  parameter Data.NLBData d "Reference data";
  parameter Integer N = 10 "Propellant duct cells per cavity";
  parameter String fileName = Modelica.Utilities.Files.loadResource(
    "modelica://NLBDynamics/Resources/h2_props.txt") "H2 property tables";

  final parameter Real share_c = d.Fcav/d.nCav "Power share of the central cavity";
  final parameter Real share_p = (1 - share_c)/(d.nCav - 1) "Power share of one peripheral cavity";
  final parameter Integer nP = d.nCav - 1 "Number of peripheral cavities";
  final parameter Real f_shared = d.BeO.f_heat + d.graphite.f_heat + d.endCaps.f_heat
    + d.shell.f_heat + d.tieRods.f_heat "Power fraction deposited in the shared solids";

  // ---- Boundary inputs ----
  Modelica.Blocks.Interfaces.RealInput rho_ext(start=0) "External reactivity (Delta k/k)"
    annotation (Placement(transformation(extent={{-140,60},{-100,100}})));
  Modelica.Blocks.Interfaces.RealInput mdotP_frac(start=1) "Propellant flow / nominal"
    annotation (Placement(transformation(extent={{-140,20},{-100,60}})));
  Modelica.Blocks.Interfaces.RealInput T_P_in(unit="K", start=2250) "Propellant inlet temperature"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealInput mdotNe_frac(start=1) "Neon buffer flow / nominal"
    annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
  Modelica.Blocks.Interfaces.RealInput coolMod_frac(start=1) "Moderator coolant flow / nominal"
    annotation (Placement(transformation(extent={{-140,-100},{-100,-60}})));
  Modelica.Blocks.Interfaces.RealInput coolSec_frac(start=1) "Secondary coolant flow / nominal"
    annotation (Placement(transformation(extent={{-140,-140},{-100,-100}})));

  // ---- Neutronics and feedback ----
  Neutronics.PointKinetics kin(
    P0=d.P0, beta=d.beta, a_i=d.a_i, lambda_i=d.lambda_i, Lambda=d.Lambda,
    lambda_F=d.mdot_F/d.M_F0, tau_loop=d.tau_loop);
  replaceable Reactivity.LinearReactivity reac(
    alpha_TF=d.alpha_TF, alpha_TP=d.alpha_TP, alpha_rhoP=d.alpha_rhoP,
    alpha_MF=d.alpha_MF, alpha_RF=d.alpha_RF, a1_BeO=d.a1_BeO, a2_BeO=d.a2_BeO,
    alpha_Gr=d.alpha_Gr) constrainedby Reactivity.PartialReactivity;

  // ---- Cavity groups ----
  CavityGroup central(d=d, N=N, share=share_c, f_shared=f_shared, fileName=fileName,
    P=kin.P, mdotP_frac=mdotP_frac, T_P_in=T_P_in, mdotNe_frac=mdotNe_frac,
    coolSec_frac=coolSec_frac) "Central cavity";
  CavityGroup peripheral(d=d, N=N, share=share_p, f_shared=f_shared, fileName=fileName,
    P=kin.P, mdotP_frac=mdotP_frac, T_P_in=T_P_in, mdotNe_frac=mdotNe_frac,
    coolSec_frac=coolSec_frac) "One representative peripheral cavity (x6)";

  // ---- Shared solids ----
  Components.CooledSolid beo(d=d.BeO, mdot_c=coolMod_frac*d.BeO.mdot_c) "BeO, all 7 cells";
  Components.CooledSolid graphite(d=d.graphite, mdot_c=coolMod_frac*d.graphite.mdot_c);
  Components.CooledSolid endCaps(d=d.endCaps, mdot_c=coolMod_frac*d.endCaps.mdot_c);
  Components.CooledSolid shell(d=d.shell, mdot_c=coolSec_frac*d.shell.mdot_c);
  Components.CooledSolid tieRods(d=d.tieRods, mdot_c=coolSec_frac*d.tieRods.mdot_c);
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow qBeO, qGr, qEC, qShell, qTR;

  // ---- Diagnostics ----
  SI.Power P = kin.P "Fission power";
  SI.Temperature T_out_c = central.T_out "Central duct outlet temperature";
  SI.Temperature T_out_p = peripheral.T_out "Peripheral duct outlet temperature";
  SI.Power Q_balance "Fission power minus all heat sinks (should be ~0 at steady state)";

equation
  // Feedback inputs: cavity-weighted means (1 central + nP peripheral)
  reac.T_F = (central.fuel.T_F + nP*peripheral.fuel.T_F)/d.nCav;
  reac.dT_P = (central.dT_P + nP*peripheral.dT_P)/d.nCav;
  reac.r_rhoP = (central.r_rhoP + nP*peripheral.r_rhoP)/d.nCav;
  reac.r_MF = (central.r_MF + nP*peripheral.r_MF)/d.nCav;
  reac.r_RF = (central.r_RF + nP*peripheral.r_RF)/d.nCav;
  reac.T_BeO = beo.T;
  reac.T_Gr = graphite.T;
  reac.rho_ext = rho_ext;
  kin.rho = reac.rho;

  // Shared solids: nuclear heating as fractions of the total power
  qBeO.Q_flow = d.BeO.f_heat*kin.P;
  qGr.Q_flow = d.graphite.f_heat*kin.P;
  qEC.Q_flow = d.endCaps.f_heat*kin.P;
  qShell.Q_flow = d.shell.f_heat*kin.P;
  qTR.Q_flow = d.tieRods.f_heat*kin.P;
  connect(qBeO.port, beo.port);
  connect(qGr.port, graphite.port);
  connect(qEC.port, endCaps.port);
  connect(qShell.port, shell.port);
  connect(qTR.port, tieRods.port);

  // Global energy check: heat leaving through coolants, propellant and neon,
  // plus storage, must equal the fission power
  Q_balance = kin.P
    - (central.cell[N].h - central.h_inlet)*central.cell[N].mdot
    - nP*(peripheral.cell[N].h - peripheral.h_inlet)*peripheral.cell[N].mdot
    - (central.neon.mdot*d.cp_Ne*(central.neon.T - central.T_Ne_in)
       + nP*peripheral.neon.mdot*d.cp_Ne*(peripheral.neon.T - peripheral.T_Ne_in))
    - (beo.Q_c + graphite.Q_c + endCaps.Q_c + shell.Q_c + tieRods.Q_c)
    - (central.tw.Q_c + central.liner.Q_c + central.nozzle.Q_c)
    - nP*(peripheral.tw.Q_c + peripheral.liner.Q_c + peripheral.nozzle.Q_c)
    - (central.fuel.mdot_out*central.fuel.cv*(central.fuel.T_F - d.T_F_inj)
       + nP*peripheral.fuel.mdot_out*peripheral.fuel.cv*(peripheral.fuel.T_F - d.T_F_inj));

  annotation (Documentation(info="<html>
<p>Level 1 reactor-only network. All coolants, the neon and the fuel are boundary
flows; the duct pressure is fixed. Calibrated at initialisation:</p>
<ul>
<li>kin.rho0: reactivity bias for a steady reference (= beta - beta_circ);</li>
<li>per group: eta, hA_Ne, T_Ne_in, reference cell states;</li>
<li>per solid: hA from its reference temperature.</li>
</ul>
<p>Q_balance includes storage terms only implicitly: it is exactly zero at the
steady state and shows the stored-energy rate during transients.</p>
</html>"));
end Reactor;
