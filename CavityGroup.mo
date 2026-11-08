within NLBDynamics;
model CavityGroup
  "One representative cavity of a group (central or peripheral): fuel, neon, wall, duct, liner, nozzle"
  import SI = Modelica.Units.SI;
  import Modelica.Constants.pi;

  parameter Data.NLBData d "Reference data";
  parameter Integer N = 10 "Number of propellant duct cells";
  parameter Real share "Fraction of the total power produced by ONE cavity of this group";
  parameter Real f_shared "Fraction of total power deposited in the shared solids";
  parameter String fileName "H2 property table file";

  // ---- Inputs from the reactor (bound by the parent) ----
  input SI.Power P "Total fission power";
  input Real mdotP_frac "Propellant flow / nominal";
  input SI.Temperature T_P_in "Propellant inlet temperature";
  input Real mdotNe_frac "Neon buffer flow / nominal";
  input Real coolSec_frac "Secondary coolant flow / nominal (wall, liner, nozzle)";

  parameter Real eta = d.eta "Fuel-side radiation factor";

  // ---- Calibrated at initialisation (Turner's method) ----
  parameter SI.ThermalConductance hA_Ne(fixed=false, start=1e3)
    "Fuel-to-neon conductance: closes the fuel energy balance at T_F0";
  parameter SI.Temperature T_Ne_in(fixed=false, start=2500)
    "Neon inlet temperature that gives the reference buffer temperature";

  // ---- Geometry ----
  final parameter SI.Length dz = d.L/N;
  final parameter SI.Radius r_tw = d.R3 + d.d_tw "Outer radius of the transparent wall";
  final parameter SI.Radius r_edge[N + 1] = {d.r_in + (d.r_out - d.r_in)*(k - 1)/N for k in 1:N + 1}
    "Liner-side duct radius at the cell edges";
  final parameter SI.Volume V_cell[N] = {pi*dz/3*(r_edge[k]^2 + r_edge[k]*r_edge[k + 1] + r_edge[k + 1]^2)
    - pi*r_tw^2*dz for k in 1:N};
  final parameter SI.Area A_rad = 2*pi*d.R2*dz "Fuel radiating area facing one cell";
  final parameter Real f_out = d.f_H2dir + d.TW.f_heat + d.liner.f_heat + d.nozzle.f_heat + f_shared
    "Fraction of power NOT deposited in the fuel plasma";

  // ---- Reference states for the reactivity feedback (found at initialisation) ----
  parameter SI.Temperature T_ref[N](each fixed=false, each start=4000);
  parameter SI.Density rho_ref[N](each fixed=false, each start=2.5);

  SI.Power P_cav "Power of one cavity of this group";
  SI.SpecificEnthalpy h_inlet "Inlet enthalpy from (p, T_P_in)";
  SI.Temperature T_out "Duct outlet temperature";
  SI.TemperatureDifference dT_P "Mean cell temperature change";
  Real r_rhoP "Mean relative cell density change";
  Real r_MF "Relative fuel mass change";
  Real r_RF "Relative fuel radius change (0 at Level 1: no fuel EOS yet)";
  SI.HeatFlowRate Q_duct "Net radiation absorbed by the duct";

  Components.FuelPlasma fuel(
    M_F0=d.M_F0/d.nCav, cv=d.cv_F, T0=d.T_F0, T_inj=d.T_F_inj,
    lambda_F=d.mdot_F/d.M_F0, Q_fis=(1 - f_out)*P_cav);
  Components.NeonBuffer neon(M=d.M_Ne, cp=d.cp_Ne, T0=d.T_Ne0,
    mdot=mdotNe_frac*d.mdot_Ne, T_in=T_Ne_in);
  Components.Conductance convNe(G=hA_Ne) "Fuel-to-neon convection";
  Components.RadiationLink rad[N](each A=A_rad, each eps=d.eps_eff, each f_avg=d.f_avg, each eta=eta);
  Components.DuctCell cell[N](
    V=V_cell, each fileName=fileName,
    h_start={31.06e6 + (230.6e6 - 31.06e6)*(k - 0.5)/N for k in 1:N},
    each p=d.p0, each mdot=mdotP_frac*d.mdot_P);
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow qH2[N] "Direct nuclear heating of the propellant";
  Components.CooledSolid tw(d=d.TW, mdot_c=coolSec_frac*d.TW.mdot_c) "Transparent wall";
  Components.CooledSolid liner(d=d.liner, mdot_c=coolSec_frac*d.liner.mdot_c) "Liner tubes";
  Components.CooledSolid nozzle(d=d.nozzle, mdot_c=coolSec_frac*d.nozzle.mdot_c) "Nozzle wall";
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow qTW, qLiner, qNozzle;

  Modelica.Blocks.Tables.CombiTable2Ds tabHin(tableOnFile=true, tableName="h_pT",
    fileName=fileName, smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

equation
  P_cav = share*P;

  // ---- Propellant duct: upwind chain ----
  tabHin.u1 = d.p0;
  tabHin.u2 = T_P_in;
  h_inlet = tabHin.y;
  cell[1].h_in = h_inlet;
  for k in 2:N loop
    cell[k].h_in = cell[k - 1].h;
  end for;
  for k in 1:N loop
    connect(fuel.port, rad[k].port_a);
    connect(rad[k].port_b, cell[k].port);
    qH2[k].Q_flow = d.f_H2dir*P_cav/N;
    connect(qH2[k].port, cell[k].port);
  end for;

  // ---- Neon buffer ----
  connect(fuel.port, convNe.port_a);
  connect(convNe.port_b, neon.port);

  // ---- Cooled structures (prescribed deposition at Level 1) ----
  qTW.Q_flow = d.TW.f_heat*P_cav;
  qLiner.Q_flow = d.liner.f_heat*P_cav;
  qNozzle.Q_flow = d.nozzle.f_heat*P_cav;
  connect(qTW.port, tw.port);
  connect(qLiner.port, liner.port);
  connect(qNozzle.port, nozzle.port);

  // ---- Outputs for the reactivity model ----
  T_out = cell[N].T;
  dT_P = sum(cell[k].T - T_ref[k] for k in 1:N)/N;
  r_rhoP = sum(cell[k].rho/rho_ref[k] - 1 for k in 1:N)/N;
  r_MF = fuel.M_F/fuel.M_F0 - 1;
  r_RF = 0;
  Q_duct = sum(rad.Q);

initial equation
  // Fuel: steady, and the duct takes f_duct of the cavity power. With eta fixed,
  // this sets the group's fuel temperature (central ~28 400 K, peripheral ~26 100 K)
  // and the fuel-to-neon conductance closes the fuel energy balance
  der(fuel.U_F) = 0;
  Q_duct = d.f_duct*P_cav;
  // Neon: reference temperature, steady
  neon.T = d.T_Ne0;
  der(neon.T) = 0;
  // Duct: steady profile
  der(cell.h) = zeros(N);
  // References for the feedback
  T_ref = cell.T;
  rho_ref = cell.rho;
end CavityGroup;
