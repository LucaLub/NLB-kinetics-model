within NLBDynamics;
package Data "Reference parameter sets"
  extends Modelica.Icons.RecordsPackage;

  record SolidData "Lumped solid node cooled by a boundary coolant flow"
    extends Modelica.Icons.Record;
    parameter Modelica.Units.SI.Mass m "Solid mass";
    parameter Modelica.Units.SI.SpecificHeatCapacity cp "Solid specific heat (constant at Level 1)";
    parameter Modelica.Units.SI.Temperature T0 "Reference (steady-state) solid temperature";
    parameter Real f_heat "Heat deposited in the node, fraction of the relevant power (total or per cavity)";
    parameter Modelica.Units.SI.MassFlowRate mdot_c "Nominal coolant flow through the node";
    parameter Modelica.Units.SI.SpecificHeatCapacity cp_c "Coolant specific heat";
    parameter Modelica.Units.SI.Temperature Tc_in "Coolant inlet temperature (boundary condition)";
    parameter Modelica.Units.SI.Mass m_c "Coolant hold-up in the node";
  end SolidData;

  record NLBData "Reference data of the NLB reactor (Level 1)"
    extends Modelica.Icons.Record;
    import SI = Modelica.Units.SI;

    // Source codes: [J] Jan's MSc draft, [U] UARL report (to obtain),
    // [O] Jan's OpenMC model (to extract), [E] engineering estimate (replace)

    // ---- Engine ----
    parameter SI.Power P0 = 4.6e9 "Total fission power [J]";
    parameter Integer nCav = 7 "Number of cavities [J]";
    parameter Real Fcav = 1.286 "Central-cavity power factor [J]";
    parameter SI.Pressure p0 = 500*101325 "Cavity / duct pressure [J]";

    // ---- Neutronics ----
    parameter Real beta = 283e-5 "Static effective delayed neutron fraction [J]";
    parameter Real a_i[6] = {0.086, 0.299, 0.252, 0.278, 0.051, 0.034}
      "Relative group abundances, Keepin U-233 thermal [placeholder: JEFF-4.0 / O]";
    parameter SI.DecayConstant lambda_i[6] = {0.0126, 0.0337, 0.139, 0.325, 1.13, 2.50}
      "Group decay constants, Keepin U-233 thermal [placeholder: JEFF-4.0 / O]";
    parameter SI.Time Lambda = 1010e-6 "Effective generation time [J]";
    parameter SI.Mass M_F0 = 15.0 "Total fuel mass in the 7 cavities [J]";
    parameter SI.MassFlowRate mdot_F = 3.5 "Total fuel circulation flow [J]";
    parameter SI.Time tau_loop = 2.0 "External fuel loop delay for precursor re-entry [E; U]";

    // ---- Power partition (fractions of total fission power) ----
    parameter Real f_duct = 0.8735 "Thermal radiation absorbed by the propellant [J]";
    parameter Real f_H2dir = 7.5/4600 "Direct nuclear heating of the propellant [J]";
    parameter Real f_TW = 0.002 "Transparent wall (UV + nuclear) [E; U]";
    parameter Real f_L = 0.005 "Liner tubes (radiation leak + nuclear) [E; U]";
    parameter Real f_Nz = 0.001 "Nozzle wall [E; U]";
    // shared solids: see SolidData.f_heat below

    // ---- Geometry (per cavity) ----
    parameter SI.Length L = 1.8288 "Active length [J]";
    parameter SI.Radius R2 = 0.2075 "Edge-of-fuel radius [J]";
    parameter SI.Radius R3 = 0.2445 "Transparent wall inner radius [J]";
    parameter SI.Length d_tw = 0.127e-2 + 2*0.0127e-2 "Transparent wall tube outer diameter [J]";
    parameter SI.Radius r_in = 0.2777 "Duct liner radius at inlet [J]";
    parameter SI.Radius r_out = 0.4023 "Duct liner radius at outlet [J]";

    // ---- Radiation ----
    parameter Real eps_eff = 0.85 "Effective emissivity fuel-to-propellant [J]";
    parameter Real f_avg = 3 "Volume-average to effective radiating temperature ratio [J]";
    parameter Real eta = 1.0
      "Fuel-side radiation factor (Turner's eta); 1 = Jan's radiative balance, the fuel temperature then follows";

    // ---- Fuel plasma ----
    parameter SI.Temperature T_F0 = 26000
      "Start value of the fuel temperature only; the steady value is solved per cavity group [J, layer mean]";
    parameter SI.SpecificHeatCapacity cv_F = 300 "Effective plasma heat capacity incl. ionisation [E; Saha at L3]";
    parameter SI.Temperature T_F_inj = 1500 "Temperature of injected fuel [E; U]";

    // ---- Neon buffer (per cavity) ----
    parameter SI.MassFlowRate mdot_Ne = 86.7/7 "Buffer neon flow per cavity [J]";
    parameter SI.SpecificHeatCapacity cp_Ne = 1030 "Neon specific heat";
    parameter SI.Mass M_Ne = 6.0 "Neon hold-up per cavity [E: buffer volume x mean density]";
    parameter SI.Temperature T_Ne0 = 4722.2 "Mean buffer temperature [J]";

    // ---- Propellant (per cavity) ----
    parameter SI.MassFlowRate mdot_P = 19.3/7 "Propellant flow per cavity [J]";
    parameter SI.Temperature T_P_in = 2250 "Duct inlet temperature [J]";

    // ---- Cavity-group solids (per cavity; f_heat = fraction of cavity power) ----
    parameter SolidData TW(m=2.8, cp=1250, T0=1698.3, f_heat=f_TW,
      mdot_c=19.3/7, cp_c=14800, Tc_in=850, m_c=0.05)
      "Transparent wall silica + secondary H2 [m: E from 1310 tori; T0: J; coolant: U]";
    parameter SolidData liner(m=16, cp=3100, T0=1100, f_heat=f_L,
      mdot_c=19.3/7, cp_c=14800, Tc_in=700, m_c=0.1)
      "Beryllium liner tubes + secondary H2 [E; U]";
    parameter SolidData nozzle(m=20, cp=500, T0=750, f_heat=f_Nz,
      mdot_c=19.3/7, cp_c=14500, Tc_in=400, m_c=0.1)
      "Nozzle wall + secondary H2 [E; U]";

    // ---- Shared solids (f_heat = fraction of total power) ----
    parameter SolidData BeO(m=8800, cp=2050, T0=1698.3, f_heat=216.9/4600,
      mdot_c=19.3, cp_c=16000, Tc_in=700, m_c=2)
      "BeO moderator, all 7 cells [m: E from geometry, O; heat: J; coolant inlet: U]";
    parameter SolidData graphite(m=7000, cp=2050, T0=2392.2, f_heat=58.7/4600,
      mdot_c=19.3, cp_c=17000, Tc_in=1700, m_c=2)
      "Graphite moderator [m: E, O; heat: J; coolant: U]";
    parameter SolidData endCaps(m=3000, cp=2050, T0=2000, f_heat=66.8/4600,
      mdot_c=19.3, cp_c=16000, Tc_in=1500, m_c=1)
      "End caps + end-wall liners [m, T0: E, O; heat: J; coolant: U]";
    parameter SolidData shell(m=2000, cp=1100, T0=700, f_heat=17.6/4600,
      mdot_c=19.3, cp_c=14500, Tc_in=350, m_c=2)
      "Pressure shell [m: E, O; heat: J; T0 < 800 K: J; coolant: U]";
    parameter SolidData tieRods(m=100, cp=3100, T0=700, f_heat=0.001,
      mdot_c=19.3, cp_c=14500, Tc_in=450, m_c=0.5)
      "Tie rods + flow divider [E; U]";

    // ---- Reactivity coefficients (Jan, ch. 5) ----
    parameter Real alpha_TF = 0.023 "Fuel temperature [pcm/K]";
    parameter Real alpha_TP = -0.24 "Propellant temperature, all cells shifted [pcm/K]";
    parameter Real alpha_rhoP = -8.5 "Propellant density [pcm/%]";
    parameter Real alpha_MF = 395 "Fuel density at fixed radius = fuel mass [pcm/%]";
    parameter Real alpha_RF = 40.8 "Fuel cloud radius at fixed mass [pcm/%]";
    parameter Real a1_BeO = -0.021 "BeO linear term [pcm/K, E: fit to Jan Tab. 5.11/5.12]";
    parameter Real a2_BeO = 5.1e-4 "BeO quadratic term [pcm/K2, E: fit to Jan Tab. 5.11/5.12]";
    parameter Real alpha_Gr = 0 "Graphite temperature [pcm/K; not computed, O]";
  end NLBData;
end Data;
