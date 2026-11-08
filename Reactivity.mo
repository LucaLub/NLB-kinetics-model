within NLBDynamics;
package Reactivity "Reactivity feedback closures (replaceable)"
  extends Modelica.Icons.Package;

  partial model PartialReactivity "Interface of every reactivity model"
    import SI = Modelica.Units.SI;
    input SI.Temperature T_F "Mean volume-averaged fuel temperature (cavity-weighted)";
    input SI.TemperatureDifference dT_P "Mean propellant cell temperature change from reference";
    input Real r_rhoP "Mean relative propellant density change (rho/rho0 - 1)";
    input Real r_MF "Relative fuel mass change (M/M0 - 1)";
    input Real r_RF "Relative fuel cloud radius change (R/R0 - 1)";
    input SI.Temperature T_BeO "BeO temperature";
    input SI.Temperature T_Gr "Graphite temperature";
    input Real rho_ext "External reactivity (Delta k/k)";
    output Real rho "Feedback + external reactivity (Delta k/k)";
  end PartialReactivity;

  model LinearReactivity
    "Linear coefficients from Jan's OpenMC sweeps, quadratic BeO term (Levels 0-2)"
    extends PartialReactivity;
    import SI = Modelica.Units.SI;

    parameter Real alpha_TF "Fuel temperature [pcm/K]";
    parameter Real alpha_TP "Propellant temperature [pcm/K]";
    parameter Real alpha_rhoP "Propellant density [pcm/%]";
    parameter Real alpha_MF "Fuel mass at fixed radius [pcm/%]";
    parameter Real alpha_RF "Fuel radius at fixed mass [pcm/%]";
    parameter Real a1_BeO "BeO linear term [pcm/K]";
    parameter Real a2_BeO "BeO quadratic term [pcm/K2]";
    parameter Real alpha_Gr "Graphite temperature [pcm/K]";

    parameter SI.Temperature T_F_ref(fixed=false) "Reference fuel temperature";
    parameter SI.Temperature T_BeO_ref(fixed=false) "Reference BeO temperature";
    parameter SI.Temperature T_Gr_ref(fixed=false) "Reference graphite temperature";

    // Components in Delta k/k, for plotting the feedback breakdown
    Real rho_TF, rho_TP, rho_rhoP, rho_MF, rho_RF, rho_BeO, rho_Gr, rho_fb;

  protected
    constant Real pcm = 1e-5;

  equation
    rho_TF = alpha_TF*pcm*(T_F - T_F_ref);
    rho_TP = alpha_TP*pcm*dT_P;
    rho_rhoP = alpha_rhoP*pcm*100*r_rhoP;
    rho_MF = alpha_MF*pcm*100*r_MF;
    rho_RF = alpha_RF*pcm*100*r_RF;
    rho_BeO = pcm*(a1_BeO*(T_BeO - T_BeO_ref) + a2_BeO*(T_BeO - T_BeO_ref)^2);
    rho_Gr = alpha_Gr*pcm*(T_Gr - T_Gr_ref);
    rho_fb = rho_TF + rho_TP + rho_rhoP + rho_MF + rho_RF + rho_BeO + rho_Gr;
    rho = rho_fb + rho_ext;

  initial equation
    T_F_ref = T_F;
    T_BeO_ref = T_BeO;
    T_Gr_ref = T_Gr;

    annotation (Documentation(info="<html>
<p>Jan's coefficients are per % for densities and radius, per K for temperatures.
The fuel-density coefficient was computed at fixed cloud radius, so it acts on
the fuel mass; the radius coefficient at fixed mass. Check in Jan's code that
the fuel-temperature sweep held the density fixed, otherwise the density effect
is counted twice.</p>
<p>BeO: the reference sits at a reactivity minimum, so a quadratic is used
(a1, a2 fitted to Jan's -598 K and +602 K points; replace with a proper fit).</p>
</html>"));
  end LinearReactivity;
end Reactivity;
