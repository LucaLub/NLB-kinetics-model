within NLBDynamics;
package Neutronics "Point kinetics with circulating fuel"
  extends Modelica.Icons.Package;

  model PointKinetics
    "Point kinetics, 6 precursor groups, drift loss and re-entry after the external loop delay"
    import SI = Modelica.Units.SI;

    parameter SI.Power P0 "Initial (steady) power";
    parameter Real beta "Static effective delayed neutron fraction";
    parameter Real a_i[6] "Relative group abundances";
    parameter SI.DecayConstant lambda_i[6] "Group decay constants";
    parameter SI.Time Lambda "Effective generation time";
    parameter SI.DecayConstant lambda_F "Fuel drift constant (fuel flow / fuel mass)";
    parameter SI.Time tau_loop(min=0) "External loop delay before precursors re-enter";
    final parameter Real beta_i[6] = a_i*beta "Group delayed fractions";

    parameter Real rho0(fixed=false)
      "Reactivity bias making the reference state steady (found at initialisation)";

    Modelica.Blocks.Interfaces.RealInput rho "Feedback + external reactivity (Delta k/k)"
      annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
    Modelica.Blocks.Interfaces.RealOutput P(start=P0, unit="W") "Fission power"
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));

    final parameter Real beta_circ = sum(beta_i[i]*lambda_i[i]/(lambda_i[i]
      + lambda_F*(1 - exp(-lambda_i[i]*tau_loop))) for i in 1:6)
      "Effective delayed fraction with drift and re-entry (one circulating-fuel dollar)";

    Real C[6](each unit="W") "Precursor concentrations, power-equivalent units";
    Real rho_tot "Total reactivity seen by the kinetics";

  equation
    rho_tot = rho + rho0;
    der(P) = (rho_tot - beta)/Lambda*P + lambda_i*C;
    for i in 1:6 loop
      der(C[i]) = beta_i[i]/Lambda*P - (lambda_i[i] + lambda_F)*C[i]
        + lambda_F*exp(-lambda_i[i]*tau_loop)*delay(C[i], tau_loop, tau_loop);
    end for;

  initial equation
    P = P0;
    der(P) = 0;
    der(C) = zeros(6);

    annotation (Documentation(info="<html>
<p>dP/dt = (rho - beta)/Lambda P + sum(lambda_i C_i)</p>
<p>dC_i/dt = beta_i/Lambda P - (lambda_i + lambda_F) C_i + lambda_F exp(-lambda_i tau_loop) C_i(t - tau_loop)</p>
<p>The static beta is used in the power equation; the circulation loss comes from
the drift terms only, so Jan's beta_circ must not be used as well.
rho0 is the reactivity bias that makes the reference state steady. With fuel
circulation it equals beta - beta_circ: the reactivity the static core must
carry to make up for precursors decaying outside it.</p>
</html>"));
  end PointKinetics;
end Neutronics;
