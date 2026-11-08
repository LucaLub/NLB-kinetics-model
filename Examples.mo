within NLBDynamics;
package Examples "Level 1 transients (Turner's set where Level 1 allows it)"
  extends Modelica.Icons.ExamplesPackage;

  model BaseCase "Reference steady state: nothing perturbed, check that nothing drifts"
    extends Modelica.Icons.Example;
    Reactor reactor annotation (Placement(transformation(extent={{20,-20},{60,20}})));
    Modelica.Blocks.Sources.Step rhoExt(height=0, offset=0, startTime=1)
      annotation (Placement(transformation(extent={{-80,70},{-60,90}})));
    Modelica.Blocks.Sources.Step mdotP(height=0, offset=1, startTime=1)
      annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Modelica.Blocks.Sources.Step TPin(height=0, offset=reactor.d.T_P_in, startTime=1)
      annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
    Modelica.Blocks.Sources.Step mdotNe(height=0, offset=1, startTime=1)
      annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
    Modelica.Blocks.Sources.Step coolMod(height=0, offset=1, startTime=1)
      annotation (Placement(transformation(extent={{-80,-50},{-60,-30}})));
    Modelica.Blocks.Sources.Step coolSec(height=0, offset=1, startTime=1)
      annotation (Placement(transformation(extent={{-80,-80},{-60,-60}})));
  equation
    connect(rhoExt.y, reactor.rho_ext);
    connect(mdotP.y, reactor.mdotP_frac);
    connect(TPin.y, reactor.T_P_in);
    connect(mdotNe.y, reactor.mdotNe_frac);
    connect(coolMod.y, reactor.coolMod_frac);
    connect(coolSec.y, reactor.coolSec_frac);
    annotation (experiment(StopTime=10, Tolerance=1e-6, Interval=1e-3));
  end BaseCase;

  model ReactivityStep "Turner Fig. 7: +0.1 % (100 pcm) step at t = 1 s"
    extends BaseCase(rhoExt(height=1e-3));
    annotation (experiment(StopTime=10, Tolerance=1e-6, Interval=1e-3));
  end ReactivityStep;

  model DollarStep "Step of one circulating-fuel dollar (beta_circ) at t = 1 s"
    extends BaseCase(rhoExt(height=reactor.kin.beta_circ));
    annotation (experiment(StopTime=10, Tolerance=1e-6, Interval=1e-4));
  end DollarStep;

  model PropellantFlowLoss
    "Turner Fig. 10: -10 % propellant flow at t = 1 s (Level 1: duct mass storage neglected)"
    extends BaseCase(mdotP(height=-0.1));
    annotation (experiment(StopTime=5, Tolerance=1e-6, Interval=1e-3));
  end PropellantFlowLoss;

  model PropellantInletTemperature "Inlet temperature -100 K at t = 1 s"
    extends BaseCase(TPin(height=-100));
    annotation (experiment(StopTime=5, Tolerance=1e-6, Interval=1e-3));
  end PropellantInletTemperature;

  model ModeratorCoolantLoss "Turner Fig. 15: total loss of moderator coolant at t = 1 s"
    extends BaseCase(coolMod(height=-1));
    annotation (experiment(StopTime=20, Tolerance=1e-6, Interval=1e-2));
  end ModeratorCoolantLoss;
end Examples;
