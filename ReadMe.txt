================================================================================
README - NLBDynamics (Modelica package, version 0.1.0)
================================================================================

WHAT THIS CODE DOES
--------------------
NLBDynamics is a reactor-only lumped dynamics model of the 7-cavity
Nuclear Light Bulb (NLB) gas-core reactor (UARL reference design, as
modelled in Jan Piaskowy's MSc draft, Aug 2026), written in Modelica.

The model contains:
. point kinetics with 6 precursor groups, precursor drift out of the core
  with the circulating fuel and re-entry after the external loop delay;
. per cavity group (central x1, peripheral x6, one representative cavity
  each): one fuel plasma node (e = c_v T, fixed mass), a flow-through neon
  buffer, grey radiation from the fuel to a 10-cell propellant duct with
  equilibrium H2/H properties from tables, direct nuclear heating of the
  propellant, and three cooled solids (transparent wall, liner, nozzle
  wall);
. shared solids: BeO and graphite moderators, end caps, pressure shell,
  tie rods;
. linear reactivity feedback (Jan's coefficients, quadratic for BeO).
The duct and cavity pressure is fixed at p0 = 500 atm. All flows (fuel,
neon, propellant, coolants) are imposed boundary flows.

The steady reference state is calibrated at initialisation (Turner's
method): the unknown parameters are declared with fixed=false and are
found by the initial equations together with the steady state:
  - kin.rho0: reactivity bias that makes the reference state steady
    (= beta - beta_circ);
  - per cavity group: the fuel temperature that makes the duct absorb
    f_duct of the cavity power (central ~28 400 K, peripheral ~26 100 K),
    the fuel-to-neon conductance hA_Ne that closes the fuel energy
    balance, the neon inlet temperature T_Ne_in that gives the reference
    buffer temperature, the reference cell temperatures and densities;
  - per solid: the solid-to-coolant conductance hA from its reference
    temperature;
  - the reference temperatures of the feedback model.

Running an example produces a result file (OpenModelica: <model>_res.mat)
with every variable of the model. The main ones are:
  - reactor.P                              fission power
  - reactor.central.fuel.T_F (peripheral)  fuel temperature
  - reactor.T_out_c, reactor.T_out_p       duct outlet temperatures
  - reactor.central.cell[k].T              axial propellant profile
  - reactor.central.neon.T                 neon buffer temperature
  - reactor.reac.rho_TF ... rho_Gr, rho_fb feedback breakdown (Delta k/k)
  - reactor.beo.T, reactor.graphite.T      shared solids
  - reactor.Q_balance                      global energy check (see below)

WHAT'S IN THE FILE SET, AND HOW IT FITS TOGETHER
---------------------------------------------------
  - package.mo / package.order -> package definition (uses Modelica
    Standard Library 4.0.0) and order of the sub-packages.
  - Data.mo        -> the reference data record NLBData and the SolidData
    record of the cooled solids. Every value carries its source code:
    [J] Jan's MSc draft, [U] UARL report (to obtain), [O] Jan's OpenMC
    model (to extract), [E] engineering estimate (to replace).
  - Neutronics.mo  -> PointKinetics: static beta in the power equation,
    drift loss lambda_F C_i and re-entry through delay(), beta_circ (one
    circulating-fuel dollar) and the calibrated bias rho0.
  - Reactivity.mo  -> PartialReactivity (interface) and LinearReactivity
    (coefficients per K or per %, quadratic BeO term, feedback breakdown).
    It is declared replaceable in Reactor.
  - Components.mo  -> the nodes and links of the network: FuelPlasma,
    NeonBuffer, CooledSolid (solid node + coolant node), DuctCell,
    RadiationLink, Conductance.
  - CavityGroup.mo -> one representative cavity: fuel, neon, 10 duct
    cells, radiation links, cavity solids, and the initial equations of
    the cavity calibration.
  - Reactor.mo     -> the whole reactor: 6 boundary inputs, kinetics,
    reactivity, the central and peripheral groups, shared solids, the
    cavity-weighted feedback inputs and the energy check Q_balance.
  - ReactorFMU.mo  -> Reactor with its main results as outputs, to be
    exported as an FMU.
  - Examples.mo    -> the transients (see below), ready to simulate.
  - Resources/h2_props.txt -> equilibrium H2/H tables (T_ph, rho_ph,
    h_pT are used), generated with python/make_h2_tables.py (Cantera).
    The package finds the file through
    modelica://NLBDynamics/Resources/, so it must stay inside
    NLBDynamics/.

WHAT YOU NEED TO PROVIDE / CAN CHANGE
--------------------------------------
  - Examples (Examples.mo), all with a step at t = 1 s:
      BaseCase                    nothing perturbed (drift check), 10 s
      ReactivityStep              +100 pcm (Turner Fig. 7), 10 s
      DollarStep                  +1 circulating-fuel dollar, 3 s
                                  (output interval 0.1 ms)
      PropellantFlowLoss          -10 % propellant flow (Fig. 10), 5 s
      PropellantInletTemperature  inlet temperature -100 K, 5 s
      ModeratorCoolantLoss        moderator coolant -> 0 (Fig. 15), 20 s
    New cases: extend BaseCase and change the height of one of the 6 Step
    sources (rhoExt, mdotP, TPin, mdotNe, coolMod, coolSec).
  - Reactor data: record Data.NLBData (power, pressure, delayed neutron
    data, fuel inventory and loop delay, power partition, geometry,
    radiation factors, fuel and neon properties, propellant flow and
    inlet temperature, the 8 cooled solids, reactivity coefficients).
  - Reactor parameters: N (duct cells), fileName (table path), and the
    replaceable reactivity model reac.
  - Solver settings: in each example's experiment annotation (stop time,
    tolerance 1e-6, output interval).

CAPABILITIES
-------------
  - Acausal, component-based model: each physical part is a reusable
    component and the equations are written as balances; the tool does
    the causality and the steady-state calibration.
  - Steady state and calibration solved together at initialisation, with
    no hard-coded reference temperatures.
  - Variable-step implicit solver (OpenModelica default DASSL) with error
    control.
  - Circulating-fuel kinetics with precursor re-entry after the loop
    delay, and the circulating-fuel dollar beta_circ as an output.
  - Feedback breakdown per effect, for plotting.
  - Global energy check Q_balance: fission power minus everything leaving
    through the propellant, neon, coolants and fuel extraction. Exactly
    zero at the steady state; during transients it is the rate of energy
    stored in the reactor.
  - FMU export (ReactorFMU, FMI 2.0, Model Exchange and Co-Simulation):
      inputs   rho_ext, mdotP_frac, T_P_in, mdotNe_frac, coolMod_frac,
               coolSec_frac (nominal 0, 1, 2250 K, 1, 1, 1)
      outputs  y_P, y_T_F_central, y_T_F_peripheral, y_T_out_central,
               y_T_out_peripheral, y_T_BeO, y_rho_fb
    so the reactor can be coupled to a control system or an engine model.

LIMITATIONS - PLEASE READ BEFORE TRUSTING THE RESULTS
--------------------------------------------------------
  - Fixed pressure and imposed flows: the duct and cavity stay at p0, the
    duct cells have no mass storage (a flow change acts at once along the
    whole duct), there is no nozzle, and no pressure load on the
    transparent wall is computed.
  - Fuel: constant mass (injection = extraction = lambda_F M_F), fixed
    cloud volume and radius, and a constant effective heat capacity
    c_v = 300 J/(kg K) [E] in place of a plasma equation of state. The
    injected fuel enters with c_v T_F_inj (T_F_inj = 1500 K [E]). As the
    fuel mass and radius do not change, the fuel-mass and fuel-radius
    feedback terms are always zero.
  - Neon: fixed hold-up M_Ne = 6 kg [E], one well-mixed node, c_p only.
  - Radiation: grey, with an effective radiating temperature T_F / f_avg
    (f_avg = 3), constant emissivity and a fixed radiating area
    (2 pi R2 dz per cell). The fuel-to-neon conductance is a calibrated
    constant.
  - Cooled solids: the solid-to-coolant conductance hA is constant; it
    does NOT change with the coolant flow. In ModeratorCoolantLoss the
    heat still reaches the (now stagnant) coolant node at the nominal hA,
    which makes the solid heat-up slower than with a flow-dependent hA.
  - Deposition fractions and the power split between cavities are fixed.
  - Kinetics: constant generation time Lambda and constant drift constant
    lambda_F. The loop delay tau_loop = 2 s is an estimate [E; U]; with it,
    the precursors re-enter the core during the transients, so the
    results after t = 3 s depend directly on this value. The delayed
    neutron data are labelled as placeholders (Keepin U-233 thermal, to be
    replaced by JEFF-4.0 / OpenMC data).
  - Reactivity: linear coefficients, separable effects, graphite
    coefficient 0 (not computed). The BeO quadratic is fitted to only two
    points (-598 K / +602 K). As noted in LinearReactivity, check that
    Jan's fuel-temperature sweep held the fuel density fixed, otherwise
    the density effect is counted twice.
  - H2 tables: they cover 1500-20000 K. Outside the grid the CombiTable2Ds
    extrapolate without warning.
  - Q_balance contains the storage terms only implicitly: it is not zero
    during transients, and is not a conservation error by itself.
  - Some documentation strings are out of date: Reactor lists eta as
    calibrated at initialisation, and RadiationLink says eta is set by the
    cavity, but eta is a fixed parameter (d.eta = 1); with eta fixed, the
    initialisation sets the fuel temperature instead.
  - FMU: h2_props.txt must be inside the FMU resources folder; otherwise
    set fileName to the absolute path of the table on the target machine.

REQUIREMENTS TO RUN
---------------------
  - A Modelica tool with the Modelica Standard Library 4.0.0 (developed
    with OpenModelica / OMEdit).
  - NLBDynamics/Resources/h2_props.txt.
  - OMEdit: File > Open Model/Library File, select NLBDynamics/package.mo,
    then open an example in NLBDynamics.Examples and press Simulate. The
    experiment settings of the example are used.
  - FMU export: select ReactorFMU, File > Export > FMU, or the script
    command buildModelFMU(NLBDynamics.ReactorFMU, version="2.0",
    fmuType="me_cs").

IMPROVEMENTS
--------------
  - Scale the solid-to-coolant conductance with the coolant flow (e.g.
    hA ~ flow^0.8, Dittus-Boelter), so that coolant-loss transients are
    realistic.
  - Replace the [E] estimates (loop delay, deposition fractions, solid
    masses, neon hold-up, plasma c_v) with UARL or OpenMC data, starting
    with tau_loop, which drives the kinetics after 2 s.
  - Replace the placeholder delayed neutron data with JEFF-4.0 values for
    U-233.
  - Model fuel plasma behaviour with EOS (Saha balance eq., …)
  - Add the explicit storage terms to Q_balance, so that it is zero at
    all times and becomes a true conservation check.
  - Model in- and out-flows for fuel and neon.
  - Thrust and Isp calculations (isentropic flows from choked nozzle).
  - Better feedbacks with OpenMC sweeps.

================================================================================