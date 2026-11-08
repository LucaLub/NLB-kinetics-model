within;
package NLBDynamics "Lumped dynamics model of the Nuclear Light Bulb gas-core reactor"
  extends Modelica.Icons.Package;

  annotation (
    uses(Modelica(version="4.0.0")),
    version="0.1.0",
    Documentation(info="<html>
<p>Reactor-only lumped dynamics model of the 7-cavity Nuclear Light Bulb
(UARL reference design, as modelled in Piaskowy's MSc draft, Aug 2026).</p>
<p><b>Level 1</b>: point kinetics with precursor drift and loop re-entry,
one fuel plasma node per cavity, neon buffer, transparent wall, 10-cell
propellant duct with equilibrium H2 tables, liner, nozzle wall, shared
moderators and structures. Boundary flows and pressure are fixed; the
reactivity feedback is linear (quadratic for BeO).</p>
<p>Later levels redeclare closures (reactivity, fuel EOS, flows) without
changing the network.</p>
<p>Property tables: run <code>python/make_h2_tables.py</code> to regenerate
<code>Resources/h2_props.txt</code>.</p>
</html>"));
end NLBDynamics;
