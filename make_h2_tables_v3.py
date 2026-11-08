"""Generate equilibrium H2/H property tables for the NLB dynamics model.

The tables are written in the Modelica CombiTable2Ds text format (one file,
several named tables). Pressure is the first table input (rows), and the
second input is specific enthalpy (forward tables) or temperature (inverse
table).

Thermochemistry: Jan Piaskowy's mechanism (Appendix A of the MSc draft),
Redlich-Kwong phase, chemical equilibrium (Damkoehler >> 1 in the duct, so
equilibrium is valid there).

Forward tables, inputs (p [Pa], h [J/kg]):
    T_ph       temperature                      K
    rho_ph     density                          kg/m3
    YH2_ph     H2 mass fraction                 -
    M_ph       mean molar mass                  kg/mol
    cpeq_ph    equilibrium cp = dh/dT|p         J/(kg K)
    drdp_ph    d rho / d p at constant h        kg/(m3 Pa)
    drdh_ph    d rho / d h at constant p        kg2/(m3 J)
    gamf_ph    frozen cp/cv (for the nozzle)    -

Inverse table, inputs (p [Pa], T [K]):
    h_pT       specific enthalpy                J/kg

Usage:
    python make_h2_tables.py   # writes ../NLBDynamics/Resources/h2_props.txt
"""
from pathlib import Path

import cantera as ct
import numpy as np

HERE = Path(__file__).resolve().parent
MECH = HERE / "h2_tables.yaml"
PHASE = "h2_dissociation_RK"
OUT = HERE.parent / "NLBDynamics" / "Resources" / "h2_props.txt"

ATM = 101325.0
P_GRID = np.linspace(300.0, 700.0, 9) * ATM   # transients stay well inside
T_MIN, T_MAX = 1500.0, 20000.0
T_GRID = np.linspace(T_MIN, T_MAX, 186)        # 100 K steps
N_H = 321


def gas():
    return ct.Solution(str(MECH), PHASE)


def state_hp(g, p, h):
    # Start from a mid-range state: a large jump straight to (h, p) can make
    # the Redlich-Kwong HP solver fail
    g.TPY = 4000.0, p, "H2:0.9, H:0.1"
    g.HP = h, p
    g.equilibrate("HP")
    return g


def h_of_pT(g, p, T):
    g.TP = T, p
    g.equilibrate("TP")
    return g.enthalpy_mass


def forward_props(g, p, h):
    """All forward properties at (p, h), derivatives by central differences."""
    state_hp(g, p, h)
    T = g.T
    rho = g.density
    yh2 = g.Y[g.species_index("H2")]
    M = g.mean_molecular_weight / 1000.0
    gamf = g.cp_mass / g.cv_mass

    dT = 1.0e-3 * T
    cpeq = (h_of_pT(g, p, T + dT) - h_of_pT(g, p, T - dT)) / (2 * dT)

    dp = 1.0e-4 * p
    rp = state_hp(g, p + dp, h).density
    rm = state_hp(g, p - dp, h).density
    drdp = (rp - rm) / (2 * dp)

    dh = 1.0e-4 * abs(h) + 1.0e3
    rp = state_hp(g, p, h + dh).density
    rm = state_hp(g, p, h - dh).density
    drdh = (rp - rm) / (2 * dh)
    return dict(T_ph=T, rho_ph=rho, YH2_ph=yh2, M_ph=M, cpeq_ph=cpeq,
                drdp_ph=drdp, drdh_ph=drdh, gamf_ph=gamf)


def write_table(f, name, u1, u2, values):
    n1, n2 = len(u1), len(u2)
    f.write(f"double {name}({n1 + 1},{n2 + 1})\n")
    f.write("0 " + " ".join(f"{x:.10e}" for x in u2) + "\n")
    for i, x in enumerate(u1):
        f.write(f"{x:.10e} " + " ".join(f"{v:.10e}" for v in values[i]) + "\n")


def main():
    g = gas()

    # Inverse table h(p, T)
    h_pT = np.array([[h_of_pT(g, p, T) for T in T_GRID] for p in P_GRID])

    # Enthalpy range common to all pressures, so every (p, h) node lies
    # between T_MIN and T_MAX
    h_lo = h_pT[:, 0].max()
    h_hi = h_pT[:, -1].min()
    H_GRID = np.linspace(h_lo, h_hi, N_H)

    names = ["T_ph", "rho_ph", "YH2_ph", "M_ph", "cpeq_ph",
             "drdp_ph", "drdh_ph", "gamf_ph"]
    fwd = {n: np.zeros((len(P_GRID), N_H)) for n in names}
    for i, p in enumerate(P_GRID):
        for j, h in enumerate(H_GRID):
            props = forward_props(g, p, h)
            for n in names:
                fwd[n][i, j] = props[n]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        f.write("#1\n")
        f.write("# Equilibrium H2/H properties, Redlich-Kwong, Cantera "
                f"{ct.__version__}\n")
        f.write(f"# p grid {P_GRID[0]/ATM:.0f}-{P_GRID[-1]/ATM:.0f} atm, "
                f"T range {T_MIN:.0f}-{T_MAX:.0f} K\n")
        for n in names:
            write_table(f, n, P_GRID, H_GRID, fwd[n])
        write_table(f, "h_pT", P_GRID, T_GRID, h_pT)

    # Short report at the reference duct pressure (500 atm)
    i500 = int(np.argmin(abs(P_GRID - 500 * ATM)))
    print(f"wrote {OUT}")
    print(f"h range {h_lo/1e6:.2f} to {h_hi/1e6:.2f} MJ/kg")
    for T in (2250.0, 4000.0, 6633.0, 8000.0):
        h = np.interp(T, T_GRID, h_pT[i500])
        pr = forward_props(g, P_GRID[i500], h)
        print(f"500 atm, T={T:6.0f} K: h={h/1e6:7.2f} MJ/kg rho={pr['rho_ph']:.3f} "
              f"Y_H2={pr['YH2_ph']:.3f} M={pr['M_ph']*1e3:.3f} g/mol "
              f"cp_eq={pr['cpeq_ph']/1e3:6.1f} kJ/kgK gamma_f={pr['gamf_ph']:.3f}")


if __name__ == "__main__":
    main()
