import math
import sys

R_FARADAY2 = 8.931e-10          
UNITS = {"m²/s": 1.0, "cm²/s": 1e4, "mm²/s": 1e6,
         "ft²/h": 38750.0775, "ft²/s": 10.7639104, "in²/s": 1550.0031}
P_ATM = {"atm": 1.0, "bar": 0.986923, "kPa": 0.00986923, "psi": 0.0680460, "mmHg": 1 / 760}

# ───────────────────────── BASE DE DATOS ─────────────────────────
# M [g/mol], sig [Å], eps [K] (eps/k), mu [debye], Tb [K], Vb [cm³/mol, Le Bas], vf (vol. difusión FSG)
GAS = {
    "Acetileno":dict(M=26.04,  sig=3.281, eps=444.0,  mu=0.0,  Tb=189.3,  Vb=40.7, vf=24.0),
    "Aire":     dict(M=28.97,  sig=3.711, eps=78.6,   mu=0.0,  Tb=78.7,   Vb=29.9, vf=20.1),
    "Argon":    dict(M=39.948, sig=3.542, eps=93.3,   mu=0.0,  Tb=87.3,   Vb=29.0, vf=16.1),
    "Arsina":   dict(M=77.95,  sig=4.000, eps=400.0,  mu=0.0,  Tb=187.3,  Vb=55.0, vf=40.0),
    "Benceno":  dict(M=78.11,  sig=5.349, eps=412.3,  mu=0.0,  Tb=353.2,  Vb=96.5, vf=90.96),
    "Bromo":    dict(M=79.90,  sig=4.400, eps=245.0,  mu=0.0,  Tb=332.0,  Vb=69.0, vf=55.0),
    "iso-Butano": dict(M=58.12,  sig=4.000, eps=200.0,  mu=0.0,  Tb=272.7,  Vb=91.0, vf=60.0),
    "Dioxido de carbono":      dict(M=44.010, sig=3.941, eps=195.2,  mu=0.0,  Tb=194.7,  Vb=34.0, vf=26.9),
    "Disulfuro de carbono":      dict(M=76.14,  sig=4.360, eps=298.0,  mu=0.0,  Tb=319.0,  Vb=58.0, vf=45.0),
    "Monoxido de carbono":       dict(M=28.010, sig=3.690, eps=91.7,   mu=0.1,  Tb=81.7,   Vb=30.7, vf=18.9),
    "Tetracolruro de carbono": dict(M=153.82, sig=4.730, eps=398.0,  mu=0.0,  Tb=349.0,  Vb=97.0, vf=85.0),
    "Sulfuro de carbonilo": dict(M=76.14,  sig=4.360, eps=298.0,  mu=0.0,  Tb=319.0,  Vb=58.0, vf=45.0),
    "Cloro":    dict(M=35.45,  sig=3.600, eps=178.0,  mu=0.0,  Tb=239.1,  Vb=56.0, vf=40.0),
    "Cloroformo":  dict(M=119.38, sig=4.810, eps=330.0,  mu=1.04, Tb=334.0,  Vb=80.0, vf=70.0),
    "Cianogeno":  dict(M=52.04,  sig=3.800, eps=200.0,  mu=0.0,  Tb=233.0,  Vb=50.0, vf=40.0),
    "Ciclohexano":  dict(M=84.16,  sig=5.090, eps=362.0,  mu=0.0,  Tb=353.8,  Vb=108.0, vf=100.0),
    "Etano":     dict(M=30.070, sig=4.443, eps=237.0,  mu=0.0,  Tb=184.5,  Vb=58.0, vf=40.0),
    "Etanol":   dict(M=46.07,  sig=4.530, eps=362.6,  mu=1.69, Tb=351.4,  Vb=62.2, vf=51.77),
    "Etileno":  dict(M=28.054, sig=4.443, eps=237.0,  mu=0.0,  Tb=169.4,  Vb=56.0, vf=40.0),
    "Fluoreno":    dict(M=18.998, sig=3.150, eps=82.0,   mu=0.0,  Tb=85.0,   Vb=38.0, vf=20.0),
    "Helio":       dict(M=4.003,  sig=2.551, eps=10.22,  mu=0.0,  Tb=4.22,   Vb=32.5, vf=2.88),
    "n-Hexano":    dict(M=86.18,  sig=4.290, eps=230.0,  mu=0.0,  Tb=341.9,  Vb=140.6, vf=120.0),
    "Hidrogeno":       dict(M=2.016,  sig=2.827, eps=59.7,   mu=0.0,  Tb=20.4,   Vb=14.3, vf=7.07),
    "Nitrogeno":       dict(M=28.013, sig=3.798, eps=71.4,   mu=0.0,  Tb=77.3,   Vb=31.2, vf=17.9),
    "Oxigeno":       dict(M=31.999, sig=3.467, eps=106.7,  mu=0.0,  Tb=90.2,   Vb=25.6, vf=16.6),
    "Metano":      dict(M=16.043, sig=3.758, eps=148.6,  mu=0.0,  Tb=111.7,  Vb=37.7, vf=24.42),
    "Agua (vapor)":      dict(M=18.015, sig=2.641, eps=809.1,  mu=1.85, Tb=373.15, Vb=18.9, vf=12.7),
    "Amoniaco":      dict(M=17.03,  sig=2.900, eps=558.3,  mu=1.47, Tb=239.7,  Vb=25.8, vf=14.9),
    "Dioxido de azufre":      dict(M=64.06,  sig=4.112, eps=335.4,  mu=1.63, Tb=263.1,  Vb=44.8, vf=41.1),
    "Sulfuro de hidrogeno":      dict(M=34.08,  sig=3.623, eps=301.1,  mu=0.97, Tb=212.8,  Vb=32.9, vf=21.0),
    "Metanol":  dict(M=32.04,  sig=3.626, eps=481.8,  mu=1.70, Tb=337.8,  Vb=42.5, vf=29.9),
    "Cloruro de hidrogeno":  dict(M=36.46,  sig=3.282, eps=205.0,  mu=1.09,  Tb=188.2,  Vb=22.4, vf=15.0),
         "Cloruro de metilo":  dict(M=50.49,  sig=3.758, eps=305.0,  mu=1.90, Tb=238.6,  Vb=56.0, vf=40.0),
    "Acetona":  dict(M=58.08,  sig=4.600, eps=560.2,  mu=2.88, Tb=329.4,  Vb=77.5, vf=67.7),
}
#                    [g/mol],     [K],  [cm³/mol,Tb], phi (asociación), mu [cP a 25 °C], rho [g/cm³], sigma [dina/cm], par = paracoro
LIQ = {
    "Agua":     dict(M=18.015, Tb=373.15, Vb=18.9,  phi=2.6, mu=0.8903, rho=0.997, sigma=72.0, par=53.0,  water=True),
    "Anilina":   dict(M=93.13,  Tb=457.0,  Vb=101.0, phi=1.0, mu=3.5,    rho=1.021, sigma=38.0, par=200.0),
    "Metanol":  dict(M=32.04,  Tb=337.8,  Vb=42.5,  phi=1.9, mu=0.544,  rho=0.787, sigma=22.1, par=88.8),
    "Etanol":   dict(M=46.07,  Tb=351.4,  Vb=62.2,  phi=1.5, mu=1.074,  rho=0.785, sigma=22.0, par=127.0),
    "n-butanol":  dict(M=58.12,  Tb=272.7,  Vb=91.0,  phi=1.0, mu=0.294,  rho=0.620, sigma=17.9, par=145.0),
    "Glicerina":  dict(M=92.09,  Tb=563.0,  Vb=73.0,  phi=3.0, mu=1.412,  rho=1.261, sigma=63.4, par=290.0),
    "Acetona":  dict(M=58.08,  Tb=329.4,  Vb=77.5,  phi=1.0, mu=0.306,  rho=0.785, sigma=23.3, par=161.5),
    "Benceno":  dict(M=78.11,  Tb=353.2,  Vb=96.5,  phi=1.0, mu=0.604,  rho=0.874, sigma=28.2, par=205.3),
    "Tolueno":  dict(M=92.14,  Tb=383.8,  Vb=118.2, phi=1.0, mu=0.560,  rho=0.862, sigma=27.9, par=245.5),
    "n-Hexano": dict(M=86.18,  Tb=341.9,  Vb=140.6, phi=1.0, mu=0.294,  rho=0.655, sigma=17.9, par=271.0),
    "Ác. acético": dict(M=60.05, Tb=391.1, Vb=63.0, phi=1.0, mu=1.13,   rho=1.049, sigma=27.6, par=131.2),
    "Amoniaco":  dict(M=17.03,  Tb=239.7,  Vb=25.8,  phi=1.0, mu=0.240,  rho=0.681, sigma=20.0, par=37.7),
    "Freon-12":  dict(M=120.91, Tb=237.0,  Vb=90.0,  phi=1.0, mu=0.138,  rho=1.188, sigma=11.8, par=90.0),
    "Cloroformo": dict(M=119.38, Tb=334.0,  Vb=80.0,  phi=1.0, mu=0.563,  rho=1.489, sigma=27.1, par=160.0),
    "Acetato de etilo": dict(M=88.11, Tb=351.5, Vb=97.0, phi=1.0, mu=0.426,  rho=0.897, sigma=23.3, par=132.0),
    "Piridina": dict(M=79.10, Tb=388.0, Vb=88.0, phi=1.0, mu=0.248,  rho=0.981, sigma=20.8, par=150.0),
    
    # solutos gaseosos disueltos (solo como soluto A)
    "O2 (disuelto)":  dict(M=31.999, Tb=90.2,  Vb=25.6, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "CO2 (disuelto)": dict(M=44.01,  Tb=194.7, Vb=34.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "N2 (disuelto)":  dict(M=28.013, Tb=77.3,  Vb=31.2, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Acetona (disuelto)": dict(M=58.08, Tb=329.4, Vb=77.5, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Acido acetico (disuelto)": dict(M=60.05, Tb=391.1, Vb=63.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Acido benzoico (disuelto)": dict(M=122.12, Tb=523.0, Vb=122.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Bromobenceno (disuelto)": dict(M=157.01, Tb=332.0, Vb=69.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Ciclohexano (disuelto)": dict(M=84.16, Tb=353.8, Vb=108.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Acido formico (disuelto)": dict(M=46.03, Tb=338.0, Vb=40.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "n-heptano (disuelto)": dict(M=100.20, Tb=372.0, Vb=160.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Metil etil cetona (MEK) (disuelto)": dict(M=72.11, Tb=329.0, Vb=80.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Naftaleno (disuelto)": dict(M=128.17, Tb=353.0, Vb=100.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Tolueno (disuelto)": dict(M=92.14, Tb=383.8, Vb=118.2, phi=1.0, mu=None, rho=None, sigma=None, par=None),
    "Cloruro de vinilo (disuelto)": dict(M=62.50, Tb=350.0, Vb=60.0, phi=1.0, mu=None, rho=None, sigma=None, par=None),
}
# conductancias iónicas límite a 25 °C [cm²·S/eq] y carga |n|
IONES_CAT = {"H+": (349.8, 1), "Li+": (38.7, 1), "Na+": (50.1, 1), "K+": (73.5, 1),
             "NH4+": (73.5, 1), "Ca2+": (59.5, 2), "Mg2+": (53.1, 2), "Cu2+": (73.5, 2), "Zn2+": (53.0, 2), "Al3+": (51.0, 3)}
IONES_AN = {"Cl-": (76.3, 1), "Br-": (78.1, 1), "NO3-": (71.4, 1), "OH-": (198.0, 1), "C2O4 2-": (138.0, 2),
            "CH3COO-": (40.9, 1), "HCO3-": (44.5, 1), "SO4 2-": (80.0, 2)}

def mu_agua(T):
    return 0.02414 * 10 ** (247.8 / (T - 140.0))

def mu_de(sp, T):
    if sp.get("water"):
        return mu_agua(T), None
    if sp["mu"] is None:
        raise ValueError("Esta especie no tiene viscosidad: úsela solo como soluto.")
    adv = None
    if abs(T - 298.15) > 5:
        adv = ("La viscosidad de la base de datos es a 25 °C; a otra T ingrese el valor "
               "real de μ en 'Parámetros'.")
    return sp["mu"], adv


# ───────────────────────── FASE GASEOSA ─────────────────────────
def omega_D(Ts):
    """Integral de colisión de difusión (Neufeld)."""
    return (1.06036 / Ts ** 0.15610 + 0.19300 / math.exp(0.47635 * Ts)
            + 1.03587 / math.exp(1.52996 * Ts) + 1.76474 / math.exp(3.89411 * Ts))


def _lj(a, b):
    """Combina parámetros de Lennard-Jones (σ, ε/k) de dos especies."""
    return (a["sig"] + b["sig"]) / 2, math.sqrt(a["eps"] * b["eps"])


def gas_hbs(a, b, T, P):
    """Chapman-Enskog / Hirschfelder-Bird-Spotz. T [K], P [atm]."""
    s, e = _lj(a, b)
    return 0.001858 * T ** 1.5 * math.sqrt(1 / a["M"] + 1 / b["M"]) / (P * s ** 2 * omega_D(T / e))


def gas_wilke_lee(a, b, T, P):
    """Wilke-Lee. T [K], P [atm] (internamente en bar)."""
    s, e = _lj(a, b)
    Mab = 2 / (1 / a["M"] + 1 / b["M"])
    Pbar = P * 1.01325
    return (3.03 - 0.98 / math.sqrt(Mab)) * 1e-3 * T ** 1.5 / (Pbar * math.sqrt(Mab) * s ** 2 * omega_D(T / e))


def gas_fsg(a, b, T, P):
    """Fuller-Schettler-Giddings. T [K], P [atm]."""
    return (1e-3 * T ** 1.75 * math.sqrt(1 / a["M"] + 1 / b["M"])
            / (P * (a["vf"] ** (1 / 3) + b["vf"] ** (1 / 3)) ** 2))


def _stockmayer(g):
    d = 1.94e3 * g["mu"] ** 2 / (g["Vb"] * g["Tb"])
    return d, 1.18 * (1 + 1.3 * d ** 2) * g["Tb"], (1.585 * g["Vb"] / (1 + 1.3 * d ** 2)) ** (1 / 3)


def gas_brokaw(a, b, T, P):
    """Brokaw (gases polares). T [K], P [atm]."""
    da, ea, sa = _stockmayer(a)
    db, eb, sb = _stockmayer(b)
    dab, eab, sab = math.sqrt(da * db), math.sqrt(ea * eb), math.sqrt(sa * sb)
    Ts = T / eab
    om = omega_D(Ts) + 0.19 * dab ** 2 / Ts
    return 0.001858 * T ** 1.5 * math.sqrt(1 / a["M"] + 1 / b["M"]) / (P * sab ** 2 * om)


def gas_wilke_mezcla(D_Aj, y_j):
    """D_A,m de A difundiendo en mezcla estancada. y_j: fracciones de los demás componentes."""
    s = sum(y_j)
    return 1.0 / sum((y / s) / d for y, d in zip(y_j, D_Aj))


def gas_escalar(D1, T1, P1, T2, P2, eps_AB=None):
    """Escala D con T y P. Con eps_AB usa T^1.5·Ω(T1)/Ω(T2); si no, T^1.75 (FSG)."""
    if eps_AB is None:
        return D1 * (P1 / P2) * (T2 / T1) ** 1.75
    return D1 * (P1 / P2) * (T2 / T1) ** 1.5 * omega_D(T1 / eps_AB) / omega_D(T2 / eps_AB)


GAS_MODELOS = {"Wilke–Lee": gas_wilke_lee, "Hirschfelder–Bird–Spotz": gas_hbs,
               "Fuller–Schettler–Giddings (FSG)": gas_fsg, "Brokaw (polares)": gas_brokaw}


def gas_avisos(modelo, species, T, P):
    adv = []
    if P > 10:
        adv.append("P > 10 atm: los modelos son para gases a baja presión (comportamiento ideal).")
    pol = [n for n, g in species.items() if g["mu"] > 0.4]
    if pol and "Brokaw" not in modelo:
        adv.append(f"Momento dipolar elevado en {', '.join(pol)}: el modelo no polar puede ser "
                   "impreciso; use Brokaw.")
    if "Brokaw" in modelo and not pol:
        adv.append("Ningún componente es polar: Brokaw se reduce a un modelo no polar; "
                   "Wilke–Lee/FSG son suficientes.")
    if T < 200 or T > 1500:
        adv.append("T fuera del rango usual de validez (≈200–1500 K).")
    return adv


# ───────────────────────── FASE LÍQUIDA ─────────────────────────
def liq_wilke_chang(VA, MB, phi, T, muB):
    """Wilke-Chang. VA [cm³/mol a Tb], muB [cP]. -> cm²/s"""
    return 7.4e-8 * math.sqrt(phi * MB) * T / (muB * VA ** 0.6)


def liq_hayduk_minhas_acuoso(VA, T, muB):
    """Hayduk-Minhas, soluto en agua. -> cm²/s"""
    eps = 9.58 / VA - 1.12
    return 1.25e-8 * (VA ** -0.19 - 0.292) * T ** 1.52 * muB ** eps


def liq_hayduk_minhas_organico(VB, VA, T, muB, sigB):
    """Hayduk-Minhas, solvente orgánico. B = solvente, A = soluto. -> cm²/s"""
    return 1.55e-8 * VB ** 0.217 * T ** 1.29 * sigB ** 0.125 / (muB ** 0.92 * VA ** 0.45)


def _alpha(xA, A_margules):
    return 1 - 2 * A_margules * xA * (1 - xA)       # 1 + dlnγA/dlnxA (Margules 1 parámetro)


def liq_vignes(D0AB, D0BA, xA, alpha=1.0):
    return D0AB ** (1 - xA) * D0BA ** xA * alpha


def liq_leffler_cullinan(D0AB, D0BA, xA, muA, muB, mu_m, alpha=1.0):
    return (D0AB * muB) ** (1 - xA) * (D0BA * muA) ** xA * alpha / mu_m


def liq_stokes_einstein(T, mu_cP, r_m):
    """D = kT/(6πμr). r en m. -> cm²/s"""
    return 1.380649e-23 * T / (6 * math.pi * mu_cP * 1e-3 * r_m) * 1e4


def liq_polson(T, mu_cP, MA):
    """Polson (proteínas globulares, M>1000). -> cm²/s"""
    return 9.40e-15 * T / (mu_cP * 1e-3 * MA ** (1 / 3)) * 1e4


def liq_nernst_haskell(T, lam_p, n_p, lam_n, n_n):
    """Electrolito fuerte diluido. λ en cm²·S/eq. -> cm²/s"""
    return R_FARADAY2 * T * (1 / n_p + 1 / n_n) / (1 / lam_p + 1 / lam_n)


def liq_perkins_geankoplis(D_Aj, x_j, mu_j, mu_m):
    """D_AM·μ_M^0.8 = Σ x_j' D°_Aj μ_j^0.8  (x_j' fracciones libres de soluto)."""
    s = sum(x_j)
    return sum((x / s) * d * m ** 0.8 for x, d, m in zip(x_j, D_Aj, mu_j)) / mu_m ** 0.8


def liq_escalar(D1, T1, mu1, T2, mu2):
    """D·μ/T = constante (Stokes-Einstein)."""
    return D1 * (T2 / T1) * (mu1 / mu2)


def liq_tyn_calus(D1, T1, T2, Tc, n):
    """Tyn–Calus (Welty): D2/D1 = ((Tc − T1)/(Tc − T2))^n. Tc del solvente [K]."""
    return D1 * ((Tc - T1) / (Tc - T2)) ** n


# ───────────────────────── VALIDACIÓN ─────────────────────────
def correr_validacion():
    """Compara la app con las respuestas de clase (Answers). D en cm²/s (1 m²/s = 1e4 cm²/s).
    Devuelve (descripción, calculado, referencia, error %, estado)."""
    T = 298.15
    TC_AGUA, N_AGUA = 647.3, 6.0
    # --- Sistema A: etanol en agua, x = 0.30, 25 °C (H-M + Vignes, Γ = 0.45)
    D12 = liq_hayduk_minhas_acuoso(LIQ["Etanol"]["Vb"], T, mu_agua(T))
    D21 = liq_hayduk_minhas_organico(LIQ["Etanol"]["Vb"], LIQ["Agua"]["Vb"], T,
                                     LIQ["Etanol"]["mu"], LIQ["Etanol"]["sigma"])
    DA = liq_vignes(D12, D21, 0.30, 0.45)
    # --- Sistema B: BSA, 37 °C (Polson, M = 66500)
    TB = 310.15
    DB = liq_polson(TB, mu_agua(TB), 66500.0)
    # --- Sistema C: acetona diluida, 15 °C (Wilke-Chang)
    TCc = 288.15
    DC = liq_wilke_chang(LIQ["Acetona (disuelto)"]["Vb"], LIQ["Agua"]["M"], LIQ["Agua"]["phi"],
                         TCc, mu_agua(TCc))
    T50 = 323.15
    casos = [
        ("D°12 etanol en agua, H-M acuoso, 25 °C", D12, 1.335e-5),
        ("D°21 agua en etanol, H-M orgánico, 25 °C", D21, 2.168e-5),
        ("Sistema A: D a x=0.30, 25 °C", DA, 6.950e-6),
        ("Sistema B: BSA, 37 °C", DB, 1.035e-6),
        ("Sistema C: acetona, 15 °C", DC, 9.686e-6),
        ("Sistema A a 50 °C (Tyn–Calus, n=6)", liq_tyn_calus(DA, T, T50, TC_AGUA, N_AGUA), 1.086e-5),
        ("Sistema B a 50 °C (Tyn–Calus, n=6)", liq_tyn_calus(DB, TB, T50, TC_AGUA, N_AGUA), 1.311e-6),
        ("Sistema C a 50 °C (Tyn–Calus, n=6)", liq_tyn_calus(DC, TCc, T50, TC_AGUA, N_AGUA), 1.792e-5),
    ]
    out = []
    for d, c, r in casos:
        e = (c - r) / r * 100
        out.append((d, c, r, e, "OK" if abs(e) < 5 else "REVISAR"))
    return out


# ───────────────────────── INTERFAZ STREAMLIT ─────────────────────────
def main():
    import streamlit as st

    st.set_page_config(page_title="Difusión molar", layout="wide")
    st.title("Estimación de coeficientes de difusión molar")

    # ---- entradas globales
    sb = st.sidebar
    sb.header("Condiciones de operación")
    tu = sb.selectbox("Unidad de T", ["K", "°C", "°F"], 1)
    tv = sb.number_input(f"Temperatura [{tu}]", value=25.0 if tu == "°C" else (298.15 if tu == "K" else 77.0))
    T = tv if tu == "K" else (tv + 273.15 if tu == "°C" else (tv - 32) * 5 / 9 + 273.15)
    pu = sb.selectbox("Unidad de P", list(P_ATM), 0)
    P = sb.number_input(f"Presión [{pu}]", value=1.0, min_value=1e-6) * P_ATM[pu]
    sb.caption(f"T = {T:.2f} K | P = {P:.4f} atm")
    unidad = sb.selectbox("Unidad de salida", list(UNITS), 1)

    def editor(nombre, sp, campos, key):
        sp = dict(sp)
        with st.expander(f"Parámetros de {nombre} (editables)"):
            cols = st.columns(len(campos))
            for c, (k, lab) in zip(cols, campos):
                if sp.get(k) is not None:
                    sp[k] = c.number_input(lab, value=float(sp[k]), format="%.4f", key=f"{key}_{nombre}_{k}")
        return sp

    def mostrar(D, avisos):
        st.success(f"**D = {D:.4e} cm²/s**   →   {D * 1e-4 * UNITS[unidad]:.4e} {unidad}")
        st.table({u: [f"{D * 1e-4 * f:.4e}"] for u, f in UNITS.items()})
        for a in avisos:
            st.warning(a)

    GC = [("M", "M [g/mol]"), ("sig", "σ [Å]"), ("eps", "ε/k [K]"), ("mu", "μp [D]"),
          ("Tb", "Tb [K]"), ("Vb", "Vb [cm³/mol]"), ("vf", "Σv FSG")]
    LC = [("M", "M [g/mol]"), ("Tb", "Tb [K]"), ("Vb", "V_b [cm³/mol]"), ("phi", "φ"),
          ("mu", "μ [cP]"), ("par", "Paracoro"), ("sigma", "σ [dina/cm]")]

    tabs = st.tabs(["Cálculo", "Escalado T y P"])

    # ================= CÁLCULO =================
    with tabs[0]:
        estado = st.radio("Estado físico", ["Gas", "Líquido"], horizontal=True)
        ncomp = st.radio("Número de componentes", ["Binario", "Ternario"], horizontal=True)

        # ---------- GAS ----------
        if estado == "Gas":
            modelo = st.selectbox("Modelo", list(GAS_MODELOS))
            nm = list(GAS)
            A = st.selectbox("Componente A (soluto)", nm, nm.index("Dioxido de carbono"))
            B = st.selectbox("Componente B", nm, nm.index("Aire"))
            comps = [A, B]
            if ncomp == "Ternario":
                C = st.selectbox("Componente C", nm, nm.index("Nitrogeno"))
                comps.append(C)
            if len(set(comps)) < len(comps):
                st.error("Los componentes deben ser distintos."); st.stop()
            sp = {n: editor(n, GAS[n], GC, "g") for n in comps}
            f = GAS_MODELOS[modelo]
            avisos = gas_avisos(modelo, sp, T, P)
            D_AB = f(sp[A], sp[B], T, P)
            if ncomp == "Binario":
                mostrar(D_AB, avisos)
            else:
                st.subheader("Composición molar (A difunde en la mezcla B + C estancada)")
                c1, c2 = st.columns(2)
                yA = c1.number_input("y_A", 0.0, 1.0, 0.2, 0.01)
                yB = c2.number_input("y_B", 0.0, 1.0 - yA, min(0.5, 1 - yA), 0.01)
                yC = 1 - yA - yB
                st.caption(f"y_C = {yC:.3f}")
                if yB + yC <= 0:
                    st.error("Se requiere y_B + y_C > 0."); st.stop()
                D_AC = f(sp[A], sp[C], T, P)
                st.write(f"D_AB = {D_AB:.4e} cm²/s | D_AC = {D_AC:.4e} cm²/s")
                mostrar(gas_wilke_mezcla([D_AB, D_AC], [yB, yC]), avisos)

        # ---------- LÍQUIDO ----------
        else:
            nm = list(LIQ)
            solv = [n for n in nm if LIQ[n]["mu"] is not None]
            if ncomp == "Binario":
                modelo = st.selectbox("Modelo", [
                    "Wilke–Chang", "Hayduk–Minhas", "Vignes / Leffler–Cullinan (concentradas)",
                    "Stokes–Einstein", "Polson (proteínas)", "Nernst–Haskell (electrolitos)"])
            else:
                modelo = "Perkins–Geankoplis"
                st.info("Ternario líquido: soluto diluido en mezcla de solventes (Perkins–Geankoplis).")
            avisos = []

            if modelo == "Wilke–Chang":
                A = st.selectbox("Soluto A", nm, nm.index("CO2 (disuelto)"))
                B = st.selectbox("Solvente B", solv, 0)
                xA = st.number_input("x_A", 0.0, 1.0, 0.01, 0.005, format="%.4f")
                a, b = editor(A, LIQ[A], LC, "l"), editor(B, LIQ[B], LC, "l")
                muB, w = mu_de(b, T); muB = b["mu"] if not b.get("water") else muB
                avisos += [w] if w else []
                if xA > 0.1: avisos.append("x_A > 0.1: el modelo es para solutos diluidos.")
                mostrar(liq_wilke_chang(a["Vb"], b["M"], b["phi"], T, muB), avisos)

            elif modelo == "Hayduk–Minhas":
                tipo = st.radio("Caso", ["Acuoso (solvente = agua)", "Orgánico no acuoso"], horizontal=True)
                A = st.selectbox("Soluto A", nm, nm.index("CO2 (disuelto)"))
                a = editor(A, LIQ[A], LC, "l")
                if tipo.startswith("Acuoso"):
                    b = editor("Agua", LIQ["Agua"], LC, "l")
                    muB, w = mu_de(b, T)
                    avisos += [w] if w else []
                    mostrar(liq_hayduk_minhas_acuoso(a["Vb"], T, muB), avisos)
                else:
                    B = st.selectbox("Solvente orgánico B", [s for s in solv if s != "Agua"], 3)
                    b = editor(B, LIQ[B], LC, "l")
                    muB, w = mu_de(b, T)
                    avisos += [w] if w else []
                    mostrar(liq_hayduk_minhas_organico(b["Vb"], a["Vb"], T, muB, b["sigma"]), avisos)

            elif modelo.startswith("Vignes"):
                A = st.selectbox("Componente A", solv, solv.index("Etanol"))
                B = st.selectbox("Componente B", solv, 0)
                if A == B: st.error("A y B deben ser distintos."); st.stop()
                xA = st.slider("x_A", 0.0, 1.0, 0.5, 0.01)
                a, b = editor(A, LIQ[A], LC, "l"), editor(B, LIQ[B], LC, "l")
                muA, wa = mu_de(a, T); muB, wb = mu_de(b, T)
                avisos += [w for w in (wa, wb) if w]
                D0AB = liq_wilke_chang(a["Vb"], b["M"], b["phi"], T, muB)   # A infinitamente diluido en B
                D0BA = liq_wilke_chang(b["Vb"], a["M"], a["phi"], T, muA)
                c1, c2 = st.columns(2)
                D0AB = c1.number_input("D°_AB [cm²/s] (Wilke–Chang por defecto)", value=D0AB, format="%.4e")
                D0BA = c2.number_input("D°_BA [cm²/s] (Wilke–Chang por defecto)", value=D0BA, format="%.4e")
                mum = st.number_input("μ de la mezcla [cP] (por defecto media log.)",
                                      value=math.exp(xA * math.log(muA) + (1 - xA) * math.log(muB)), format="%.4f")
                Am = st.number_input("Parámetro de Margules A (0 = mezcla ideal, α = 1)", value=0.0)
                al = _alpha(xA, Am)
                st.caption(f"α = 1 + dlnγ_A/dlnx_A = {al:.4f}")
                if al <= 0: avisos.append("α ≤ 0: mezcla inestable (dentro de la zona de inmiscibilidad).")
                if Am == 0: avisos.append("Mezcla supuesta ideal (α = 1); si no lo es, ingrese Margules/γ.")
                mostrar(liq_vignes(D0AB, D0BA, xA, al), avisos)
                st.write(f"Leffler–Cullinan: **{liq_leffler_cullinan(D0AB, D0BA, xA, muA, muB, mum, al):.4e} cm²/s**")

            elif modelo == "Stokes–Einstein":
                muB = st.number_input("Viscosidad del solvente [cP]", value=mu_agua(T), format="%.4f")
                modo = st.radio("Tamaño", ["Radio [nm]", "Estimar con M y ρ"], horizontal=True)
                if modo.startswith("Radio"):
                    r = st.number_input("Radio hidrodinámico [nm]", value=2.0, min_value=0.01) * 1e-9
                    M = st.number_input("M [g/mol] (solo para validar)", value=5000.0)
                else:
                    M = st.number_input("M [g/mol]", value=5000.0)
                    rho = st.number_input("ρ [g/cm³]", value=1.35)
                    r = (3 * M / (4 * math.pi * 6.02214076e23 * rho)) ** (1 / 3) * 1e-2
                    st.caption(f"r = {r * 1e9:.3f} nm")
                if M < 1000: avisos.append("M < 1000 g/mol: el modelo es para macromoléculas esféricas grandes.")
                if 1e9 * r < 0.5: avisos.append("Radio pequeño: el continuo no es válido (use Wilke–Chang).")
                mostrar(liq_stokes_einstein(T, muB, r), avisos)

            elif modelo.startswith("Polson"):
                MA = st.number_input("M de la proteína [g/mol]", value=66500.0)
                muB = st.number_input("Viscosidad del solvente [cP]", value=mu_agua(T), format="%.4f")
                if MA <= 1000: avisos.append("Polson requiere M_A > 1000 g/mol.")
                avisos.append("Válido para proteínas globulares diluidas.")
                mostrar(liq_polson(T, muB, MA), avisos)

            elif modelo.startswith("Nernst"):
                cat = st.selectbox("Catión", list(IONES_CAT)); an = st.selectbox("Anión", list(IONES_AN))
                lp, np_ = IONES_CAT[cat]; ln_, nn = IONES_AN[an]
                lp = st.number_input("λ+ [cm²·S/eq]", value=lp); ln_ = st.number_input("λ− [cm²·S/eq]", value=ln_)
                c = st.number_input("Concentración [mol/L]", value=0.01, format="%.4f")
                if c > 0.1: avisos.append("c > 0.1 M: Nernst–Haskell es para dilución infinita/electrolito diluido.")
                if abs(T - 298.15) > 5: avisos.append("Los λ listados son a 25 °C; corrija λ± a la T de trabajo.")
                mostrar(liq_nernst_haskell(T, lp, np_, ln_, nn), avisos)

            else:  # Perkins–Geankoplis
                A = st.selectbox("Soluto A", nm, nm.index("CO2 (disuelto)"))
                B = st.selectbox("Solvente B", solv, solv.index("Agua"))
                C = st.selectbox("Solvente C", solv, solv.index("Etanol"))
                if len({B, C}) < 2: st.error("B y C deben ser distintos."); st.stop()
                xB = st.slider("x_B (fracción de B en la mezcla libre de soluto)", 0.0, 1.0, 0.5, 0.01)
                xA = st.number_input("x_A (soluto en la solución real)", 0.0, 1.0, 0.01, 0.005, format="%.4f")
                a = editor(A, LIQ[A], LC, "l")
                b, c_ = editor(B, LIQ[B], LC, "l"), editor(C, LIQ[C], LC, "l")
                (muB, wb), (muC, wc) = mu_de(b, T), mu_de(c_, T)
                avisos += [w for w in (wb, wc) if w]
                D_AB = liq_wilke_chang(a["Vb"], b["M"], b["phi"], T, muB)
                D_AC = liq_wilke_chang(a["Vb"], c_["M"], c_["phi"], T, muC)
                mum = st.number_input("μ de la mezcla de solventes [cP] (por defecto media log.)",
                                      value=math.exp(xB * math.log(muB) + (1 - xB) * math.log(muC)), format="%.4f")
                if xA > 0.1: avisos.append("x_A > 0.1: Perkins–Geankoplis es para soluto diluido.")
                avisos.append("D°_Aj estimados con Wilke–Chang; la media log. de μ es solo una aproximación.")
                st.write(f"D°_AB = {D_AB:.4e} | D°_AC = {D_AC:.4e} cm²/s")
                mostrar(liq_perkins_geankoplis([D_AB, D_AC], [xB, 1 - xB], [muB, muC], mum), avisos)

    # ================= ESCALADO =================
    with tabs[1]:
        st.subheader("Corrección de D con T y P a partir de un valor de referencia")
        est = st.radio("Fase", ["Gas", "Líquido"], horizontal=True, key="esc")
        D1 = st.number_input("D de referencia [cm²/s]", value=0.16, format="%.4e")
        c1, c2 = st.columns(2)
        T1 = c1.number_input("T1 [K]", value=298.15); T2 = c2.number_input("T2 [K]", value=T)
        if est == "Gas":
            P1 = c1.number_input("P1 [atm]", value=1.0); P2 = c2.number_input("P2 [atm]", value=P)
            met = st.radio("Dependencia con T", ["T^1.75 (tipo FSG)", "T^1.5·Ω(T1)/Ω(T2) (Chapman–Enskog)"])
            eps = None
            if "Ω" in met:
                nm = list(GAS); a1 = st.selectbox("A", nm, nm.index("Dioxido de carbono"), key="ea")
                b1 = st.selectbox("B", nm, nm.index("Aire"), key="eb")
                eps = _lj(GAS[a1], GAS[b1])[1]
            st.success(f"D(T2,P2) = {gas_escalar(D1, T1, P1, T2, P2, eps):.4e} cm²/s")
            if max(P1, P2) > 10: st.warning("P > 10 atm: la relación D ∝ 1/P deja de ser válida.")
        else:
            met_l = st.radio("Método", ["Tyn–Calus / Welty: ((Tc−T1)/(Tc−T2))^n",
                                        "Stokes–Einstein: D·μ/T = cte"])
            if met_l.startswith("Tyn"):
                Tc = c1.number_input("Tc del solvente [K] (agua = 647.3)", value=647.3)
                n_ = c2.number_input("n (según ΔHvap; agua = 6)", value=6.0)
                if T1 >= Tc or T2 >= Tc:
                    st.error("T1 y T2 deben ser menores que Tc.")
                else:
                    st.success(f"D(T2) = {liq_tyn_calus(D1, T1, T2, Tc, n_):.4e} cm²/s")
            else:
                m1 = c1.number_input("μ1 [cP] a T1", value=0.8903, format="%.4f")
                m2 = c2.number_input("μ2 [cP] a T2", value=mu_agua(T), format="%.4f")
                st.success(f"D(T2) = {liq_escalar(D1, T1, m1, T2, m2):.4e} cm²/s")
            st.info("En líquidos el efecto de P es despreciable a presiones moderadas.")


if __name__ == "__main__":
    if "--test" in sys.argv:
        for d, c, r, e, s in correr_validacion():
            print(f"{d:55s} calc={c:.4e} ref={r:.4e} err={e:+.1f}%  {s}")
    else:
        main()
