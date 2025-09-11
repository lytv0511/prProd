# forces.py
import numpy as np
import math

def hat(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 0 else v

# Surfaces
class Plane:
    """Finite square plane (center r0, unit normal n, half-size 'size')."""
    def __init__(self, r0, n, name="plane", size=5.0):
        self.r0 = np.array(r0, dtype=float)
        self.n = hat(np.array(n, dtype=float))
        self.name = name
        self.size = float(size)

    def gap(self, r, R):
        """Signed gap: positive if ball center is outside contact (above), <=0 if penetrating/contact.
           If the contact point projects outside the finite square, return +inf (no collision).
        """
        r_rel = r - self.r0
        # plane-local tangent basis
        if abs(self.n[0]) < 1e-6 and abs(self.n[1]) < 1e-6:
            v1 = np.array([1.0, 0.0, 0.0])
        else:
            v1 = np.cross(self.n, [0.0, 0.0, 1.0])
            v1 = v1 / np.linalg.norm(v1)
        v2 = np.cross(self.n, v1)
        v2 = v2 / np.linalg.norm(v2)
        proj1 = np.dot(r_rel, v1)
        proj2 = np.dot(r_rel, v2)
        if abs(proj1) > self.size or abs(proj2) > self.size:
            return float('inf')
        # distance from plane along normal minus ball radius
        return np.dot(r_rel, self.n) - R

    def normal(self, r):
        return self.n

# Ball and aerodynamics
class Ball:
    def __init__(self, m=0.045, R=0.021, I_type="solid"):
        self.m = float(m)
        self.R = float(R)
        if I_type == "solid":
            self.kappa = 2.0/5.0
        elif I_type == "shell":
            self.kappa = 2.0/3.0
        else:
            self.kappa = 2.0/5.0
        self.I = self.kappa * self.m * self.R**2
        self.V = 4.0/3.0 * math.pi * self.R**3

class Aero:
    def __init__(self, rho=1.225, Cd=0.3, Cl0=0.2, spin_cap=0.5, Cm=0.02, wind=(0,0,0)):
        self.rho = float(rho)
        self.Cd = float(Cd)
        self.Cl0 = float(Cl0)
        self.spin_cap = float(spin_cap)
        self.Cm = float(Cm)
        self.wind = np.array(wind, dtype=float)

# Forces and torques
def forces_and_torque(r, v, w, ball, aero, g_vec=np.array([0.0, 0.0, -9.81])):
    """Return F (3,) and tau (3,) acting on ball. z-axis is vertical here (consistent with main).
       Note: main uses z as vertical, so gravity is -9.81 in z.
    """
    m, R, Vball = ball.m, ball.R, ball.V
    v_rel = v - aero.wind
    Vmag = np.linalg.norm(v_rel)

    # gravity (acting on negative z), buoyancy (acting on positive z)
    Fg = m * g_vec
    Fb = -aero.rho * Vball * g_vec

    # drag (quadratic)
    Fd = np.zeros(3)
    if Vmag > 1e-8:
        A = math.pi * R**2
        Fd = -0.5 * aero.rho * aero.Cd * A * Vmag * v_rel

    # Magnus (lift) - direction ~ cross(v, w)
    Fm = np.zeros(3)
    if Vmag > 1e-8 and np.linalg.norm(w) > 1e-8:
        S = (R * np.linalg.norm(w)) / Vmag
        Cl = np.clip(aero.Cl0 * S, -aero.spin_cap, aero.spin_cap)
        dir_vec = np.cross(v_rel, w)
        nd = np.linalg.norm(dir_vec)
        if nd > 0:
            e_dir = dir_vec / nd
            Fm = 0.5 * aero.rho * Cl * math.pi * R**2 * Vmag**2 * e_dir

    # aerodynamic spin damping torque (simple model)
    tau_air = -0.5 * aero.rho * aero.Cm * math.pi * R**3 * (Vmag**2) * w

    F = Fg + Fb + Fd + Fm
    return F, tau_air

# Integrator (RK4)
def rk4_step(state, dt, deriv):
    r, v, w = state
    def add(s, k, a):
        return (s[0] + a*k[0], s[1] + a*k[1], s[2] + a*k[2])
    k1 = deriv(state)
    k2 = deriv(add(state, k1, dt/2.0))
    k3 = deriv(add(state, k2, dt/2.0))
    k4 = deriv(add(state, k3, dt))
    dr = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) * (dt/6.0)
    dv = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) * (dt/6.0)
    dw = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) * (dt/6.0)
    return (r + dr, v + dv, w + dw)

# Collision impulse (sphere vs massive surface)
def collide(v_in, w_in, n, ball, e_n=0.5, mu=0.2):
    """Return v_out, w_out after instantaneous contact with surface having unit normal n.
       Uses normal restitution e_n and Coulomb friction mu. The contact point is at -R*n.
    """
    n = hat(n)
    m, R, kappa, I = ball.m, ball.R, ball.kappa, ball.I

    # velocity of contact point (center velocity + omega x r_c; r_c = -R n)
    r_c = -R * n
    v_c = v_in + np.cross(w_in, r_c)   # linear + rotational
    # decompose into normal and tangential parts
    u_n = np.dot(v_c, n)
    u_t_vec = v_c - u_n * n
    u_t = np.linalg.norm(u_t_vec)

    # normal impulse magnitude using restitution
    Jn = - (1 + e_n) * m * u_n

    # effective mass for tangential impulse (for sphere vs fixed plane)
    m_eff = m / (1.0 + kappa)   # derivation: contact impulse partitioning
    Jt_stick = - m_eff * u_t_vec
    Jt_max = mu * abs(Jn)
    Jt_mag = np.linalg.norm(Jt_stick)

    if Jt_mag <= Jt_max:
        Jt = Jt_stick
    else:
        if Jt_mag > 0:
            Jt = - Jt_max * (Jt_stick / Jt_mag)
        else:
            Jt = np.zeros(3)

    # total impulse
    J = Jn * n + Jt

    # update linear and angular velocities
    v_out = v_in + J / m
    # angular impulse: delta L = r_c x J  -> delta w = (r_c x J) / I
    dw = np.cross(r_c, J) / I
    w_out = w_in + dw

    return v_out, w_out

class SphereSurface:
    """Spherical surface obstacle (center c, radius Rs)."""
    def __init__(self, c, Rs, name="sphere"):
        self.c = np.array(c, dtype=float)
        self.Rs = float(Rs)
        self.name = name

    def gap(self, r, R):
        """Signed gap: distance between ball center and sphere surface minus (Rs + R)."""
        d = np.linalg.norm(r - self.c)
        return d - (self.Rs + R)

    def normal(self, r):
        """Normal points outward from sphere center to contact point."""
        diff = r - self.c
        nrm = np.linalg.norm(diff)
        if nrm < 1e-12:
            return np.array([0.0, 0.0, 1.0])  # fallback
        return diff / nrm