import numpy as np
import matplotlib.pyplot as plt
from forces import Ball, Aero, Plane, SphereSurface
from bouncesim import Simulator

def rotate_vector(vec, axis, angle):
    """
    Rotate a vector 'vec' around 'axis' by 'angle' radians using Rodrigues' formula.
    """
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)
    vec = np.asarray(vec)
    return (vec * np.cos(angle) +
            np.cross(axis, vec) * np.sin(angle) +
            axis * np.dot(axis, vec) * (1 - np.cos(angle)))

def perturb_normal(normal, max_angle_deg, rng):
    """
    Perturb a normal vector by a random angle up to max_angle_deg.
    """
    normal = np.asarray(normal)
    normal = normal / np.linalg.norm(normal)
    # sample angle from 0 to max_angle_deg (in radians)
    max_angle_rad = np.deg2rad(max_angle_deg)
    theta = rng.uniform(0, max_angle_rad)
    phi = rng.uniform(0, 2 * np.pi)
    # Find a random perpendicular axis
    if np.allclose(normal, [0, 0, 1]):
        perp = np.array([1, 0, 0])
    else:
        perp = np.cross(normal, [0, 0, 1])
        perp /= np.linalg.norm(perp)
    # rotate normal by theta around perp, then by phi around normal
    perturbed = rotate_vector(normal, perp, theta)
    perturbed = rotate_vector(perturbed, normal, phi)
    return perturbed

def build_perturbed_surfaces(surfaces, max_angle_deg, rng):
    """
    Given a list of surfaces, return a new list with perturbed normals.
    Only planar surfaces are perturbed (for SphereSurface, you could implement a similar logic if needed).
    """
    perturbed = []
    for surf in surfaces:
        if isinstance(surf, Plane):
            new_normal = perturb_normal(surf.n, max_angle_deg, rng)
            perturbed.append(Plane(surf.r0, new_normal))
        elif isinstance(surf, SphereSurface):
            # For spheres, maybe perturb the center by a small amount if desired
            perturbed.append(SphereSurface(surf.c, surf.Rs))
        else:
            perturbed.append(surf)
    return perturbed

def make_sim_factory(ball, aero, surfaces, e_n, mu, dt, t_max, max_bounces):
    """
    Returns a function that creates a Simulator with the given parameters, rebuilding surfaces and optionally perturbing aero.
    """
    def factory(perturb_deg=0, rng=None, perturb_aero=False):
        perturbed_surfs = build_perturbed_surfaces(surfaces, perturb_deg, rng) if perturb_deg > 0 else surfaces
        used_aero = aero
        # Optionally, you could perturb aero parameters here if desired and perturb_aero is True.
        return Simulator(
            ball=ball,
            aero=used_aero,
            surfaces=perturbed_surfs,
            e_n=e_n,
            mu=mu,
            dt=dt,
            t_max=t_max,
            max_bounces=max_bounces
        )
    return factory

def _check_required_sequence(bounces, required_sequence):
    """
    Check if the sequence of surface names in bounces matches required_sequence.
    Each bounce is a tuple (time, name, pos), so use name.
    """
    if required_sequence is None:
        return True
    bounce_names = [b[1] for b in bounces]
    return bounce_names[:len(required_sequence)] == list(required_sequence)

def _get_surface_position_from_bounces(bounces, surface_name):
    """
    Return the position of the first bounce on the specified surface_name, or None.
    Each bounce is a tuple (time, name, pos).
    """
    for b in bounces:
        if b[1] == surface_name:
            return b[2]
    return None

def run_monte_carlo(r0, v0, w0, sim_factory, num_samples, perturb_deg=0, required_sequence=None, surface_index=None, seed=None):
    """
    Run a Monte Carlo simulation, perturbing surfaces for each sample.
    Returns a list of positions (np.array) for the specified surface_index.
    """
    rng = np.random.default_rng(seed)
    positions = []
    base_sim = sim_factory()
    base_surfaces = base_sim.surfaces
    for i in range(num_samples):
        perturbed_surfaces = build_perturbed_surfaces(base_surfaces, perturb_deg, rng)
        sim = Simulator(ball=base_sim.ball,
                        aero=base_sim.aero,
                        surfaces=perturbed_surfaces,
                        e_n=base_sim.e_n,
                        mu=base_sim.mu,
                        dt=base_sim.dt,
                        t_max=base_sim.t_max,
                        max_bounces=base_sim.max_bounces)
        _, _, _, bounces, _ = sim.integrate(r0, v0, w0)
        if not _check_required_sequence(bounces, required_sequence):
            continue
        pos = _get_surface_position_from_bounces(bounces, surface_index)
        if pos is not None:
            positions.append(np.array(pos))
    return np.array(positions)

def plot_dispersion(positions, title='Dispersion', show=True, ax=None):
    """
    Plot the 2D scatter of impact positions.
    """
    if positions is None or len(positions) == 0:
        print("No positions to plot.")
        return
    if ax is None:
        fig, ax = plt.subplots()
    ax.scatter(positions[:,0], positions[:,1], alpha=0.6)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(title)
    ax.axis('equal')
    if show:
        plt.show()

if __name__ == "__main__":
    # Example setup (same as in main.py)
    ball = Ball(m=0.045, R=0.021, I_type="solid")
    aero = Aero(rho=1.225, Cd=0.28, Cl0=0.2, spin_cap=0.7, Cm=0.02, wind=(0.0,0.0,0.0))

    # Surfaces (copy from main.py)
    plane1 = Plane(r0=(12.0, -1.5, 4.0), n=(-0.9, 0.435, 0.0), name="plane1", size=4.0)
    plane2 = Plane(r0=(6.0, 0.5, 5.5), n=(0.866, 0.5, 0.0), name="plane2", size=4.0)
    sphere1 = SphereSurface(c=(10.0, 3.0, 4.0), Rs=1.0, name="sphere1")

    surfaces = [plane1, plane2, sphere1]

    # Initial state
    start_pos = np.array([4.0, 10.0, 6.0])
    target_pos = np.array([14.0, -2.0, 5.5])
    direction = (target_pos - start_pos) / np.linalg.norm(target_pos - start_pos)
    v0 = direction * 40.0
    w0 = np.array([0.0, 0.0, 0.0])

    # Make sim factory
    factory = make_sim_factory(ball, aero, surfaces,
                               e_n=0.45, mu=0.35,
                               dt=1e-4, t_max=6.0, max_bounces=20)

    # Run Monte Carlo
    hits = run_monte_carlo(start_pos, v0, w0, factory, num_samples=10,
                           perturb_deg=0.5,
                           required_sequence=["plane1", "plane2", "sphere1"],
                           surface_index="sphere1",
                           seed=42)

    print(f"Monte Carlo finished, {len(hits)} successful runs")

    # Plot
    plot_dispersion(hits, title="Monte Carlo dispersion on sphere1")