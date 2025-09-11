# main.py
import numpy as np
import matplotlib.pyplot as plt
from forces import Ball, Aero, Plane, SphereSurface
from bouncesim import Simulator

def plot_sphere(ax, c, Rs, color='green', alpha=0.3, resolution=30):
    u = np.linspace(0, 2*np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = c[0] + Rs * np.outer(np.cos(u), np.sin(v))
    y = c[1] + Rs * np.outer(np.sin(u), np.sin(v))
    z = c[2] + Rs * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, shade=True)

def plot_plane(ax, r0, n, size=5.0, color='blue', alpha=0.3):
    n = np.array(n, dtype=float)
    # pick a vector not parallel to n to construct basis
    if abs(n[0]) < 1e-6 and abs(n[1]) < 1e-6:
        v1 = np.array([1.0, 0.0, 0.0])
    else:
        v1 = np.cross(n, [0.0, 0.0, 1.0])
    v1 = v1 / np.linalg.norm(v1)
    v2 = np.cross(n, v1)
    v2 = v2 / np.linalg.norm(v2)

    grid = np.linspace(-size, size, 20)
    a, b = np.meshgrid(grid, grid)
    pts = r0.reshape(3,1,1) + v1.reshape(3,1,1)*a + v2.reshape(3,1,1)*b
    ax.plot_surface(pts[0], pts[1], pts[2], color=color, alpha=alpha, linewidth=0, shade=True)

if __name__ == "__main__":
    # Parameters
    z_start = 10.0 # starting height (m)
    plane_size = 4.0 # half-width of plane (m)
    speed = 40.0 # initial speed in m/s

    # Ball and aero
    ball = Ball(m=0.045, R=0.021, I_type="solid") # 45 g (golf ball)
    aero = Aero(rho=1.225, Cd=0.28, Cl0=0.2, spin_cap=0.7, Cm=0.02, wind=(0.0,0.0,0.0))

    # Initial state (for ball)
    start_pos = np.array([4.0, 10.0, 6.0])
    
    target_pos = np.array([14.0, -2.0, 5.5]) # tweak for position of ball target direction

    direction = target_pos - start_pos
    direction = direction / np.linalg.norm(direction)
    v0 = direction * speed
    w0 = np.array([0.0, 0.0, 0.0]) # initial spin

    # Raw wooden planes
    plane1 = Plane(r0=(12.0, -1.5, 4.0), n=(-0.9, 0.435, 0.0), name="plane1", size=plane_size)
    # plane2 after bouncing from plane1
    plane2 = Plane(r0=(6.0, 0.5, 5.5), n=(0.866, 0.5, 0.0), name="plane2", size=plane_size)

    # Define the spherical obstacle
    sphere1 = SphereSurface(c=(5.0, 8.0, 0.0), Rs=2.0, name="sphere1")

    # Add it to surfaces
    surfaces = [plane1, plane2, sphere1]

    # Run sim (no changes needed!)
    sim = Simulator(ball=ball, aero=aero, surfaces=surfaces, e_n=0.45, mu=0.35,
                    dt=1e-4, t_max=6.0, max_bounces=20)

    # Plot it

    # wood is relatively inelastic: e_n ~ 0.45, moderate friction mu ~ 0.35
    sim = Simulator(ball=ball, aero=aero, surfaces=surfaces, e_n=0.45, mu=0.35,
                    dt=1e-4, t_max=6.0, max_bounces=20)

    times, traj, seg_ids, bounces, target_hits = sim.integrate(start_pos, v0, w0)

    # Debug code for proximity with plane2
    for i, r in enumerate(traj):
        dist = abs(np.dot(r - plane2.r0, plane2.n))
        # if dist < 0.05:
        #     print(f"Near plane2 at t={times[i]:.4f}s, r={r}, normal-dist={dist:.4f}")

    # Result in text
    required_order = ["plane1", "plane2", "sphere1"]
    bounce_names = [name for (_, name, _) in bounces]
    failed = False
    last_idx = -1
    for req in required_order:
        if req not in bounce_names:
            print(f"FAILED: Ball did not hit {req}")
            failed = True
            break
        idx = bounce_names.index(req)
        if idx <= last_idx:
            print(f"FAILED: Ball hit {req} out of order")
            failed = True
            break
        last_idx = idx
    if not failed:
        print("SUCCESS: Required planes were hit in order.")

    # Print bounces
    last_printed = None
    for tb, name, pos in bounces:
        if name == last_printed:
            continue
        print(f"Bounce on {name} at t={tb:.4f}s, r={pos}")
        last_printed = name

    # Trajectory and surfaces
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111, projection='3d')
    plot_sphere(ax, sphere1.c, sphere1.Rs, color='green', alpha=0.4)
    ax.plot(traj[:,0], traj[:,1], traj[:,2], lw=2, label='trajectory')
    ax.scatter(start_pos[0], start_pos[1], start_pos[2], color='red', s=60, label='start')

    # Bounces for planes
    if bounces:
        B = np.array([b[2] for b in bounces])
        ax.scatter(B[:,0], B[:,1], B[:,2], s=50, marker='o', label='bounces')

    # Planes
    plot_plane(ax, plane1.r0, plane1.n, size=plane1.size, color='saddlebrown', alpha=0.5)
    plot_plane(ax, plane2.r0, plane2.n, size=plane2.size, color='peru', alpha=0.5)

    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]'); ax.set_zlabel('z [m]')
    ax.set_title('Trickshot trajectory (golf ball vs wood planes)')
    ax.set_xlim(0, 18)
    ax.set_ylim(-8, 8)
    ax.set_zlim(0, z_start + 1.0)
    ax.legend()
    plt.tight_layout()
    plt.show()