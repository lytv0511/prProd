# bouncesim.py
import numpy as np
from forces import forces_and_torque, rk4_step, collide

class Simulator:
    def __init__(self, ball, aero, surfaces, target_plane=None,
                 e_n=0.82, mu=0.35, dt=1e-4, t_max=5.0, max_bounces=50):
        self.ball = ball
        self.aero = aero
        self.surfaces = surfaces
        self.target_plane = target_plane
        self.e_n = float(e_n)
        self.mu = float(mu)
        self.dt = float(dt)
        self.t_max = float(t_max)
        self.max_bounces = int(max_bounces)

    def integrate(self, r0, v0, w0):
        r = np.array(r0, dtype=float)
        v = np.array(v0, dtype=float)
        w = np.array(w0, dtype=float)
        ball = self.ball
        aero = self.aero

        times = [0.0]
        traj = [r.copy()]
        seg_ids = [0]
        bounces = []
        target_hits = []

        t = 0.0
        seg = 0
        bcount = 0

        def deriv(state):
            r_, v_, w_ = state
            F, tau = forces_and_torque(r_, v_, w_, ball, aero)
            a = F / ball.m
            alpha = tau / ball.I
            return (v_, a, alpha)

        def crossed_target(r_a, r_b):
            """Return (t_hit, r_hit) if target_plane is crossed between r_a and r_b else None.
               Use robust interpolation similar to earliest_collision.
            """
            if self.target_plane is None:
                return None
            n = self.target_plane.n
            a = np.dot(r_a - self.target_plane.r0, n)
            b = np.dot(r_b - self.target_plane.r0, n)
            # require a sign change or touching
            if not (a > 0 and b <= 0 or a < 0 and b >= 0 or abs(a) <= 1e-12 or abs(b) <= 1e-12):
                return None
            denom = (a - b)
            if abs(denom) < 1e-12:
                s = 0.5  # fallback to midpoint
            else:
                s = a / denom
            s = np.clip(s, 0.0, 1.0)
            r_hit = r_a + s * (r_b - r_a)
            return t + s * self.dt, r_hit

        def earliest_collision(r_pre, r_post):
            """Return (alpha, surf) for earliest collision in [0,1] between r_pre->r_post.
               Robust to tiny denominators and to surfaces returning inf from gap().
            """
            earliest = None
            surf_hit = None
            # displacement for interpolation checks
            dr = r_post - r_pre
            for s in self.surfaces:
                try:
                    g0 = s.gap(r_pre, ball.R)
                    g1 = s.gap(r_post, ball.R)
                except Exception:
                    # if gap computation fails for any reason, skip surface
                    continue

                # skip surfaces where point projects outside finite plane (gap returns inf)
                if not np.isfinite(g0) or not np.isfinite(g1):
                    continue

                # We only consider a crossing from positive (outside) to <=0 (contact/penetration).
                # This avoids detecting trivial grazing or already-interpenetrating states.
                if not (g0 > 0 and g1 <= 0):
                    continue

                denom = (g0 - g1)
                if abs(denom) < 1e-12:
                    # both nearly equal -> fallback to midpoint for alpha
                    a = 0.5
                else:
                    a = g0 / denom

                # clip to [0,1]; if after clipping the alpha is on boundary, still allow but sanity-check
                if not np.isfinite(a):
                    continue
                a = float(a)
                if a < -1e-6 or a > 1.0 + 1e-6:
                    # outside step (shouldn't happen) -> skip
                    continue
                a = np.clip(a, 0.0, 1.0)

                # compute candidate collision point and sanity-check it lies close to surface
                r_cand = r_pre + a * dr
                # final robust gap check at r_cand; allow small negative penetration tolerance
                g_cand = s.gap(r_cand, ball.R)
                if not np.isfinite(g_cand):
                    continue
                # Accept candidate if it's touching or slightly penetrating (tolerance)
                if g_cand <= 1e-6:
                    if earliest is None or a < earliest:
                        earliest = a
                        surf_hit = s
            return earliest, surf_hit

        # main loop
        while t < self.t_max and bcount < self.max_bounces:
            # RK4 prediction for full dt
            r_next, v_next, w_next = rk4_step((r, v, w), self.dt, deriv)

            # target crossing check
            thit = crossed_target(r, r_next)
            if thit is not None:
                target_hits.append(thit)

            # collision detection
            alpha, surf = earliest_collision(r, r_next)
            if surf is None:
                # accept full step
                t += self.dt
                r, v, w = r_next, v_next, w_next
                times.append(t); traj.append(r.copy()); seg_ids.append(seg)
                continue

            # Interpolate state at collision
            # alpha is guaranteed finite and in [0,1]
            r_c = r + alpha * (r_next - r)
            v_c = v + alpha * (v_next - v)
            w_c = w + alpha * (w_next - w)

            # Sanity: ensure no NaNs
            if np.isnan(r_c).any() or np.isnan(v_c).any() or np.isnan(w_c).any():
                # Skip this collision if interpolation produced NaNs (shouldn't happen)
                # Accept full step to try to recover numerically
                t += self.dt
                r, v, w = r_next, v_next, w_next
                times.append(t); traj.append(r.copy()); seg_ids.append(seg)
                continue

            # perform collision impulse using surface normal at contact
            n = surf.normal(r_c)
            v_plus, w_plus = collide(v_c, w_c, n, ball, e_n=self.e_n, mu=self.mu)

            # log bounce
            tb = t + alpha * self.dt
            bcount += 1
            bounces.append((tb, surf.name, r_c.copy()))

            # advance remainder of the step after collision
            rem_dt = (1.0 - alpha) * self.dt
            if rem_dt > 0:
                # integrate remainder with post-collision state
                r_post, v_post, w_post = rk4_step((r_c, v_plus, w_plus), rem_dt, deriv)
                t += self.dt
                r, v, w = r_post, v_post, w_post
            else:
                # collision exactly at end of step
                t += alpha * self.dt
                r, v, w = r_c, v_plus, w_plus

            seg += 1
            times.append(t); traj.append(r.copy()); seg_ids.append(seg)

        return np.array(times), np.array(traj), np.array(seg_ids), bounces, target_hits