Physics Debate Prep: Trick Shot Credibility with Superball, Wooden Planes, and Bowling Ball

1. Core Question
	•	Investigate credibility of trick shots involving:
	•	Ball: Superball
	•	Surfaces: 2 wooden planar surfaces + 1 non-planar surface (bowling ball)
	•	Focus is on error propagation and credibility analysis.
	•	Goal: critique claims (e.g., “this bounce is impossible due to friction limits”, not just probability scores).

⸻

2. Conceptual Approach
	•	Treat the system as a composition of maps:
	•	Free flight → Collision → Free flight → Collision → …
	•	Each stage is a deterministic map with uncertainties in parameters (e.g., restitution, friction, spin).
	•	State vector representation:
\mathbf{x}(t) =
\begin{bmatrix}
\mathbf{r}(t) \\
\mathbf{v}(t) \\
\boldsymbol{\omega}(t)
\end{bmatrix}
	•	\mathbf{r}: position
	•	\mathbf{v}: velocity
	•	\boldsymbol{\omega}: angular velocity

⸻

3. Collision Physics (Impulse Model)
	•	Contact point relative to center: \mathbf{r}_c
	•	Relative velocity:
\mathbf{v}_\text{rel}^- = \mathbf{v}^- + \boldsymbol{\omega}^- \times \mathbf{r}_c
	•	Split into normal & tangential:
v_n^- = \mathbf{n}\cdot\mathbf{v}_\text{rel}^-, \quad
\mathbf{v}t^- = \mathbf{v}\text{rel}^- - v_n^- \mathbf{n}

Impulse equations:
	•	Normal impulse:
J_n = -\frac{(1+e) v_n^-}{\tfrac{1}{m} + \mathbf{n}\cdot(\mathbf{r}_c \times I^{-1}(\mathbf{r}_c \times \mathbf{n}))}
	•	Tangential impulse:
	•	If sticking possible: solve for zero tangential velocity.
	•	Otherwise slip:
\mathbf{J}_t = -\mu J_n \frac{\mathbf{v}_t^-}{\|\mathbf{v}_t^-\|}

⸻

4. Backpropagation via Jacobians
	•	Each map has Jacobian:
J_\mathcal{C} = \frac{\partial \mathbf{x}^+}{\partial \mathbf{x}^-}
	•	Total chain (flight + collisions):
J_\text{total} = J_{F_2} J_{\mathcal{C}2} J{F_1} J_{\mathcal{C}1} J{F_0}
	•	This is mathematically identical to neural network backpropagation.

⸻

5. Error / Covariance Propagation
	•	If pre-impact uncertainty is \Sigma^-:
\Sigma^+ \approx J \Sigma^- J^\top + \Sigma_\text{model}
	•	For the whole sequence:
\Sigma_f = J_\text{total}\Sigma_0 J_\text{total}^\top + \sum_i J_{\text{post-i}}\Sigma_{\text{param},i}J_{\text{post-i}}^\top
	•	This lets you build confidence ellipsoids for possible final states.

⸻

6. Impossibility & Constraint Checks
	•	Friction cone violation: required |\mathbf{J}_t| > \mu J_n → impossible.
	•	Energy increase: kinetic energy grows without input → impossible.
	•	Spin bound: required ω exceeds physically plausible spin.
	•	Geometry violation: implied contact outside surface area.

Debate phrasing example:

“Reconstructing the observed track requires a tangential impulse 3.4× larger than the maximum static friction impulse measured on the wooden plane — therefore the trajectory is physically infeasible under the documented surfaces.”

⸻

7. Practical Implementation Plan
	1.	Calibrate parameters
	•	Mass, radius, inertia.
	•	Coefficient of restitution (drop test).
	•	Friction coefficient μ (inclined plane test).
	•	Plane geometry.
	2.	Record experiment
	•	High-speed video (≥1000 fps if possible).
	•	Use markers for angular velocity.
	3.	Forward model
	•	Implement rigid impulse collision maps.
	4.	Error model
	•	Build covariance for measurement & parameter uncertainty.
	5.	Inverse / feasibility test
	•	Solve for required μ, e, ω → check plausibility.
	6.	Monte Carlo
	•	Simulate 10k–100k random trials.
	•	Estimate match probability.

⸻

8. Worked Example: Sphere-Plane Collision
	•	Inertia: I = \tfrac{2}{5}mR^2.
	•	Effective mass for normal impulse:
m_\text{eff} = \frac{m}{1+\kappa}
(where \kappa is rotational coupling factor).
	•	Normal impulse:
J_n = -(1+e)m_\text{eff}v_n^-

⸻

9. Sensitivity Ranking
	•	Dominant uncertainties:
	•	Impact angle.
	•	Spin.
	•	μ (friction coefficient).
	•	Less important: mass, radius (often cancel out).

⸻

10. Practical High-School Lab Feasibility
	•	Yes, feasible with careful design:
	•	Consumer high-fps cameras (240–960 fps).
	•	Simple drop rigs and inclined planes.
	•	Free tools: Tracker, Python (NumPy/SciPy).

Challenges:
	•	Measuring spin reliably.
	•	Camera calibration error.
	•	Nonlinear stick/slip behavior.

⸻

11. Recommended Workflow
	1.	Measure m, R, e, μ.
	2.	Record high-fps video (track positions + spin).
	3.	Build forward simulator in Python.
	4.	Use Jacobians for sensitivity analysis.
	5.	Run Monte Carlo for final credibility claims.
	6.	Report results as impossibility or implausibility with measured bounds.

⸻

12. Debate-Ready Statements
	•	“Given measured μ ∈ [0.30, 0.40] and e ∈ [0.62, 0.66], reproducing the observed exit velocity would require μ ≥ 1.2 — outside the measured range.”
	•	“Monte Carlo with 100k samples gave 0 matches; 95% upper bound on probability < 3×10⁻⁵.”
	•	“Required spin is >3× measured max; therefore implausible.”

⸻

13. Final Verdict
	•	Your approach (Jacobian backpropagation) is conceptually correct.
	•	Use it for sensitivity analysis.
	•	For credibility conclusions: Monte Carlo simulations + physical impossibility checks are more robust in a high-school lab context.
	•	With careful parameter measurement and conservative assumptions, you can make strong, defensible claims.