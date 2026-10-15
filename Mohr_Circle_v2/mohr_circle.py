"""Interactive plane-stress element and Mohr's circle.

Usage in a notebook::

    from mohr_circle import mohr_circle_app
    mohr_circle_app("hibbler")

Supported sign conventions (see ``CONVENTIONS``):

- ``"hibbler"``          : Hibbler transformation, tau axis pointing downward
- ``"hibbler_opposite"`` : Hibbler transformation, tau axis pointing upward
- ``"timoshenko"``       : Timoshenko transformation, tau axis pointing upward
"""

import matplotlib.pyplot as plt
import numpy as np
from ipywidgets import FloatSlider, interact

# shear_sign: sign of the 2*tau*sin*cos term in the sigma_x' equation
# tau_down:   True -> positive tau is plotted downward on the Mohr's circle
CONVENTIONS = {
    "hibbler": {"shear_sign": +1, "tau_down": True},
    "hibbler_opposite": {"shear_sign": +1, "tau_down": False},
    "timoshenko": {"shear_sign": -1, "tau_down": False},
}

_SQUARE = np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [-0.5, -0.5]])
_ARROW = dict(angles="xy", scale_units="xy", scale=1)


def transform_stress(sigma_x, sigma_y, tau_xy, angle, convention="hibbler"):
    """Return (sigma_x', sigma_y', tau_x'y') for an element rotated by `angle` degrees."""
    k = CONVENTIONS[convention]["shear_sign"]
    t = np.radians(angle)
    s, c = np.sin(t), np.cos(t)
    sx = sigma_x * c**2 + sigma_y * s**2 + k * 2 * tau_xy * s * c
    sy = sigma_x * s**2 + sigma_y * c**2 - k * 2 * tau_xy * s * c
    txy = -k * (sigma_x - sigma_y) * s * c + tau_xy * (c**2 - s**2)
    return sx, sy, txy


def mohr_center_radius(sigma_x, sigma_y, tau_xy):
    """Return (center, radius) of the Mohr's circle."""
    return (sigma_x + sigma_y) / 2, np.hypot((sigma_x - sigma_y) / 2, tau_xy)


def _face_arrows(R, midpoint, normal_dir, shear_dir, sigma, tau):
    """Start points and vectors of the normal and shear arrows on one face."""
    n = R @ (normal_dir * (1 if sigma >= 0 else -1))
    v = R @ (shear_dir * (1 if tau >= 0 else -1))
    mid = R @ midpoint
    n_start = mid if sigma >= 0 else mid - n
    s_start = (2 * n_start + n) / 3 - 0.5 * v
    return n_start, n, s_start, v


def _draw_element(ax, angle, sx, sy, txy):
    t = np.radians(angle)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    square = _SQUARE @ R.T
    ax.plot(square[:, 0], square[:, 1], color="tab:blue", lw=2, label="Stress Element")

    faces = [
        ("Right", [0.5, 0], [1, 0], [0, 1], sx, "r", "g"),
        ("Top", [0, 0.5], [0, 1], [1, 0], sy, "tab:orange", "tab:green"),
    ]
    for name, mid, nd, sd, sigma, c_n, c_s in faces:
        n_start, n, s_start, v = _face_arrows(R, np.array(mid), np.array(nd), np.array(sd), sigma, txy)
        ax.quiver(*n_start, *n, color=c_n, label=f"{name} Normal Stress", **_ARROW)
        ax.quiver(*s_start, *v, color=c_s, label=f"{name} Shear Stress", **_ARROW)
        # Label at the arrow head for tension, at the tail for compression
        n_label = n_start + n if sigma >= 0 else n_start
        ax.text(*(n_label + 0.1), f"{abs(sigma):.1f}", color="k", fontsize=12)
        ax.text(*(s_start + v + 0.1), f"{abs(txy):.1f}", color="k", fontsize=12)

    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.axvline(0, color="gray", lw=0.5, ls="--")
    ax.set_title("Stress Element")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper right")


def _draw_mohr(ax, angle, center, radius, sx, sy, txy, tau_down):
    phi = np.linspace(0, 2 * np.pi, 200)
    ax.plot(center + radius * np.cos(phi), radius * np.sin(phi), color="tab:blue")
    ax.scatter([sx], [txy], color="r", label="Right Stress Point")
    ax.scatter([sy], [-txy], color="tab:orange", label="Top Stress Point")
    ax.scatter([center], [0], color="black")
    ax.plot([sx, sy], [txy, -txy], color="gray", ls="--")

    ax.text(center + 2, 2, f"({center:.1f}, 0)", color="black", fontsize=10)
    ax.text(sx + 2, txy + 2, f"({sx:.1f}, {txy:.1f})", color="r", fontsize=10)
    ax.text(sy + 2, -txy + 2, f"({sy:.1f}, {-txy:.1f})", color="tab:orange", fontsize=10)
    ax.text(sx + 2, txy - 16, f"Angle: {2 * angle:.1f}°", color="tab:blue", fontsize=12)

    lim = radius + 10
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.axvline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Normal Stress (σ)")
    ax.set_ylabel("Shear Stress (τ)")
    ax.set_xlim(center - lim, center + lim)
    ax.set_ylim((lim, -lim) if tau_down else (-lim, lim))
    ax.set_title("Mohr's Circle")
    ax.legend(loc="upper right")
    ax.set_aspect("equal", adjustable="box")


def plot_stress_and_mohr_circle(angle, sigma_x, sigma_y, tau_xy, convention="hibbler"):
    """Plot the rotated stress element (left) and its Mohr's circle (right)."""
    sx, sy, txy = transform_stress(sigma_x, sigma_y, tau_xy, angle, convention)
    center, radius = mohr_center_radius(sigma_x, sigma_y, tau_xy)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    _draw_element(axes[0], angle, sx, sy, txy)
    _draw_mohr(axes[1], angle, center, radius, sx, sy, txy,
               CONVENTIONS[convention]["tau_down"])
    plt.tight_layout()
    plt.show()


def mohr_circle_app(convention="hibbler", sigma_x=200, sigma_y=-200, tau_xy=-300,
                    stress_range=300, step=5):
    """Show the interactive stress element / Mohr's circle with sliders."""
    if convention not in CONVENTIONS:
        raise ValueError(f"convention must be one of {list(CONVENTIONS)}")

    def slider(value, label):
        return FloatSlider(value=value, min=-stress_range, max=stress_range,
                           step=step, description=label)

    interact(
        lambda angle, sigma_x, sigma_y, tau_xy: plot_stress_and_mohr_circle(
            angle, sigma_x, sigma_y, tau_xy, convention),
        angle=FloatSlider(value=0, min=-180, max=180, step=0.5, description="Angle:"),
        sigma_x=slider(sigma_x, "σx:"),
        sigma_y=slider(sigma_y, "σy:"),
        tau_xy=slider(tau_xy, "τxy:"),
    )
