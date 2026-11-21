"""Mohr's circle for second moments of areas (Section 10-3.1).

Companion module for the notebooks in this folder. Point names follow
Figs. 10-22 and 10-23 of Riley, Sturges & Morris, *Statics and Mechanics
of Materials*:

    O  origin                         C  centre, OC = (Ix + Iy)/2
    A  (Ix,  Ixy)   -> x-axis         B  (Iy, -Ixy)   -> y-axis
    D  (Imax, 0)    -> u-axis         E  (Imin, 0)    -> v-axis
    F  (Ix', Ix'y') -> x'-axis        G  (Iy', -Ix'y')-> y'-axis

Sign convention (same as the text): Ix' is plotted to the right, Ix'y' is
plotted UPWARD. The axes are rotated by θ counterclockwise; on the circle
the point moves 2θ counterclockwise. 2θp is the angle from CD to CA
(Eq. c: "ACD ... clockwise"), so the u-axis lies θp clockwise from x.

Usage in a notebook::

    from mohr_circle_area import mohr_circle_area_app
    mohr_circle_area_app(Ix=30, Iy=10, Ixy=8)
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Rectangle

# ---------------------------------------------------------------------------
# Equations
# ---------------------------------------------------------------------------


def transform_inertia(Ix, Iy, Ixy, theta_deg):
    """(I_x', I_y', I_x'y') for axes rotated θ CCW (Eqs. 10-15, 10-19, d, e)."""
    t2 = 2 * np.radians(theta_deg)
    avg, half = (Ix + Iy) / 2, (Ix - Iy) / 2
    Ixp = avg + half * np.cos(t2) - Ixy * np.sin(t2)
    Iyp = avg - half * np.cos(t2) + Ixy * np.sin(t2)
    Ixyp = Ixy * np.cos(t2) + half * np.sin(t2)
    return Ixp, Iyp, Ixyp


def mohr_geometry(Ix, Iy, Ixy):
    """Everything that defines the circle of Fig. 10-22.

    Returns a dict with
      OC, R            centre distance and radius
      Imax, Imin       Eqs. (a), (b)  (points D and E)
      two_theta_p      angle from CD to CA in degrees, CCW positive (Fig. 10-22b);
                       Eq. (c) gives its size: tan(ACD) = 2Ixy / (Ix - Iy)
      theta_u          CCW angle from the x-axis to the u-axis (= -two_theta_p / 2)
    """
    OC = (Ix + Iy) / 2
    half = (Ix - Iy) / 2
    R = float(np.hypot(half, Ixy))
    two_tp = float(np.degrees(np.arctan2(Ixy, half))) if R > 1e-12 else 0.0
    theta_u = -two_tp / 2 + 0.0
    if theta_u <= -90 + 1e-9:   # Ixy = 0, Ix < Iy: point u along +y rather than -y
        theta_u += 180.0
    return dict(OC=OC, R=R, Imax=OC + R, Imin=OC - R,
                two_theta_p=two_tp, theta_u=theta_u)


def is_feasible(Ix, Iy, Ixy):
    """A real area needs Imin >= 0, i.e. Ixy^2 <= Ix*Iy (Cauchy-Schwarz)."""
    return Ix >= 0 and Iy >= 0 and Ixy**2 <= Ix * Iy + 1e-12


# ---------------------------------------------------------------------------
# Composite sections (rectangles given as lower-left x0, y0, width b, height h)
# ---------------------------------------------------------------------------


class Section:
    """Area made of rectangles; second moments about centroidal x, y axes."""

    def __init__(self, name, rects, unit="mm"):
        self.name, self.unit = name, unit
        self.rects = [tuple(map(float, r)) for r in rects]
        a = np.array([b * h for _, _, b, h in self.rects])
        cx = np.array([x0 + b / 2 for x0, _, b, _ in self.rects])
        cy = np.array([y0 + h / 2 for _, y0, _, h in self.rects])
        self.A = a.sum()
        self.xc, self.yc = (a * cx).sum() / self.A, (a * cy).sum() / self.A
        dx, dy = cx - self.xc, cy - self.yc
        self.Ix = sum(b * h**3 / 12 for _, _, b, h in self.rects) + (a * dy**2).sum()
        self.Iy = sum(h * b**3 / 12 for _, _, b, h in self.rects) + (a * dx**2).sum()
        self.Ixy = (a * dx * dy).sum()

    @property
    def inertia(self):
        return self.Ix, self.Iy, self.Ixy

    def table(self, scale=1e6):
        """Parallel-axis table as text (part, A, dx, dy, Ix, Iy, Ixy)."""
        lines = [f"{'part':6}{'A':>8}{'dx':>9}{'dy':>9}{'Ix':>10}{'Iy':>10}{'Ixy':>10}"
                 f"   [mm, mm^2, 10^{int(round(np.log10(scale)))} mm^4]"]
        for i, (x0, y0, b, h) in enumerate(self.rects, 1):
            a = b * h
            dx, dy = x0 + b / 2 - self.xc, y0 + h / 2 - self.yc
            lines.append(f"{i:<6}{a:8.0f}{dx:9.2f}{dy:9.2f}"
                         f"{(b * h**3 / 12 + a * dy**2) / scale:10.3f}"
                         f"{(h * b**3 / 12 + a * dx**2) / scale:10.3f}{a * dx * dy / scale:10.3f}")
        lines.append(f"{'total':6}{self.A:8.0f}{'':18}{self.Ix / scale:10.3f}{self.Iy / scale:10.3f}{self.Ixy / scale:10.3f}")
        lines.append(f"centroid: x̄ = {self.xc:.2f} mm, ȳ = {self.yc:.2f} mm")
        return "\n".join(lines)


SECTIONS = {
    # Example Problem 10-13 (Fig. 10-24): unequal-leg angle 80 x 60 x 10 mm
    # A = vertical leg 10 x 80, B = horizontal leg 50 x 10
    "Example 10-13": Section("Example 10-13: angle 80 x 60 x 10 mm", [(0, 0, 10, 80), (10, 0, 50, 10)]),
    # Ix > Iy, Ixy > 0: the same situation as Fig. 10-22
    "Z-section": Section("Z-section (Ixy > 0)", [(80, 180, 100, 20), (80, 20, 20, 160), (0, 0, 100, 20)]),
    # Ixy < 0: A lies below the axis, u-axis is counterclockwise from x
    "L-angle": Section("L-angle 150 x 100 x 20 (Ixy < 0)", [(0, 0, 20, 150), (20, 0, 80, 20)]),
    # Ix < Iy: A lies left of C, D is on the other side
    "L-angle, long leg horizontal": Section("L-angle 100 x 150 x 20 (Ix < Iy)", [(0, 0, 20, 100), (20, 0, 130, 20)]),
    # Axis of symmetry: Ixy = 0, A and B on the horizontal axis
    "T-section": Section("T-section (symmetric about y)", [(0, 130, 120, 20), (50, 0, 20, 130)]),
    # Rectangle b = 60, h = 120: two axes of symmetry, Ixy = 0, A = D and B = E
    "Rectangle": Section("Rectangle 60 x 120 (b x h)", [(0, 0, 60, 120)]),
    # Ix = Iy, Ixy = 0: the circle shrinks to a point
    "Square": Section("Square 100 x 100", [(0, 0, 100, 100)]),
}


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _f(v):
    """Format to 2 decimals without printing -0.00."""
    return f"{(0.0 if abs(v) < 0.005 else v):.2f}"


CIRCLE = "#e0245e"   # the textbook draws the circle in red/pink
C_A, C_B, C_P, C_F, C_G = "black", "black", "tab:green", "tab:blue", "tab:purple"


def _label(ax, P, text, C, R, color="k", size=11, push=0.09):
    """Write `text` next to point P, pushed radially away from the centre C."""
    v = np.array(P) - np.array([C, 0.0])
    n = np.linalg.norm(v)
    d = v / n if n > 1e-12 else np.array([0.7, 0.7])
    q = np.array(P) + push * max(R, 1e-9) * d
    if push < 0:   # inward label: lift it off the diameter line
        below = abs(P[1]) < 1e-6 * max(R, 1e-9)   # point on the I_x' axis: go below
        q = q + np.array([0.0, (-0.07 if below else 0.07) * max(R, 1e-9)])
        d = np.array([-d[0], -1.0 if below else 1.0])
    ax.text(*q, text, color=color, fontsize=size,
            ha="left" if d[0] > 0.2 else ("right" if d[0] < -0.2 else "center"),
            va="bottom" if d[1] > 0.2 else ("top" if d[1] < -0.2 else "center"))


def _arc(ax, center, r, a0, a1, color, text=None, size=10, text_r=1.25):
    a = np.radians(np.linspace(a0, a1, 60))
    ax.plot(center[0] + r * np.cos(a), center[1] + r * np.sin(a), color=color, lw=1.2)
    if abs(a1 - a0) > 1:   # small arrow head at the end of the arc
        ax.annotate("", xy=(center[0] + r * np.cos(a[-1]), center[1] + r * np.sin(a[-1])),
                    xytext=(center[0] + r * np.cos(a[-4]), center[1] + r * np.sin(a[-4])),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))
    if text:
        m = a[len(a) // 2]
        ax.text(center[0] + text_r * r * np.cos(m), center[1] + text_r * r * np.sin(m), text,
                color=color, fontsize=size, ha="center", va="center")


def _dim(ax, x0, x1, y, text, color="dimgray"):
    """Horizontal dimension line from x0 to x1 at height y (like Fig. 10-22a)."""
    if abs(x1 - x0) < 1e-12:
        return
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<|-|>", color=color, lw=1))
    ax.text((x0 + x1) / 2, y, text, color=color, fontsize=9, ha="center", va="top")


def _draw_point_circle(ax, C, unit, dims):
    """Ix = Iy and Ixy = 0: R = 0, so Mohr's circle shrinks to one point."""
    h = 0.35 * max(C, 1e-9)                    # plot height, for nice proportions
    xmax = 1.3 * C
    ax.annotate("", xy=(xmax, 0), xytext=(-0.15 * C, 0), arrowprops=dict(arrowstyle="-|>", color="k", lw=1))
    ax.annotate("", xy=(0, h), xytext=(0, -h), arrowprops=dict(arrowstyle="-|>", color="k", lw=1))
    ax.text(xmax, 0.03 * h, "$I_{x'}$", fontsize=12, ha="right", va="bottom")
    ax.text(0.02 * C, h, "$I_{x'y'}$", fontsize=12, ha="left", va="top")
    ax.text(0, 0, "O ", fontsize=11, ha="right", va="top")
    ax.plot(C, 0, "o", color=CIRCLE, ms=11, zorder=3)
    ax.annotate(f"A = B = C = D = E = F = G\n({C:.2f}, 0)", xy=(C, 0), xytext=(C, 0.45 * h),
                ha="center", fontsize=11, arrowprops=dict(arrowstyle="->", color="gray"))
    if dims:
        _dim(ax, 0, C, -0.25 * h, f"OC = (Ix+Iy)/2 = {C:.2f}")
    ax.text(0.5 * C, -0.55 * h,
            "R = 0: the circle shrinks to a single point.\n"
            f"Every axis has $I_{{x'}}$ = {C:.2f} and $I_{{x'y'}}$ = 0,\n"
            "so every axis through O is a principal axis.",
            ha="center", va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="#fff4f7", ec=CIRCLE))
    ax.set_xlim(-0.2 * C, xmax)
    ax.set_ylim(-1.05 * h, 1.05 * h)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Mohr's circle for second moments  [{unit}]", fontsize=12)


def draw_mohr(ax, Ix, Iy, Ixy, theta=None, steps=4, dims=True, unit="$10^6$ mm$^4$"):
    """Draw Mohr's circle for second moments in the style of Fig. 10-22.

    steps : 1 -> point A;  2 -> + point B;  3 -> + diameter AB and centre C;
            4 -> + circle, radius R, D, E and angle 2θp (complete Fig. 10-22a).
    theta : if given, also draw F, G and the angle 2θ (Fig. 10-22b).
    """
    g = mohr_geometry(Ix, Iy, Ixy)
    C, R = g["OC"], g["R"]
    if steps >= 2 and R <= 1e-9 * max(abs(C), 1.0):
        _draw_point_circle(ax, C, unit, dims)
        return
    Rv = max(R, 0.06 * max(C, 1e-9))           # visual size when R is ~0

    # axes through O
    xmin = min(0.0, C - R) - 0.18 * (C + R + Rv)
    xmax = C + R + 0.35 * Rv + 0.05 * (C + R)
    ylim = 1.45 * Rv + (0.35 * Rv if dims else 0)
    ax.annotate("", xy=(xmax, 0), xytext=(xmin, 0), arrowprops=dict(arrowstyle="-|>", color="k", lw=1))
    ax.annotate("", xy=(0, 1.3 * Rv), xytext=(0, -ylim), arrowprops=dict(arrowstyle="-|>", color="k", lw=1))
    ax.text(xmax, 0.04 * Rv, "$I_{x'}$", fontsize=12, ha="right", va="bottom")
    ax.text(0.02 * (xmax - xmin), 1.3 * Rv, "$I_{x'y'}$", fontsize=12, ha="left", va="top")
    ax.text(0, 0, "O ", fontsize=11, ha="right", va="top")

    A, B = (Ix, Ixy), (Iy, -Ixy)
    ax.plot(*A, "o", color=C_A, ms=6)
    on_axis = abs(Ixy) <= 1e-9 * max(Rv, 1.0)   # Ixy = 0: A and B are the ends D, E
    if on_axis:   # label above the axis so it does not sit on O or the D/E labels
        ax.text(Ix, 0.08 * Rv, f"A ({_f(Ix)}, 0)", fontsize=11, ha="center", va="bottom")
    else:
        _label(ax, A, f"A ({_f(Ix)}, {_f(Ixy)})", C, Rv)
    if dims:   # vertical I_xy dimension at A
        ax.plot([Ix, Ix], [0, Ixy], ":", color="dimgray", lw=1)
    if steps >= 2:
        ax.plot(*B, "o", color=C_B, ms=6)
        if on_axis:
            ax.text(Iy, 0.08 * Rv, f"B ({_f(Iy)}, 0)", fontsize=11, ha="center", va="bottom")
        else:
            _label(ax, B, f"B ({_f(Iy)}, {_f(-Ixy)})", C, Rv)
    if steps >= 3:
        ax.plot([A[0], B[0]], [A[1], B[1]], color="k", lw=1.2)
        ax.plot(C, 0, "k|", ms=12, mew=1.5)
        ax.text(C, -0.06 * Rv, "C", fontsize=11, ha="left", va="top")
    if steps >= 4:
        phi = np.linspace(0, 2 * np.pi, 300)
        ax.plot(C + R * np.cos(phi), R * np.sin(phi), color=CIRCLE, lw=2)
        ax.plot([C + R, C - R], [0, 0], "s", color=C_P, ms=6)
        ax.text(C + R + 0.04 * Rv, -0.05 * Rv, f"D\n$I_{{max}}$={C + R:.2f}", color=C_P,
                fontsize=10, ha="left", va="top")
        ax.text(C - R - 0.04 * Rv, -0.05 * Rv, f"E\n$I_{{min}}$={C - R:.2f}", color=C_P,
                fontsize=10, ha="right", va="top")
        # radius label on CA
        if R > 1e-12:   # label R beside the radius CA, offset perpendicular to it
            u = np.array([Ix - C, Ixy]) / R
            p = np.array([C, 0]) + 0.62 * R * u + 0.09 * R * np.array([-u[1], u[0]])
            ax.text(*p, f"R={R:.2f}", fontsize=10, ha="center", va="center")
        if R > 1e-12 and 0.05 < abs(g["two_theta_p"]) < 179.95:   # no arc when A is on the axis
            tp2 = g["two_theta_p"]
            _arc(ax, (C, 0), 0.26 * R, 0, tp2, C_P, f"2θp={abs(tp2):.1f}°", text_r=1.8)
        if dims:
            yd = -1.18 * Rv
            _dim(ax, 0, C, yd, f"(Ix+Iy)/2 = {C:.2f}")
            _dim(ax, C, Ix, yd - 0.22 * Rv, f"(Ix−Iy)/2 = {(Ix - Iy) / 2:.2f}")

    if theta is not None and steps >= 4:
        Ixp, Iyp, Ixyp = transform_inertia(Ix, Iy, Ixy, theta)
        Ixyp = Ixyp + 0.0
        F, G = (Ixp, Ixyp), (Iyp, -Ixyp)
        ax.plot([F[0], G[0]], [F[1], G[1]], color=C_F, lw=1.6)
        ax.plot(*F, "o", color=C_F, ms=7)
        ax.plot(*G, "o", color=C_G, ms=7)
        named = {"A": A, "B": B, "D": (C + R, 0.0), "E": (C - R, 0.0)}
        for P, nm, col in [(F, "F", C_F), (G, "G", C_G)]:
            same = [k for k, Q in named.items() if np.hypot(P[0] - Q[0], P[1] - Q[1]) < 0.02 * Rv]
            if same:   # coincides with an existing point: short label on the inside
                _label(ax, P, f"{nm} = {same[0]}", C, Rv, color=col, push=-0.2)
            else:
                _label(ax, P, f"{nm} ({_f(P[0])}, {_f(P[1])})", C, Rv, color=col)
        ax.plot([Ixp, Ixp], [0, Ixyp], ":", color=C_F, lw=1)
        ax.plot([C, Ixp], [Ixyp, Ixyp], ":", color=C_F, lw=1)
        if R > 1e-12 and abs(theta) > 0.05:
            tp2 = g["two_theta_p"]
            _arc(ax, (C, 0), 0.5 * R, tp2, tp2 + 2 * theta, C_F, f"2θ={2 * theta:.1f}°", text_r=1.22)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-ylim, 1.45 * Rv)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Mohr's circle for second moments  [{unit}]", fontsize=12)


def _blob():
    p = np.linspace(0, 2 * np.pi, 200)
    r = 1.0 + 0.12 * np.cos(3 * p + 0.8) + 0.08 * np.sin(2 * p) + 0.05 * np.cos(5 * p)
    return np.c_[1.15 * r * np.cos(p), r * np.sin(p)]


def draw_axes(ax, Ix, Iy, Ixy, theta=0.0, section=None):
    """The area with its x, y, x', y', u, v axes, labelled with the matching
    point on Mohr's circle (Fig. 10-23)."""
    g = mohr_geometry(Ix, Iy, Ixy)
    if section is None:
        ax.add_patch(Polygon(_blob(), closed=True, fc="#9fb6c8", ec="none", alpha=0.8))
        L = 0.95
        ax.set_title("Area and axes (Fig. 10-23)", fontsize=12)
    else:
        for x0, y0, b, h in section.rects:
            ax.add_patch(Rectangle((x0 - section.xc, y0 - section.yc), b, h,
                                   fc="#9fb6c8", ec="#38607e", alpha=0.8))
        ext = max(max(abs(x0 - section.xc), abs(x0 + b - section.xc),
                      abs(y0 - section.yc), abs(y0 + h - section.yc))
                  for x0, y0, b, h in section.rects)
        L = 0.85 * ext + 15
        ax.set_title(section.name, fontsize=12)

    axes_list = [(0.0, "x", "A", "k", "-", 1.0), (90.0, "y", "B", "k", "-", 1.0)]
    if g["R"] > 1e-12:
        axes_list += [(g["theta_u"], "u", "D", C_P, "--", 1.4), (g["theta_u"] + 90, "v", "E", C_P, "--", 1.4)]
    if theta is not None:
        axes_list += [(theta, "x′", "F", C_F, "-", 1.4), (theta + 90, "y′", "G", C_G, "-", 1.4)]

    groups = []   # axes pointing the same way share one label, e.g. "y = x′ (B = F)"
    for ang, name, pt, color, ls, lw in axes_list:
        d = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        ax.plot([-0.25 * L * d[0], L * d[0]], [-0.25 * L * d[1], L * d[1]], ls=ls, color=color, lw=lw)
        ax.annotate("", xy=L * d, xytext=0.9 * L * d, arrowprops=dict(arrowstyle="-|>", color=color))
        for grp in groups:
            if abs((ang - grp[0] + 180) % 360 - 180) < 0.5:
                grp[1].append(name); grp[2].append(pt); grp[3] = color
                break
        else:
            groups.append([ang, [name], [pt], color])
    for ang, names, pts, color in groups:
        d = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        ax.text(*(1.12 * L * d), f"{' = '.join(names)} ({' = '.join(pts)})", color=color, fontsize=11,
                ha="center", va="center", fontweight="bold")

    if g["R"] > 1e-12 and abs(g["theta_u"]) > 0.05:
        _arc(ax, (0, 0), 0.42 * L, 0, g["theta_u"], C_P, f"θp={abs(g['theta_u']):.1f}°", text_r=1.3)
    if theta is not None and abs(theta) > 0.05:
        _arc(ax, (0, 0), 0.25 * L, 0, theta, C_F, f"θ={theta:.1f}°", text_r=1.45)
    ax.plot(0, 0, "ko", ms=4)
    ax.text(-0.04 * L, -0.04 * L, "O", ha="right", va="top", fontsize=11)
    ax.set_xlim(-1.35 * L, 1.35 * L)
    ax.set_ylim(-1.35 * L, 1.35 * L)
    ax.set_aspect("equal")
    ax.axis("off")


def report(Ix, Iy, Ixy, theta=None):
    """Text summary: the textbook equations evaluated."""
    g = mohr_geometry(Ix, Iy, Ixy)
    if g["R"] < 1e-9:
        direction = "R = 0: the circle is a single point, every axis is a principal axis"
    elif abs(g["theta_u"]) < 1e-9:
        direction = "A = D: the x-axis is already the u-axis (Ixy = 0)"
    elif abs(abs(g["theta_u"]) - 90) < 1e-9:
        direction = "A = E, B = D: the y-axis is the u-axis and the x-axis is the v-axis (Ixy = 0, Ix < Iy)"
    else:
        direction = (f"angle from CD to CA = 2θp = {abs(g['two_theta_p']):.2f}°  ->  u-axis is "
                     f"{abs(g['theta_u']):.2f}° {'clockwise' if g['theta_u'] < 0 else 'counterclockwise'} from x")
    s = (f"OC = (Ix+Iy)/2 = {g['OC']:.3f}     R = {g['R']:.3f}\n"
         f"Imax = OC + R = {g['Imax']:.3f}   (D, u-axis)      Eq. (a)\n"
         f"Imin = OC − R = {g['Imin']:.3f}   (E, v-axis)      Eq. (b)\n"
         f"{direction}     Eq. (c)")
    if theta is not None:
        a, b, c = transform_inertia(Ix, Iy, Ixy, theta)
        s += (f"\nθ = {theta:.1f}°:  Ix' = {a:.3f},  Iy' = {b:.3f},  Ix'y' = {c + 0.0:.3f}   Eqs. (d), (e)")
    if not is_feasible(Ix, Iy, Ixy):
        s += ("\n\n⚠ Imin < 0: no real area has these values. Second moments are always positive, "
              "which requires Ixy² ≤ Ix·Iy (the whole circle must lie right of O).")
    return s


def plot_mohr_area(Ix, Iy, Ixy, theta=None, section=None, steps=4, dims=True,
                   unit="$10^6$ mm$^4$", figsize=(15, 7)):
    """Area/axes (left) and Mohr's circle (right), side by side."""
    fig, axes = plt.subplots(1, 2, figsize=figsize, gridspec_kw=dict(width_ratios=[1, 1.35]))
    draw_axes(axes[0], Ix, Iy, Ixy, theta, section)
    draw_mohr(axes[1], Ix, Iy, Ixy, theta, steps=steps, dims=dims, unit=unit)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Interactive apps
# ---------------------------------------------------------------------------


def auto_scale(section):
    """10^3 mm^4 for small sections, 10^6 mm^4 for large ones."""
    return 1e3 if max(section.Ix, section.Iy) < 1e7 else 1e6


def unit_label(scale):
    return f"$10^{int(round(np.log10(scale)))}$ mm$^4$"


STEP_NAMES = {1: "1: point A (x-axis)", 2: "2: + point B (y-axis)",
              3: "3: + diameter AB, centre C", 4: "4: + circle, D, E, F, G"}


def mohr_circle_area_app(Ix=808.0, Iy=388.0, Ixy=-323.0, value_range=1000.0, step=1.0,
                         unit="$10^3$ mm$^4$"):
    """Sliders for the construction step, θ, Ix, Iy and Ixy.

    Defaults are Example Problem 10-13 (values in 10^3 mm^4).
    """
    from ipywidgets import Dropdown, FloatSlider, interact

    def _show(steps, theta, Ix, Iy, Ixy):
        plot_mohr_area(Ix, Iy, Ixy, theta if steps >= 4 else None, steps=steps, unit=unit)
        print(report(Ix, Iy, Ixy, theta if steps >= 4 else None))

    interact(_show,
             steps=Dropdown(options=[(v, k) for k, v in STEP_NAMES.items()], value=4,
                            description="Construction:"),
             theta=FloatSlider(value=0, min=-180, max=180, step=0.5, description="θ (deg):"),
             Ix=FloatSlider(value=Ix, min=0, max=value_range, step=step, description="Ix:"),
             Iy=FloatSlider(value=Iy, min=0, max=value_range, step=step, description="Iy:"),
             Ixy=FloatSlider(value=Ixy, min=-value_range, max=value_range, step=step, description="Ixy:"))


def _show_section(s, theta):
    """Parallel-axis table, area + Mohr's circle, and report for one Section."""
    sc = auto_scale(s)
    I = [v / sc for v in s.inertia]
    print(s.table(scale=sc), "\n")
    plot_mohr_area(*I, theta, section=s, unit=unit_label(sc))
    print(report(*I, theta))


def custom_section_app(section):
    """θ slider for one user-defined Section (no dropdown)."""
    from ipywidgets import FloatSlider, interact

    if not isinstance(section, Section):
        raise TypeError("pass a Section, e.g. mc.Section('name', [(x0, y0, b, h), ...])")
    interact(lambda theta: _show_section(section, theta),
             theta=FloatSlider(value=0, min=-180, max=180, step=0.5, description="θ (deg):",
                               continuous_update=False))


def section_app(section="Example 10-13"):
    """Dropdown of sections in SECTIONS and a θ slider (values about the centroid)."""
    from ipywidgets import Dropdown, FloatSlider, interact

    def _show(section, theta):
        _show_section(SECTIONS[section], theta)

    interact(_show,
             section=Dropdown(options=list(SECTIONS), value=section, description="Section:"),
             theta=FloatSlider(value=0, min=-180, max=180, step=0.5, description="θ (deg):",
                               continuous_update=False))


def rectangle_app(b=60.0, h=120.0, size_max=300.0, step=5.0, scale=1e3):
    """Sliders for the width b, height h (mm) and θ of a rectangular cross-section.

    Second moments are about the centroid, in units of `scale` mm^4 (default 10^3 mm^4).
    """
    from ipywidgets import FloatSlider, interact

    def _show(b, h, theta):
        s = Section(f"Rectangle {b:g} x {h:g} mm (b x h)", [(0, 0, b, h)])
        Ix, Iy = b * h**3 / 12 / scale, h * b**3 / 12 / scale
        print(f"Ix  = b h^3 / 12 = {b:g} x {h:g}^3 / 12 = {Ix:,.2f}")
        print(f"Iy  = h b^3 / 12 = {h:g} x {b:g}^3 / 12 = {Iy:,.2f}")
        print("Ixy = 0   (x and y are axes of symmetry)")
        print(f"[units: 10^{int(round(np.log10(scale)))} mm^4]\n")
        plot_mohr_area(Ix, Iy, 0.0, theta, section=s, unit=unit_label(scale))
        print(report(Ix, Iy, 0.0, theta))

    size = dict(min=step, max=size_max, step=step, continuous_update=False)
    interact(_show,
             b=FloatSlider(value=b, description="b (mm):", **size),
             h=FloatSlider(value=h, description="h (mm):", **size),
             theta=FloatSlider(value=0, min=-180, max=180, step=0.5, description="θ (deg):",
                               continuous_update=False))
