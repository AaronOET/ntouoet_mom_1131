# Mohr's Circle for Second Moments of Areas

Interactive notebooks for Section 10-3.1 *Mohr's Circle for Second Moments of Areas* (Riley, Sturges & Morris, *Statics and Mechanics of Materials*, Ch. 10, Figs. 10-22 to 10-25, Example Problem 10-13).

All calculation and plotting code is in [mohr_circle_area.py](mohr_circle_area.py). The notebooks import it, in the same way the [Mohr_Circle_v2](../Mohr_Circle_v2/README.md) notebooks import `mohr_circle.py`. For the analytical derivation of the transformation equations, see [Principal_Second_Moments](../Principal_Second_Moments/README.md).

## Notebooks

| # | Notebook | Content | Launch |
| --- | --- | --- | --- |
| 1 | [MohrCircle_SecondMoments_1_Example.ipynb](MohrCircle_SecondMoments_1_Example.ipynb) | Sliders for $I_x$, $I_y$, $I_{xy}$, $\theta$ and a construction-step selector (A → B → diameter and C → circle), starting from Example Problem 10-13 ($I_x = 808$, $I_y = 388$, $I_{xy} = -323$ $\times 10^3$ mm$^4$; $I_{max} = 983$, $I_{min} = 213$, $\theta_p = 28.5^\circ$ counterclockwise) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_for_Second_Moments_of_Areas%2FMohrCircle_SecondMoments_1_Example.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_for_Second_Moments_of_Areas/MohrCircle_SecondMoments_1_Example.ipynb) |
| 2 | [MohrCircle_SecondMoments_2_Sections.ipynb](MohrCircle_SecondMoments_2_Sections.ipynb) | Dropdown of real cross-sections (Example 10-13 angle, Z-section, L-angles, T-section, rectangle, square) and a worked rectangle example with parallel-axis tables and a $\theta$ slider | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_for_Second_Moments_of_Areas%2FMohrCircle_SecondMoments_2_Sections.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_for_Second_Moments_of_Areas/MohrCircle_SecondMoments_2_Sections.ipynb) |
| 3 | [MohrCircle_SecondMoments_3_Custom_Section.ipynb](MohrCircle_SecondMoments_3_Custom_Section.ipynb) | Define your own section from rectangles and explore its Mohr's circle with a $\theta$ slider (custom section only, no dropdown) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_for_Second_Moments_of_Areas%2FMohrCircle_SecondMoments_3_Custom_Section.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_for_Second_Moments_of_Areas/MohrCircle_SecondMoments_3_Custom_Section.ipynb) |
| 4 | [MohrCircle_SecondMoments_4_Rectangle.ipynb](MohrCircle_SecondMoments_4_Rectangle.ipynb) | Rectangle with sliders for width $b$, height $h$ and $\theta$: hand calculation of $I_x = bh^3/12$, $I_y = hb^3/12$, then Mohr's circle. Tall vs. wide vs. square | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_for_Second_Moments_of_Areas%2FMohrCircle_SecondMoments_4_Rectangle.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_for_Second_Moments_of_Areas/MohrCircle_SecondMoments_4_Rectangle.ipynb) |

## Points on the circle (Figs. 10-22, 10-23)

| Point | Coordinates | Axis on the area |
| --- | --- | --- |
| O | origin | |
| A | $(I_x,\ I_{xy})$ | $x$ |
| B | $(I_y,\ -I_{xy})$ | $y$ |
| C | $\left(\frac{I_x+I_y}{2},\ 0\right)$ | centre |
| D | $(I_{max},\ 0)$ | $u$ |
| E | $(I_{min},\ 0)$ | $v$ |
| F | $(I_{x'},\ I_{x'y'})$ | $x'$ |
| G | $(I_{y'},\ -I_{x'y'})$ | $y'$ |

In each figure, the axes drawn on the area are tagged with their point on the circle, e.g. `x′ (F)`.

## Sign convention

The same as the textbook:

- Rectangular second moments $I_{x'}$ are plotted along the horizontal axis. They are always positive, so the whole circle lies right of O.
- Mixed second moments $I_{x'y'}$ are plotted **upward** when positive.
- Rotating the axes by $\theta$ **counterclockwise** on the area moves F by $2\theta$ **counterclockwise** on the circle.
- $2\theta_p = \angle ACD$. When A is above the axis (as in Fig. 10-22), D is reached clockwise from A, so the $u$-axis is $\theta_p$ clockwise from $x$.

$$
OC = \frac{I_x+I_y}{2},\qquad R = \left[\left(\frac{I_x-I_y}{2}\right)^2 + I_{xy}^2\right]^{1/2},\qquad I_{max,\,min} = OC \pm R
$$

$$
I_{x'} = OC + R\cos(2\theta_p + 2\theta),\qquad I_{x'y'} = R\sin(2\theta_p + 2\theta)
$$

## Built-in sections

Rectangles `(x0, y0, width, height)` in mm, properties about the centroid:

| Key | Case |
| --- | --- |
| `"Example 10-13"` | Unequal-leg angle 80 × 60 × 10 mm (Fig. 10-24): $I_{xy} < 0$, $u$ at 28.5° counterclockwise from $x$ |
| `"Z-section"` | $I_x > I_y$, $I_{xy} > 0$ (same as Fig. 10-22) |
| `"L-angle"` | $I_{xy} < 0$: A below the axis, $u$ counterclockwise from $x$ |
| `"L-angle, long leg horizontal"` | $I_x < I_y$: A left of C |
| `"T-section"` | Axis of symmetry: $I_{xy} = 0$, A = D |
| `"Rectangle"` | 60 × 120 mm: two axes of symmetry, $I_{xy} = 0$, A = D and B = E |
| `"Square"` | $I_x = I_y$, $I_{xy} = 0$: circle shrinks to a point |

## Using the module directly

```python
import mohr_circle_area as mc

mc.plot_mohr_area(808, 388, -323, theta=28.5, unit="$10^3$ mm$^4$")   # area + circle
print(mc.report(808, 388, -323, theta=28.5))  # Eqs. (a)-(e) evaluated
mc.mohr_circle_area_app()                     # sliders, Example 10-13 values

sec = mc.Section("My angle", [(0, 0, 15, 120), (15, 0, 65, 15)])
print(sec.table(scale=1e3))                   # parallel-axis table in 10^3 mm^4
mc.custom_section_app(sec)                    # θ slider for this section only
mc.rectangle_app(b=60, h=120)                 # b, h and θ sliders for a rectangle
mc.section_app("Example 10-13")               # dropdown of the built-in sections
```

## Running

- **Online:** click a Binder or Colab badge above. On Colab each notebook downloads `mohr_circle_area.py` from this repository automatically.
- **Locally:** install `numpy`, `matplotlib`, `ipywidgets` and `sympy` then run all cells.
