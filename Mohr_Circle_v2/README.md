# Mohr's Circle

Interactive notebooks for plane-stress transformation. Each notebook shows the rotated stress element (left) next to its Mohr's circle (right), with sliders for the rotation angle $\theta$, $\sigma_x$, $\sigma_y$ and $\tau_{xy}$.

The three notebooks differ only in the sign convention they use. All plotting code lives in [mohr_circle.py](mohr_circle.py); each notebook just calls `mohr_circle_app(...)` with a different convention name.

## Notebooks

| Notebook | Convention | Transformation | Positive $\tau$ plotted | Launch |
| --- | --- | --- | --- | --- |
| [MohrCircle_Hibbler.ipynb](MohrCircle_Hibbler.ipynb) | `"hibbler"` | Hibbler | Downward | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_v2%2FMohrCircle_Hibbler.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_v2/MohrCircle_Hibbler.ipynb) |
| [MohrCircle_Hibbler_Opposite.ipynb](MohrCircle_Hibbler_Opposite.ipynb) | `"hibbler_opposite"` | Hibbler | Upward | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_v2%2FMohrCircle_Hibbler_Opposite.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_v2/MohrCircle_Hibbler_Opposite.ipynb) |
| [MohrCircle_Timoshenko.ipynb](MohrCircle_Timoshenko.ipynb) | `"timoshenko"` | Timoshenko | Upward | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AaronOET/ntouoet_mom_1131/main?labpath=Mohr_Circle_v2%2FMohrCircle_Timoshenko.ipynb) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AaronOET/ntouoet_mom_1131/blob/main/Mohr_Circle_v2/MohrCircle_Timoshenko.ipynb) |

### 1. Hibbler — `MohrCircle_Hibbler.ipynb`

Follows R. C. Hibbler, *Mechanics of Materials*. The stress transformation is

$$
\sigma_{x'} = \sigma_x\cos^2\theta + \sigma_y\sin^2\theta + 2\tau_{xy}\sin\theta\cos\theta
$$

$$
\tau_{x'y'} = -(\sigma_x - \sigma_y)\sin\theta\cos\theta + \tau_{xy}(\cos^2\theta - \sin^2\theta)
$$

and the $\tau$ axis of the Mohr's circle points **downward**, as in the textbook. With this choice, rotating the element by $\theta$ counterclockwise moves the stress point $2\theta$ counterclockwise on the circle. Default input: $\sigma_x = 200$, $\sigma_y = -200$, $\tau_{xy} = -300$.

### 2. Hibbler, opposite axis — `MohrCircle_Hibbler_Opposite.ipynb`

Uses the same Hibbler transformation equations, but plots positive $\tau$ **upward** (the usual math orientation). Compare it with notebook 1 to see that flipping the $\tau$ axis mirrors the circle and reverses the apparent direction of rotation, while the stress values themselves are unchanged. Default input: $\sigma_x = 200$, $\sigma_y = -200$, $\tau_{xy} = 300$.

### 3. Timoshenko — `MohrCircle_Timoshenko.ipynb`

Uses the Timoshenko shear sign convention, which flips the sign of the shear terms:

$$
\sigma_{x'} = \sigma_x\cos^2\theta + \sigma_y\sin^2\theta - 2\tau_{xy}\sin\theta\cos\theta
$$

$$
\tau_{x'y'} = (\sigma_x - \sigma_y)\sin\theta\cos\theta + \tau_{xy}(\cos^2\theta - \sin^2\theta)
$$

Positive $\tau$ is plotted **upward**. Default input: $\sigma_x = 200$, $\sigma_y = -200$, $\tau_{xy} = -300$.

## What the plots show

- **Stress element**: the square rotated by $\theta$, with normal-stress arrows (red / orange) and shear-stress arrows (green) on the right and top faces, labelled with their magnitudes.
- **Mohr's circle**: center $C = (\sigma_x + \sigma_y)/2$ and radius $R = \sqrt{\left(\frac{\sigma_x - \sigma_y}{2}\right)^2 + \tau_{xy}^2}$. The red point is the right face $(\sigma_{x'}, \tau_{x'y'})$, the orange point is the top face $(\sigma_{y'}, -\tau_{x'y'})$, and the dashed diameter joins them. The $2\theta$ angle is printed on the plot.

## Running

- **Online**: click a Binder or Colab badge above. On Colab the notebook downloads `mohr_circle.py` from this repository automatically.
- **Locally**: install `numpy`, `matplotlib` and `ipywidgets`, open any notebook in Jupyter, and run all cells.

You can also use the module directly:

```python
from mohr_circle import mohr_circle_app

mohr_circle_app("timoshenko", sigma_x=100, sigma_y=50, tau_xy=40, stress_range=500)
```

`stress_range` sets the slider limits ($\pm$ value, default 300) and `step` sets the slider step (default 5).
