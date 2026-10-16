"""Fig. 4: regulation signals (delta_H*, delta_tp*) vs pi1-pi4, all simulations (pi1 panels: main ensemble)."""
import style as S
from scatter_matrix import matrix

S.apply()
matrix(S.load(), ["delta_H*", "delta_tp*"], "fig04_scatter_regulation_vs_pi")
