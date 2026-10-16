"""Fig. 5: process metrics (delta_leakage*, delta_storage*, t_r*) vs pi1-pi4, all simulations (pi1 panels: main
ensemble)."""
import style as S
from scatter_matrix import matrix

S.apply()
matrix(S.load(), ["delta_leakage*", "delta_storage*", "t_r*"], "fig05_scatter_processes_vs_pi")
