# -*- coding: utf-8 -*-
"""Differenciális térenergia ΔW(L) = W_fő(<L) − W_ref(<L) (I = 1 A), 2026-10-21. V19: https://doi.org/10.5281/zenodo.23497128
Radiális Gauss–Legendre panelek, közös szögpontok (ugyanaz a Fibonacci-gömb mindkét középpontnál, ugyanazon r-en), interpoláció nélkül.
  L ≤ L_J (0.95 m; a referencia-gömbben nincs tekercs):  ΔW = W_fő,be(L) − W_ref,be(L)           (belső forma)
  L_J < L < 2.6 m:                                         ΔW = W_ref,ki(L) − W_fő,ki(L)           (A út: panelek L→3 m + r = 3/t farok)
  L ≥ 2.6 m:                                               ΔW = W_ref,ki(L) − W_fő,ki(L)           (B út: r = L/t, t ∈ (0,1], Gauss–Legendre; nincs levágás)
A külső formában nem szerepel teljes energia, így két függetlenül becsült teljes energia különbsége nem jelenik meg."""
import math
import numpy as np
import fields as F
MU0 = F.MU0
XG8, WG8 = np.polynomial.legendre.leggauss(8)
def fib(n):
    i = np.arange(n) + 0.5; ph = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)], -1)
class Integrator:
    def __init__(self, cs, a, ell, N, pref, nlo=6000, nhi=24000, href=1, nt=32, hi_range=(0.25, 2.4), LJ=0.95):
        self.cs, self.a, self.ell, self.N, self.pref = cs, a, ell, N, np.asarray(pref, float)
        self.Ulo, self.Uhi = fib(nlo), fib(nhi); self.hi = hi_range; self.LJ = LJ
        self.XT, self.WT = np.polynomial.legendre.leggauss(nt); self.T = 0.5 * (self.XT + 1); self.WTt = 0.5 * self.WT
        self.E = np.unique(np.concatenate([np.linspace(0, 0.25, 12 * href + 1), np.linspace(0.25, 0.75, 50 * href + 1),
                                           np.linspace(0.75, 2.6, 92 * href + 1), np.linspace(2.6, 3.0, 8 * href + 1)]))
        self.cache = {}
    def shell(self, ci, r):                       # 4πr²·<B²>/(2µ0) a ci középpontú r sugarú gömbön
        key = (ci, float(r))
        if key not in self.cache:
            c = np.zeros(3) if ci == 0 else self.pref
            U = self.Uhi if self.hi[0] <= r <= self.hi[1] else self.Ulo
            B = F.coils_B(c + r * U, self.cs, self.a, self.ell, self.N * 1.0, [1] * len(self.cs))
            self.cache[key] = 4 * np.pi * r * r * float(np.mean(np.sum(B * B, -1))) / (2 * MU0)
        return self.cache[key]
    def gl(self, ci, p, q):
        if q <= p: return 0.0
        rr = 0.5 * (q - p) * XG8 + 0.5 * (q + p)
        return 0.5 * (q - p) * sum(w * self.shell(ci, r) for r, w in zip(rr, WG8))
    def tail(self, ci, R0):                       # ∫_R0^∞, r = R0/t
        return float(sum(w * self.shell(ci, R0 / t) * R0 / t ** 2 for t, w in zip(self.T, self.WTt)))
    def w_in(self, ci, L):                        # ∫_0^L (L ≤ 3 m)
        k = np.searchsorted(self.E, L, side="right") - 1
        return sum(self.gl(ci, self.E[j], self.E[j + 1]) for j in range(k)) + self.gl(ci, self.E[k], L)
    def w_out_A(self, ci, L, Rt=3.0):             # ∫_L^Rt panelek + farok
        if L >= Rt:                               # panelek L → L + 0.5 m, majd farok
            e = np.linspace(L, L + 0.5, 11)
            return sum(self.gl(ci, p, q) for p, q in zip(e[:-1], e[1:])) + self.tail(ci, L + 0.5)
        k = np.searchsorted(self.E, L, side="right")
        return self.gl(ci, L, self.E[k]) + sum(self.gl(ci, self.E[j], self.E[j + 1]) for j in range(k, len(self.E) - 1)) + self.tail(ci, Rt)
    def w_out_B(self, ci, L): return self.tail(ci, L)
    def w_tot(self, ci): return self.w_in(ci, 3.0) + self.tail(ci, 3.0)
    def dW(self, L):
        if L <= self.LJ: return self.w_in(0, L) - self.w_in(1, L)
        if L < 2.6: return self.w_out_A(1, L) - self.w_out_A(0, L)
        return self.w_out_B(1, L) - self.w_out_B(0, L)
    def dW_A(self, L, Rt=3.0): return self.w_out_A(1, L, Rt) - self.w_out_A(0, L, Rt)
    def dW_B(self, L): return self.w_out_B(1, L) - self.w_out_B(0, L)
    def dW_inner(self, L): return self.w_in(0, L) - self.w_in(1, L)
