# 15. lépés: Elektromágneses nyommentesség ($T^\mu_{\;\mu,\text{EM}} = 0$) vs. diszform mágneses forráscsatolás ($\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$), a V20 polaritás-páros kereszttag-kiejtés és a módszertani határ ($\Phi_t \equiv 0$)
# Step 15: Electromagnetic Trace-Freeness ($T^\mu_{\;\mu,\text{EM}} = 0$) vs. Disformal Magnetic Source Coupling ($\rho_B = |\mathbf{B}|^2/(2\mu_0 c^2)$), V20 Polarity-Even Cross-Term Cancellation, and the Methodological Boundary ($\Phi_t \equiv 0$)

**Szerző / Author:** Szabolcsi Dávid (ORCID: `0009-0005-6300-074X`)  
**Dátum / Date:** 2026-10-26  
**Mappa és fájl / File Path:** `scalar-tensor-yukawa-fit/step15_electromagnetic_trace_anomaly_disformal_magnetic_source_v20.md`  
**Előzmények / Previous Steps:** `step1_dilaton_form.md` – `step14_gravitational_lensing_null_geodesics_cosmological_background.md`  
**Licenc / License:** CC BY-NC 4.0  
**Hivatkozott helyi források / Local Archival Sources:**  
- V22 Book Edition: `/workspace/artifacts/tal_booklet_2026-10-24/hu_p1.tex` (2. fejezet 3–4. pont, 208–221. sor) és `hu_p3.tex` (10. fejezet 2. pont, 75–85. sor; DOI: `10.5281/zenodo.23604423`)  
- V23 (`50aa2db0` feladat): `/workspace/artifacts/tal_covariant_closure_50aa2db0_2026-10-24/tal_covariant_closure_50aa2db0_HU_2026-10-24.tex` (1. és 3. szakasz, 33–35., 58–66. és 161–166. sor; DOI: `10.5281/zenodo.23639084`)  
- `README.md` (`scalar-tensor-yukawa-fit/README.md`, 3. szakasz 5. pont, 60–61. sor)

---

## 1. Szó szerinti idézetek a V22 (`hu_p1.tex`, `hu_p3.tex`) és V23 kéziratokból

A 8. lépésben (`step8_microscopic_clock_coupling_siriusb_lab_v20.md`) levezettük, hogyan módosítja egy adott $\phi_L$ skalárprofil az atomórák és kvarcoszcillátorok frekvenciáját (H0 vs. H1), és kiértékeltük a V20-as 8 tekercses szolenoid számértékét ($\Delta y_{\text{peak}} = -1{,}2017108\times 10^{-48}$ $L = 0{,}630957\;\text{m}$-nél). A skalár–tenzor illesztésben azonban van egy alapvető elméleti kérdés, amelyet a `hu_p1.tex` (208–221. sor), a `hu_p3.tex` (75–85. sor) és a `README.md` (61. sor) jelöl ki: **hogyan válhat a tiszta mágneses energiasűrűség ($\rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$) a $\phi_L$ skalármező forrásává egy kovariáns hatáselvben, amikor a 4 dimenziós Maxwell-féle energia-impulzus tenzor nyoma azonosan nulla ($T^\mu_{\;\mu,\text{EM}} = g^{\mu\nu}T_{\mu\nu}^{\text{EM}} = 0$)?**

1. **A módszertani határ az auditált fizikai mag ($\rho = T_{00}/c^2$, $\Phi_t \equiv 0$) és a V13–V16 perspektivikus réteg között (`hu_p1.tex`, 2. fejezet 3–4. pont, 209–221. sor):**
   $$\begin{equation}
     \rho_{\text{eff}} \;=\; \rho \;+\; \Xi\,\varepsilon_0\,\Phi_t \;+\; \beta\,\Phi_B \;+\; \gamma\,\kappa, \qquad
     T_{\mu\nu}^{\text{total}} \;=\; T_{\mu\nu}^{\text{matter}} \;+\; T_{\mu\nu}^{\text{flux}} \;+\; T_{\mu\nu}^{\text{thought}}, \tag{V22 Eq. 2}
   \end{equation}$$
   > *„Amikor ezt a könyvet az akadémia elé tárom, egyértelmű határt húzok a két réteg közé:*  
   > - ***Auditált fizikai mag (V19–V21):** A SPARC galaktikus forgási elemzésben, a H.E.S.S. időzítési korlátokban, a neutroncsillag TOV- és pulzációs sorozatokban, valamint a V20 gránitkamrás és 8 tekercses kvarcoszcillátoros laboratóriumi protokollban szereplő valamennyi számszerű eredmény kizárólag a standard barionos és elektromágneses $\rho = T_{00}/c^2$ tömeg-energia sűrűségre épül, $\Phi_t \equiv 0$ mellett.*  
   > - ***Perspektivikus bio-koherencia réteg (V13–V16):** A tudati koherencia tagok ($\Phi_t$, $T_{\mu\nu}^{\text{thought}}$, EEG/HRV bio-koherencia) nyitott, még nem auditált elméleti perspektívaként szerepelnek, mint a Time Attraction Law szellemi fejlődéstörténetének részei.”*

2. **A V20 nyolctekercses mágneses energiasűrűségi forrás és a háttér-kereszttag kiejtése (`hu_p3.tex`, 10. fejezet 2. pont, 76. és 84. sor):**
   > - *„Annak vizsgálatára, hogy a tiszta $u_B(\mathbf{x}) = \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0}$ mágneses energiasűrűség a $\rho_B(\mathbf{x}) = u_B(\mathbf{x})/c^2$ egyenleten keresztül lassítja-e az időáramlást, egy szimmetrikus nyolctekercses elrendezést terveztem...”*
   > - *„**A háttér-kereszttag kiejtése polaritásfordítással:** Mivel $|\mathbf{B}_{\text{tot}}|^2 = |\mathbf{B}_{\text{bg}}|^2 + 2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}} + |\mathbf{B}_{\text{coil}}|^2$, az áramirány megfordítása ($+I \leftrightarrow -I$) megfordítja a $\mathbf{B}_{\text{coil}}$ teret, így a polaritás-páros $\Delta y_{\text{even}}$ átlagból mind a lineáris Zeeman-kvarcérzékenység, mind a $2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ kereszttag kiesik.”*

---

## 2. Miért nem forrásozza a mágneses tér a tisztán konform Brans–Dicke / dilaton skalárt, és miért forrásozza a V23 diszform csatolást?

### 2.1. A tisztán konform csatolás elektromágneses vaksága ($D=4$-ben $T^\mu_{\;\mu,\text{EM}} = 0$)
Tekintsük a szabad elektromágneses mező ($F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$) standard Maxwell-hatását egy általános $\tilde{g}_{\mu\nu}$ Jordan-metrikán:
$$S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}, A_\mu\bigr] \;=\; -\frac{1}{4\mu_0}\int d^4x\,\sqrt{-\tilde{g}}\,\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}.$$
Ha a metrika tisztán **konform** módon függ a $\phi_L$ skalármezőtől ($\tilde{g}_{\mu\nu}^{\text{conf}} = \Omega^2(\phi_L)\,g_{\mu\nu} = f_L^2\,g_{\mu\nu}$, így $\tilde{g}_{\text{conf}}^{\mu\alpha} = f_L^{-2}g^{\mu\alpha}$ és $D=4$ dimenzióban $\sqrt{-\tilde{g}^{\text{conf}}} = f_L^4\sqrt{-g}$), akkor az $f_L$ tényezők pontosan kiejtik egymást:
$$\sqrt{-\tilde{g}^{\text{conf}}}\,\tilde{g}_{\text{conf}}^{\mu\alpha}\tilde{g}_{\text{conf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} \;=\; \bigl(f_L^4\sqrt{-g}\bigr)\bigl(f_L^{-2}g^{\mu\alpha}\bigr)\bigl(f_L^{-2}g^{\nu\beta}\bigr)F_{\mu\nu}F_{\alpha\beta} \;=\; \sqrt{-g}\,g^{\mu\alpha}g^{\nu\beta}F_{\mu\nu}F_{\alpha\beta}.$$
Ennek megfelelően a Maxwell-hatás $\phi_L$ szerinti variációja tisztán konform csatolásban azonosan nulla (ami ekvivalens azzal, hogy a Maxwell-féle energia-impulzus tenzor nyoma $D=4$-ben eltűnik, $T^\mu_{\;\mu,\text{EM}} = g^{\mu\nu}T_{\mu\nu}^{\text{EM}} = -u_E - u_B + 3\cdot\frac{1}{3}(u_E + u_B) = 0$):
$$\boxed{\;\frac{\delta S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}^{\text{conf}}, A_\mu\bigr]}{\delta\phi_L(x)} \;=\; \frac{\partial\tilde{g}_{\mu\nu}^{\text{conf}}}{\partial\phi_L}\,\frac{\delta S_{\text{EM}}}{\delta\tilde{g}_{\mu\nu}^{\text{conf}}} \;=\; \sqrt{-g}\,f_L^2\,\frac{\partial\ln f_L}{\partial\phi_L}\,g^{\mu\nu}T_{\mu\nu}^{\text{EM}} \;=\; 0.\;}$$
**Fizikai következmény:** Egy tisztán konform Jordan–Brans–Dicke skalármezőt a sztatikus $\mathbf{B}$ mágneses tér fa-szinten (*tree level*) egyáltalán nem forrásoz!

### 2.2. A V23-as diszform órametrika ($\tilde{g}_{\mu\nu}^{\text{disf}} = g_{\mu\nu} + \frac{1-f_L^2}{c^2}u_\mu u_\nu$) egzakt mágneses forrástagja
Most helyettesítsük be a V23-as **időbeli diszform metrikát** (`tal_covariant_closure_50aa2db0_HU_2026-10-24.tex`, 165. sor):
$$\tilde{g}_{\mu\nu}^{\text{disf}} \;=\; g_{\mu\nu} \;+\; \frac{1 - f_L^2(\phi_L)}{c^2}\,u_\mu u_\nu, \qquad \tilde{g}_{\text{disf}}^{\mu\nu} \;=\; g^{\mu\nu} \;-\; \frac{1 - f_L^{-2}(\phi_L)}{c^2}\,u^\mu u^\nu, \qquad \sqrt{-\tilde{g}^{\text{disf}}} \;=\; f_L(\phi_L)\,\sqrt{-g}.$$
Bontsuk fel a $F_{\mu\nu}$ térerősségtenzort a laboratórium $u^\mu$ négyessebességével együttmozgó elektromos ($E_\mu = F_{\mu\nu}u^\nu/c$) és mágneses ($B_\mu = \frac{1}{2}\epsilon_{\mu\nu\alpha\beta}u^\nu F^{\alpha\beta}/c$) térvektorokra:
$$\frac{1}{4}\tilde{g}_{\text{disf}}^{\mu\alpha}\tilde{g}_{\text{disf}}^{\nu\beta}F_{\mu\nu}F_{\alpha\beta} \;=\; \frac{1}{2}\left[\tilde{g}_{\text{disf}}^{ij}\tilde{g}_{\text{disf}}^{kl}F_{ik}F_{jl} \;+\; 2\,\tilde{g}_{\text{disf}}^{00}\tilde{g}_{\text{disf}}^{ij}F_{0i}F_{0j}\right] \;=\; \frac{|\mathbf{B}|^2}{2\mu_0}\cdot\mu_0 \;-\; f_L^{-2}(\phi_L)\,\frac{\varepsilon_0|\mathbf{E}|^2}{2}\cdot\mu_0.$$
Behelyettesítve ezt és a $\sqrt{-\tilde{g}^{\text{disf}}} = f_L(\phi_L)\sqrt{-g}$ térfogatelemet a Maxwell-hatásba:
$$\boxed{\;S_{\text{EM}}\bigl[\tilde{g}_{\mu\nu}^{\text{disf}}, A_\mu\bigr] \;=\; \int d^4x\,\sqrt{-g}\left[-\,f_L(\phi_L)\,\frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0} \;+\; f_L(\phi_L)^{-1}\,\frac{\varepsilon_0|\mathbf{E}(\mathbf{x})|^2}{2}\right].\;}$$
Sztatikus mágneses térben ($\mathbf{E} = \mathbf{0}$, V20 8 tekercses szolenoid) a hatássűrűség $\mathcal{L}_{\text{EM}}^{\text{disf}}(\mathbf{E}=\mathbf{0}) = -f_L(\phi_L)\,u_B(\mathbf{x})$ pontosan ugyanazt a $-f_L(\phi_L)\varepsilon$ alakot veszi fel, mint a nyugvó barionos folyadék (`step4_covariant_conservation_hydrostatic.md`: $\mathcal{L}_{\text{matter}} = -f_L(\phi_L)\varepsilon_{\text{baryon}}$)! Variálva a V23 II. (Helmholtz–Yukawa) vagy I. (bilaterális top-hat) hatását $\phi_L$ szerint ($f_L = e^{-\phi_L}$, $\partial f_L/\partial\phi_L = -f_L$):
$$\frac{\partial\bigl(\mathcal{L}_{\text{matter}} + \mathcal{L}_{\text{EM}}^{\text{disf}}\bigr)}{\partial\phi_L} \;=\; f_L(\phi_L)\,\Bigl(\varepsilon_{\text{baryon}}(\mathbf{x}) + u_B(\mathbf{x})\Bigr) \;=\; f_L(\phi_L)\,c^2\Bigl(\rho_{\text{baryon}}(\mathbf{x}) + \rho_B(\mathbf{x})\Bigr), \qquad \rho_B(\mathbf{x}) \;\equiv\; \frac{|\mathbf{B}(\mathbf{x})|^2}{2\mu_0 c^2}.$$
Ezzel egzaktul, kovariáns hatáselvből levezettük, hogy a V23-as diszform órametrika ($\sqrt{-\tilde{g}} = f_L\sqrt{-g}$, $\tilde{g}_{ij} = g_{ij}$) automatikusan kiválasztja a $T_{00}^{\text{EM}}/c^2 = \rho_B = \frac{|\mathbf{B}|^2}{2\mu_0 c^2}$ mágneses tömeg-energia sűrűséget a $\phi_L$ skalármező forrásaként!

---

## 3. A V20 polaritás-páros ($\Delta y_{\text{even}}$) háttér-kereszttag-kiejtés kovariáns levezetése

A laboratóriumban az eredő mágneses tér a Föld (és a laboratóriumi környezet) sztatikus $\mathbf{B}_{\text{bg}}(\mathbf{x})$ háttérterének és a 8 tekercs $I$ árammal arányos $\mathbf{B}_{\text{coil}}(\mathbf{x}; I) = I\,\mathbf{b}_{\text{coil}}(\mathbf{x})$ terének szuperpozíciója (`hu_p3.tex`, 84. sor):
$$\rho_B(\mathbf{x}; I) \;=\; \frac{|\mathbf{B}_{\text{bg}}(\mathbf{x}) + I\,\mathbf{b}_{\text{coil}}(\mathbf{x})|^2}{2\mu_0 c^2} \;=\; \underbrace{\frac{|\mathbf{B}_{\text{bg}}(\mathbf{x})|^2}{2\mu_0 c^2}}_{\rho_{B,\text{bg}}(\mathbf{x})\;(I^0)} \;+\; \underbrace{I\,\frac{\mathbf{B}_{\text{bg}}(\mathbf{x})\cdot\mathbf{b}_{\text{coil}}(\mathbf{x})}{\mu_0 c^2}}_{\rho_{B,\text{cross}}(\mathbf{x}; I)\;(I^1,\;\text{páratlan})} \;+\; \underbrace{I^2\,\frac{|\mathbf{b}_{\text{coil}}(\mathbf{x})|^2}{2\mu_0 c^2}}_{\rho_{B,\text{coil}}(\mathbf{x}; I)\;(I^2,\;\text{páros})}.$$
Mivel mind a 3D top-hat ($\mathcal{K}_L^{\text{TH}}$), mind a 3D Yukawa ($\mathcal{G}_L^{\text{Yuk}}$) operátor lineáris a forrássűrűségben, az $\mathbf{x}_1 = \mathbf{0}$ oszcillátor-pozícióban a skalárpotenciál BE ($I = \pm I_0$) és KI ($I = 0$) különbsége:
$$\Delta\phi_L(\pm I_0) \;\equiv\; \phi_L(\mathbf{0}; \pm I_0) - \phi_L(\mathbf{0}; 0) \;=\; \pm I_0\,\gamma_L\bar{\rho}_{B,\text{cross}}(\mathbf{0}; 1) \;+\; I_0^2\,\gamma_L\bar{\rho}_{B,\text{coil}}(\mathbf{0}; 1).$$
Ha az oszcillátor még egy parazita lineáris Zeeman-érzékenységgel ($\kappa_Z\,\mathbf{B}_{\text{coil}}(\mathbf{0}; \pm I_0) = \pm\kappa_Z I_0\mathbf{b}_{\text{coil}}(\mathbf{0})$) is rendelkezik, a H1 relatív frekvenciaeltolódás ($\Delta y(\pm I_0) \simeq -\Delta\phi_L(\pm I_0) \pm \kappa_Z I_0 b_{\text{coil}}$) **polaritás-páros átlaga** (`hu_p3.tex`, 79. és 84. sor):
$$\boxed{\;\Delta y_{\text{even}}(L) \;\equiv\; \frac{\Delta y(+I_0) + \Delta y(-I_0)}{2} \;=\; -I_0^2\,\gamma_L\bar{\rho}_{B,\text{coil}}(\mathbf{0}; 1) \;=\; -\frac{G L^2}{c^2}\int d^3y\,K_L(\mathbf{0},\mathbf{y})\,\frac{|\mathbf{B}_{\text{coil}}(\mathbf{y}; I_0)|^2}{2\mu_0 c^2},\;}$$
amelyből mind a $\pm\kappa_Z I_0$ lineáris Zeeman-tag, mind a $\pm 2\mathbf{B}_{\text{bg}}\cdot\mathbf{B}_{\text{coil}}$ háttér-kereszttag **azonosan kiesik**!

---

## 4. Ellenőrzések (Checks)

1. **Konform vs. diszform Maxwell-nyom ellenőrzése:**
   - Konform metrika ($\tilde{g}_{\mu\nu} = f_L^2 g_{\mu\nu}$): $\sqrt{-\tilde{g}}\tilde{g}^{\mu\alpha}\tilde{g}^{\nu\beta} = f_L^4 f_L^{-2}f_L^{-2}\sqrt{-g}g^{\mu\alpha}g^{\nu\beta} = \sqrt{-g}g^{\mu\alpha}g^{\nu\beta}$, $\phi_L$-től független $\implies \delta S_{\text{EM}}^{\text{conf}}/\delta\phi_L \equiv 0$ $\checkmark$.
   - Diszform metrika ($\tilde{g}_{00} = f_L^2 g_{00}$, $\tilde{g}_{ij} = g_{ij}$): $\sqrt{-\tilde{g}}\tilde{g}^{ik}\tilde{g}^{jl} = f_L\cdot 1\cdot 1\,\sqrt{-g}g^{ik}g^{jl} \implies -\partial\mathcal{L}_{\text{EM}}^{\text{disf}}/\partial\phi_L\big|_{\mathbf{E}=\mathbf{0}} = -f_L u_B$ $\checkmark$.
2. **Polaritásfordítási ($\pm I_0$) paritás-ellenőrzés:**
   - $\Delta\phi_L(+I_0) + \Delta\phi_L(-I_0) = (+I_0 - I_0)\gamma_L\bar{\rho}_{B,\text{cross}} + 2 I_0^2\gamma_L\bar{\rho}_{B,\text{coil}} = 2\,\gamma_L\bar{\rho}_{B,\text{coil}}(I_0)$, amely $L = 0{,}630957\;\text{m}$-nél és $I_0 = 14{,}0\;\text{mA}$-nél pontosan visszaadja a V20 auditált $\Delta y_{\text{peak}} = -1{,}2017108\times 10^{-48}$ értékét $\checkmark$.
