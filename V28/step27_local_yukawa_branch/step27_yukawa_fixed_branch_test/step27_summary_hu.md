# 27. lépés (a V28 része): a rögzített Yukawa-ág lokális gravitációs tesztje

Sorozat: Time Attraction Law — Scalar-Tensor Yukawa Fit, V28. A V28 a V27 26. lépésének kozmológiai csomagját, az interaktív szimulátort és ezt a 27. lépést foglalja magában. V28 DOI még nincs. Korábbi DOI-k: V24 10.5281/zenodo.23680715, V25 10.5281/zenodo.23699685, V26 10.5281/zenodo.23727585, V27 10.5281/zenodo.23761574.

Hatókör: ez a vizsgálat kizárólag a TAL skalár–tenzor Yukawa-illesztés **rögzített, árnyékolatlan, univerzális lokális Yukawa-ágára** vonatkozik (L = 10 km, α_Y = 1, β_Yuk = 1/(4π), C = 0 attraktor). Más TAL-változatokra az eredmény nem vihető át automatikusan.

## 1. Mit jósol a csomag

Konvenció: V(r) = −G_∞ m₁m₂/r · (1 + α e^{−r/λ}). Laboratóriumi mérésnél (r ≪ λ) G_lab = G_∞(1+α).

- A hatás a `hubble_desi_bao_test/tal_cosmology_common.py` 50–64. sorában található. A kinetikus tag c⁴/(4πGβ_Yuk), az anyag pedig univerzálisan csatolódik az A²(φ)g metrikához, ahol A = e^{α_Y φ}. Pontforrás tere: φ = −α_Y β_Yuk GM e^{−r/L}/(c²r). A próbatömegre ható potenciál ebből c²α_Y φ = −α_Y²β_Yuk GM e^{−r/L}/r.
  - Ezért **α = q = α_Y²β_Yuk = 1/(4π) = 0,0795775** és **λ = L = 10 km**.
  - Az előjel pozitív, vagyis vonzó.
- Ugyanez az eredmény a csomag három másik helyén is szerepel:
  - `step26_md_templates.py` 103. sor: G_lab ≈ G_*(1+q), „r ≪ L, árnyékolatlan kis testek” esetén.
  - `step5_coupling_dictionary_ppn_screening_en.md` 54., 73. és 101. sor: α_std = α_Y/(4π), valamint G_eff(r) = G[1+α_std(1+r/L)e^{−r/L}].
  - `run_step26_revision_controls.py` 244. sor: G_eff/G_lab = 1/(1+q) = 0,9262883.
- Árnyékolás: **nincs**.
  - A `step26_md_templates.py` 106. sora és a `run_step26_revision_controls.py` 253. sora szerint „no screening branch is introduced”.
  - A 11. lépés kaméleon/szimmetron mechanizmusa (`step11_..._en.md`, 18. sor) külön nemlineáris kiterjesztési hipotézis, a rögzített ágnak nem része.
  - Az 5. lépés 100. sorában a „screened” szó csak az e^{−r/L} lecsengést jelenti.
- Összetételfüggés: a csatolás univerzális, nincs ekvivalenciaelv-sértő jel. A teszt így a fordított négyzetes törvény próbájára szűkül.

CODATA G: minden laboratóriumi G-mérés r ≪ 10 km mellett történik, ezért mindegyik ugyanazt a G_∞(1+α) értéket méri. A mérések szórása így nem hordoz jelet. Jel csak az 1–10 km-es skálán jelenik meg.


## 2. Korlátok és modellmegfigyelhetők (λ = 10 km)

A **határhoz viszonyított arány**: R = 0,0795775 / α_lim(10 km).
- Csak a fenti konvencióban és a megadott feltevések mellett értelmezhető.
- R > 1 azt jelenti, hogy a modellpont az *adott rekord* megengedett tartományán kívül esik.
- R **nem** szignifikancia, σ-ban nem adjuk meg.
- **A táblázat minden R-értéke feltételes**: egyik nevező sem ellenőrzött, közvetlenül idézhető elsődleges korlát 10 km-nél.

Adatcímkék:
- **másodlagos forrásból vett adat (nem közvetlenül ellenőrzött elsődleges adat)**: egy későbbi cikk idézi az eredetit (Zumberge-bemenetek March és mtsai. 2019 alapján).
- **elsődleges absztrakt szövege**: az eredeti cikk absztraktjából; a cikk szövegét nem értük el.
- **elsődleges jelentés**: magát az eredeti dokumentumot megszereztük (BREN, OSTI 6190822).

Oszlopok: **Modell (a)** a rögzített modellpont, α = 0,0796 λ = 10 km-nél. **Modell (b)** a kísérlet saját megfigyelhetője az adott elrendezésben számolva.

| Forrás | Adatcímke | Megfigyelhető | Feltevések / CL | α_lim 10 km-nél | Modell (a) | Modell (b) | R (feltételes) |
|---|---|---|---|---|---|---|---|
| Zumberge és mtsai., PRL 67, 3051 (1991), doi:10.1103/PhysRevLett.67.3051; bemenetek March és mtsai. alapján, PRD 100, 042002 (2019), arXiv:1904.12789 | másodlagos forrásból vett adat (nem közvetlenül ellenőrzött elsődleges adat) | A reziduum meredeksége a vízoszlopban, 500–4800 m: −0,060 ± 0,172 mGal/km (March 74. egy.) | α > 0, univerzális, árnyékolatlan; idealizált rétegmodell; idézett hiba ×1,96 feltételezett Gauss 1σ-olvasat mellett (egyoldali 97,5% csak e feltevés alatt) | Publikált 0,002 csak 1–5000 m-re; **10 km-re nincs érték**. Újraszámolva (pontos inverzió): **0,00226** (≈0,0023); lineáris közelítés 0,00225 | 0,0796 | Meredekség **+9,1 mGal/km** | **≈35** (35,27; lineáris közelítéssel 35,35; kerekített 0,0023-mal 34,6) |
| Jekeli, Eckhardt, Romaides, PRL 64, 1204 (1990) | elsődleges absztrakt szövege | Kb. 600 m-es torony vs. felfelé folytatás | Előjel, λ-tartomány, CL nincs megadva | **Nem alkalmazzuk.** Absztraktbeli érték α < 0,001, λ-tartomány ismeretlen; 10 km-es érték nincs (csak szöveges hivatkozás) | 0,0796 | nem számoltuk | nem alkalmazzuk |
| Thomas és mtsai., PRL 63, 1902 (1989) | elsődleges absztrakt szövege | 465 m-en (−60 ± 95) µGal | Idealizált féltér, ρ = 2000–2670 (feltételezett), terep elhagyva; a ± 1σ-nak olvasva | Újraszámolva: 0,0049–0,0065 | 0,0796 | −2,8 … −3,7 mGal | ≈12–16 |
| Romaides és mtsai., PRD 50, 3608 (1994) | elsődleges absztrakt szövege | 493 m-en (−33 ± 30) µGal | Mint fent | Újraszámolva: 0,0017–0,0023 | 0,0796 | −3,0 … −4,0 mGal | ≈35–47 |
| Romaides és mtsai., PRD 55, 4532 (1997) | elsődleges absztrakt szövege | 568 m-en (32 ± 32) µGal | Mint fent | Újraszámolva: 0,00050–0,00066 | 0,0796 | −3,4 … −4,6 mGal | ≈120–160 |
| Thomas, Kasameyer és mtsai., LLNL-jelentés, OSTI 6190822 (előzetes BREN-elemzés) | elsődleges jelentés | „minden magasságban 93 µGal-on belül”, bizonytalanság legfeljebb 95 µGal (1σ, mintavétel); előjel nincs megadva | Óvatos választás: μ = −93 µGal 454 m-en; ρ feltételezett | Újraszámolva: 0,0056–0,0076 (pozitív előjel esetén 0,0019–0,0025); a kb. 100 µGal-os szisztematika nincs benne (lásd 2c) | 0,0796 | −2,7 … −3,7 mGal | ≈11–14 |
| Cruz, Harrison, Speake, Niebauer, McHugh, Keyser, Faller, Mäkinen, Beruff, J. Geophys. Res. 96(B12), 20073 (1991), doi:10.1029/91JB01756 | elsődleges absztrakt szövege | 300 m-es Erie-torony; rms eltérés 10 µGal | \|α\|-ra; CL nincs megadva | **\|α\| < 0,001 λ = 1000 m-nél**. 10 km-nél ≈0,00088 csak **feltételes becslés** (lásd lent) | 0,0796 | nem számoltuk | 80 (1 km-nél); 10 km-nél nem használjuk |
| Kontextus: Adelberger, Heckel, Nelson, Annu. Rev. Nucl. Part. Sci. 53, 77 (2003), 4. ábra | digitalizált összefoglaló | — | 95% | ≈4×10⁻⁴, digitalizálva; ≈ ±0,1 dex | 0,0796 | nem számoltuk | csak kontextus |
| Kontextus: Fischbach & Talmadge, arXiv:hep-ph/9606249 (1996), 1. ábra | digitalizált összefoglaló | — | 2σ; a tengely −α, vagyis **csak α < 0** | ≈4×10⁻⁴ (3–5×10⁻⁴), digitalizálva; ≈ ±0,15 dex | 0,0796 | nem számoltuk | nem használjuk (ellentétes előjel) |

**Erie (Cruz és mtsai. 1991), 10 km-es érték.** A 0,00088 feltételes becslés, nem korlát. Egy 1 km-es határ nem vihető át 10 km-re pusztán a Yukawa-exponenciál alapján. Ehhez a torony mérési különbségének részletes válaszfüggvénye kellene, a forrásgeometriával és a G-normálással együtt, ez pedig nincs meg (a kiadói szöveg HTTP 403-at adott). Skálázásunk feltevései: (i) az 1 km-es határt a toronytető reziduuma szabja meg h = 300 m-en; (ii) sík, homogén sűrűségű féltér; (iii) a felszíni tér Yukawa-részét a folytatás állandónak veszi; (iv) mindkét λ-nál azonos G-normálás. Így a faktor S(1 km)/S(10 km), ahol S(λ) = λ(1 − e^{−h/λ}), azaz 259,2/295,5 = 0,877.

**A korlátok G-normálása (1. bírálati kör, S01).** A mért G_lab minden próba-α-nál rögzített, ezért G_∞(α) = G_lab/(1+α). Tornyok: Δg(α) = −Kα/(1+α), K = 2πG_lab ρλ(1 − e^{−h/λ}), a pontos gyök α_lim = E/(K − E), ahol E = −(μ − 1,96σ). Óceán: s(α) = Sα/(1+α), α_lim = él/(S − él). Az előző változat az α-keresés közben a q modellpontnál rögzítette G_∞-t; toronyhatárai 7–8%-kal túl nagyok voltak (pl. BREN ρ = 2000: 0,00810 → 0,00756). A modell-reziduumok α = q-nál változatlanok. A kódban lévő visszahelyettesítési teszt igazolja, hogy minden határpont visszaadja a megengedett reziduum szélét (relatív eltérés ≤ 1×10⁻¹⁰).

Az újraszámolás módszere (képletek, forrásbeli egyenletszámok, ismert hiányok): `step27_rederivation_method_hu.md`.

## 2a. A döntő korlát és státusza

- **Döntő rekord: az óceáni újraszámolás**, α_lim(10 km) ≈ 0,0023, **R ≈ 35 (feltételes)**. R-státusz: *feltételes: a nevező, α_lim ≈ 0,0023, ellenőrizetlen újraszámolás*. A bemenetek másodlagos forrásból vett adatok (nem közvetlenül ellenőrzött elsődleges adatok).
- Érzékenység olyan választásokra, amelyeket nem tudtunk ellenőrizni:
  - ha a ±0,172-t intervallumhatárnak olvassuk, ahogy March és mtsai. a (75)–(76) egyenletben: α_lim ≈ 0,00091, R ≈ 87;
  - ha a rétegintegrálra March (81) egyenletét használjuk a nyomtatott formában: α_lim ≈ 0,0035, R ≈ 23 (pontosan 0,00350, R = 22,7). Ez **külön forgatókönyv**, nem statisztikai hibasáv. Saját közvetlen rétegintegrálunk és a nyers numerikus integrál eltér a nyomtatott (81)-től; az eltérés puszta előjel- vagy G-konvencióváltással nem oldható fel, a szerzői szándékot vagy helyesbítést nem ellenőriztük.
- Az óceáni 0,0009–0,0035 (R ≈ 23–87) értékek **forgatókönyvek** (hibaértelmezés, rétegképlet), nem egyetlen kombinált konfidenciaintervallum. R soha nem σ-szignifikancia.
- Az óceáni α_lim a pontos inverzió, él/(S − él) = 0,00225606 (R = 35,27); a lineáris él/S = 0,00225098 (R = 35,35) csak közelítésként szerepel.
- Nem találtunk olyan elsődleges forrást, amely α > 0-ra pontosan λ = 10 km-nél számszerű határt adna.

## 2b. Hiányzó elsődleges források

A felhasznált adatoktól külön jelölve. Az itt használt Zumberge-értékek **March és mtsai. 2019 révén átvett másodlagos forrásbeli adatok**, nem közvetlenül ellenőrzött elsődleges adatok. A legfontosabb hiányzó elsődleges források: **Zumberge és mtsai. 1991** (nem értük el), **Cruz és mtsai. 1991 cikkszövege**, **Jekeli és mtsai. 1990 cikkszövege**.

- Zumberge és mtsai., PRL 67, 3051 (1991): **nem értük el**; az óceáni bemenetek March és mtsai. 2019-ből származnak.
- Cruz és mtsai., JGR 96, 20073 (1991): csak absztrakt (kiadói PDF: HTTP 403).
- Jekeli és mtsai., PRL 64, 1204 (1990): csak absztrakt.
- Thomas és mtsai., PRL 63, 1902 (1989): csak absztrakt; az előzetes LLNL-jelentés megvan.
- Romaides és mtsai. 1994 és 1997: csak absztrakt.
- Nem értük el: Stacey és mtsai. RMP 59, 157 (1987); Moore és mtsai. PRD 38, 1023 (1988); Ander és mtsai. PRL 62, 985 (1989); Bartlett & Tew PRL 63, 1531 (1989); Achilli és mtsai. 1997; Rapp 1987; Fischbach & Talmadge 1999 (könyv).

## 2c. Nyitott kérdések (a 27. lépésben nincs újrakonstruálva)

- **Óceán:** az eredeti bizonytalanság, a mérési szelvények súlyozása és a nyomásból számolt mélység modellfüggése nincs újrakonstruálva. Az analitikus rétegkontroll csak az idealizált válasz helyességét igazolja.
- **Tornyok:** a tényleges felszíni mintavétel, a terep, a mélységi sűrűség és a Poisson-folytatás együttes Yukawa-válasza nincs végigszámolva.
- **BREN:** a jelentés kb. 100 µGal-os szisztematikája külön nuisance-kezelést igényel; nem adható automatikusan független Gauss-szórásként a 95 µGal-hoz. Csak illusztráció: egy korlátos, ±100 µGal-os nuisance a legrosszabb szélén a ρ = 2000 féltér α-szélét 0,00756-ról 0,01029-re módosítaná (ρ = 2670: 0,00565 → 0,00769). Ez nem új adatillesztés és nem empirikus korlát.
- **Az R-tartományok** a különböző hibaértelmezések, rétegképletek és feltételezett sűrűségek forgatókönyvei, nem egyetlen kombinált konfidenciaintervallum; R egyik esetben sem σ-szignifikancia.
- **March (81):** külön forgatókönyvként marad; az eltérés előjel- vagy G-konvencióváltással nem oldható fel; a szerzői szándékot nem ellenőriztük.
- **λ < 1 cm:** csak újratesztelendő, nem vizsgált jelölt (lásd 3. szakasz); a 27. lépés nem vizsgálta, hogy ez a tartomány megengedett-e.

## 3. Döntés

**A saját újraszármaztatás erős feszültséget jelez a rögzített, árnyékolatlan Yukawa-ággal; a kizárás megerősítéséhez a mérési válasz és a bizonytalanságkezelés ellenőrzése szükséges.**

Az „erős feszültség” jelentése: **az idealizált réteg-/féltérmodellekben, a feltételezett hibaértelmezés mellett nagy numerikus eltérés adódik.** Ez nem konfidenciaszinttel rendelkező empirikus kizárás: a tényleges mérési válasz, a háttérkorrekciók és a szisztematikus bizonytalanságok nincsenek újrakonstruálva.

Hatókör: kizárólag a rögzített, árnyékolatlan, univerzális lokális Yukawa-ág (α = 0,0796, λ = 10 km); nem minden TAL-változat, és nem a 26. lépés újranyitása. R a **határhoz viszonyított arány**, és **FELTÉTELES** marad; minden alkalmazott feltételes R 1 fölött van (≈11–160, forgatókönyvek, nem konfidenciaintervallum). R arány, nem σ-érték.

Nem vizsgált jelöltek (mindegyik új ágat jelent, amelyet újra kell tesztelni; egyik sem igazolt):
- földi és laboratóriumi forrásokra ható árnyékolás (a csomagban nincs ilyen);
- α ≲ 10⁻³ λ = 10 km-nél, vagyis α_Y²β_Yuk legalább kb. 80-szoros csökkentése;
- más λ (például λ < 1 cm) változatlan α-val: **újratesztelendő jelölt**, **nem** igazolt menekülő tartomány; a 27. lépés nem vizsgálta, hogy ez a tartomány megengedett-e.

Az 1. bírálati kör változásai: `step27_review1_changes_hu.md`.

Futtatás: `python run_step27.py` → `step27_results.json` (szigorú JSON, benne a `back_substitution_tests` mezővel), valamint az ábrák: `step27_alpha_lambda_{en,hu}.png` és `step27_observables_{en,hu}.png`.
