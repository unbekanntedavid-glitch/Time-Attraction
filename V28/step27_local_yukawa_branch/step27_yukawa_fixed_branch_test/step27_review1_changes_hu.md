# 27. lépés (a V28 része): az 1. bírálati kör változásai

Alap: a `step27_V28_review_bundle.zip` külső bírálata (2026-10-29), S01–S04 pontok és az 5. szakasz nyitott kérdései. V28 DOI még nincs. A hatókör változatlan: kizárólag a rögzített, árnyékolatlan, univerzális lokális Yukawa-ág (α = q = 0,0795775, λ = 10 km); ez nem a 26. lépés újranyitása.

## S01: G-kalibráció az α-határkeresésben (`run_step27.py`)
- Új `g_inf_of(α) = G_lab/(1+α)`; minden próba-α ezt használja. A q-nál rögzített globális `G_INF` már csak α = q-nál szerepel.
- Tornyok: Δg(α) = −Kα/(1+α), K = 2πG_lab ρλ(1 − e^{−h/λ})·10⁸ µGal (`tower_K_uGal`); pontos gyök α_lim = E/(K − E), E = −(μ − 1,96σ).
- A `tower_model_residual_uGal(1.0, …)` most G_N/2-vel számol, a dokumentált G(α)-val összhangban.
- Óceán: az elsődleges érték a pontos inverzió, α_lim = él/(S − él); a lineáris él/S közelítésként szerepel. S az α → 0 érzékenység G_lab-bal.
- A March-végponti húrok G_lab-bal: Yukawa 32,22081 → **34,78485** mGal/km; Bouguer **86,97917** (változatlan).
- Új `back_substitution_tests()`: minden határpontnak vissza kell adnia a megengedett reziduum szélét; eredmény a JSON `back_substitution_tests` mezőjében (tűrés 1×10⁻⁹ relatív). **Mind a 11 teszt megfelelt**; legnagyobb relatív eltérés 8,8×10⁻¹¹ (óceán, véges differenciás S); tornyok ≤ 3,5×10⁻¹⁶.
- A modell-reziduumok α = q-nál változatlanok (minden toronyra és sűrűségre 0,0 eltérés; az óceáni meredekség +9,07449725 mGal/km, változatlan).

| Korlát | Régi α_lim | Új α_lim (pontos) | Régi R | Új R |
|---|--:|--:|--:|--:|
| óceán (elsődleges) | 0,00225098 (lineáris) | **0,00225606** | 35,352 | **35,273** |
| óceán, lineáris közelítés | — | 0,00225098 | — | 35,352 |
| óceán, March ±1× intervallum | 0,00090977 | 0,00091060 | 87,47 | 87,39 |
| óceán, March (81) forgatókönyv | 0,0034647 | 0,0035012 (lineáris 0,0034890) | 22,97 | 22,73 |
| thomas1989 / 2000 | 0,00697470 | 0,00650259 | 11,409 | 12,238 |
| thomas1989 / 2670 | 0,00522449 | 0,00486292 | 15,232 | 16,364 |
| romaides1994_493m / 2000 | 0,00245635 | 0,00228048 | 32,397 | 34,895 |
| romaides1994_493m / 2670 | 0,00183996 | 0,00170725 | 43,250 | 46,612 |
| romaides1997_568m / 2000 | 0,00071609 | 0,00066375 | 111,127 | 119,891 |
| romaides1997_568m / 2670 | 0,00053640 | 0,00049711 | 148,355 | 160,081 |
| bren / 2000 | 0,00809680 | 0,00755665 | 9,828 | 10,531 |
| bren / 2670 | 0,00606502 | 0,00564970 | 13,121 | 14,085 |
| bren, pozitív előjelű változat / 2000; 2670 | 0,00270; 0,00202 | 0,00251; 0,00188 | 29,4; 39,3 | 31,7; 42,4 |

A (81) forgatókönyv azért változott, mert a March-végponti érzékenység most G_lab-bal számol (S81 = 79,98 → 79,43 mGal/km), és pontos inverziót használunk. Az alkalmazott R-tartomány: 9,83–148,4 → **10,53–160,1** (forgatókönyvek, nem konfidenciaintervallum).

## S02: Jekeli 1990
- JSON: `alpha_lim_10km_published` = null; új mezők: `abstract_alpha_value` = 0,001, `lambda_status` = "unknown", `applied` = false; `ratio_to_bound_R` = null (korábban 79,58).
- A 10 km-es ábrapontot és a q/0,001 tagot eltávolítottuk a 10 km-es R-összesítésből (`verdict.min_R/max_R`). Szöveges, nem alkalmazott hivatkozásként marad.
- Erie (Cruz 1991): az 1 km-es érték és feltételes 10 km-es becslése (0,00088) külön marad, a döntésen kívül.

## S03: a regionális trendről szóló mondat (módszertan A3, HU és EN)
- Régi: „Regionális trend: csak a lineáris meredekség-illesztés távolítja el.”
- Új: „A regionális korrekció a forrás g_c/δg előfeldolgozásából származik; ezt itt nem rekonstruáljuk. A lineáris illesztés a megmaradó meredekséget becsüli, és egy maradó lineáris regionális hibát nem választ szét a Yukawa-jeltől.” Hasonló állítás máshol nem volt.

## S04: ábrafeliratok
- Óceáni sáv: „mért meredekség-sáv (95%)” → „idézett hiba ×1.96; feltételezett Gauss 1σ” (EN: "quoted error ×1.96; assumed Gaussian 1σ").
- Saját óceáni korlátgörbe: FELTÉTELES jelölés. Modellgörbék: „idealizált réteg-/féltérmodell” / „idealizált féltér”. A torony-panel címe ugyanezt a hibaszöveget használja. Az Adelberger-görbe 95%-os címkéje változatlan.

## 5. szakasz, nyitott kérdések (összefoglalók 2c, módszertan A8/B6/C, JSON)
- Óceán: az eredeti bizonytalanság, a szelvények súlyozása és a nyomásból számolt mélység nincs újrakonstruálva (JSON `open_questions`).
- Tornyok: a tényleges felszíni mintavétel, terep, mélységi sűrűség és a Poisson-folytatás válasza nincs végigszámolva.
- BREN kb. 100 µGal-os szisztematika: külön nuisance-kezelés kell; egy illusztratív, korlátos ±100 µGal-os nuisance a ρ = 2000 szélt 0,00756 → 0,01029-re tolja (ρ = 2670: 0,00565 → 0,00769), csak illusztráció.
- λ < 1 cm: a szövegben és a JSON-ban (`verdict.untested_candidates_to_retest`, a korábbi `what_would_relieve_the_tension` helyett) újratesztelendő, nem vizsgált jelölt, nem igazolt menekülő tartomány.
- Az R-tartományok forgatókönyvek, nem kombinált konfidenciaintervallum, soha nem σ.
- March (81): külön forgatókönyv; az eltérés előjel- vagy G-konvencióváltással nem oldható fel; a szerzői szándékot nem ellenőriztük.

## Döntési megfogalmazás
- Megmaradt: „A saját újraszármaztatás erős feszültséget jelez a rögzített, árnyékolatlan Yukawa-ággal; a kizárás megerősítéséhez a mérési válasz és a bizonytalanságkezelés ellenőrzése szükséges.”
- Kiegészítés: az erős feszültség csak abban az értelemben áll, hogy az idealizált réteg-/féltérmodellekben, a feltételezett hibaértelmezés mellett nagy numerikus eltérés adódik; nem konfidenciaszinttel rendelkező empirikus kizárás. R = határhoz viszonyított arány, FELTÉTELES marad. A Zumberge-adatok March 2019 révén átvett másodlagos forrásbeli adatok. A hiányzó elsődleges források külön szerepelnek (Zumberge 1991, Cruz 1991 cikkszöveg, Jekeli 1990 cikkszöveg).
