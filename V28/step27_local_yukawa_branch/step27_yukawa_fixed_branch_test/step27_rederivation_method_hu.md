# 27. lépés (a V28 része): az újraszámolás módszere – óceán és BREN

Sorozat: Time Attraction Law — Scalar-Tensor Yukawa Fit, V28 (V28 DOI még nincs; korábbi DOI-k: V24 10.5281/zenodo.23680715, V25 10.5281/zenodo.23699685, V26 10.5281/zenodo.23727585, V27 10.5281/zenodo.23761574).

Ez a dokumentum azt mutatja meg, hogyan lesz egy publikált gravitációs reziduumból α-határ λ = 10 km-nél, a bíráló által kért két újraszámolásra: az óceáni kísérletre (Zumberge és mtsai. 1991, bemenetek March és mtsai. 2019 alapján) és a BREN-toronyra (LLNL-jelentés, OSTI 6190822). Minden itteni α_lim **saját, ellenőrizetlen újraszámolás**. Minden rá épülő R **feltételes**.

Konvenció: V(r) = −G_∞ m₁m₂/r · (1 + α e^{−r/λ}). Laboratóriumi mérésnél (r ≪ λ) G_lab = G_N = G_∞(1+α), vagyis **G_∞ = G_N/(1+α)**. Modellpont: α = 0,0795775, λ = 10 km, α > 0, univerzális, árnyékolatlan.

Kódtérkép (`run_step27.py`):

| Lépés | Függvény | Feladat |
|---|---|---|
| 1 | `yukawa_slab_field_infinite_layers` | Vízszintes rétegek Yukawa-tere (térfogati integrál) |
| 2 | `ocean_model_residual_slope` | Modell-reziduum a kísérlet saját definíciójával, és illesztett meredeksége |
| 3 | `ocean_alpha_lim_from_residual` | Az az α, amelynél a modell meredeksége eléri a mért intervallum szélét (pontos inverzió: él/(S − él)) |
| 4 | `march_eq81_F_difference` | Az 1. lépés ellenőrzése March (81) egyenlete és numerikus nyers integrál ellen |
| 5 | `tower_model_residual_uGal` | Féltér-toronyreziduum (mért − folytatott előrejelzés) |
| 6 | `tower_alpha_lim_from_residual` | Az az α, amelynél a modell-reziduum eléri μ − kσ-t (pontos gyök: E/(K − E)) |
| 7 | `cruz_conditional_scaling_factor` | Feltételes 1 km → 10 km skálázás az Erie-toronyra (nem korlát) |
| 8 | `back_substitution_tests` | Ellenőrzi, hogy minden határpont visszaadja a megengedett reziduum szélét (megfelelt/nem felelt meg a JSON-ban) |

G-kezelés minden α-keresésben (1. bírálati kör, S01): a mért G_lab rögzített, ezért minden próba-α G_∞(α) = G_lab/(1+α)-t használ; az α → 0 érzékenységekben G_lab szerepel.

Minden közbülső szám a `step27_results.json` → `review_trace` mezőben található.

## A. Óceán: Zumberge és mtsai. 1991, March és mtsai. 2019 alapján

**Adatstátusz: másodlagos forrásból vett adat (nem közvetlenül ellenőrzött elsődleges adat).** Magát a Zumberge-cikket (PRL 67, 3051) **nem értük el**. Minden alábbi bemenet March és mtsai., Phys. Rev. D 100, 042002 (2019), arXiv:1904.12789v2 idézetéből származik. Az oldalszámok erre a PDF-re vonatkoznak.

### A1. A reziduum definíciója (mit vonnak le, melyik G-vel)
- March 11. o., (72) egy.: korrigált gravitációs különbségek, g_c(z_N) − g_c(z_N,1). A nyers gravimétermérést korrigálják a műszerdriftre, valamint (i) időbeli hatásokra (Eötvös-hatás, földi és óceáni árapály, a tengeralattjáró függőleges gyorsulása) és (ii) gravitációs zavarokra δg (izosztatikusan kompenzált helyi tengerfenék-topográfia, regionális sűrűség-inhomogenitások az óceán alatt). Ezek mértékére és bizonytalanságára March Zumbergére hivatkozik.
- March 11. o., (73) egy.: **Δg_obs(z_N) = g_c(z_N) − g_c(z_N,1) − β_N (z_N − z_N,1)**. Ez a **G = G_N (laboratóriumi G)** értékkel számolt newtoni érték feletti többlet. β_N a newtoni függőleges gradiens a tengervízben (March 9. o.: β_N = 2,226×10⁻⁶ s⁻² az oceanográfiai használatban; 10. o., (67) egy.: c₀ = β_N + 4πG_Nρ_w).
- A mért szelvények tartománya z_N,1 = 500 m … z_N,2 = 4800 m (March 11. o.).

### A2. Közölt érték és bizonytalanság
- March 11. o., (74) egy.: a négy szelvény átlagos Δg_obs-meredekségére illesztett érték **dΔg_obs/dz_N = −0,060 ± 0,172 mGal/km**.
- **A ±0,172 statisztikai jelentését March nem adja meg, Zumbergénél pedig nem tudtuk ellenőrizni.** March a (75)–(76) egyenletben merev intervallumként használja: −0,232 ≤ meredekség ≤ +0,112 mGal/km (±1 × 0,172).
- **Saját feltevésünk:** a ±0,172 1σ Gauss-bizonytalanság.

### A3. Geofizikai háttér
- March 7. o., (46) egy.: víz 0–z_w, z_w = 5000 m; óceáni kéreg z_w–z_c, z_c − z_w = 7000 m, ρ_c = 2860 kg/m³; köpeny z_c alatt, ρ_m = 3400 kg/m³.
- March 12. o., (79) egy. és szöveg: a tengervíz sűrűsége 1023,6 → 1050,5 kg/m³; lineáris profil átlaga ρ̄_w = 1037,05 kg/m³.
- A rétegek végtelen vízszintes lapok. March ΔY és ΔG oldalirányú tagjait elhanyagoljuk (12. o., a (81) előtt). A 36 m-es üledékréteg hatása elhanyagolható (7. o.).
- Terep, izosztázia és regionális hatások csak Zumberge g_c-be épített korrekcióin keresztül szerepelnek (A1). Saját terepmodellt nem adunk hozzá.
- Szabad levegő: a szabadlevegő-rész a pályákból ismert GM_⊕-ből jön, az adatban és a referenciában azonos, így kiesik.
- A regionális korrekció a forrás g_c/δg előfeldolgozásából származik; ezt itt nem rekonstruáljuk. A lineáris illesztés a megmaradó meredekséget becsüli, és egy maradó lineáris regionális hibát nem választ szét a Yukawa-jeltől.

### A4. Képlet: reziduum → α
Modell-reziduum, a (73) egyenlettel azonos definícióval (March 63–64. egyenletei ugyanezt a szerkezetet adják):

  Δg_m(z) = −4π (G_∞ − G_N) ρ̄_w (z − z₁)  +  2π G_∞ α [F(z) − F(z₁)],   ahol G_∞ = G_N/(1+α)

  F(z) = ∫₀^∞ ρ(z′) sign(z′ − z) e^{−|z′−z|/λ} dz′   (térfogati integrál végtelen lapokra; March 47. és 66. egy., F = 𝓕)

- Az első tag a G-renormálás: az adat referenciája G_N-t használ, a valódi newtoni Bouguer-típusú gradiens viszont G_∞-t. Kis α-ra ez +4πG_N ρ̄_w α (z − z₁).
- A modell meredeksége, s(α), Δg_m legkisebb négyzetes meredeksége a z ∈ [500, 4800] m tartományon.
- Rögzített G_lab mellett mindkét tag pontosan α/(1+α) szerint skálázik: **s(α) = S·α/(1+α)**.
- α → 0 érzékenység: S = **123,108 mGal/km egységnyi α-ra**. A 0–5000 m-es végponti húrok G_lab-bal: **86,97917** (G-renormálás, Bouguer-típusú tag) és **34,78485** (Yukawa-térfogati integrál). Ezek húrértékek, nem az illesztett S közvetlen összetevői. (Az előző változat 32,2-t adott a Yukawa-húrra, mert G_∞-t q-nál rögzítette.)
- A modellpontban s(q) = **+9,07 mGal/km** (= S·q/(1+q)).
- Él: s_obs + k·σ, k = 1,95996 mellett 0,27711 mGal/km.
- **Pontos inverzió (elsődleges): α_lim = él/(S − él) = 0,00225606**, R = 35,27 (kerekített 0,0023-mal 34,6).
- Lineáris közelítés, csak közelítésként: α_lim ≈ él/S = 0,00225098, R = 35,35 (0,23%-os eltérés; a kerekített ≈0,0023, ≈35 nem változik). **R feltételes.**

### A5. Az α_lim statisztikai jelentése
Ha σ = 0,172 1σ Gauss-hiba és k = 1,96, a határ a meredekség központi 95%-os intervallumának felső széle. α > 0 pozitív meredekséget jósol, ezért α_lim **egyoldali 97,5%-os felső határ**. Egyoldali 95%-ra (k = 1,645) az érték 0,00181 lenne. Ez a jelentés csak akkor igaz, ha a hiba Gauss-eloszlású, 1σ, és statisztikai hiba dominálja. Ezek egyike sincs ellenőrizve.

### A6. Saját feltevéseink, amelyeket a forrás nem tesz
1. Tiszta, univerzális, árnyékolatlan Yukawa. March NMC „extra erejét” (nyomás–mélység tag) elhagyjuk.
2. A ±0,172-t 1σ Gauss-hibának olvassuk (March ±1× intervallumként használja).
3. A köpeny-féltér a végtelenig tart. March (7. o.) ezzel szemben azt feltételezi, hogy a Yukawa-erő a felső köpeny alatt elnyomott, ezt szerinte „utólag ellenőrizni kell”.
4. A meredekséget a saját modellünkre illesztjük 500–4800 m-en, egyenletes rácson. Zumberge a mért szelvényekre illesztett.
5. A mélységet közvetlenül z-nek vesszük. March (10. o.) szerint a nyomás–mélység átváltás hatása α ~ 10⁻³-nál cm-es nagyságrendű; α ≈ 0,08-ra ezt nem értékeltük újra.
6. A λ = 10 km kívül esik Zumberge publikált határának 1–5000 m-es tartományán. A kiterjesztés csak a fenti rétegmodellen alapul.

### A7. Ellenőrzés (a kódból)
- Saját rétegintegrálunk: F(z_w) − F(0) = 4,147×10⁶ kg/m². Független numerikus nyers integrál: 4,148×10⁶ (egyezés 10⁻⁴-nél jobb).
- **March (81) egyenlete a nyomtatott formában** (12. o.) −9,006×10⁵ kg/m²-t ad. A különbség, 5,048×10⁶ kg/m², pontosan λ(1 − e^{−z_w/λ})(2ρ_c − ρ_w − ρ_m). Közvetlen rétegintegrálunk összhangban van March (46)–(47) bemeneti profiljával; a nyomtatott (48)/(50)/(81) redukált kifejezések nem ugyanazt az algebrai rétegeredményt adják. Az eltérés **puszta előjel- vagy G-konvencióváltással nem oldható fel**; a szerzői szándékot vagy helyesbítést **nem ellenőriztük**.
- **Külön forgatókönyv** (nem statisztikai hibasáv): ha a (81) egyenletet nyomtatott formájában használnánk, a végponti húr-érzékenység 86,98 − 7,55 = 79,4 mGal/km per α, így α_lim = 0,00350 (pontos inverzió; lineárisan 0,00349), R = 22,7.
- March saját ±1× intervallumával a határ 0,112, így α_lim = 0,00091 (pontos), R = 87,4.
- Ellenőrzés λ ≲ 100 m-nél: a modell a ×1,96-os szélnél α_lim ≈ 0,0032-t ad, ez nagyságrendben egyezik Zumberge publikált 0,002-es értékével (March 6. o. idézete).
- Összességében az óceáni α_lim 10 km-nél **≈0,0009–0,0035, R ≈ 23–87**, attól függően, melyik ellenőrizhetetlen választást tesszük. Ezek **forgatókönyvek**, nem egyetlen kombinált konfidenciaintervallum; R soha nem σ-szignifikancia.

### A8. Nyitott kérdések (nincs újrakonstruálva)
- A meredekség eredeti bizonytalansága, a mérési szelvények súlyozása és a nyomásból számolt mélység modellfüggése nincs újrakonstruálva.
- Az analitikus rétegkontroll csak az idealizált válasz helyességét igazolja, nem a kísérlet tényleges válaszát.

## B. BREN-torony: LLNL-jelentés, OSTI 6190822

**Adatstátusz: elsődleges dokumentum, előzetes eredmények** (Thomas, Kasameyer és mtsai., „A test of Newton's Law of Gravity using the BREN Tower, Nevada”, LLNL konferenciajelentés, https://www.osti.gov/biblio/6190822; 16 oldal). A folyóirat-változatból, Thomas és mtsai., PRL 63, 1902 (1989), csak az absztrakt volt elérhető.

### B1. A reziduum definíciója
- Jelentés 6–7. o. (4–7. egy.): a felszíni gravitációs anomáliát, Δg-t, Poisson-integrállal folytatják felfelé, gömbi közelítésben, minden felszíni adatot azonos magasságúnak tekintve. A tornyon kapott előrejelzést a mért Δg(r)-különbségekkel vetik össze. A centrifugális és az árapályhatást minden adatból eltávolítják (7. o.).
- Reziduum = a Δg magasság szerinti változásának mért és előrejelzett értéke közti különbség. **A vezető tagban nincs laboratóriumi G**: az előrejelzés a mért gravitáció harmonikus folytatása.

### B2. Közölt érték és bizonytalanság
- Jelentés 11. o. (4. szakasz): „Our predictions are within 93 µGal of observations at all heights on the tower, and the uncertainties at each height are as large as 95 µGal. The error bars in Figure 2 indicate the single standard deviation uncertainty from sampling errors.”
- Az absztrakt (3. o.) szerint az eredmény „93 ± 95 µGal-on belül” egyezik a newtoni hipotézissel a torony tetején.
- **A 93 µGal előjele nincs megadva.** Óvatosan **μ = −93 µGal** értéket használunk (az α > 0 jellel azonos előjel), **σ = 95 µGal** mellett (1σ, csak mintavételi hiba), a toronytetőn, h = 454 m-en.

### B3. Geofizikai háttér
- Jelentés 11. o.: 91 állomás 2,5 km-en belül, valamint 60 000 USGS/NOAA-pont 300 km-ig. Gyűrűnként szektorátlagokat képeznek, és a bizonytalanság becslése előtt egy „regionális” síkfelületet vonnak le (13. o.).
- Jelentés 12. o.: a gömbi közelítés „potentially our most significant source of systematic error”. Egy terepmodell (sűrűség 1,8 g/cm³) szerint **kb. 100 µGal-lal alulbecsüli Δg(r₁) − Δg(r₀)-t**. Ez a szisztematikus hiba **nincs** benne a 95 µGal-ban.
- Saját modellünk a terepet egyáltalán nem veszi figyelembe (sík féltér).

### B4. Képlet: reziduum → α
Homogén ρ sűrűségű féltér (feltételezett 2000 vagy 2670 kg/m³; nem a jelentésből). A felszíni tér Yukawa-részét, 2πG_∞αρλ-t, a folytatás állandónak veszi, a valódi tér viszont h magasságban 2πG_∞αρλe^{−h/λ}:

  Δg(h) = −2π G_∞(α) α ρ λ (1 − e^{−h/λ}) = −K α/(1+α),   G_∞(α) = G_lab/(1+α),   K = 2π G_lab ρ λ (1 − e^{−h/λ})

A G_lab minden próba-α-nál rögzített (1. bírálati kör, S01). h = 454 m-en és λ = 10 km-nél K = 37 226 µGal (ρ = 2000), illetve 49 697 µGal (ρ = 2670); ez az α → 0 érzékenység. (α = 1-nél a reziduum −K/2, vagyis G_N/2-vel számolva.)

k = 1,96 mellett az él μ − kσ = −279,2 µGal, E = 279,2 µGal, a pontos gyök **α_lim = E/(K − E) = 0,00756 (ρ = 2000), illetve 0,00565 (ρ = 2670), így R ≈ 10,5–14,1 (feltételes)**. Az előző változat G_∞-t q-nál rögzítette, és 0,0081/0,0061-et adott. Ha az előjel +93 µGal volna, az eredmény 0,0025 (ρ = 2000), illetve 0,0019 (ρ = 2670) lenne, R ≈ 32–42.

### B5. Statisztikai jelentés
Egyoldali 97,5%-os felső határ, de csak akkor, ha σ Gauss-eloszlású 1σ hiba és a szisztematikus hibák elhanyagolhatók. Ez utóbbinak ellentmond a jelentés saját, kb. 100 µGal-os gömbi közelítési szisztematikája.

### B6. Saját feltevéseink, amelyeket a forrás nem tesz
1. Tiszta univerzális Yukawa; sík féltér; homogén, feltételezett sűrűség; a terepet elhagyjuk.
2. A 93 µGal reziduum előjele (óvatos választás).
3. Csak a toronytető pontját használjuk; a 2. ábra magassági profilját nem digitalizáltuk.
4. A kb. 100 µGal-os szisztematikus hibát nem adjuk hozzá. Ez külön nuisance-kezelést igényel, és nem adható automatikusan független Gauss-szórásként a 95 µGal-hoz. Csak illusztráció: egy korlátos, ±100 µGal-os nuisance a legrosszabb szélén a ρ = 2000 szélt 0,00756-ról 0,01029-re módosítaná (ρ = 2670: 0,00565 → 0,00769); ez nem új adatillesztés és nem empirikus korlát.
5. A felszíni Yukawa-teret a felmért terület egészén egyenletesnek vesszük, holott az adatok 300 km-ig, vagyis λ-nál jóval messzebbre terjednek.

## C. Ismert hiányok
- Zumberge 1991-et nem értük el, így a ±0,172 jelentése, a korrekciók és bizonytalanságaik ellenőrizetlenek.
- Saját rétegintegrálunk és March nyomtatott (81) egyenlete közti eltérés tisztázatlan.
- Egyik újraszámolás sem tartalmazza a kísérlet saját válaszát egy Yukawa-térre a valódi geometriában (véges felmérési terület, terep, a referenciamodell felépítése). Ez az a „mérési válasz”, amelyet minden kizárási állítás előtt ellenőrizni kell.
- Tornyok: a tényleges felszíni mintavétel, a terep, a mélységi sűrűség és a Poisson-folytatás (a jelentés (7) egyenlete) együttes Yukawa-válasza nincs végigszámolva. A tornyok sűrűsége feltételezett; a publikált magasság szerinti reziduumokat nem használtuk.
- Minden határpont átmegy a visszahelyettesítési teszten (visszaadja a megengedett reziduum szélét, relatív eltérés ≤ 1×10⁻¹⁰); ez csak az inverziót ellenőrzi, a fizikai modellt nem.
- Jekeli és mtsai. 1990: az absztraktbeli α < 0,001 λ-tartománya ismeretlen. Csak szöveges hivatkozás, nem alkalmazzuk (nincs az ábrán, nincs a 10 km-es R-összesítésben).
- Az Erie-torony (Cruz 1991) 1 km-es határa a kísérlet válaszfüggvénye nélkül nem vihető át 10 km-re, ezért a 0,00088 csak feltételes becslés (`cruz_conditional_scaling_factor`).
- Hiányzó elsődleges források (külön jelölve; az itt használt Zumberge-értékek March és mtsai. 2019 révén átvett másodlagos forrásbeli adatok): Zumberge és mtsai. 1991; a cikkszöveg Cruz és mtsai. 1991, Jekeli és mtsai. 1990, Thomas és mtsai. 1989, Romaides és mtsai. 1994/1997 esetén; nem értük el: Stacey és mtsai. 1987, Moore és mtsai. 1988, Ander és mtsai. 1989, Bartlett & Tew 1989, Achilli és mtsai. 1997, Rapp 1987, Fischbach & Talmadge 1999 (könyv).
