# Harmadik változat – célzott független ellenőrzés

**Dátum:** 2026-10-21. **Bemenet:** `final_report (3).zip`, `report_hu (3).pdf`. **Kért ellenőrzés:** a korábbi L≈3 m-es integrálillesztési hiba és a háttérmágneses kereszttag kezelése. Az eredeti fájlokat nem módosítottam; publikációs felületet nem érintettem. A Bothwell-gradiensfeladat ettől külön marad.

## Ítélet

**A két kért javítást a megadott modell és az új szöveges lehatárolás mellett elfogadom.**

1. Az L≈3 m-es korábbi nagy numerikus eltérés megszűnt. Az új külső integrál más szögkvadratúrával végzett ellenőrzése megerősíti a javítást.
2. A jelentés most kimondja a háttértér 2B_bg·B_coil kereszttagját, az ideális ±I átlagolás feltételét és azt, hogy a közölt jel az izolált, polaritásra páros tekercs-B² tag.

A fő jel továbbra is **körülbelül −1,20×10⁻⁴⁸**. Ez elfogadható számított nagyságrend a választott ideális áramlapmodellben, **nem TAL-kimutatás és nem mért készülékérzékenység**.

Három kisebb, illetve kísérleti értelmezési pontosítás marad: a távoli tartomány „0 szórása” nem hibabecslés; az L=0,95 m-nél nyíltan közölt kis illesztési hiba megmarad; az azonos összteljesítményű hőkontroll nem garantál automatikusan azonos hőmérsékletmezőt. Ezeket alább elkülönítem a már kijavított fő problémáktól.

## 1. Mit ellenőriztem most ténylegesen?

- **18/18 SHA-256 bejegyzés egyezik.** A külön feltöltött PDF bájtról bájtra azonos a ZIP-beli jelentéssel.
- Átolvastam az új jelentést, az új `energy_integrator.py` modult és a `final_fixed.py` változásait. A PDF táblázatait a nyomtatott oldalaképen is ellenőriztem.
- **Célzottan újrafuttattam** az archivált integrátort a fő csúcspontban, valamint a külső A/B útvonalakat a releváns L-eknél.
- Külön köráram-képlettel ellenőriztem a teret; Descartes-koordinátás Sobol-kubatúrával a fő jelet; más szög- és radiális kvadratúrával a távoli energiakülönbséget.

**Nem** állítok új teljes `final_fixed.py`-futtatást: a teljes négyváltozatos konvergencia-, statikus Monte Carlo- és ábragenerálást ebben a körben nem ismételtem meg. A kért két javítás és a fő csúcs célzott ellenőrzése megtörtént.

Archív SHA-256: `160981d3adb435b4fa719ca1a9d9a04d4e889a27f9d4c7bb0381e5c057ad2ea8`.

## 2. Fő jel: reprodukálható és külön is ellenőrizhető

Az új, szerzői mérési jegyzőkönyv 291,7 mm-es menetátmérőt ad meg; a kód most valóban ezt használja, nem az ellenállásból származtatja a modell sugarát. A 300 mm-es méretet peremátmérőként határozza meg. A geometriai eredet korábbi tisztázatlansága ezzel dokumentumszinten rendezett. A jegyzőkönyvet szerző által szolgáltatott mérési adatként kezelem, nem általam hitelesített helyszíni mérésként.

| Mennyiség | Mostani ellenőrzés |
|---|---:|
| B a fő órán, 14 mA | **11,632900 µT** |
| B a referencián | **0,430539 µT** |
| Független köráram-teszt legnagyobb relatív eltérése, 5 pont | **1,04×10⁻¹⁴** |
| Archivált integrátor célzott újrafuttatása, L=0,630957 m | **−1,2017108×10⁻⁴⁸** |
| Eltérés az archivált fő eredményhez képest | **2,2×10⁻¹⁶ relatív**, kerekítési szint |
| Független térfogati Sobol-integrál ugyanott | **−1,2008435×10⁻⁴⁸** |
| Független magnitúdó eltérése a csomagtól | **−0,0722%** |
| Sobol-keverésekből becsült standard hiba | **3,04×10⁻⁵²** |

A Sobol-ellenőrzés a fő gömbhöz 4×2¹⁹ pontot, a sima referencia-gömbhöz 4×2¹² pontot használ. A megadott standard hiba nem a teljes modell bizonytalansága. A tekercsgeometria, áram, valós tekercselés és környezet hibáját nem foglalja magában. A rácsmaximumot nem azonosítom pontos folytonos maximumhellyel.

A mért D alapján számolt ellenállás **1788,160 Ω**, eltérése a megadott 1789 Ω-tól **−0,0470%**: a jelentés 0,05%-os összevetése helyes. Öt azonos leolvasás mintaszórása valóban nulla; ez nem jelenti a méret teljes mérési bizonytalanságának nulláját.

## 3. Az L≈3 m-es javítás elfogadható

A régi, külön teljes energiákból származó **0,01268 J/A²** eltolás már nem kerül bele a 3 m körüli külső ágba. Az új A út panelintegrált és végtelenbe menő farkat, a B út közvetlen r=L/t helyettesítést használ. A 2,6–3,01 m-es tartományban a célzott újrafuttatás a két utat körülbelül **10⁻¹⁵ relatív pontossággal** azonosnak adja.

Ezen túl a külső energiák különbségét külön **Gauss–Legendre poláris és radiális kvadratúrával, egyenletes azimutszögekkel** számoltam. Ez nem az archivált Fibonacci-szögátlag. A végtelen tartományt ugyanúgy r=L/t képezi le; nincs véges levágás vagy radiális interpoláció. A mezőfüggvény közös az archívummal, de öt ponton külön köráram-integrál ellenőrizte. A külön integrál durva/finom rácsa az itt táblázott közeli pontokon legfeljebb néhány 10⁻⁶ relatív eltérést ad; ez konvergenciateszt, nem abszolút hibabizonyítás.

| L (m) | Új archivált integrátor ΔW (J/A²) | Külön integrál (J/A²) | Relatív eltérés |
|---:|---:|---:|---:|
| 2,6 | 0,00713721108 | **0,00713720352** | +0,000106% |
| 2,818383 | 0,00369084579 | **0,00369083957** | +0,000169% |
| 3,00 | 0,00233131494 | **0,00233130964** | +0,000227% |
| 3,162 | 0,00162105726 | **0,00162105265** | +0,000284% |
| 10 | 2,77763231×10⁻⁶ | **2,77754429×10⁻⁶** | +0,00317% |
| 31,62 | 8,37811251×10⁻⁹ | **8,37717664×10⁻⁹** | +0,0112% |
| 100 | 2,63649938×10⁻¹¹ | **2,63556282×10⁻¹¹** | +0,0355% |

A korábbi több száz százalékos átmeneti és több százalékos távoli eltérések így megszűntek az ellenőrzött pontokon. A 100 m-nél látható kis maradék szögkvadratúra-finomítással tovább csökken: 6000/12000/24000 távoli szögpontnál rendre **2,6364994 / 2,6357023 / 2,6355308 ×10⁻¹¹ J/A²**. A radiális t-kvadratúra itt 32, majd 64 pontos. Nem állítok minden L-re vagy tetszőleges távolságra tanúsított hibakorlátot.

### Két numerikus minősítés marad

**A. A távoli „0 szórás” nem bizonyít nulla hibát.** A konvergenciatáblázat távoli sorában valóban azonos mind a négy érték. A négy futás azonban a magas felbontású `nhi` szögpontszámot és a belső panelezést változtatja; a távoli `nlo=6000`, `nt=32` beállítás mindegyikben ugyanaz. Emiatt ez a sor nem teszteli az adott tartomány felbontásfüggését.

**Javítás – magas bizonyosság, kisebb jelentési pontosítás:** a távoli sor „0” értékeihez ezt a megjegyzést kell hozzáadni:

> „A négy futás ebben a tartományban azonos nlo és nt beállítást használ. A nulla mintaszórás ezért nem távoli hibabecslés; a távoli konvergenciát nlo és nt külön változtatásával kell vizsgálni.”

**B. Az új, kisebb illesztési hiba közölt és valós.** A jelentés már feltünteti, hogy L=0,95 m-nél a belső/külső váltás nem pontosan folytonos. Az archivált értékek **2,17738947 és 2,17597607 J/A²**, különbségük **0,00141340 J/A²**, azaz a belső érték **0,0649%-a**. Ez a teljesenergia-becslések eltérésével algebrailag egyezik; nem a régi 3 m-es nagy hiba fennmaradása. Ezt most archivált adatokból és az integrálok azonosságából ellenőriztem, nem egy új teljesenergia-integrálással.

**Minősítés:** az L=0,95 m körüli kis hibát a görbe numerikus korlátjaként meg kell tartani, nem szabad az egész görbét gépi pontossággal folytonosnak nevezni. A dokumentum ezt lényegében már megteszi, ezért nem kezelem új, a fő következtetést megdöntő hibaként.

## 4. A háttértér kezelése most megfelelően lehatárolt

**Hely:** §4, polaritásváltás. A hozzáadott szöveg:

> A tiszta tekercs-B² tag páros a polaritásra. Állandó háttértérben az energiasűrűség-változásnak 2B_bg·B_coil kereszttagja is van; ezt az azonos nagyságú ±I állapotok OFF-hoz viszonyított átlaga ideális esetben törli. A közölt −1.20×10⁻⁴⁸ érték erre az izolált páros tagra vonatkozik, nem minden lehetséges egyszerű ON/OFF mérés teljes előrejelzése.

**Ezt a javítást elfogadom.** A gyenge jelre használt lineáris frekvenciaeltolódási közelítésben, változatlan háttértér és azonos nagyságú ±I áram mellett:

\[
\frac{1}{2}\{[|\mathbf B_{bg}+\mathbf b|^2-|\mathbf B_{bg}|^2]
+[|\mathbf B_{bg}-\mathbf b|^2-|\mathbf B_{bg}|^2]\}=|\mathbf b|^2.
\]

A kód továbbra is a tiszta tekercs-B² tagot számolja. Ez most egyezik a szövegben kifejezetten lehatárolt számítási céllal. Nem igazol tényleges háttértér-mérést, ±I szimmetriát vagy teljes egyszerű ON/OFF előrejelzést; ezekhez kísérleti adatok kellenek. A §5 első sorában is érdemes az „izolált páros tag” megnevezést használni, hogy a feltétel a végkövetkeztetésben se vesszen el.

## 5. Egy fennmaradt kísérleti túlállítás: hőkompenzáció

**Hely:** §4, 5. pont. Az azonos hőterhelésű kontroll bevezetése helyes irány, de a következő következtetés túl erős:

> A kikapcsolt szakaszokban ugyanakkora hőterhelést (2.81 W) kell leadni egy mágneses teret nem keltő (bifiláris) ellenálláson, a tekercsek közelében. Így a be/ki különbség nem tartalmaz hőmérsékleti lépést.

Azonos összteljesítmény nem garantál azonos térbeli hőleadást vagy azonos hővezetési és sugárzási csatolást az órákhoz. A beállás kivárása az időfüggő tranziens mérséklésére alkalmas, de a két beállt állapot hőmérséklete továbbra is eltérhet. A bifiláris kialakítás maradéktere is ellenőrizendő, nem feltétlenül pontosan nulla.

**Javítás – magas bizonyosság, a kísérleti bizonyító erőt érinti:** az idézet második mondatát ezzel kell felváltani:

> „Az azonos összteljesítményű, kis maradékterű kontroll célja a kapcsolással korrelált termikus változás csökkentése. A be- és kikapcsolt állapotok hőmérsékleti azonosságát és a kontroll maradékterét méréssel ellenőrizni kell; a maradék eltérést a frekvenciabecslés bizonytalanságába kell foglalni.”

Ez nem a hőkontroll elutasítása, hanem a garantált törlés helyett az ellenőrizhető maradék meghatározása.

## 6. Egyéb javítások és a mérési állítás határa

**Elfogadott javítások:**

- az offset nem zajküszöb; a főóra-egyenértékű mágneses együtthatók nem mért felső korlátok;
- a geometriai jegyzőkönyv és a perem/menetfelület különválasztása megjelent;
- az Sr-érték rendszeres bizonytalanságként, nem univerzális zajküszöbként szerepel, a rendkívüli fehér-FM idő pedig algebrai extrapolációként;
- a házhőmérséklet adatlap-összevetése feltételes lett;
- a nagyáramú számításból nem vezet le biztonságos impulzushosszt: a kért figyelmeztetés bekerült.

**Nincs új fázis- vagy frekvenciaidősor, ON/OFF/±I állapotnapló és kapcsolt termikus mérés.** A méretjegyzőkönyv nem pótolja ezeket. Emiatt a valós érzékenység továbbra sem állapítható meg. A 65 perces, valódi 2π elektromos fázisciklus átszámítása változatlanul **0,256410 mHz** és **2,56410×10⁻¹¹ átlagos relatív offset**.

A §4 kapcsolási ütemét továbbra is a beat-periódushoz köti. Ez nem önmagában érvénytelenné tevő választás, de a 65 percet nem az offset teszi optimálissá; zaj- és termikus beállásvizsgálattal kell indokolni. A szükséges következő bizonyíték az időbélyeges fázis és az állapotcímkék alapján számított be/ki becslés és bizonytalansága, nem egy újabb számított beat-arány.

**Végső minősítés:** a kért két fő javítás megtörtént; a számított ≈−1,20×10⁻⁴⁸-as tiszta páros térenergiajel számszerűen alátámasztott. A maradék numerikus és kontroll-megfogalmazási korlátokkal ez használható feltételes tervezési előrejelzés, de nem hitelesített mérési eredmény.

## Reprodukció és kísérőadat

A `revision3_checks.json` tartalmazza a 18 integritási tesztet, az öt külön tértesztet, a célzott főcsúcs-újrafuttatást, a két Sobol-pontot, kilenc külső ellenőrző integrált, a távoli felbontásváltoztatást és a közölt 0,95 m-es illesztési maradék számtanát. Ezek számított ellenőrzések; nem órás nyersadatok.

A helyi ellenőrző kód:

```bash
python3 user/workspace/co-scientist-outputs/ocxo-revision3-review/check_revision3.py
```

Python 3, NumPy és SciPy szükséges. A kód a ZIP-et olvassa, a futtatáshoz szükséges két modult scratch területre másolja, és az eredmény-JSON-t újragenerálja. Az eredeti csomagot, a kéziratot és a publikációs felületeket nem írja felül.
