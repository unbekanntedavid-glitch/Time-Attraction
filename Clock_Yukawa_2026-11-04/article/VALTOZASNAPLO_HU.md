# A közlemény változásnaplója (2026-11-04)

Új közlemény, első változat ebben a készletben. A számítási modulok nem változtak: a három r4 archívum bájtra azonosan
szerepel (SHA-256 a `../SHA256SUMS.txt` fájlban).

1. A fejezetterv írás előtt rögzítve: 1 bevezetés és kategóriák; 2 hatás, anyagcsatolás, jelölések; 3 óra és gyorsulás
   egyetlen forrásból (egy potenciál, laboratóriumi kalibráció, megfigyelhetők, deklarált geometria); 4 numerikus módszer,
   ellenőrzések, reprodukálhatóság; 5 eredmények P1–P3-ra (S megvalósíthatósági mutatóval); 6 Yukawa-korlátok (BREN-lánc,
   szub-cm lánc, arányok); 7 következtetés; függelékek: A BREN-operátor, B digitalizálás és korlátkonstrukciók,
   C nem univerzális kiegészítés; adatok és kód elérhetősége; hivatkozások.
2. Állítás–forrás rendszer: a `src/claims_spec.py` minden számhoz megadja az archívumot, a fájlt, a kulcsot vagy sort és a
   képletet; a `src/make_values.py` írja az EN/HU értékmakrókat, a táblázatokat és a `CLAIMS_SOURCES.csv/.md` fájlt (295 állítás).
3. A szöveg csak értékmakrókat használ; a táblázatok ugyanezekből a makrókból készülnek.
4. Az ábrák a modulábrák PNG-másolatai (nincs szövegréteg, így az EN/HU PDF-számösszevetést nem zavarják a feliratok);
   az óramodul ábrája mindkét változatban angol feliratú, ezt a magyar ábraaláírás jelzi.
5. Írás közben: a számjegykereső által talált kézzel beírt számokat értékmakróra cseréltük (δy GR-hez viszonyított
   nagysága, a konvergencia-konfigurációk csonkolási sugarai, a kvadratúrarács 1 cm-es alsó határa, most 10^quadLo m
   alakban); az „R = 10 nem 10σ” mondatot számok nélkül fogalmaztuk át.
6. A nem használt sorok kikerültek az állítás–forrás táblázatból (P2/P3 nyers oszlopai, két mértékegység nélküli
   operátor-önellenőrzés, a három PASS-jelző; a PASS-eredményeket a `check_claims.py` 10. ellenőrzése közvetlenül a naplókban nézi).
7. Dekádonkénti pontszám: a kvadratúrarács és a csak korlátos rács saját értékével szerepel.
8. A DOI-k törhetetlen dobozban, így DOI-n belüli sortörés nem változtathatja meg a PDF-számösszevetést.
9. A hivatkozások a modulcsomagokban megnevezett forrásokra, a DOI-k a négy TAL-rekordra korlátozódnak.
10. Az ellenőrzési mód magyar megfogalmazása külön mondatban, a rögzített kifejezéssel.
