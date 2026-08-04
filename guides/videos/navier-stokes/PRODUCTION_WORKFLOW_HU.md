# Manuális videó-AI gyártási munkafolyamat

Ez nem one-shot generálás. Hat külön tiszta klipet készítesz, mindegyiket
ellenőrzöd, csak az elfogadott snitt után lépsz tovább, majd a ritmust és a
magyar információs réteget vágóprogramban építed rá.

## Ajánlott útvonal

```text
Kulcsképek
  ↓
Flow / Veo 3.1 Fast: 4 próba snittenként
  ↓ vizuális ellenőrzés
Flow / Veo 3.1 Quality: 2 döntő változat
  ↓ csak ha a kamera hibás
Higgsfield / Kling 3.0: problémás snitt újragenerálása
  ↓
Vágás + időrámpák + megállítások
  ↓
Magyar overlay + narráció + SRT
  ↓
Telefonos végellenőrzés
```

Az OpenArt külön, teljes alternatív útvonal és gyors kompozíciós labor. Nem kell
ugyanazt a snittet mindhárom szolgáltatásban legenerálnod.

## Nincs szükséges plugin

A teljes film elkészíthető a Flow normál promptmezőjével, a Video → Frames
funkcióval és egy hagyományos vágóval. Nem használunk Agentet, Storyboard
Studiót, Smart Shotot, AI Directort, Motion Controlt, Castot, Soul ID-t vagy
külső plugint. Ezek említése a platform runbookokban kizárólag azt jelzi, hogy
ne válaszd őket.

## 0. Előkészítés

Nyisd meg egyszerre:

- ezt a fájlt;
- a választott platform `RUNBOOK_HU.md` fájlját;
- a platform `PROMPTS.md` vagy `READY_TO_PASTE.md` fájlját;
- az `assets/` könyvtárat;
- a `REVIEW_SHEET.csv` táblát.

Minden AI-kimenetet tölts le azonnal. A javasolt elnevezés:

```text
S01_flow_fast_take01.mp4
S01_flow_quality_take01.mp4
S02_higgsfield_kling_take03.mp4
```

Tedd őket a `renders/` megfelelő alkönyvtárába. A szolgáltatás projektneve vagy
előnézete nem helyettesíti a helyi letöltést.

## 1. Snittciklus — ezt ismételd S01–S06 esetén

### 1.1 Bemenet

1. Ellenőrizd a snitthez rendelt start frame-et a platform runbookjában.
2. Állítsd be a 9:16 képarányt és a megadott időtartamot.
3. Másold be a snitt teljes promptját a megjelölt mezőbe.
4. Negatív promptmező esetén a negatív blokk oda kerüljön; ha nincs ilyen
   mező, maradjon a fő prompt végén.
5. Magyar feliratot, képletet vagy UI-t ne kérj a videogenerátortól.

### 1.2 Próbagenerálás

1. Készíts 3–4 változatot ugyanazzal a bemenettel.
2. Ne változtass egyszerre modellen, prompton és képen: különben nem tudod,
   melyik módosítás javított.
3. Töltsd le a változatokat, majd nézd meg őket egyszer normál sebességgel és
   egyszer 0,5× sebességgel.

### 1.3 Kötelező ellenőrzési kapu

Csak akkor fogadd el a snittet, ha minden rá vonatkozó pont teljesül:

- a kocka merev, nem olvad és nem törik;
- pontosan egy korall örvényfonál marad;
- a türkiz áramlás folyadéknak, nem füstnek vagy textilnek látszik;
- nincs generált szöveg, vízjel, ember vagy pszeudo-képlet;
- a kamera a kért végkompozícióban áll meg;
- az elején és végén van legalább 8–12 tiszta vágási képkocka;
- a következő snitthez szükséges tárgy és képkivágás megmarad.

Írd az eredményt a `REVIEW_SHEET.csv` másolatába. Ha egy feltétel hibás, a
snitt nem mehet a Quality körbe.

### 1.4 Egyetlen hibát javíts

| Hiba | Következő lépés |
| --- | --- |
| A korall fonál megsokszorozódik | Ugyanaz a start frame; erősítsd: `one single continuous coral filament, never duplicate` |
| A kocka deformálódik | Erősítsd: `rigid unchanged glass cube`; csökkentsd a kamera sebességét |
| Füstnek látszik | Erősítsd: `coherent heavy liquid, not smoke, not vapor` |
| A kamera célt téveszt | Ugyanaz a modell és kép; csak a kamera mondatát pontosítsd |
| Túl nagy a mozgás | Lassabb kamerautasítás vagy Higgsfield-retake |
| Megjelenik valamilyen felirat | Add hozzá a negatív blokk elejéhez: `absolutely no visible text or symbols` |
| S06-ban változik a 16 jel | Ne bízd a számlálást a mozgó modellre: használd a generált mozgást csak a pullbackhez, majd vágj a fix `keyframe-06-frontier.png` képre |

Legfeljebb két célzott újrapróbálkozást csinálj ugyanazon a modellen. Ha a
kamera a harmadik körben sem stabil, vidd át a snittet Higgsfield/Klingre.

### 1.5 Final generálás

1. Tartsd meg az elfogadott start frame-et, promptot és – ha elérhető – seedet.
2. Válts a platform jobb minőségű módjára.
3. Készíts két döntő változatot.
4. Ismételd meg a teljes ellenőrzési kaput.
5. Válaszd ki a győztes fájlt, és csak ezután lépj a következő snittre.

## 2. Snittenkénti külön ellenőrzés

| Snitt | A legfontosabb elfogadási feltétel |
| --- | --- |
| S01 | A felső 38% nyugodt és sötét; a főcím olvasható lesz rajta |
| S02 | A kamera áttöri a képsíkot, de az üveg nem törik; a korall kanyar center-right helyen megáll |
| S03 | A makró kanyar végig fókuszban marad; nincs ember vagy portré |
| S04 | Ez a leggyorsabb kamera; négy külön tesztpulzus és stabil héjgeometria látszik |
| S05 | Ez lassabb és pontosabb; az amber scan valóban keresztezi a korall fonalat |
| S06 | Ez a leglassabb; nyitott türkiz út és hosszú stabil utolsó képkocka marad |

## 3. Vágási folyamat

1. Rakd be az elfogadott klipeket a `shot-manifest.csv` sorrendjében.
2. Vágd őket pontosan 29,70 másodpercre.
3. Alkalmazd az `EDIT_PLAN.md` sebesség- és megállítási pontjait.
4. Az ellenpéldák, `ELVETVE` bélyegek, pointerek és T+ óra mind szerkesztői
   grafikák; ne próbáld a generált klipbe újrageneráltatni őket.
5. S06 végén vágj vagy oldj át a fix `keyframe-06-frontier.png` képre, hogy a
   16 jel pontos maradjon.
6. S07 ehhez a fix képhez adott 4,5 másodperces kitartás, 72%-os sötét fátyollal.

## 4. Szöveg és hang

1. Másold a képernyőszöveget az `EDITOR_OVERLAYS_HU.md` fájlból.
2. Importáld a `captions/navier-stokes-hu.srt` fájlt.
3. Vedd fel vagy szintetizáld a `narration/hu.md` szövegét.
4. Célozd a 29,2–29,7 másodperces narrációt.
5. S07 alatt kapcsold ki a normál feliratot, mert a teljes oldalas kártya
   tartalmazza az eredményt és a becslést.

## 5. Végső elfogadás

Futtasd végig a `DELIVERY_QC.md` listát, majd nézd meg:

1. telefonon, hanggal;
2. telefonon, némítva;
3. 0,75× sebességen;
4. egyszer úgy, hogy csak a feliratokat olvasod.

Csak ezután exportáld a mastert és a social változatot.
