# Google Flow — kattintásról kattintásra

Ez az elsődleges útvonal. Csak a Flow normál promptmezőjét és a beépített
**Frames** képből-videó funkciót használja.

## Amit nem kell megnyitnod

- Agent;
- Storyboard Studio vagy jelenetépítő;
- Ingredients;
- karakter-, avatar- vagy hangreferencia;
- Extend vagy automatikus történetgenerálás.

## Első beállítás

1. Nyisd meg a https://flow.google címet és jelentkezz be.
2. Kattints új projekt létrehozására.
3. Projekt neve: `Navier Stokes — Lemma Stress Test — 9x16`.
4. Keresd meg az alsó vagy középső nagy promptmezőt.
5. Ha az **Agent** kapcsoló aktív, kapcsold ki.
6. Kattints a promptmezőnél látható modellnévre. Alapból ez lehet egy
   képmodell, például Nano Banana Pro.
7. A megnyíló választóban kattints: **Video**.
8. A videó bemeneti típusánál válaszd: **Frames**.
9. A `+ Add start frame` helyre húzd be az adott snitt PNG-jét.
10. Az `+ Add end frame` mezőt hagyd üresen.
11. Válaszd ki:
    - aspect ratio: `9:16`;
    - model: `Veo 3.1 Fast`;
    - outputs: `4`;
    - generation length: a lenti táblázat szerinti legközelebbi érték.
12. Hangnál válassz ambient/no speech beállítást, vagy kapcsold ki, ha van
    külön sound kapcsoló.

Ha a gombfelirat kissé eltér, a logika ugyanaz: **standard prompt → Video →
Frames → Add start frame**.

## Snittenként ezt csináld

| Sorrend | Snitt | Start frame | Hossz | Másolandó prompt |
| ---: | --- | --- | ---: | --- |
| 1 | S01 | `../assets/keyframe-01-world.png` | 5 s | `READY_TO_PASTE.md` → S01 teljes blokk |
| 2 | S02 | `../assets/keyframe-01-world.png` | 6 s | `READY_TO_PASTE.md` → S02 teljes blokk |
| 3 | S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | `READY_TO_PASTE.md` → S03 teljes blokk |
| 4 | S04 | `../assets/keyframe-01-world.png` | 8 s | `READY_TO_PASTE.md` → S04 teljes blokk |
| 5 | S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | `READY_TO_PASTE.md` → S05 teljes blokk |
| 6 | S06 | `../assets/keyframe-06-frontier.png` | 8 s | `READY_TO_PASTE.md` → S06 teljes blokk |

Példa S01-re:

1. Húzd a `keyframe-01-world.png` képet a `+ Add start frame` mezőre.
2. Nyisd meg a `READY_TO_PASTE.md` fájlt.
3. Keresd meg: `S01 — 5 s`.
4. A három backtick közötti teljes angol blokkot jelöld ki.
5. Másold a Flow fő promptmezőjébe.
6. Ellenőrizd újra: Video, Frames, 9:16, Veo 3.1 Fast, 4 outputs, 5 s.
7. Kattints **Generate**. A felület egyes nyelveken hibásan `Generate Image`
   szöveget is mutathat; a kiválasztott Video mód a fontos.
8. Várd meg mind a négy kimenetet.
9. Töltsd le őket a `../renders/flow/drafts/` könyvtárba
   `S01_flow_fast_take01.mp4`–`take04.mp4` néven.
10. Futtasd végig a `../PRODUCTION_WORKFLOW_HU.md` ellenőrzési kapuját.

Ugyanezt ismételd a táblázat következő sorával. Minden új snittnél töröld az
előző start frame-et és promptot, majd töltsd be az újat.

## Ha megvan a jó Fast változat

1. Ne változtasd meg a start frame-et és a promptot.
2. Kattints újra a modellnévre.
3. Válts `Veo 3.1 Quality` modellre.
4. Outputs: `2`.
5. A hossz és a 9:16 maradjon változatlan.
6. Kattints Generate.
7. Töltsd le a két eredményt a `../renders/flow/finals/` könyvtárba.
8. Nézd meg mindkettőt normál és 0,5× sebességen.
9. A kiválasztott fájl nevének végére írd: `_SELECTED.mp4`.

## Mikor ne lépj tovább?

- Ha a kocka deformálódik.
- Ha egynél több piros/korall fonál jelenik meg.
- Ha olvashatatlan AI-szöveg vagy képlet jelenik meg.
- Ha a kamera nem a kívánt végpontban áll meg.
- Ha nincs tiszta első és utolsó 8–12 képkocka.

Ilyenkor ne nyiss Storyboard Studiót vagy Agentet. A
`../PRODUCTION_WORKFLOW_HU.md` hibatáblája alapján módosíts egyetlen mondatot,
majd generálj újra.
