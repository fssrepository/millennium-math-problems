# Higgsfield — kattintásról kattintásra, alap módban

Ezt elsősorban S02, S05 vagy S06 kamerahibájának javítására használd. A normál
képből-videó generálás elég; a kamera leírása benne van a promptban.

## Amit nem használunk

- AI Director vagy automatikus storyboard;
- Cast vagy Soul ID;
- Motion Control és referencia-mozgás;
- beszélő karakter, lip sync vagy automatikus narráció;
- többjelenetes generálás.

Ha a Higgsfield a normál videópanelt `Cinema Studio` név alatt mutatja, beléphetsz
abba, de csak a **Video / Image to Video** panelt használd. A director- és
karaktereszközök nem szükségesek.

## Első beállítás

1. Nyisd meg a https://higgsfield.ai/ai/video címet és jelentkezz be.
2. Nyisd meg a **Video** vagy **Image to Video** felületet.
3. Modell: `Kling 3.0`.
4. Aspect ratio: `9:16`.
5. Audio/Sound: off.
6. Output count: `3`, ha a csomagod engedi; különben generálj háromszor azonos
   beállítással.

## Minden snittnél

1. Töltsd fel a lenti táblázatban jelzett PNG-t az image/start-frame mezőbe.
2. Nyisd meg a `PROMPTS.md` fájlt.
3. A fő promptmezőbe másold:
   - először a `Shared continuity` fenced blokkot;
   - közvetlenül utána a kiválasztott `S0X` fenced blokkot.
4. Ha van negatív promptmező, oda másold a `Shared negative` blokkot. Ha nincs,
   illeszd a fő prompt végére.
5. Állítsd be a lenti időtartamot.
6. Ha a normál panel külön kínál lens/camera mezőt, használd a táblázat értékét.
   Ha nem kínál, ne keress plugint: a prompt már tartalmazza ugyanezt.
7. Kattints Generate.
8. Töltsd le az eredményt a `../renders/higgsfield/` könyvtárba.
9. Ellenőrizd a `../PRODUCTION_WORKFLOW_HU.md` szerint.

| Snitt | Start image | Hossz | Kamera/lencse, ha külön mező van |
| --- | --- | ---: | --- |
| S01 | `../assets/keyframe-01-world.png` | 5 s | 50 mm, slow dolly in |
| S02 | `../assets/keyframe-01-world.png` | 5–6 s | 24 mm, push in / tracking |
| S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | 85 mm macro, 30° orbit |
| S04 | `../assets/keyframe-01-world.png` | 8 s | 35 mm, fast fly-through |
| S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | 100 mm macro, locked |
| S06 | `../assets/keyframe-06-frontier.png` | 8 s | 24 mm, slow dolly out |

## Javítási szabály

Egy újrapróbálkozásnál csak a hibás kamera-mondatot módosítsd. Ne adj hozzá
AI Directort, storyboardot vagy Motion Controlt. A cél tiszta, egyetlen kamera-
snitt; a pontos ritmust később a vágóprogram készíti.
