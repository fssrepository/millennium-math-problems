# OpenArt — kattintásról kattintásra, alap módban

Ez alternatív teljes gyártási útvonal. A normál **Image to Video** felületet
használja.

## Amit nem használunk

- One Click Story;
- Smart Shot;
- OpenArt Director;
- Motion Control és mozgásreferencia-videó;
- AI Character, lip sync vagy beépített narrátor.

## Első beállítás

1. Nyisd meg a https://openart.ai/ai-video-generator/ címet és jelentkezz be.
2. Nyisd meg a **Create Video** felületet.
3. Válaszd az **Image to Video** módot. Ha csak `Start Frame` felirat látszik,
   az is megfelelő.
4. Modell: a normál `Kling 3.0` image-to-video változat. Ne a `Motion Control`
   változatot válaszd, mert ehhez a filmhez nincs emberi mozgásreferencia.
5. Aspect ratio: `9:16`.
6. Output count: `4` a próbákhoz.
7. Audio/voice: kikapcsolva vagy ambient only.

## Minden snittnél

1. Töltsd fel a start képet a lenti táblázat alapján.
2. Nyisd meg a `PROMPTS.md` fájlt.
3. A fő **Prompt** mezőbe ebben a sorrendben másold:
   - `Shared continuity` teljes fenced blokk;
   - a kiválasztott `S0X` fenced blokk.
4. Ha van külön **Negative prompt** mező, oda másold a `Shared negative`
   blokkot. Ha nincs, illeszd a fő prompt végére.
5. Állítsd be a snitthosszt.
6. Kattints **Generate**.
7. Töltsd le mind a négy eredményt a `../renders/openart/` könyvtárba.
8. Ellenőrizd őket a `../PRODUCTION_WORKFLOW_HU.md` szerint.

| Snitt | Start image | Hossz |
| --- | --- | ---: |
| S01 | `../assets/keyframe-01-world.png` | 5 s |
| S02 | `../assets/keyframe-01-world.png` | 6 s vagy a legközelebbi elérhető |
| S03 | `../assets/keyframe-03-checkpoint.png` | 5 s |
| S04 | `../assets/keyframe-01-world.png` | 8–10 s |
| S05 | `../assets/keyframe-03-checkpoint.png` | 8–10 s |
| S06 | `../assets/keyframe-06-frontier.png` | 8–10 s |

## Jó változat véglegesítése

1. Ha a felület mutat seedet, jegyezd fel a `../REVIEW_SHEET.csv` fájlba.
2. Tartsd ugyanazt a képet, promptot és seedet.
3. Állítsd a kimenetszámot 2-re és válaszd a jobb minőséget, ha van külön
   quality kapcsoló.
4. Generálj újra.
5. A győztes fájl végére írd `_SELECTED.mp4`.

Ha két célzott újrapróbálkozás után sem stabil a kamera, ne nyiss Smart Shotot;
vidd át csak ezt a snittet a Higgsfield alap image-to-video útvonalára.
