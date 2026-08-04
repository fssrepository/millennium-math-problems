# EZT NYISD MEG ELŐSZÖR

Az egész munkafolyamat alapváltozata Google Flow-ban készül. Nem kell hozzá
plugin, Agent vagy Storyboard Studio.

Nyisd meg: https://flow.google

## Mit hol találsz?

| Amit keresel | Fájl / könyvtár |
| --- | --- |
| Feltöltendő kezdőképek | `assets/` |
| Pontos Flow-gombok | `flow/RUNBOOK_HU.md` |
| Egyben bemásolható Flow-promptok | `flow/READY_TO_PASTE.md` |
| Mit ellenőrizz minden generálás után | `PRODUCTION_WORKFLOW_HU.md` |
| Hova írd, melyik take jó | `REVIEW_SHEET.csv` |
| Klipek pontos sorrendje és hossza | `shot-manifest.csv` |
| Gyorsítások és megállítások | `EDIT_PLAN.md` |
| Képernyőre kerülő magyar szöveg | `EDITOR_OVERLAYS_HU.md` |
| Feliratfájl | `captions/navier-stokes-hu.srt` |
| Narráció és TTS-stílus | `narration/hu.md` |
| Letöltött klipek helye | `renders/` |

## A teljes folyamat hét mondatban

1. Nyisd meg a `flow/RUNBOOK_HU.md` fájlt és a Google Flow weboldalát.
2. A Flow-ban válaszd: normál promptmező → Video → Frames → 9:16.
3. Töltsd fel az adott snittnél megnevezett képet az `assets/` könyvtárból.
4. Másold be egyben a megfelelő blokkot a `flow/READY_TO_PASTE.md` fájlból.
5. Generálj négy próbát, töltsd le őket, és ellenőrizd a
   `PRODUCTION_WORKFLOW_HU.md` listájával.
6. A jó snittből készíts két Quality változatot, majd ugyanezt ismételd S01-től
   S06-ig.
7. A hat elfogadott klipet rakd össze az `EDIT_PLAN.md` szerint, majd add rá az
   `EDITOR_OVERLAYS_HU.md` szövegeit, az SRT-t és a narrációt.

Ha egy kamera két javítás után sem sikerül, csak azt az egy snittet vidd át a
`higgsfield/RUNBOOK_HU.md` folyamatába. Az OpenArt teljes alternatív útvonala az
`openart/RUNBOOK_HU.md` fájlban van.

Alternatív felületek:

- OpenArt video: https://openart.ai/ai-video-generator/
- Higgsfield video: https://higgsfield.ai/ai/video
