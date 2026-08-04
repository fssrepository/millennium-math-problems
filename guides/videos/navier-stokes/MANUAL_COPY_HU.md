# Rövid magyar kézikönyvszöveg

Ez a blokk a hosszabb leírás vagy a kézikönyv videó előtti bevezetőjeként
használható.

## Mi a kérdés?

A háromdimenziós Navier–Stokes-egyenlet a folyadék mozgását írja le. Az örvény
egyre finomabb szerkezeteket nyújthat ki, miközben a viszkozitás kisimítja őket.
Az egymillió dolláros kérdés leegyszerűsítve: sima kezdetből marad-e mindig sima
az áramlás, vagy kialakulhat végtelenül durva pont?

## Mit csinált a kétnapos labor?

Nem a teljes problémát „próbálta ki”. Köztes lemmajelölteket támadott
ellenpéldákkal: skálázási tesztekkel, Fourier-triádokkal és egyre nagyobb,
egzakt gradiensű optimalizálásokkal. Ha egy állítás elbukott, a rendszer
mentette az akadályt, szűkítette a következő állítást, majd folytatta. Így egy
rossz út percek vagy órák alatt kieshet, mielőtt heteket vinne el.

## Mi lett az eredmény?

33 óra 28 perc alatt 16 lemmajelöltet vetettünk el, és egy részleges
far-tail-becslést bizonyítottunk. A döntő, vágásfüggetlen L4/PNT-12 lépés továbbra
is nyitott; a véges számítási teszt túlélése nem bizonyítás. A következő PNT-12
döntési kör becslése körülbelül 1–3 nap célzott munka. A teljes Clay-megoldás
ideje nem becsülhető őszintén.

## Mihez hasonlít?

Módszertani rokonság: Protasék adjunkt-alapú szélsőáramlás-keresése és Tao
kvantitatív a priori korlát-szemlélete. Ebben a laborban az adverszárius keresés
elsősorban a köztes lemmákat stresszteszteli. Ez összehasonlítás, nem állítás
arról, hogy a munka a két szerző bizonyítását kombinálná.

## Kötelező képaláírás

> A videó folyadékképei koncepcióvizualizációk. A C++ labor tényleges kimenete
> kód, mentett állapot, tanúsítvány és véges Fourier/Galerkin-adat volt.
