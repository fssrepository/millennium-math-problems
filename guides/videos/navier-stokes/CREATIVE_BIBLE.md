# Creative bible

## Concept

The recurring hero is not a scientist or mascot. It is one recognizable fluid
system:

- **The flow** — luminous teal streamlines inside a glass cube. It represents
  the evolving 3D velocity field.
- **The dangerous filament** — one braided coral-red vortex filament. It
  represents concentrated stretching and must keep the same material and color
  in every shot.
- **Viscosity** — cyan smoothing visible as fine ripples losing contrast. It is
  a behavior, not a character or magic force.
- **The falsifier** — one restrained amber scan plane. It appears only at a
  checkpoint and never behaves like a laser weapon.
- **The ledger** — small archived coral crystals on an amber curve. Exactly 16
  marks represent the rejected candidates F000–F015.
- **The frontier** — an open cyan gap beyond nested transparent spectral
  shells. It signals an unresolved continuation, not a win.

No human character is needed. If a platform asks for “character references,”
use `assets/keyframe-01-world.png` as the recurring subject/world reference.

## World

Dark near-black scientific space, faint spectral grid, tiny particles, glass
and refractive water. The look is a premium science documentary: physically
plausible, clean and rigorous, with a small amount of cinematic haze. It must
not drift into fantasy, underwater ocean footage, lightning, fire or creature
imagery.

## Palette

| Role | Hex | Use |
| --- | --- | --- |
| Background | `#030A10` | Near-black navy |
| Main flow | `#18C7E8` | Teal/cyan streamlines and open frontier |
| Vortex / rejected | `#FF725E` | Coral filament and failure marks |
| Checkpoint | `#F2A23A` | Amber scan and ledger line |
| Text | `#F4F8FA` | Primary captions |
| Secondary text | `#A8BBC5` | Explanatory copy |

## Typography and graphics

- Typeface: `Inter Tight SemiBold` or `Manrope SemiBold`; use one family only.
- Main title: 88–104 px; sentence case, maximum three lines.
- Mid-film cards: 50–58 px; maximum two lines.
- Subtitles: 48–54 px with a 70% black rounded backing.
- Status words may be uppercase. Do not use all-uppercase paragraphs.
- Pointer: 2 px off-white line, 8 px round endpoint, no cartoon arrowhead.
- Cards: `#07131C` at 84% opacity, 22 px radius, 1 px cyan or coral edge.

## Camera grammar

- S01: 50 mm-equivalent, slow dolly-in.
- S02: 24 mm-equivalent, fast push-through with strong parallax.
- S03/S05: 85–100 mm macro, shallow but readable depth of field.
- S04: 28–35 mm fast forward flight through nested shells.
- S06: 24–28 mm deliberately slow pullback, centered geometry.
- S07: near-still wide frame, no triumphal tilt-up.
- Use one fast movement followed by a clear arrest. Never let the whole film
  wobble continuously.

## Motion rules

- Fluid motion is continuous, heavy and coherent; not smoke, cloth or hair.
- The cube remains rigid and intact.
- The coral filament may stretch, fold and thin, but must never split into
  many unrelated worms.
- Hard freezes, frame holds, stamps, counters and the three restart beats are
  created in post.
- Camera energy decreases across S04 → S05 → S06 to make the increasing search
  depth legible while the T+ research clock remains linear.
- Generate no text inside the video model. Current image/video systems may
  misspell copy and drift labels between frames.

## Master visual prompt

Paste this above any shot prompt when the platform does not accept an image
reference:

```text
Premium photorealistic scientific visualization in a dark near-black research
space: one transparent glass cube containing a physically plausible
three-dimensional incompressible fluid flow made of coherent luminous
turquoise-blue streamlines, and one recognizable thin braided coral-red vortex
filament being stretched and folded. Cool cyan rim light, restrained warm-amber
checkpoint accents, subtle spectral grid and particles, physically based 3D
render, science-documentary rigor, strong phone-readable depth, 9:16 portrait.
Preserve the exact glass material, teal flow, coral filament identity, palette
and lighting across every shot.
```

## Shared negative prompt

```text
No text, letters, numbers, equations, captions, logos, watermark, UI panels,
people, hands, scientist, ocean wave, smoke-only fluid, fire, lightning,
tentacles, creature, duplicated cube, broken glass, impossible splash,
explosion, trophy, checkmark, victory symbolism, fantasy magic, camera shake,
flicker, morphing materials, random color changes, extra red filaments.
```

## Scientific guardrail

The visuals are explanatory metaphors. Add `CONCEPT VISUALIZATION` as a small
editor overlay from S02 onward. Never call the generated fluid plate “the
simulation result.” The actual run produced code, certificates, failed-lemma
records and finite Fourier/Galerkin evidence.
