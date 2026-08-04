#!/usr/bin/env python3
"""Build the public Millennium Math Problems overview and Navier--Stokes brief.

The visual source is self-contained and intentionally follows the image-backed
presentation pattern used by the MyScoutee public pitches: Pillow renders each
page, ReportLab creates the release PDF, and headless LibreOffice packages the
same pages into a predictable cross-platform PPTX.

Run with the system Python so the LibreOffice UNO bridge is available:

    /usr/bin/python3 guides/pitches/build.py
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Callable, Sequence

from PIL import Image, ImageDraw, ImageFilter, ImageFont
from reportlab.pdfgen import canvas as pdf_canvas


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ASSETS = HERE / "assets"
OVERVIEW_OUT = HERE / "output"
NAVIER_OUT = HERE / "output"

VERSION = "1.0.0"
LANGUAGE = "EN"
RELEASE_DATE = "AUG 2026"
W, H = 1920, 1080

OVERVIEW_STEM = f"Millennium_Math_Problems_Project_Overview_v{VERSION}_{LANGUAGE}"
NAVIER_STEM = f"Navier_Stokes_Brief_v{VERSION}_{LANGUAGE}"

FONT_REG = "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf"
FONT_MED = "/usr/share/fonts/truetype/ubuntu/Ubuntu-M.ttf"
FONT_BOLD = "/usr/share/fonts/truetype/ubuntu/Ubuntu-B.ttf"
FONT_MONO = "/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf"
FONT_MONO_BOLD = "/usr/share/fonts/truetype/ubuntu/UbuntuMono-B.ttf"
FONT_MATH = "/usr/share/fonts/truetype/dejavu/DejaVuMathTeXGyre.ttf"

BG = "#07101C"
BG_2 = "#10213A"
CARD = "#12253E"
CARD_2 = "#17304C"
TEXT = "#F4F8FC"
MUTED = "#B6C5D6"
MUTED_2 = "#7E96AD"
CYAN = "#38D6E8"
BLUE = "#4E83FF"
TEAL = "#2DD4A8"
AMBER = "#F6C85F"
ORANGE = "#FF8C42"
RED = "#FF6B78"
PURPLE = "#9A7BFF"
WHITE = "#FFFFFF"
INK = "#08121E"
LINE = "#2B4865"

RESAMPLE_LANCZOS = getattr(Image, "Resampling", Image).LANCZOS


def rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[index : index + 2], 16) for index in (0, 2, 4))


def rgba(value: str, alpha: int = 255) -> tuple[int, int, int, int]:
    return (*rgb(value), alpha)


def font(size: int, weight: str = "regular", *, mono: bool = False) -> ImageFont.FreeTypeFont:
    if mono:
        path = FONT_MONO_BOLD if weight == "bold" else FONT_MONO
    elif weight == "bold":
        path = FONT_BOLD
    elif weight == "medium":
        path = FONT_MED
    else:
        path = FONT_REG
    return ImageFont.truetype(path, size=size)


def wrap_lines(value: str, fnt: ImageFont.FreeTypeFont, max_width: int) -> list[str]:
    lines: list[str] = []
    for paragraph in value.split("\n"):
        if not paragraph:
            lines.append("")
            continue
        words = paragraph.split()
        line = words[0]
        for word in words[1:]:
            candidate = f"{line} {word}"
            if fnt.getlength(candidate) <= max_width:
                line = candidate
            else:
                lines.append(line)
                line = word
        lines.append(line)
    return lines


def text(
    image: Image.Image,
    xy: tuple[int, int],
    value: str,
    *,
    size: int,
    color: str = TEXT,
    weight: str = "regular",
    max_width: int | None = None,
    spacing: int | None = None,
    mono: bool = False,
    anchor: str | None = None,
    align: str = "left",
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(image)
    fnt = font(size, weight, mono=mono)
    rendered = value
    if max_width is not None:
        rendered = "\n".join(wrap_lines(value, fnt, max_width))
    if spacing is None:
        spacing = max(4, round(size * 0.22))
    draw.multiline_text(
        xy,
        rendered,
        font=fnt,
        fill=rgb(color),
        spacing=spacing,
        anchor=anchor,
        align=align,
    )
    return draw.multiline_textbbox(
        xy,
        rendered,
        font=fnt,
        spacing=spacing,
        anchor=anchor,
        align=align,
    )


def vertical_gradient(top: str = BG, bottom: str = BG_2) -> Image.Image:
    image = Image.new("RGBA", (W, H), rgba(top))
    draw = ImageDraw.Draw(image)
    top_rgb, bottom_rgb = rgb(top), rgb(bottom)
    for y in range(H):
        ratio = y / max(1, H - 1)
        color = tuple(round(top_rgb[i] * (1 - ratio) + bottom_rgb[i] * ratio) for i in range(3))
        draw.line((0, y, W, y), fill=(*color, 255))
    return image


def add_glow(
    image: Image.Image,
    box: tuple[int, int, int, int],
    color: str,
    *,
    blur: int = 140,
    alpha: int = 70,
) -> None:
    layer = Image.new("RGBA", image.size, (0, 0, 0, 0))
    ImageDraw.Draw(layer).ellipse(box, fill=rgba(color, alpha))
    image.alpha_composite(layer.filter(ImageFilter.GaussianBlur(blur)))


def rounded_card(
    image: Image.Image,
    box: tuple[int, int, int, int],
    *,
    fill: str = CARD,
    outline: str | None = LINE,
    radius: int = 24,
    shadow: bool = True,
    alpha: int = 248,
) -> None:
    x1, y1, x2, y2 = box
    if shadow:
        layer = Image.new("RGBA", image.size, (0, 0, 0, 0))
        ImageDraw.Draw(layer).rounded_rectangle(
            (x1 + 4, y1 + 10, x2 + 4, y2 + 10), radius=radius, fill=(0, 0, 0, 100)
        )
        image.alpha_composite(layer.filter(ImageFilter.GaussianBlur(16)))
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    ImageDraw.Draw(overlay).rounded_rectangle(
        box,
        radius=radius,
        fill=rgba(fill, alpha),
        outline=rgba(outline, 230) if outline else None,
        width=2,
    )
    image.alpha_composite(overlay)


def pill(
    image: Image.Image,
    xy: tuple[int, int],
    label: str,
    *,
    fill: str,
    color: str = WHITE,
    size: int = 15,
    height: int = 34,
    pad_x: int = 15,
    outline: str | None = None,
) -> tuple[int, int, int, int]:
    fnt = font(size, "bold")
    width = math.ceil(fnt.getlength(label)) + 2 * pad_x
    x, y = xy
    ImageDraw.Draw(image).rounded_rectangle(
        (x, y, x + width, y + height),
        radius=height // 2,
        fill=rgb(fill),
        outline=rgb(outline) if outline else None,
        width=2,
    )
    ImageDraw.Draw(image).text(
        (x + width / 2, y + height / 2 - 1), label, font=fnt, fill=rgb(color), anchor="mm"
    )
    return (x, y, x + width, y + height)


def section_header(
    image: Image.Image,
    index: str,
    category: str,
    title_value: str,
    subtitle: str,
    *,
    status: str | None = None,
) -> None:
    text(image, (88, 48), f"{index}  /  {category.upper()}", size=18, color=TEAL, weight="bold")
    if status:
        pill(image, (1605, 45), status, fill=TEAL, color=INK)
    text(image, (88, 88), title_value, size=48, weight="bold", max_width=1600)
    text(image, (90, 158), subtitle, size=21, color=MUTED, max_width=1690)


def footer(image: Image.Image, page: int, total: int, label: str, *, note: str | None = None) -> None:
    draw = ImageDraw.Draw(image)
    draw.line((88, 1031, 1832, 1031), fill=rgb(LINE), width=2)
    text(
        image,
        (88, 1043),
        f"MILLENNIUM MATH PROBLEMS  •  {label.upper()}  •  V{VERSION}  •  {RELEASE_DATE}",
        size=13,
        color=MUTED_2,
        weight="medium",
    )
    if note:
        text(image, (960, 1043), note, size=13, color=MUTED_2, anchor="ma", align="center")
    text(image, (1832, 1043), f"{page} / {total}", size=13, color=MUTED_2, weight="bold", anchor="ra")


def paste_cover(
    image: Image.Image,
    source: Path,
    box: tuple[int, int, int, int],
    *,
    radius: int = 28,
    darken: int = 0,
    crop: str = "cover",
    outline: str | None = LINE,
) -> None:
    x1, y1, x2, y2 = box
    target_w, target_h = x2 - x1, y2 - y1
    original = Image.open(source).convert("RGBA")
    scale = min(target_w / original.width, target_h / original.height) if crop == "contain" else max(
        target_w / original.width, target_h / original.height
    )
    resized = original.resize(
        (max(1, round(original.width * scale)), max(1, round(original.height * scale))),
        RESAMPLE_LANCZOS,
    )
    canvas = Image.new("RGBA", (target_w, target_h), (0, 0, 0, 0))
    canvas.alpha_composite(resized, ((target_w - resized.width) // 2, (target_h - resized.height) // 2))
    if darken:
        canvas.alpha_composite(Image.new("RGBA", canvas.size, (0, 0, 0, darken)))
    mask = Image.new("L", (target_w, target_h), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, target_w, target_h), radius=radius, fill=255)
    image.paste(canvas, (x1, y1), mask)
    if outline:
        ImageDraw.Draw(image).rounded_rectangle(box, radius=radius, outline=rgb(outline), width=2)


def arrow(
    image: Image.Image,
    start: tuple[int, int],
    end: tuple[int, int],
    *,
    color: str = CYAN,
    width: int = 5,
    head: int = 14,
) -> None:
    draw = ImageDraw.Draw(image)
    draw.line((*start, *end), fill=rgb(color), width=width)
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    left = (
        end[0] - head * math.cos(angle) + head * 0.65 * math.sin(angle),
        end[1] - head * math.sin(angle) - head * 0.65 * math.cos(angle),
    )
    right = (
        end[0] - head * math.cos(angle) - head * 0.65 * math.sin(angle),
        end[1] - head * math.sin(angle) + head * 0.65 * math.cos(angle),
    )
    draw.polygon((end, left, right), fill=rgb(color))


def bullet_list(
    image: Image.Image,
    items: Sequence[str],
    *,
    x: int,
    y: int,
    width: int,
    size: int = 18,
    gap: int = 12,
    bullet_color: str = TEAL,
) -> int:
    cursor = y
    for item in items:
        ImageDraw.Draw(image).ellipse((x, cursor + 8, x + 9, cursor + 17), fill=rgb(bullet_color))
        bounds = text(image, (x + 24, cursor), item, size=size, color=TEXT, max_width=width - 24)
        cursor = bounds[3] + gap
    return cursor


def simple_table(
    image: Image.Image,
    box: tuple[int, int, int, int],
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    widths: Sequence[int],
    *,
    row_height: int,
    header_height: int = 58,
    body_size: int = 16,
    highlight_rows: set[int] | None = None,
    status_column: int | None = None,
) -> None:
    x1, y1, x2, y2 = box
    if sum(widths) != x2 - x1:
        raise ValueError("table widths must equal table width")
    highlight_rows = highlight_rows or set()
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle(box, radius=20, fill=rgb(CARD), outline=rgb(LINE), width=2)
    draw.rounded_rectangle((x1, y1, x2, y1 + header_height), radius=20, fill=rgb("#1B3855"))
    draw.rectangle((x1, y1 + header_height - 20, x2, y1 + header_height), fill=rgb("#1B3855"))
    cursor_x = x1
    for index, (header, width) in enumerate(zip(headers, widths)):
        text(image, (cursor_x + 16, y1 + 18), header.upper(), size=14, color=CYAN, weight="bold", max_width=width - 32)
        cursor_x += width
        if index < len(headers) - 1:
            draw.line((cursor_x, y1, cursor_x, y2), fill=rgb(LINE), width=2)
    cursor_y = y1 + header_height
    status_palette = {
        "ACTIVE": (TEAL, INK),
        "PLANNED": (BLUE, WHITE),
        "NOT STARTED": (MUTED_2, WHITE),
        "SOLVED EXTERNALLY": (AMBER, INK),
        "PROVED": (TEAL, INK),
        "OPEN": (ORANGE, INK),
        "FALSIFIED": (RED, WHITE),
        "VERIFIED": (CYAN, INK),
    }
    for row_index, row in enumerate(rows):
        row_top = cursor_y + row_index * row_height
        row_bottom = row_top + row_height
        if row_index in highlight_rows:
            draw.rectangle((x1 + 2, row_top, x2 - 2, row_bottom), fill=rgb("#163B49"))
        if row_index > 0:
            draw.line((x1, row_top, x2, row_top), fill=rgb(LINE), width=1)
        cursor_x = x1
        for column_index, (cell, width) in enumerate(zip(row, widths)):
            if status_column == column_index and cell in status_palette:
                fill, color = status_palette[cell]
                pill(
                    image,
                    (cursor_x + 14, row_top + (row_height - 30) // 2),
                    cell,
                    fill=fill,
                    color=color,
                    size=12,
                    height=30,
                    pad_x=11,
                )
            else:
                text(
                    image,
                    (cursor_x + 16, row_top + 14),
                    cell,
                    size=body_size,
                    color=TEXT if column_index == 0 else MUTED,
                    weight="bold" if column_index == 0 else "regular",
                    max_width=width - 32,
                    spacing=4,
                )
            cursor_x += width


def overview_slide_1() -> Image.Image:
    image = vertical_gradient("#040B14", "#10213A")
    paste_cover(image, ASSETS / "navier-stokes-flow.png", (790, 0, 1920, 1080), radius=0, darken=58, outline=None)
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    gradient = ImageDraw.Draw(overlay)
    for x in range(1100):
        alpha = max(0, 255 - round(x / 1100 * 255))
        gradient.line((x, 0, x, H), fill=(4, 11, 20, alpha))
    image.alpha_composite(overlay)
    add_glow(image, (-300, -150, 850, 850), BLUE, blur=185, alpha=58)

    pill(image, (92, 66), "PUBLIC RESEARCH GUIDE", fill=TEAL, color=INK, size=15)
    pill(image, (360, 66), "HUMAN-DIRECTED • AI-EXECUTED", fill=PURPLE, color=WHITE, size=15)
    text(image, (92, 150), "Fast falsification.\nShorter research loops.", size=64, weight="bold", max_width=770)
    text(
        image,
        (96, 330),
        "Under human direction, Codex turns proof gaps into C++ tests, attacks weak lemma routes, reads the failures and carries the surviving structure into the next candidate.",
        size=25,
        color=MUTED,
        max_width=690,
        spacing=8,
    )
    rounded_card(image, (92, 520, 720, 775), fill="#112B3D", outline=TEAL, radius=26)
    text(image, (122, 548), "PRIMARY ADVANTAGE", size=15, color=TEAL, weight="bold")
    text(image, (122, 590), "Bad directions fail quickly.", size=31, weight="bold")
    bullet_list(
        image,
        [
            "Exact numerical objectives stress a proposed lemma.",
            "Counterexamples and failure certificates become reusable data.",
            "AI uses each failure to refine the next proof obligation.",
        ],
        x=124,
        y=646,
        width=550,
        size=18,
        gap=9,
        bullet_color=CYAN,
    )
    text(image, (96, 832), "CURRENT ACTIVE LAB", size=14, color=AMBER, weight="bold")
    text(image, (96, 866), "3D Navier–Stokes existence and smoothness", size=28, weight="bold", max_width=690)
    text(image, (96, 922), "Research code and finite evidence — not a Clay solution.", size=18, color=MUTED)
    footer(image, 1, 5, "project overview", note="AI executes the loop; experts still verify the mathematics")
    return image


def overview_slide_2() -> Image.Image:
    image = vertical_gradient("#06111D", "#10253A")
    add_glow(image, (1200, -250, 2200, 700), PURPLE, blur=190, alpha=40)
    section_header(
        image,
        "01",
        "Problem portfolio",
        "Seven problems, one active laboratory",
        "The repository is organized by Millennium Prize Problem. Navier–Stokes is active; the other unsolved tracks are explicit future scope.",
    )
    headers = ("Problem", "Area", "Clay status", "Repository")
    rows = (
        ("Birch–Swinnerton-Dyer", "Number theory", "NOT STARTED", "Planned"),
        ("Hodge conjecture", "Algebraic geometry", "NOT STARTED", "Planned"),
        ("Navier–Stokes", "PDE / fluid dynamics", "ACTIVE", "C++ lemma lab"),
        ("P vs NP", "Complexity theory", "NOT STARTED", "Planned"),
        ("Poincaré conjecture", "Topology", "SOLVED EXTERNALLY", "No lab needed"),
        ("Riemann hypothesis", "Number theory", "NOT STARTED", "Planned"),
        ("Yang–Mills mass gap", "Mathematical physics", "NOT STARTED", "Planned"),
    )
    simple_table(
        image,
        (88, 235, 1832, 929),
        headers,
        rows,
        (545, 420, 375, 404),
        row_height=90,
        header_height=64,
        body_size=17,
        highlight_rows={2},
        status_column=2,
    )
    text(
        image,
        (92, 958),
        "Official problem definitions and current solved/unsolved classification: claymath.org/millennium-problems",
        size=16,
        color=MUTED_2,
    )
    footer(image, 2, 5, "project overview", note="Only Navier–Stokes currently contains implementation and certificates")
    return image


def flow_box(
    image: Image.Image,
    box: tuple[int, int, int, int],
    label: str,
    title_value: str,
    body: str,
    *,
    accent: str,
    status: str,
) -> None:
    rounded_card(image, box, fill=CARD, outline=accent, radius=22, shadow=False)
    x1, y1, x2, _ = box
    text(image, (x1 + 22, y1 + 18), label, size=13, color=accent, weight="bold")
    text(image, (x1 + 22, y1 + 52), title_value, size=22, weight="bold", max_width=x2 - x1 - 44)
    text(image, (x1 + 22, y1 + 96), body, size=16, color=MUTED, max_width=x2 - x1 - 44, spacing=5)
    pill(
        image,
        (x1 + 22, box[3] - 45),
        status,
        fill=TEAL if status == "IMPLEMENTED" else BLUE,
        color=INK if status == "IMPLEMENTED" else WHITE,
        size=11,
        height=27,
        pad_x=10,
    )


def overview_slide_3() -> Image.Image:
    image = vertical_gradient("#07101C", "#14213D")
    add_glow(image, (-220, 220, 720, 1100), CYAN, blur=190, alpha=42)
    section_header(
        image,
        "02",
        "Lemma engine",
        "AI executes the research loop",
        "Under human direction, Codex formulates candidates, writes and runs the C++, interprets failures and designs the next lemma. Every claim returns to exact tests and expert proof review.",
        status="RUN IN CAMPAIGN",
    )

    boxes = [
        ((88, 255, 386, 520), "INPUT 01", "Proof obligation", "A precise inequality, scaling law or closure target from the proof plan.", CYAN, "IMPLEMENTED"),
        ((448, 255, 746, 520), "INPUT 02", "Empirical artifacts", "TSV states, JSON certificates, gradient traces and failed-lemma records.", BLUE, "IMPLEMENTED"),
        ((808, 255, 1106, 520), "ENGINE 03", "Adversarial search", "Exact objectives, JVP/VJP checks, L-BFGS and multistart stress the statement.", PURPLE, "IMPLEMENTED"),
        ((1168, 255, 1466, 520), "SYNTHESIS 04", "AI pattern analysis", "Read recurring geometry, scaling and failure modes across machine artifacts.", AMBER, "RUN IN CAMPAIGN"),
        ((1528, 255, 1826, 520), "OUTPUT 05", "New candidate lemma", "Formulate a narrower statement with explicit assumptions, constants and tests.", TEAL, "AI-REFINED"),
    ]
    for box, label, title_value, body, accent, status in boxes:
        flow_box(image, box, label, title_value, body, accent=accent, status=status)
    for start_x in (386, 746, 1106, 1466):
        arrow(image, (start_x + 10, 387), (start_x + 52, 387), color=CYAN, width=4, head=12)

    rounded_card(image, (88, 580, 1832, 925), fill="#0E2733", outline=TEAL, radius=28)
    text(image, (120, 610), "THE HUMAN-DIRECTED, AI-EXECUTED RESEARCH LOOP", size=15, color=TEAL, weight="bold")
    loop = [
        ("1", "FORMULATE", "Write a lemma with a measurable finite analogue."),
        ("2", "ATTACK", "Search for the worst state, not a reassuring example."),
        ("3", "FALSIFY OR SURVIVE", "Store failure evidence or a bounded-looking branch."),
        ("4", "EXPLAIN", "Extract the geometry and normalization that caused the result."),
        ("5", "PROVE", "Return to cutoff-independent mathematics and peer review."),
    ]
    for index, (number, title_value, body) in enumerate(loop):
        x = 120 + index * 338
        ImageDraw.Draw(image).ellipse((x, 685, x + 52, 737), fill=rgb(TEAL if index < 3 else BLUE))
        text(image, (x + 26, 710), number, size=17, color=INK if index < 3 else WHITE, weight="bold", anchor="mm")
        text(image, (x, 758), title_value, size=17, color=CYAN, weight="bold", max_width=295)
        text(image, (x, 798), body, size=16, color=MUTED, max_width=286)
        if index < len(loop) - 1:
            arrow(image, (x + 66, 711), (x + 306, 711), color=LINE, width=3, head=10)
    text(
        image,
        (120, 890),
        "Guardrail: a finite certificate may reject a proposed route; it cannot by itself establish a uniform theorem.",
        size=17,
        color=AMBER,
        weight="medium",
    )
    footer(image, 3, 5, "project overview", note="Fast falsification is the engine’s highest-confidence benefit")
    return image


def bar(
    image: Image.Image,
    x: int,
    y: int,
    width: int,
    value: float,
    maximum: float,
    *,
    color: str,
    label: str,
    value_label: str,
) -> None:
    text(image, (x, y), label, size=15, color=MUTED, weight="medium", max_width=220)
    track_x = x + 225
    ImageDraw.Draw(image).rounded_rectangle((track_x, y + 2, track_x + width, y + 27), radius=12, fill=rgb("#243A51"))
    fill_width = max(3, round(width * value / maximum))
    ImageDraw.Draw(image).rounded_rectangle((track_x, y + 2, track_x + fill_width, y + 27), radius=12, fill=rgb(color))
    text(image, (track_x + width + 18, y + 1), value_label, size=15, color=TEXT, weight="bold", mono=True)


def overview_slide_4() -> Image.Image:
    image = vertical_gradient("#06111E", "#12243B")
    add_glow(image, (1200, 130, 2080, 1030), BLUE, blur=190, alpha=40)
    section_header(
        image,
        "03",
        "Navier–Stokes case study",
        "The engine already changes which lemmas are worth pursuing",
        "One route was rejected by exact adversarial evidence; the surviving route became a narrower joint normalization problem.",
        status="ACTIVE",
    )

    rounded_card(image, (88, 235, 906, 930), fill=CARD, outline=RED, radius=28)
    pill(image, (120, 266), "FALSIFIED ROUTE", fill=RED, color=WHITE, size=13)
    text(image, (120, 320), "PNT-13: standalone shell decorrelation", size=27, weight="bold", max_width=730)
    text(
        image,
        (120, 383),
        "The proposed exponential height-gap decay is contradicted by optimized finite states with near-perfect normalized shell correlation.",
        size=18,
        color=MUTED,
        max_width=720,
    )
    text(image, (120, 490), "NORMALIZED CORRELATION", size=14, color=RED, weight="bold")
    bar(image, 120, 536, 360, 0.985710, 1.0, color=RED, label="gap 4 • K8", value_label="0.985710")
    bar(image, 120, 588, 360, 0.998998, 1.0, color=RED, label="gap 6 • K8", value_label="0.998998")
    bar(image, 120, 640, 360, 0.994971, 1.0, color=RED, label="gap 7 • K12", value_label="0.994971")
    rounded_card(image, (120, 720, 872, 875), fill="#321F2C", outline=RED, radius=18, shadow=False)
    text(image, (146, 742), "RESEARCH TIME SAVED", size=13, color=RED, weight="bold")
    text(
        image,
        (146, 778),
        "Stop trying to prove correlation decay in isolation. Preserve the full coupled normalization instead.",
        size=20,
        weight="medium",
        max_width=680,
    )

    rounded_card(image, (946, 235, 1832, 930), fill=CARD, outline=TEAL, radius=28)
    pill(image, (978, 266), "SURVIVING TARGET", fill=TEAL, color=INK, size=13)
    text(image, (978, 320), "PNT-12: joint normalized tail bound", size=27, weight="bold", max_width=760)
    text(
        image,
        (978, 383),
        "The direct exact-gradient objective keeps the Gram row, tail mass and common normalization coupled. This route is still open.",
        size=18,
        color=MUTED,
        max_width=770,
    )
    simple_table(
        image,
        (978, 485, 1800, 758),
        ("Evidence", "Result", "Meaning"),
        (
            ("H16, K12", "4.18215e-4", "New finite stress record"),
            ("Gradient check", "1.21e-10", "Directional error"),
            ("Peak RSS", "4.88 GiB", "Down from 7.94 GiB"),
        ),
        (245, 205, 372),
        row_height=69,
        header_height=66,
        body_size=15,
    )
    rounded_card(image, (978, 790, 1800, 875), fill="#10352F", outline=TEAL, radius=18, shadow=False)
    text(image, (1004, 811), "STATUS", size=13, color=TEAL, weight="bold")
    text(image, (1120, 805), "Bounded-looking finite evidence; no uniform proof yet.", size=19, weight="medium", max_width=640)
    footer(image, 4, 5, "project overview", note="Finite values are experimental certificates, not theorem constants")
    return image


def gate(
    image: Image.Image,
    x: int,
    y: int,
    width: int,
    number: str,
    title_value: str,
    body: str,
    *,
    state: str,
) -> None:
    palettes = {
        "PROVED": (TEAL, "#0D342D"),
        "ACTIVE": (CYAN, "#123846"),
        "OPEN": (ORANGE, "#3A2A1E"),
        "LATER": (MUTED_2, CARD),
    }
    accent, fill = palettes[state]
    rounded_card(image, (x, y, x + width, y + 190), fill=fill, outline=accent, radius=22, shadow=False)
    pill(image, (x + 18, y + 16), f"{number}  {state}", fill=accent, color=INK if state in ("PROVED", "ACTIVE", "OPEN") else WHITE, size=11, height=28)
    text(image, (x + 18, y + 58), title_value, size=21, weight="bold", max_width=width - 36)
    text(image, (x + 18, y + 102), body, size=15, color=MUTED, max_width=width - 36)


def overview_slide_5() -> Image.Image:
    image = vertical_gradient("#07101C", "#172039")
    add_glow(image, (-250, 200, 720, 1150), TEAL, blur=180, alpha=38)
    section_header(
        image,
        "04",
        "Roadmap and guardrails",
        "A useful engine is honest about the distance to a proof",
        "The project has produced exact reductions, a proved partial tail lemma and strong falsification tooling. The central closure chain remains open.",
    )
    gates = [
        ("01", "Dynamic far tail", "A cutoff-independent partial lemma removes sufficiently separated shells.", "PROVED"),
        ("02", "PNT-12", "Prove the coupled tail/Gram/normalization bound uniformly in height and cutoff.", "ACTIVE"),
        ("03", "PNT-4 + RQ-11", "Close the complementary channel and the local/transition remainder.", "OPEN"),
        ("04", "L4 regularity", "Convert the local decomposition into the required trajectory estimate.", "OPEN"),
        ("05", "L5 to L6", "Complete the remaining regularity and limit-passage steps.", "LATER"),
    ]
    for index, (number, title_value, body, state) in enumerate(gates):
        x = 88 + index * 350
        gate(image, x, 258, 314, number, title_value, body, state=state)
        if index < len(gates) - 1:
            arrow(image, (x + 320, 352), (x + 344, 352), color=LINE, width=3, head=9)

    rounded_card(image, (88, 492, 1180, 918), fill=CARD, outline=BLUE, radius=28)
    text(image, (120, 522), "NEXT ENGINE MILESTONES", size=15, color=BLUE, weight="bold")
    simple_table(
        image,
        (120, 572, 1148, 852),
        ("Work item", "Purpose", "Exit condition"),
        (
            ("Height/cutoff scans", "Stress PNT-12 row by row", "Stable or growing branch classified"),
            ("AI artifact schema", "Compare JSON/TSV evidence", "Ranked candidate lemmas with provenance"),
            ("Analytic restart", "Prove joint normalization", "Cutoff-independent inequality"),
        ),
        (290, 330, 408),
        row_height=71,
        header_height=66,
        body_size=15,
    )
    text(image, (122, 873), "All machine claims remain reproducible from repository commands and certificates.", size=16, color=MUTED_2)

    rounded_card(image, (1220, 492, 1832, 918), fill="#2B1D28", outline=RED, radius=28)
    text(image, (1254, 524), "NON-CLAIM", size=15, color=RED, weight="bold")
    text(image, (1254, 574), "The Clay problem is not solved.", size=31, weight="bold", max_width=520)
    text(
        image,
        (1254, 646),
        "A proof still requires cutoff-independent analysis, the full local closure, the complementary normalization channel, and the remaining L4–L6 argument.",
        size=20,
        color=MUTED,
        max_width=520,
    )
    pill(image, (1254, 820), "RESEARCH SOFTWARE", fill=AMBER, color=INK, size=14)
    text(image, (1254, 866), "Evidence generator • falsifier • proof navigator", size=18, color=TEXT, weight="medium")
    footer(image, 5, 5, "project overview", note="Repository: ./navier-stokes  •  Brief: ./guides/problems/navier-stokes")
    return image


def brief_slide_1() -> Image.Image:
    image = vertical_gradient("#030A13", "#0E1B30")
    paste_cover(image, ASSETS / "navier-stokes-flow.png", (0, 0, W, H), radius=0, darken=52, outline=None)
    overlay = Image.new("RGBA", image.size, (0, 0, 0, 0))
    od = ImageDraw.Draw(overlay)
    for x in range(1250):
        alpha = max(0, 245 - round(x / 1250 * 245))
        od.line((x, 0, x, H), fill=(3, 10, 19, alpha))
    image.alpha_composite(overlay)
    pill(image, (90, 62), "NAVIER-STOKES BRIEF", fill=CYAN, color=INK, size=15)
    text(image, (90, 154), "Can a smooth fluid\nbecome infinitely rough?", size=63, weight="bold", max_width=840)
    text(
        image,
        (94, 344),
        "The 3D Navier–Stokes equations describe water and air. The Millennium Problem asks whether every smooth, physically reasonable starting flow stays smooth for all time—or whether a singularity can form in finite time.",
        size=25,
        color=MUTED,
        max_width=730,
        spacing=8,
    )
    rounded_card(image, (92, 610, 730, 833), fill="#102838", outline=CYAN, radius=26)
    text(image, (122, 638), "THE QUESTION IN ONE LINE", size=14, color=CYAN, weight="bold")
    text(image, (122, 684), "Smooth forever", size=30, color=TEAL, weight="bold")
    text(image, (378, 684), "or", size=24, color=MUTED, weight="medium")
    text(image, (434, 684), "finite blow-up?", size=30, color=ORANGE, weight="bold")
    text(
        image,
        (122, 750),
        "No generally accepted proof decides this in three dimensions.",
        size=19,
        color=TEXT,
        max_width=550,
    )
    text(image, (94, 894), "Our repository investigates one proof route computationally.", size=20, weight="medium")
    text(image, (94, 934), "It has not solved the Millennium Problem.", size=18, color=AMBER, weight="bold")
    footer(image, 1, 3, "Navier–Stokes brief", note="Conceptual illustration — not a simulation result")
    return image


def brief_slide_2() -> Image.Image:
    image = vertical_gradient("#06101C", "#10243C")
    section_header(
        image,
        "01",
        "Why it is difficult",
        "Stretching creates small scales; viscosity removes them",
        "The unresolved issue is whether the smoothing always wins before any derivative becomes unbounded.",
    )
    rounded_card(image, (88, 238, 1832, 385), fill="#0D2A38", outline=CYAN, radius=26)
    ImageDraw.Draw(image).text(
        (960, 298),
        "∂u/∂t + (u · ∇)u = −∇p + νΔu,     ∇ · u = 0",
        font=ImageFont.truetype(FONT_MATH, size=44),
        fill=rgba(WHITE),
        anchor="mm",
    )
    text(
        image,
        (960, 344),
        "velocity changes = transport + pressure + viscous smoothing; incompressibility forbids sources or sinks",
        size=16,
        color=MUTED,
        anchor="ma",
        align="center",
    )

    rounded_card(image, (88, 438, 770, 890), fill="#251D2A", outline=ORANGE, radius=28)
    text(image, (120, 468), "VORTEX STRETCHING", size=15, color=ORANGE, weight="bold")
    text(image, (120, 516), "Can sharpen gradients", size=31, weight="bold")
    bullet_list(
        image,
        [
            "Three-dimensional vortices can stretch and intensify.",
            "Energy may move toward smaller spatial scales.",
            "A proof must rule out uncontrolled concentration.",
        ],
        x=122,
        y=590,
        width=590,
        size=20,
        gap=17,
        bullet_color=ORANGE,
    )
    text(image, (120, 822), "NONLINEAR TERM", size=14, color=ORANGE, weight="bold")
    text(image, (320, 813), "(u · ∇)u", size=25, color=TEXT, weight="bold", mono=True)

    arrow(image, (805, 657), (1115, 657), color=AMBER, width=7, head=20)
    pill(image, (864, 608), "UNKNOWN BALANCE", fill=AMBER, color=INK, size=13)

    rounded_card(image, (1150, 438, 1832, 890), fill="#12322F", outline=TEAL, radius=28)
    text(image, (1182, 468), "VISCOSITY", size=15, color=TEAL, weight="bold")
    text(image, (1182, 516), "Smooths the flow", size=31, weight="bold")
    bullet_list(
        image,
        [
            "Diffusion damps high-frequency motion.",
            "Energy decreases in the classical energy estimate.",
            "Known estimates are not yet strong enough to close 3D regularity.",
        ],
        x=1184,
        y=590,
        width=590,
        size=20,
        gap=17,
        bullet_color=TEAL,
    )
    text(image, (1182, 822), "DIFFUSIVE TERM", size=14, color=TEAL, weight="bold")
    text(image, (1388, 813), "νΔu", size=25, color=TEXT, weight="bold", mono=True)
    text(
        image,
        (92, 935),
        "A successful proof must control every scale uniformly—not only every scale tested by a computer.",
        size=19,
        color=AMBER,
        weight="medium",
    )
    footer(image, 2, 3, "Navier–Stokes brief")
    return image


def brief_slide_3() -> Image.Image:
    image = vertical_gradient("#06101B", "#14213A")
    add_glow(image, (1250, 150, 2100, 1040), TEAL, blur=180, alpha=36)
    section_header(
        image,
        "02",
        "What we built and where we are",
        "Human-directed, AI-executed falsification—not an automated proof",
        "Codex wrote and ran the C++ engine, searched for hostile states, read failed routes and refined the next obligation while preserving restartable evidence.",
        status="NOT SOLVED",
    )

    simple_table(
        image,
        (88, 246, 1832, 735),
        ("Layer", "What the project does", "Current result", "Status"),
        (
            ("Spectral engine", "Finite Fourier / Galerkin Navier–Stokes algebra", "Reproducible C++20 kernels and certificates", "VERIFIED"),
            ("Gradient engine", "JVP/VJP, RK4 adjoint checks and L-BFGS search", "Fast adversarial state optimization", "VERIFIED"),
            ("Far-tail analysis", "Separate sufficiently distant frequency shells", "Cutoff-independent partial lemma", "PROVED"),
            ("PNT-13", "Test standalone shell decorrelation", "Near-perfect gap correlations found", "FALSIFIED"),
            ("PNT-12", "Keep tail mass, Gram row and normalization coupled", "Finite record 4.18215e-4; no theorem", "OPEN"),
        ),
        (285, 610, 525, 324),
        row_height=85,
        header_height=64,
        body_size=16,
        highlight_rows={3, 4},
        status_column=3,
    )

    rounded_card(image, (88, 775, 1155, 936), fill="#0F2E38", outline=CYAN, radius=22)
    text(image, (116, 800), "WHY THIS CAN SHORTEN RESEARCH", size=14, color=CYAN, weight="bold")
    text(
        image,
        (116, 840),
        "Codex can attack a weak lemma in seconds or minutes, store its counterexample, read the failure and redesign the next candidate around what survived.",
        size=21,
        weight="medium",
        max_width=980,
    )

    rounded_card(image, (1195, 775, 1832, 936), fill="#32252A", outline=AMBER, radius=22)
    text(image, (1224, 800), "WHAT REMAINS", size=14, color=AMBER, weight="bold")
    text(
        image,
        (1224, 840),
        "Prove PNT-12, close PNT-4 and RQ-11, then complete the L4 to L5 to L6 chain.",
        size=21,
        weight="medium",
        max_width=550,
    )
    footer(image, 3, 3, "Navier–Stokes brief", note="Full technical roadmap: navier-stokes/PROOF_PLAN.md")
    return image


def save_pdf(slide_paths: Sequence[Path], destination: Path, *, title: str, subject: str) -> None:
    width_points = 13.333333 * 72
    height_points = 7.5 * 72
    pdf = pdf_canvas.Canvas(str(destination), pagesize=(width_points, height_points))
    pdf.setTitle(title)
    pdf.setAuthor("Millennium Math Problems")
    pdf.setSubject(subject)
    pdf.setKeywords(f"Millennium Prize Problems, Navier-Stokes, lemma engine, version {VERSION}, {LANGUAGE}")
    for slide_path in slide_paths:
        pdf.drawImage(str(slide_path), 0, 0, width=width_points, height=height_points)
        pdf.showPage()
    pdf.save()


def uno_prop(name: str, value):
    from com.sun.star.beans import PropertyValue

    prop = PropertyValue()
    prop.Name = name
    prop.Value = value
    return prop


def package_pptx(slide_paths: Sequence[Path], destination: Path) -> None:
    try:
        import uno
        from com.sun.star.awt import Point, Size
    except ImportError as exc:
        raise RuntimeError("UNO is unavailable. Run with /usr/bin/python3.") from exc

    pipe_name = f"millennium_guides_{os.getpid()}_{destination.stem}"
    accept = f"pipe,name={pipe_name};urp;StarOffice.ServiceManager"
    resolve_url = f"uno:pipe,name={pipe_name};urp;StarOffice.ComponentContext"
    profile_dir = Path(tempfile.mkdtemp(prefix="millennium-guides-lo-"))
    process = subprocess.Popen(
        [
            "soffice",
            "--headless",
            "--nologo",
            "--nodefault",
            "--nofirststartwizard",
            "--norestore",
            f"-env:UserInstallation={uno.systemPathToFileUrl(str(profile_dir))}",
            f"--accept={accept}",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        local_context = uno.getComponentContext()
        resolver = local_context.ServiceManager.createInstanceWithContext(
            "com.sun.star.bridge.UnoUrlResolver", local_context
        )
        context = None
        deadline = time.time() + 20
        while time.time() < deadline:
            try:
                context = resolver.resolve(resolve_url)
                break
            except Exception:
                time.sleep(0.25)
        if context is None:
            raise RuntimeError("Could not connect to headless LibreOffice")

        service_manager = context.ServiceManager
        desktop = service_manager.createInstanceWithContext("com.sun.star.frame.Desktop", context)
        document = desktop.loadComponentFromURL("private:factory/simpress", "_blank", 0, ())
        pages = document.getDrawPages()
        while pages.getCount() > 1:
            pages.remove(pages.getByIndex(pages.getCount() - 1))

        page_width, page_height = 33867, 19050
        for index, slide_path in enumerate(slide_paths):
            page = pages.getByIndex(0) if index == 0 else pages.insertNewByIndex(index)
            while page.getCount() > 0:
                page.remove(page.getByIndex(0))
            page.Width = page_width
            page.Height = page_height
            for property_name in ("BorderLeft", "BorderRight", "BorderTop", "BorderBottom"):
                try:
                    setattr(page, property_name, 0)
                except Exception:
                    pass
            shape = document.createInstance("com.sun.star.drawing.GraphicObjectShape")
            shape.Position = Point(0, 0)
            shape.Size = Size(page_width, page_height)
            shape.GraphicURL = uno.systemPathToFileUrl(str(slide_path.resolve()))
            page.add(shape)

        temporary_odp = profile_dir / f"{destination.stem}.odp"
        document.storeAsURL(uno.systemPathToFileUrl(str(temporary_odp)), (uno_prop("FilterName", "impress8"),))
        document.storeToURL(
            uno.systemPathToFileUrl(str(destination.resolve())),
            (uno_prop("FilterName", "Impress MS PowerPoint 2007 XML"), uno_prop("Overwrite", True)),
        )
        document.close(True)
    finally:
        process.terminate()
        try:
            process.wait(timeout=8)
        except subprocess.TimeoutExpired:
            process.kill()
        shutil.rmtree(profile_dir, ignore_errors=True)


def contact_sheet(slide_paths: Sequence[Path], destination: Path) -> None:
    thumb_w, thumb_h = 576, 324
    margin = 24
    columns = 2
    rows = math.ceil(len(slide_paths) / columns)
    sheet = Image.new(
        "RGB",
        (columns * thumb_w + (columns + 1) * margin, rows * thumb_h + (rows + 1) * margin),
        rgb("#040A12"),
    )
    for index, slide_path in enumerate(slide_paths):
        slide = Image.open(slide_path).convert("RGB").resize((thumb_w, thumb_h), RESAMPLE_LANCZOS)
        x = margin + (index % columns) * (thumb_w + margin)
        y = margin + (index // columns) * (thumb_h + margin)
        sheet.paste(slide, (x, y))
    sheet.save(destination, optimize=True)


def render_deck(
    builders: Sequence[Callable[[], Image.Image]],
    output_dir: Path,
    stem: str,
    *,
    title: str,
    subject: str,
    preview_name: str,
    no_office: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{stem}-") as temporary:
        slide_dir = Path(temporary)
        slide_paths: list[Path] = []
        for index, builder in enumerate(builders, start=1):
            slide_path = slide_dir / f"{index:02d}.png"
            builder().convert("RGB").save(slide_path, optimize=True)
            slide_paths.append(slide_path)
            print(f"rendered {stem} page {index}/{len(builders)}")

        pdf_path = output_dir / f"{stem}.pdf"
        save_pdf(slide_paths, pdf_path, title=title, subject=subject)
        contact_sheet(slide_paths, ASSETS / preview_name)
        print(f"rendered {pdf_path}")
        if not no_office:
            pptx_path = output_dir / f"{stem}.pptx"
            package_pptx(slide_paths, pptx_path)
            print(f"rendered {pptx_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-office", action="store_true", help="Skip PPTX packaging")
    args = parser.parse_args()

    hero = ASSETS / "navier-stokes-flow.png"
    if not hero.is_file():
        raise FileNotFoundError(f"Missing required generated hero image: {hero}")

    render_deck(
        (overview_slide_1, overview_slide_2, overview_slide_3, overview_slide_4, overview_slide_5),
        OVERVIEW_OUT,
        OVERVIEW_STEM,
        title="Millennium Math Problems Project Overview",
        subject="Public five-page overview of the computational lemma research project",
        preview_name="project-overview-preview.png",
        no_office=args.no_office,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
