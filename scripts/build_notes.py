"""Regenerate the figures and PDF for docs/notes-for-the-exacto-author.md.

Figures are drawn from results/exacto_results.json and site/data.json, so the
document cannot quote a number the pipeline no longer produces. Run after a
fresh `python -m pipeline.build_site`:

    python scripts/build_notes.py

Needs matplotlib, markdown and wkhtmltopdf.
"""

from __future__ import annotations

import base64
import json
import pathlib
import re
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
FIGURES = DOCS / "figures"
OK, WARN, BAD, MUT = "#2f7d5d", "#b8860b", "#a83232", "#8a8a8a"


def figures(paths: dict) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 160,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "font.family": "DejaVu Sans"})

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    labels = ["Reads cover\nlocus", "Allele in\nRNA", "Exacto called\nin RNA",
              "Protein\ntranslated", "Correct residue\n(any candidate)",
              "Correct residue\n(single pick)"]
    keys = ["covered", "allele", "called", "translated", "residue", "residue_pick"]
    any_measured = False
    for row in paths["variant_funnel"]:
        by_key = {s["key"]: s for s in row["stages"]}
        ceiling = [by_key[k]["n"] for k in keys[:5]]
        if any(v is None for v in ceiling):
            continue
        colour = OK if row["method"] == "reads" else BAD
        style = {"marker": "o", "ms": 4} if row["method"] == "reads" \
            else {"marker": "s", "ms": 4, "ls": "--"}
        ax.plot(range(5), ceiling, color=colour, lw=2.0, alpha=.75, **style)

        # The last rung has two versions. "Any candidate" is a ceiling; the
        # single pick is what a caller actually gets. Drawing only the first
        # would repeat the overstatement this figure exists to correct.
        pick = by_key.get("residue_pick")
        if pick and not pick["pending"]:
            any_measured = True
            ax.plot([4, 5], [ceiling[-1], pick["n"]], color=colour, lw=2.0,
                    alpha=.45, ls=":", marker="D", ms=4)

    ax.set_xticks(range(6 if any_measured else 5))
    ax.set_xticklabels(labels[:6] if any_measured else labels[:5])
    ax.set_ylabel("mutations remaining (of 37)")
    if not any_measured:
        # Say so in the figure, not only in the caption: an axis silently
        # truncated reads as an axis that does not exist.
        ax.set_xlim(-0.3, 5.2)
        ax.axvspan(4.35, 5.2, color="#f0f0f0", zorder=0)
        ax.annotate("single-pick\nnot yet measured", xy=(4.78, ax.get_ylim()[1] * 0.55),
                    ha="center", va="center", fontsize=8, color=MUT, style="italic")
    ax.legend(handles=[
        Line2D([], [], color=OK, lw=2, marker="o", ms=4, label="reads"),
        Line2D([], [], color=BAD, lw=2, ls="--", marker="s", ms=4, label="assembly"),
        Line2D([], [], color=MUT, lw=2, ls=":", marker="D", ms=4,
               label="single pick"),
    ], frameon=False, fontsize=8, loc="upper center",
       bbox_to_anchor=(0.52, 1.0), ncol=3)
    fig.tight_layout(); fig.savefig(FIGURES / "variant-funnel.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    for path, colour in zip(paths["paths"][:2], (OK, BAD)):
        ys = [s["of_input"] * 100 for s in path["stages"]]
        ax.plot(range(len(ys)), ys, marker="o", ms=4, color=colour, lw=2,
                label=path["label"])
        ax.annotate(f"{ys[-1]:.1f}%", xy=(len(ys) - 1, ys[-1]), xytext=(4, 6),
                    textcoords="offset points", color=colour, fontsize=9,
                    weight="bold")
    ax.set_yscale("log"); ax.set_ylabel("% of input surviving (log)")
    ax.set_xticks(range(5))
    ax.set_xticklabels(["input", "assembled /\ncapped", "support\nfilter",
                        "aligned", "spliced\nfilter"])
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIGURES / "sequence-funnel.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 2.8))
    rows = [(f"{p['label']}\n{plat}", q["frameshift_fraction"] * 100,
             q["n_proteoforms"])
            for p in paths["paths"] for plat, q in p["by_platform"].items()
            if q["n_proteoforms"] and q["frameshift_fraction"] is not None]
    rows.sort(key=lambda r: -r[1])
    bars = ax.bar([r[0] for r in rows], [r[1] for r in rows],
                  color=[BAD, WARN, OK][:len(rows)], width=.55)
    for bar, (_, pct, n) in zip(bars, rows):
        ax.annotate(f"{pct:.0f}%\nn={n}",
                    xy=(bar.get_x() + bar.get_width() / 2, pct),
                    xytext=(0, 3), textcoords="offset points", ha="center",
                    fontsize=8)
    ax.set_ylabel("% of proteoforms frameshifted")
    # Room for the value labels sitting above each bar, plus title padding:
    # without both, the tallest label runs into the title.
    ax.set_ylim(0, max((r[1] for r in rows), default=1) * 1.42)
    fig.tight_layout(); fig.savefig(FIGURES / "frameshift-by-platform.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 2.8))
    prof = paths["vaf_profile"]
    bars = ax.bar([r["label"] for r in prof], [r["median_vaf"] or 0 for r in prof],
                  color=[OK, WARN, MUT], width=.55)
    for bar, r in zip(bars, prof):
        ax.annotate(f"{r['median_vaf']:.3f}\nn={r['n']}",
                    xy=(bar.get_x() + bar.get_width() / 2, r["median_vaf"] or 0),
                    xytext=(0, 3), textcoords="offset points", ha="center",
                    fontsize=8)
    ax.set_ylabel("median ONT VAF")
    ax.set_ylim(0, max((r["median_vaf"] or 0 for r in prof), default=1) * 1.42)
    fig.tight_layout(); fig.savefig(FIGURES / "vaf-profile.png"); plt.close(fig)



def orf_provenance_figure() -> None:
    """Schematic: how the reading frame gets established, best to worst.

    Not drawn from results -- it is a taxonomy, not a measurement -- but kept
    here so the document rebuilds in one command.
    """
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(4, 1, figsize=(7.2, 5.2), sharex=True)
    line = "#cfcfcf"
    tiers = [
        ("A  Anchored to annotated start", OK,
         "observed sequence reaches the annotated ATG "
         "\u2014 frame read from the annotation", 0.0, 1.0, True),
        ("B  Anchored, 5' truncated", "#5b8c6e",
         "matches a reference transcript but stops short of the start codon",
         0.22, 1.0, False),
        ("C  Stitched", WARN,
         "fragment joined to a reference transcript by sequence overlap "
         "\u2014 the join is inferred", 0.46, 0.82, False),
        ("D  De novo", BAD,
         "no annotated transcript applies \u2014 nothing external constrains "
         "the frame", 0.30, 0.72, False),
    ]
    for ax, (title, colour, note, x0, x1, reaches) in zip(axes, tiers):
        ax.add_patch(Rectangle((0.0, 0.50), 1.0, 0.15, facecolor=line,
                               edgecolor="none"))
        if title.startswith("D"):
            ax.text(0.5, 0.575, "no annotated transcript", ha="center",
                    va="center", fontsize=7.5, color="#777", style="italic")
        else:
            ax.text(0.005, 0.575, "ATG", ha="left", va="center", fontsize=7.5,
                    color="#555")
            ax.text(0.995, 0.575, "stop", ha="right", va="center", fontsize=7.5,
                    color="#555")
        ax.add_patch(Rectangle((x0, 0.20), x1 - x0, 0.19, facecolor=colour,
                               edgecolor="none", alpha=.85))
        if title.startswith("C"):
            ax.add_patch(Rectangle((0.0, 0.20), x0, 0.19, facecolor=colour,
                                   edgecolor=colour, alpha=.18, hatch="///",
                                   lw=.5))
            ax.annotate("", xy=(x0 + 0.02, 0.44), xytext=(x0 - 0.12, 0.44),
                        arrowprops={"arrowstyle": "<->", "color": WARN, "lw": 1})
            ax.text(x0 - 0.05, 0.455, "overlap", ha="center", va="bottom",
                    fontsize=7, color=WARN)
        if reaches:
            ax.plot([0.0, 0.0], [0.18, 0.67], color=OK, lw=1.2, ls=":")
        ax.text(0.0, 0.94, title, fontsize=9, weight="bold", color=colour,
                va="top")
        ax.text(0.0, 0.02, note, fontsize=7.5, color="#444", va="bottom")
        ax.set_xlim(-0.14, 1.06)
        ax.set_ylim(0, 1.0)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGURES / "orf-provenance.png", dpi=160)
    plt.close(fig)


CSS = """
@page {
  size: A4;
  margin: 20mm 18mm;
  @bottom-center { content: counter(page) " / " counter(pages); font-size: 8pt;
                   color: #888; }
}
body { font-family: "Charter","Georgia",serif; font-size: 10.5pt; line-height: 1.5;
       color: #1a1a1a; }
h1 { font-size: 20pt; margin: 0 0 .2em; line-height: 1.2; }
h2 { font-size: 14pt; margin: 1.6em 0 .4em; padding-bottom: .2em;
     border-bottom: 1px solid #ddd; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 1.2em 0 .3em; page-break-after: avoid; }
p, li { orphans: 3; widows: 3; }
code { font-family: "SF Mono","Menlo",monospace; font-size: 9pt;
       background: #f4f4f2; padding: 1px 3px; border-radius: 3px; }
pre { background: #f4f4f2; padding: 8px 10px; border-radius: 4px;
      font-size: 8.5pt; line-height: 1.4; page-break-inside: avoid; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: .8em 0; font-size: 9.5pt;
        page-break-inside: avoid; }
th, td { border-bottom: 1px solid #ddd; padding: 5px 7px; text-align: left;
         vertical-align: top; }
th { border-bottom: 1.5px solid #999; font-weight: 600; }
figure { margin: 1.1em 0; page-break-inside: avoid; text-align: center; }
figure img { max-width: 100%; }
figcaption { margin-top: .45em; font-size: 9pt; font-weight: 600; color: #333;
             text-align: left; line-height: 1.35; }
hr { border: none; border-top: 1px solid #ddd; margin: 1.6em 0; }
a { color: #1a5490; text-decoration: none; }
"""


def build_pdf() -> pathlib.Path:
    import markdown

    source = DOCS / "notes-for-the-exacto-author.md"
    body = markdown.markdown(source.read_text(),
                             extensions=["tables", "fenced_code", "attr_list"])

    def embed(match: re.Match) -> str:
        alt, path = match.group(1), match.group(2)
        data = base64.b64encode((DOCS / path).read_bytes()).decode()
        return (
            f'<figure><img alt="{alt}" src="data:image/png;base64,{data}">'
            f'<figcaption>{alt}</figcaption></figure>'
        )

    body = re.sub(r'<p><img alt="([^"]*)" src="([^"]*)"\s*/?></p>', embed, body)
    html = DOCS / "_notes.html"
    html.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Notes for the Exacto author</title>"
        f"<style>{CSS}</style></head><body>{body}</body></html>"
    )
    pdf = DOCS / "notes-for-the-exacto-author.pdf"

    # Chrome rather than wkhtmltopdf. Both embed ToUnicode maps, but
    # wkhtmltopdf's Qt WebKit lays glyphs out so that Preview cannot select
    # across them -- the text is present and extractable by pdftotext, and
    # still unusable to a reader trying to quote it.
    chrome = next(
        (c for c in (
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/usr/bin/google-chrome", "/usr/bin/chromium",
        ) if pathlib.Path(c).exists()),
        None,
    )
    if chrome is None:
        raise SystemExit("no Chrome found for PDF rendering")
    subprocess.run(
        [chrome, "--headless", "--disable-gpu", "--no-pdf-header-footer",
         f"--print-to-pdf={pdf}", html.as_uri()],
        check=True, capture_output=True, timeout=120,
    )
    html.unlink(missing_ok=True)
    return pdf


def main() -> None:
    data = ROOT / "site" / "data.json"
    if not data.exists():
        raise SystemExit("site/data.json missing — run pipeline.build_site first")
    figures(json.loads(data.read_text())["paths"])
    orf_provenance_figure()
    print(f"figures -> {FIGURES}")
    print(f"pdf     -> {build_pdf()}")


if __name__ == "__main__":
    sys.exit(main())
