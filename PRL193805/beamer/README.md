# CSU Physics Beamer

A 16:9 XeLaTeX Beamer template for physics and scientific presentations at
Central South University (中南大学物理学术报告模板).

Based on:

- [Sakulyn/CSU-Beamer](https://github.com/Sakulyn/CSU-Beamer)
- [FangWHao/THU-beamer-template](https://github.com/FangWHao/THU-beamer-template)
- the SINTEF / Sapienza Beamer theme lineage
  (original theme by Federico Zenith, SINTEF; adapted by Andrea Gasparini for Sapienza)

## Features

- 16:9 widescreen (paper 16 cm × 9 cm)
- XeLaTeX + ctex / xeCJK for full Chinese support
- Native LaTeX mathematics (`\usefonttheme[onlymath]{serif}`)
- Clean scientific slide layout: blocks, columns, side-picture slides,
  chapter slides, TikZ-ready
- CSU branding:
  - CSU logo in the headline
  - alternating footer — odd frames: page number + CSU motto (校训),
    even frames: CSU spirit (校风) + page number
  - faded CSU seal watermark on the title page
- LXGWWenKai bundled for optional Chinese titles
- PDF-first workflow: `csu-beamer.pdf` in the repository is the compiled preview

## Compile

Run XeLaTeX twice (the second pass resolves the footer and watermark positions):

```bash
xelatex -interaction=nonstopmode csu-beamer.tex
xelatex -interaction=nonstopmode csu-beamer.tex
```

## Fonts

- The bundled LXGWWenKai (`fonts/`) is used optionally for Chinese titles
  via `\titlefont`.
- Body fonts (Chinese and Latin) follow your XeLaTeX / TeX configuration;
  no special setup is required.
- Mathematics uses Beamer serif math.

## Using as a template

Treat this repository as a **read-only template**. Do not copy its `.git`
history into a research project, and do not make the target project a Git
fork of this template — each project keeps its own Git history.

Recommended: copy with `rsync`, excluding the repository metadata:

```bash
rsync -a \
  --exclude='.git' \
  --exclude='render-*' \
  CSU-Physics-Beamer/ \
  /path/to/target-project/talk/
```

(Alternatively copy the folder and delete `.git` by hand.)

Then, inside the copy:

- edit `csu-beamer.tex` (title, author, date),
- replace the contents of `sections/` with your own slides,
- add your own figures,
- keep the theme files (`beamerthemesintef.sty`, `sintefcolor.sty`),
  `fonts/` and `assets/` untouched.

To update the template later, `git pull` in your local template clone and
copy again — do not edit the template clone for a specific talk.

## License

GPL v3, inherited from the upstream template. See [LICENSE](LICENSE).
