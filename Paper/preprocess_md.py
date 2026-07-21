"""Pre-process PAPER_draft.md for pandoc -> LaTeX conversion."""
import re, sys

# Citation number -> BibTeX key mapping
CITE_MAP = {
    "1":  "Thompson2022",
    "2":  "Plimpton1995",
    "3":  "Buehler2006",
    "4":  "He2024",
    "5":  "Famprikis2019",
    "6":  "Batatia2022",
    "7":  "Erhard2024",
    "8":  "Cheng2019",
    "9":  "Bauer2012",
    "10": "Brown2011",
    "11": "Edwards2014",
    "12": "Gao2025",
    "13": "Gissinger2024a",
    "14": "Gissinger2024b",
    "15": "Maresca2020",
    "16": "Ferreira2025",
    "17": "Zhang2025",
    "18": "Aktulga2012",
    "19": "Austin2021",
    "20": "Christ2010",
    "21": "You2026",
    "22": "Shirts2008",
    "23": "Alvarado2022",
    "24": "Nijkamp2022",
    "25": "Hermann2017",
    "26": "Lahkar2023",
    "27": "Toth2025",
    "28": "Jorgensen1996",
    "29": "Beutler1994",
    "30": "Bhati2017",
    "31": "Chen2021",
    "32": "Gapsys2020",
    "33": "Joung2008",
    "34": "Xu2004",
    "35": "Xu2014",
    "36": "Bhatt2015",
    "37": "Wan2016",
    "38": "Takeuchi2012",
    "39": "Chaban2015",
    "40": "Pirani2004",
}

# Figure placeholder -> image path mapping  (relative to Paper/)
FIGURE_MAP = {
    "1": ("../tools/plots/figure1_workflow.png",       "AI-assisted LAMMPS optimization workflow."),
    "2": ("../tools/plots/figure2_omp_coverage.png",   "OpenMP acceleration coverage in LAMMPS before and after optimization."),
    "3": ("../tools/plots/figure3_build_flags.png",    "Compiler flag optimization on the 32,000-atom benchmark."),
    "4": ("../tools/plots/figure4_omp_scaling.png",    "OpenMP pair\\_style optimization."),
    "5": ("../tools/plots/figure5_summary.png",        "Combined optimization summary."),
    "6": ("../tools/plots/figure6_fep_casestudy.png",  "OPLS-AA FEP case study across three battery electrolyte solvents."),
}

def replace_citations(text):
    """Convert [1,2,3] and [1] to [@key1; @key2; @key3]."""
    def repl(m):
        nums = [n.strip() for n in m.group(1).split(",")]
        keys = [CITE_MAP.get(n) for n in nums]
        if all(keys):
            return "[@" + "; @".join(keys) + "]"
        return m.group(0)  # leave unchanged if not in map
    # Match [digits] or [digits,digits,...] not preceded by ! (images) or ( (links)
    return re.sub(r'(?<![!(])\[(\d+(?:,\s*\d+)*)\]', repl, text)

def replace_figure_blocks(text):
    """Replace **[FIGURE N]** + italic caption lines with pandoc figure syntax."""
    lines = text.split("\n")
    out = []
    i = 0
    while i < len(lines):
        m = re.match(r'\*\*\[FIGURE (\d+)\]\*\*\s*$', lines[i].strip())
        if m:
            fignum = m.group(1)
            # Skip blank lines and collect the caption line (italic)
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            caption_line = ""
            if j < len(lines) and lines[j].strip().startswith("*") and lines[j].strip().endswith("*"):
                caption_line = lines[j].strip()[1:-1]  # strip outer *
                j += 1
            path, default_cap = FIGURE_MAP.get(fignum, (None, ""))
            cap = caption_line if caption_line else default_cap
            # Escape underscores in caption for LaTeX (pandoc handles, but be safe)
            if path:
                out.append(f"![Figure {fignum}: {cap}]({path}){{width=\\textwidth}}")
                out.append("")
            i = j
        else:
            out.append(lines[i])
            i += 1
    return "\n".join(out)

def replace_table_markers(text):
    """Remove **[TABLE N]** markers (the actual markdown tables follow)."""
    text = re.sub(r'\*\*\[TABLE \d+\]\*\*\s*\n', '', text)
    return text

def replace_unicode_math(text):
    """Replace common Unicode math/chemistry chars with LaTeX equivalents."""
    # Superscripts/subscripts already handled by sequence regex before this function
    # Greek letters (outside code spans)
    greek = {
        'λ': r'$\lambda$', 'ε': r'$\varepsilon$', 'σ': r'$\sigma$',
        'ρ': r'$\rho$', 'δ': r'$\delta$', 'τ': r'$\tau$', 'α': r'$\alpha$',
        'β': r'$\beta$', 'γ': r'$\gamma$', 'ω': r'$\omega$', 'μ': r'$\mu$',
    }
    # Math symbols
    math_sym = {
        '⟨': r'$\langle$', '⟩': r'$\rangle$',
        '≈': r'$\approx$', '≥': r'$\geq$', '≤': r'$\leq$',
        '×': r'$\times$', '±': r'$\pm$',
        '→': r'$\rightarrow$', '∫': r'$\int$',
        '∂': r'$\partial$',
        '∑': r'$\sum$', '∞': r'$\infty$',
        '–': '--', '—': '---',
        'Å': r'\AA{}',
        '■': r'$\blacksquare$',
        '𝜎': r'$\sigma$',
    }
    all_map = {**greek, **math_sym}
    # Apply replacements outside of backtick spans
    def repl_outside_code(text, char, replacement):
        parts = text.split('`')
        for i in range(0, len(parts), 2):  # even indices = outside backticks
            parts[i] = parts[i].replace(char, replacement)
        return '`'.join(parts)
    for char, repl in all_map.items():
        text = repl_outside_code(text, char, repl)
    return text


def insert_supplementary_figures(text):
    """Insert actual figure includes for supplementary figures."""
    sfig_map = {
        "S1": ("../tools/plots/figure_s1_omp_scaling.png",
               "Full OMP scaling benchmarks for three newly parallelized pair styles."),
        "S2": ("../tools/plots/figure_s2_nve_energy.png",
               "NVE energy conservation traces comparing serial and OMP-4-thread execution."),
        "S3": ("../tools/plots/figure_s3_a3_backport.png",
               "A-3 code transformation (before/after pair hotloop)."),
    }
    for key, (path, cap) in sfig_map.items():
        # Replace the italic placeholder after "Supplementary Figure SN:" heading
        pattern = r'(## Supplementary Figure ' + re.escape(key) + r':.*?\n\n)(\*\[.*?\]\*\n\n)?'
        insert = f'\n\n![Figure {key}: {cap}]({path}){{width=\\textwidth}}\n\n'
        text = re.sub(pattern, lambda m: m.group(1) + insert, text, count=1)
    return text


def clean_supplementary_list(text):
    """The supplementary list at end of paper is just informational."""
    return text

def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else "../PAPER_draft.md"
    with open(infile, encoding="utf-8") as f:
        text = f.read()

    # Remove the journal submission line at end
    text = re.sub(r'\*Manuscript submitted to.*\*\s*$', '', text, flags=re.MULTILINE)

    # Handle √N patterns first (context-aware)
    text = re.sub(r'√(\d+)', r'$\\sqrt{\1}$', text)
    text = text.replace('√', r'$\surd$')

    # Handle combined super/subscript digit sequences (scientific notation)
    # Map superscript and subscript digit chars
    sup_digits = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁺':'+','⁻':'-'}
    sub_digits = {'₀':'0','₁':'1','₂':'2','₃':'3','₄':'4','₅':'5','₆':'6','₇':'7','₈':'8','₉':'9'}
    sup_class = ''.join(sup_digits.keys())
    sub_class = ''.join(sub_digits.keys())
    def convert_supscript_seq(m):
        seq = m.group(0)
        return '$^{' + ''.join(sup_digits[c] for c in seq) + '}$'
    def convert_subscript_seq(m):
        seq = m.group(0)
        return '$_{' + ''.join(sub_digits[c] for c in seq) + '}$'
    # Replace sequences of superscript chars (including sign)
    text = re.sub(f'[{re.escape(sup_class)}]+', convert_supscript_seq, text)
    text = re.sub(f'[{re.escape(sub_class)}]+', convert_subscript_seq, text)

    # Replace unicode math symbols
    text = replace_unicode_math(text)

    # Replace figure blocks with pandoc image syntax
    text = replace_figure_blocks(text)

    # Remove [TABLE N] markers (leave the italic caption that follows)
    text = replace_table_markers(text)

    # Replace citation numbers with BibTeX keys
    text = replace_citations(text)

    # Remove the leading title (in YAML front matter)
    text = re.sub(r'^# AI-Assisted.*\n', '', text, count=1)

    # Remove author/affiliation/correspondence lines (in YAML front matter)
    text = re.sub(r'^\*\*Chun-Yeol You.*\*\*\s*\n', '', text, flags=re.MULTILINE, count=1)
    text = re.sub(r'^¹ Department.*\n', '', text, flags=re.MULTILINE, count=1)
    text = re.sub(r'^\*Correspondence:.*\*\s*\n', '', text, flags=re.MULTILINE, count=1)

    # Insert supplementary figures
    text = insert_supplementary_figures(text)

    # Remove repeated --- horizontal rules (keep one)
    text = re.sub(r'(\n\n---\n\n){2,}', '\n\n---\n\n', text)

    # Prepend YAML front matter
    yaml_header = """---
title: "AI-Assisted Systematic Optimization of the LAMMPS Molecular Dynamics Simulator via Large Language Model Coding Agents"
author:
  - "Chun-Yeol You^1^"
institute:
  - "^1^Department of Physics and Chemistry, DGIST, Daegu, Republic of Korea"
date: "2026"
bibliography: refs.bib
biblio-style: unsrt
geometry: "margin=1in"
fontsize: 11pt
documentclass: article
header-includes:
  - \\usepackage{booktabs}
  - \\usepackage{longtable}
  - \\usepackage{pdflscape}
  - \\usepackage[hyphens]{url}
  - \\usepackage{amsmath}
  - \\usepackage{amssymb}
  - \\usepackage{graphicx}
  - \\usepackage{caption}
  - \\usepackage{float}
  - \\usepackage{microtype}
  - \\usepackage[hidelinks]{hyperref}
  - \\usepackage{xcolor}
  - \\usepackage{listings}
  - \\usepackage{setspace}
  - \\usepackage{array}
  - \\usepackage{tabularx}
  - \\lstset{basicstyle=\\small\\ttfamily, breaklines=true, frame=single, backgroundcolor=\\color{gray!10}, columns=flexible}
  - \\setlength{\\parskip}{6pt}
  - \\setlength{\\parindent}{0pt}
  - \\renewcommand{\\figurename}{Figure}
  - \\captionsetup{labelfont=bf, font=small}
---

"""
    text = yaml_header + text

    sys.stdout.buffer.write(text.encode("utf-8"))

if __name__ == "__main__":
    main()
