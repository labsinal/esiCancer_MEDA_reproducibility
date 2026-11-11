"""
Supplementary Tables S1–S4 as separate .xlsx files (repo/tables/).
  S1 esiCancer parameters and the values used   S2 esiTables (ESCC table, six-gene table and interactions)
  S3 configuration and outcome of every simulated condition   S4 endpoints, denominators and calculations
Inputs: repo/esitables/*.csv, output/conditions_occurrence_meda.csv, output/germline/germline_summary.csv,
output/turnover/exp2 results (values copied from output/turnover_results.txt), and the Table S4 text in
table_s4_text.py (same text as the manuscript). Rules of esiCancer described here were checked in its source code
(cellsystem.cpp, cell.cpp, mutationtable.cpp).
"""
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

import table_s4_text as v4
from paths import ROOT, OUT

TAB = ROOT / "tables"; TAB.mkdir(exist_ok=True)
ESI = ROOT / "esitables"
FONT, BOLD = Font(name="Arial", size=9), Font(name="Arial", size=9, bold=True)
TITLE, NOTE = Font(name="Arial", size=11, bold=True), Font(name="Arial", size=8, italic=True)
FLAG = Font(name="Arial", size=9, color="C00000")
HEAD_FILL = PatternFill("solid", fgColor="E7E6E6")
THIN = Side(style="thin", color="A6A6A6")
BOX = Border(top=THIN, bottom=THIN, left=THIN, right=THIN)
WRAP = Alignment(wrap_text=True, vertical="top")


def sheet(ws, title, legend, header, rows, widths, notes=(), freeze_cols=1):
    """Title, legend, header (row 4), data from row 5, notes below. Returns the first data row."""
    ws["A1"] = title; ws["A1"].font = TITLE
    ws["A2"] = legend; ws["A2"].font = FONT; ws["A2"].alignment = Alignment(wrap_text=False)
    for j, h in enumerate(header, 1):
        c = ws.cell(4, j, h); c.font = BOLD; c.fill = HEAD_FILL; c.border = BOX; c.alignment = WRAP
    for i, row in enumerate(rows, 5):
        for j, v in enumerate(row, 1):
            c = ws.cell(i, j, v); c.border = BOX; c.alignment = WRAP
            c.font = FLAG if isinstance(v, str) and ("[TO " in v) else FONT
    r = 5 + len(rows) + 1
    for n in notes:
        ws.cell(r, 1, n).font = NOTE; r += 1
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = ws.cell(5, freeze_cols + 1)
    return 5


# ------------------------------------------------------------------ Table S1
def table_s1():
    header = ["Group", "Parameter (esiCancer field)", "Definition", "Value(s) used", "Varied in", "Source / note"]
    rows = [
        ("Population", "Initial number of cells (Number of Cells)", "Cells at generation 0; all identical, no event.",
         "1,000", "–", "Default of esiCancer (ref. 17)."),
        ("Population", "Probability of division (Prob of Division)", "Probability that a cell divides in one generation.",
         "0.01", "0.02 (Figure 3G; Figure S4), 0.2 (Figure 3H; Figure S4)", "Set equal to the probability of death (CTOR)."),
        ("Population", "Probability of death (Prob of Death)", "Probability that a cell dies in one generation.",
         "0.01", "As probability of division", "Equal to the probability of division, so the expected population is constant."),
        ("Population", "Cell turnover rate (CTOR)", "Probability of division = probability of death.",
         "0.01", "0.02, 0.2", "Ref. 21. Division = death = CTOR in every condition; the outputs give about 0.01, 0.02 and 0.2 divisions per cell per generation."),
        ("Population", "Maximum divisions per cell (Max Divisions)", "Divisions after which a cell stops dividing (senescence).",
         "1,000", "10,000 in the matched-turnover runs (Figure S4)", "Ref. 21. At CTOR 0.01 a lineage divides about 48 times in 92 years, so the limit never binds."),
        ("Mutation input", "Events per division (Events Per Division)", "Number of events (mutations) drawn at each division.",
         "1 (default)", "Constant 1, 3, 5, 10 (Figure 3A, D); 2, 3, 4, 6 (six-gene model, Figure 3C, J); increasing with age (Figure 3B, E, G, H)",
         "Integer values only."),
        ("Mutation input", "Increasing schedule (modified esiCancer)",
         "Events per division = 1 + floor(generation / step); the value set in the interface is not used.",
         "Steps of 20, 10, 7 and 5 years (522 generations for the 10-year step); nominal levels 5, 10, 14, 20 (= 100 / step), "
         "which reach 5, 10, 14 and 19 events per division at 92 years", "Figure 3B, E, G, H; Figure 2C",
         "Function mutation_rate() in cellsystem.cpp (esicancer/increasing_events_per_division.patch)."),
        ("Mutation input", "Genome size (esiTable)", "Denominator of the hit probability of every gene.",
         "16,647,721 (sum of gene lengths, ESCC table); 750 (six-gene table)", "×2 and ×20 in the matched-turnover runs (Figure S4)",
         "Refs. 19, 20. Multiplying the genome size by k divides each gene's chance of being hit by k."),
        ("Mutation input", "Oncogenome size", "Sum of the event weights (PROBEVENT) of the driver records.", "58", "–", "Derived from Table S2."),
        ("Gene records", "PROBEVENT", "Relative event weight; a gene's chance of being hit by one event is its weight ÷ genome size.",
         "1 for driver and essential records; 5–45 for passenger records", "–", "Passenger weights reflect mutational target size (e.g. TTN 45, MUC16 17)."),
        ("Gene records", "DOMINANCE", "Fraction of the effect expressed when only one allele carries the event.",
         "0–1 per record", "–", "20/20 rule (ref. 3). A hit on one allele adds effect × dominance; a hit on the second allele "
         "adds effect × (1 − dominance), so both alleles together give the full effect."),
        ("Gene records", "DIVFUNC, DIVRATE", "Change of the probability of division when the event occurs (ADD = added, MULT = multiplied).",
         "Drivers: +f × 0.1; essential genes: −f × 0.1; passengers: 0 (f = alteration frequency in tumours)", "–", "Frequencies from cBioPortal ESCC datasets."),
        ("Gene records", "DEFUNC, DEATHRATE", "Change of the probability of death when the event occurs.",
         "Drivers: −f × 0.1; essential genes: +f × 0.1; passengers: 0", "–", "As DIVRATE, with opposite sign."),
        ("Gene records", "mutationFUNC, mutation", "Change of the events per division of the cell (mutation-rate modifier).",
         "+1 for TP53 and ATM; 0 for all other records", "–",
         "Not active: esiCancer reads this field only from columns named 'mutfunc' and 'mut', and the esiTable names them "
         "'mutationFUNC' and 'mutation', so no record changed the events per division in any simulation."),
        ("Gene records", "MAXDIVFUNC, MAXDIVRATE", "Change of the maximum number of divisions.", "+100 for TERT; 0 otherwise", "–", "Telomere maintenance."),
        ("Gene records", "MEFUNC, MICROENVIRONMENT", "Effect on the tumour environment.", "MULT 1 (no effect)", "–", "Tumour environment switched off."),
        ("Bounds", "Probability limits", "What happens when an update would take a probability below 0 or above 1.",
         "Upper: probabilities of division and death are capped at 1, and when their sum exceeds 1 both are rescaled so that "
         "it equals 1. Lower: no clipping; a negative value acts as 0 (the cell never divides, or never dies), and later "
         "effects are added to the stored negative value.", "–", "Requested by the reviewer (comment 8)."),
        ("Interactions", "Interaction table (DIVMOD, DEMOD, MAXDIVMOD)", "Modifiers applied when both genes of a pair carry an event.",
         "Six-gene model: every ordered pair, DIVMOD 3, DEMOD 1/3, MAXDIVMOD 1,000; none in the other models", "Figure 3C, J, K",
         "For each partner event already present in the cell, the effect of the new event on division is multiplied by DIVMOD "
         "and its effect on death by DEMOD; with DEMOD 1/3 the death effect is divided by 3. MAXDIVMOD has no effect "
         "(MAXDIVRATE is 0)."),
        ("Manual events", "Aging events", "Event applied to 100% of cells at a set age.",
         "4 events at 60, 70, 80, 90 years; each multiplies the probability of division by 0.9 and of death by 1.1",
         "Figure 3D, E, G, H; Figure 2B", "Records aging_1–aging_4 of Table S2."),
        ("Manual events", "Inherited first hit", "Gene A mutated in 100% of cells at generation 0.", "Carriers only",
         "Figure 3J, K", "Six-gene model, 3 events per division."),
        ("Time and stop", "Generation length", "Time step of the model.", "1 week (52.143 generations per year)", "–",
         "Not clinically calibrated."),
        ("Time and stop", "Number of generations (stop)", "End of follow-up.", "4,800 generations ≈ 92 years", "Truncated at 60 and 75 years in the analysis (Figure S3)", ""),
        ("Time and stop", "Number of cells (stop)", "esiTumor threshold.", "10,000 cells (10× the initial population)", "–",
         "Populations that crossed it always continued to grow."),
        ("Time and stop", "Seeds", "Random seeds; one run (individual) per seed.", "1–1,000", "1–3,000 per arm in the inherited first-hit experiment", ""),
        ("Other", "Tumour environment; asymmetric mitosis", "Optional modules of esiCancer.", "Off", "–", ""),
    ]
    wb = Workbook(); ws = wb.active; ws.title = "Table S1"
    sheet(ws, "Supplementary Table S1. Parameters of esiCancer and the values used in this study",
          "Field names as in the esiCancer interface and esiTable.",
          header, rows, [14, 30, 42, 40, 34, 44],
          ["CTOR, cell turnover rate. f, alteration frequency of the gene in ESCC tumours."], freeze_cols=2)
    wb.save(TAB / "TableS1_parameters.xlsx")


# ------------------------------------------------------------------ Table S2
COLS = ["GENES", "EVENTS", "PROBEVENT", "DOMINANCE", "DIVFUNC", "DIVRATE", "DEFUNC", "DEATHRATE",
        "mutationFUNC", "mutation", "MAXDIVFUNC", "MAXDIVRATE", "MEFUNC", "MICROENVIRONMENT"]


def gene_type(r):
    if str(r.GENES).startswith("aging"):
        return "aging"
    return "driver" if r.DIVRATE > 0 else "essential" if r.DIVRATE < 0 else "passenger"


def esitable_sheet(ws, f, title, legend):
    raw = pd.read_csv(f, encoding="utf-8-sig")
    genome = int(raw["genome size"].dropna().iloc[0])
    d = raw[COLS].copy()
    d.insert(1, "Type", d.apply(gene_type, axis=1))
    rows = [tuple(x.item() if hasattr(x, "item") else x for x in r) for r in d.itertuples(index=False)]
    sheet(ws, title, f"{legend} Genome size: {genome:,}.", ["Gene"] + list(d.columns[1:]), rows,
          [11, 10, 13, 10, 10, 8, 9, 8, 10, 12, 8, 11, 11, 8, 16])
    return len(rows)


def table_s2():
    wb = Workbook(); ws = wb.active; ws.title = "ESCC esiTable"
    n = esitable_sheet(ws, ESI / "esiTable_escc.csv",
                       "Supplementary Table S2. esiTables used in the simulations",
                       "ESCC table (Figure 3A, B, D–H; Figure S4): 179 gene-event records (173 genes) plus 4 aging records.")
    first, last = 5, 4 + n
    ws2 = wb.create_sheet("Six-gene esiTable")
    esitable_sheet(ws2, ESI / "esiTable_six_genes.csv", "Six-gene esiTable (multistage model)",
                   "Figure 3C, J, K. Aging records present but not scheduled.")
    ws3 = wb.create_sheet("Six-gene interactions")
    it = pd.read_csv(ESI / "interaction_six_genes.csv").drop(columns=["Link"])
    sheet(ws3, "Six-gene interaction table", "Modifiers applied when both genes of an ordered pair carry an event.",
          list(it.columns), [tuple(x.item() if hasattr(x, "item") else x for x in r) for r in it.itertuples(index=False)],
          [14, 14, 12, 12, 10, 10, 12],
          ["For each partner event already present in the cell, the effect of the new event on division is multiplied by DIVMOD "
           "and its effect on death by DEMOD (DEMOD 1/3 divides the death effect by 3)."])

    ws4 = wb.create_sheet("Summary")
    ref = "'ESCC esiTable'!"
    rng = lambda col: f"{ref}${col}${first}:${col}${last}"  # noqa: E731
    summary = [
        ("Records, all", f"=COUNTA({rng('A')})"),
        ("Driver records", f'=COUNTIF({rng("B")},"driver")'),
        ("Essential records", f'=COUNTIF({rng("B")},"essential")'),
        ("Passenger records", f'=COUNTIF({rng("B")},"passenger")'),
        ("Aging records", f'=COUNTIF({rng("B")},"aging")'),
        ("Gene-event records (all minus aging)", "=B5-B9"),
        ("Oncogenome size (sum of driver PROBEVENT)", f'=SUMIF({rng("B")},"driver",{rng("D")})'),
        ("Passenger PROBEVENT, minimum", f'=_xlfn.MINIFS({rng("D")},{rng("B")},"passenger")'),
        ("Passenger PROBEVENT, maximum", f'=_xlfn.MAXIFS({rng("D")},{rng("B")},"passenger")'),
        ("Records with a mutation-rate modifier (mutation ≠ 0)", f'=COUNTIF({rng("K")},"<>0")'),
    ]
    sheet(ws4, "Summary of the ESCC esiTable", "Counts computed from the 'ESCC esiTable' sheet.",
          ["Quantity", "Value"], summary, [52, 12],
          ["Type: driver, DIVRATE > 0; essential, DIVRATE < 0; passenger, DIVRATE = 0; aging, records aging_1–4.",
           "The mutation-rate modifier is +1 for TP53 and ATM in the file, but esiCancer does not read this field (see column "
           "definitions), so it had no effect in any simulation."])

    ws5 = wb.create_sheet("Column definitions")
    defs = [
        ("Gene", "Gene name; a gene can appear in more than one record when it has more than one event type."),
        ("Type", "Driver, essential, passenger or aging (added for this table; not an esiCancer field)."),
        ("EVENTS", "Type of alteration: mutation, amplification or deletion."),
        ("PROBEVENT", "Relative event weight; chance of being hit by one event = weight ÷ genome size."),
        ("DOMINANCE", "Fraction of the effect expressed when only one allele carries the event (20/20 rule)."),
        ("DIVFUNC, DIVRATE", "How the event changes the probability of division: ADD (added) or MULT (multiplied), and by how much."),
        ("DEFUNC, DEATHRATE", "Same for the probability of death."),
        ("mutationFUNC, mutation", "Intended change of the events per division of the cell. Not read by esiCancer, which "
         "expects the column names 'mutfunc' and 'mut'; it had no effect in any simulation."),
        ("MAXDIVFUNC, MAXDIVRATE", "Change of the maximum number of divisions (TERT +100)."),
        ("MEFUNC, MICROENVIRONMENT", "Effect on the tumour environment; not used (MULT 1)."),
        ("Gene_Before, Gene_After", "Ordered gene pair of the interaction table."),
        ("DIVMOD, DEMOD, MAXDIVMOD", "Modifiers applied when both genes of the pair carry an event."),
    ]
    sheet(ws5, "Column definitions", "Fields of the esiTable and interaction table.", ["Column", "Definition"], defs, [28, 90])
    wb.save(TAB / "TableS2_esiTables.xlsx")


# ------------------------------------------------------------------ Table S3
def table_s3():
    co = pd.read_csv(OUT / "conditions_occurrence_meda.csv")
    co = co[co["T"] == 92]
    series = [
        ("Figure 3A", "A_constant", "ESCC", "16,647,721", "none", "constant", "none", 0.01, 1000, "1–1000",
         {1: "constant_1", 3: "constant_3", 5: "constant_5", 10: "constant_10"}, False),
        ("Figure 3B", "B_increasing", "ESCC", "16,647,721", "none", "increasing (modified esiCancer)", "none", 0.01, 1000, "1–1000",
         {5: "linear_5", 10: "linear_10", 14: "linear_14", 20: "linear_20"}, True),
        ("Figure 3C", "C_multistage", "six-gene", "750", "six-gene", "constant", "none", 0.01, 1000, "1–1000",
         {2: "6_genes_mut_002", 3: "6_genes_mut_003", 4: "6_genes_mut_004", 6: "6_genes_mut_006"}, False),
        ("Figure 3D", "D_constant_aging", "ESCC", "16,647,721", "none", "constant", "60, 70, 80, 90 y", 0.01, 1000, "1–1000",
         {1: "mut_rate_const_001", 3: "mut_rate_const_003", 5: "mut_rate_const_005", 10: "mut_rate_const_010"}, False),
        ("Figure 3E", "E_incr_aging_ctor001", "ESCC", "16,647,721", "none", "increasing (modified esiCancer)", "60, 70, 80, 90 y", 0.01, 1000, "1–1000",
         {5: "mut_rate_PA5_mut_2", 10: "mut_rate_PA10_mut", 14: "mut_rate_PA14_mut", 20: "mut_rate_PA20_mut"}, True),
        ("Figure 3G", "G_incr_aging_ctor002", "ESCC", "16,647,721", "none", "increasing (modified esiCancer)", "60, 70, 80, 90 y", 0.02, 1000, "1–1000",
         {5: "ctor_002_5", 10: "ctor_002_10", 14: "ctor_002_14", 20: "ctor_002_20"}, True),
        ("Figure 3H", "H_incr_aging_ctor02", "ESCC", "16,647,721", "none", "increasing (modified esiCancer)", "60, 70, 80, 90 y", 0.2, 1000, "1–1000",
         {5: "ctor_02_5", 10: "ctor_02_10", 14: "ctor_02_14", 20: "ctor_02_20"}, True),
    ]
    step = {5: 20, 10: 10, 14: 7, 20: 5}
    rows = []
    for fig, key, tab, genome, inter, sched, aging, ctor, maxdiv, seeds, files, incr in series:
        for lvl, prefix in files.items():
            r = co[(co.series == key) & (co.input == lvl)].iloc[0]
            ev = f"1 + floor(age/{step[lvl]} y) (level {lvl})" if incr else str(lvl)
            rows.append([fig, prefix, tab, genome, inter, ev, aging, "none", ctor, maxdiv, seeds,
                         int(r.n_runs), int(r.n_tumour), int(r.n_extinct), None, None, round(r.meda, 1)])
    gs = pd.read_csv(OUT / "germline" / "germline_summary.csv").set_index("level")
    g, c = gs.loc["germline_003"], gs.loc["input_003_to_004"]
    rows += [
        ["Figure 3J, K", "germ_003_non_carrier", "six-gene", "750", "six-gene", "3", "none", "none", 0.01, 1000, "1–3000",
         int(g.n0), int(g.events0), 0, None, None, round(g.meda0, 1)],
        ["Figure 3J, K", "germ_003_carrier", "six-gene", "750", "six-gene", "3", "none", "gene A, 100% of cells, generation 0", 0.01, 1000,
         "1–3000", int(g.n1), int(g.events1), 0, None, None, round(g.meda1, 1)],
        ["Figure 3J, K", "germ_004_non_carrier (= 6_genes_mut_004)", "six-gene", "750", "six-gene", "4", "none", "none", 0.01, 1000,
         "1–1000", int(c.n1), int(c.events1), 0, None, None, round(c.meda1, 1)],
    ]
    # Matched-input turnover (output/turnover_results.txt, Experiment 2); CTOR 0.01 arms are Figure 3A.
    rows += [
        ["Figure S4", "ctor002_e3_g2", "ESCC, genome ×2", "33,295,442", "none", "3", "none", "none", 0.02, 10000, "1–1000", 1000, 442, 0, None, None, 42.4],
        ["Figure S4", "ctor002_e10_g2", "ESCC, genome ×2", "33,295,442", "none", "10", "none", "none", 0.02, 10000, "1–1000", 1000, 833, 0, None, None, 29.0],
        ["Figure S4", "ctor02_e3_g20", "ESCC, genome ×20", "332,954,420", "none", "3", "none", "none", 0.2, 10000, "1–1000", 1000, 58, 942, None, None, 15.3],
        ["Figure S4", "ctor02_e10_g20", "ESCC, genome ×20", "332,954,420", "none", "10", "none", "none", 0.2, 10000, "1–1000", 1000, 238, 762, None, None, 20.7],
        ["Figure 2B", "singles_with_aging_seed1–14", "ESCC", "16,647,721", "none", "20 (constant; recorded in the parameter files)",
         "60, 70, 80, 90 y", "none", 0.01, 1000, "1–14", 14, 6, 6, None, None, None],
    ]
    header = ["Figure", "Condition (file prefix)", "esiTable", "Genome size", "Interaction table", "Events per division",
              "Aging events", "Other manual events", "CTOR", "Max divisions", "Seeds", "Runs", "Tumours by 92 y",
              "Extinctions by 92 y", "Neither by 92 y (censored)", "Cumulative occurrence by 92 y", "MEDA (years)"]
    wb = Workbook(); ws = wb.active; ws.title = "Table S3"
    first = sheet(ws, "Supplementary Table S3. Configuration and outcome of every simulated condition",
                  "Common settings: 1,000 initial cells; probability of division = probability of death = CTOR; stop at "
                  "4,800 generations (92 years) or 10,000 cells; tumour environment and asymmetric mitosis off.",
                  header, rows, [11, 30, 18, 12, 11, 22, 15, 20, 7, 9, 9, 7, 9, 10, 11, 12, 9],
                  ["CTOR, cell turnover rate. MEDA, median age at detection among runs with an esiTumor (Data Analysis).",
                   "Increasing schedule: events per division = 1 + floor(age / step); levels 5, 10, 14, 20 correspond to steps of 20, 10, 7 and 5 years "
                   "(level = 100 / step; the 5-year step reaches 19 events per division at 92 years).",
                   "Figure S4: CTOR 0.01 arms are the Figure 3A runs at 3 and 10 events per division. Genome ×k divides each gene's chance of being hit by k.",
                   "Counts from output/conditions_occurrence_meda.csv (T = 92 y), output/germline/germline_summary.csv and output/turnover_results.txt.",
                   "esiCancer Linux build (ref. 17); increasing schedules with the modification in esicancer/increasing_events_per_division.patch. "
                   "esiTable of each series checked against the gene list of its outputs."],
                  freeze_cols=2)
    for i in range(first, first + len(rows)):
        ws.cell(i, 15, f"=L{i}-M{i}-N{i}")
        ws.cell(i, 16, f"=M{i}/L{i}").number_format = "0.0%"
        for j in (15, 16):
            ws.cell(i, j).font = FONT; ws.cell(i, j).border = BOX
    wb.save(TAB / "TableS3_conditions.xlsx")


# ------------------------------------------------------------------ Table S4
def table_s4():
    header = ["Item"] + v4.COLS[1:]
    wb = Workbook(); ws = wb.active; ws.title = "Table S4"
    rows = [list(r) for r in v4.ROWS]
    rows[0][4] = "2, 3 (3J, K: inherited first hit); S2–S4"
    sheet(ws, "Supplementary Table S4. Endpoints, denominators and calculations for the epidemiological examples and the simulations",
          "Answers reviewer comment 1: event definition, population at risk, denominator, observation period, age range and standardisation for every dataset.",
          header, rows, [22, 45, 45, 45, 45],
          ["BCAC, Breast Cancer Association Consortium; CTOR, cell turnover rate; GBD, Global Burden of Disease; MEDA, median age at detection; "
           "NMSC, non-melanoma skin cancer; PH, proportional hazards; PTV, protein-truncating variant."])
    for c in ws["A"][4:]:
        c.font = BOLD
    wb.save(TAB / "TableS4_endpoints.xlsx")


if __name__ == "__main__":
    table_s1(); table_s2(); table_s3(); table_s4()
    print("written to", TAB)
