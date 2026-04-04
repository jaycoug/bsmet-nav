# Biostatistical Methods Decision System

An interactive R Shiny application that guides you through the complete decision-making workflow for biostatistical analysis — from identifying your study design to selecting the correct test, interpreting results, and running analyses live.

Built as a study and reference tool for the two-semester **PHST 681 Biostatistical Methods** sequence.

---

## Features

**Scenario navigator** — Walk through a step-by-step decision tree that mirrors the analytical workflow taught in class. Each decision gate asks a question about your data and routes you to the correct method.

**Live analysis** — Enter cell counts from a 2×2 table (or stratified tables for Mantel-Haenszel) and run the analysis directly in the app. Results are computed in R and displayed instantly.

**Decision gates** — Key routing checkpoints are built into the flow:
- Sample size check (χ² vs Fisher's exact)
- Causal pathway check (should you control for this confounder?)
- MH eligibility (binary/categorical vs continuous variables)
- Homogeneity of ORs (confounding vs effect modification)
- OR interpretation (binary vs continuous vs categorical predictors)

**R code templates** — Every method includes a copy-ready R code template so you can reproduce the analysis in your own scripts or R Markdown files.

**Browse mode** — Skip the decision tree and jump directly to any method from a dropdown menu.

---

## Sections

| Tab | Status | Coverage |
|-----|--------|----------|
| Epidemiologic Studies | ✅ Active | Study design, measures of effect, contingency table tests, confounding/MH, logistic regression, survival, meta-analysis, equivalence, cross-over |
| Person-Time & Survival | 🔲 Placeholder | Incidence rates, Kaplan-Meier, log-rank, Cox PH |
| Categorical Data | 🔲 Placeholder | Binomial tests, 2×2 tables, R×C tables, McNemar's, Kappa (semester 1) |
| Multisample Inference | 🔲 Placeholder | ANOVA, Kruskal-Wallis, multiple comparisons (semester 1) |
| Regression & Correlation | 🔲 Placeholder | Linear regression, partial correlation, Fisher's z (semester 1) |
| Help | 🔲 Stub | Guided wizard for users who aren't sure where to start |

---

## Installation

### Prerequisites

- **R** (≥ 4.0)
- **RStudio** (recommended, not required)

### Required packages

```r
install.packages(c("shiny", "bslib"))
```

### Optional packages

```r
install.packages("survival")    # for Kaplan-Meier, log-rank, Cox PH
install.packages("DescTools")   # for Breslow-Day test (homogeneity of ORs)
```

---

## Usage

### From RStudio

1. Open `app.R` in RStudio
2. Click **Run App**

### From the R console

```r
shiny::runApp("path/to/biostat_app")
```

### Navigating the app

1. **Choose an entry point** — categorical outcome, person-time data, or browse directly
2. **Answer each decision gate** — the sidebar reveals the next question based on your previous answer
3. **View the method card** — description, formula, key notes, and R code template
4. **Enter your data** (where available) — fill in cell counts and click **Run analysis**
5. **Copy the R template** to use in your own R Markdown homework files

---

## Project structure

```
biostat_app/
├── app.R          # Single-file Shiny application (UI + server + methods database)
└── README.md      # This file
```

The entire app lives in a single `app.R` file. The `methods_db` list at the top is the content database — every method is defined as a named list with title, category, description, formula, notes, R code function (for live analysis), and R code template (for copy-paste). To add a new method, add an entry to `methods_db` and wire it into the navigator's `conditionalPanel` logic.

---

## Extending the app

### Adding a new method

Add an entry to `methods_db`:

```r
new_method = list(
  title = "Method Name",
  icon = "🔧",
  category = "Category Name",
  cat_color = "#HEX",
  description = "What it does.",
  formula = "The formula",
  notes = "Key caveats.",
  has_input = TRUE,          # FALSE if no data input panel
  input_type = "2x2",        # "2x2", "stratified", or "none"
  r_code = function(a, b, cc, d) {
    # Compute and return a string of results
    paste0("Result: ", round(some_value, 4))
  },
  r_template = '# Copyable R code template',
  r_packages = "base R"
)
```

### Adding a new section tab

1. Add a `tabPanel()` to the `navbarPage()` in the UI
2. Create a new navigator sidebar with `conditionalPanel()` gates
3. Wire the server logic to resolve the selected method

---

## Roadmap

- [ ] Complete Epidemiologic Studies section with all homework/quiz-relevant methods
- [ ] Build Person-Time & Survival module
- [ ] Port semester 1 content (Categorical Data, Multisample Inference, Regression & Correlation)
- [ ] Wire up the Help tab as an interactive guided wizard
- [ ] Add data upload (CSV/file) support for logistic regression and survival analysis nodes
- [ ] Add printable summary / export-to-Rmd feature

---

## License

This is a personal educational tool. Not intended for clinical or production use.

---

## Acknowledgments

Course materials from PHST 681 Biostatistical Methods (University of Louisville), based on Rosner, *Fundamentals of Biostatistics*.
