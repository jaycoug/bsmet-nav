# =============================================================================
# PHST 681 — Biostatistical Methods Decision System
# app.R — v3: Data input, Help tab, polished design
#
# To run: shiny::runApp("path/to/this/folder")
# Required: install.packages(c("shiny", "bslib", "survival"))
# Optional: install.packages("DescTools")  # for Breslow-Day test
# =============================================================================

library(shiny)
library(bslib)

# =============================================================================
# METHODS DATABASE
# =============================================================================
methods_db <- list(

  cohort = list(
    title = "Cohort (Prospective) Study", icon = "\U0001F52C",
    category = "Study Design", cat_color = "#0077b6",
    description = "Subjects classified by exposure status and followed forward in time to observe disease incidence. All three effect measures \u2014 RD, RR, and OR \u2014 are directly estimable.",
    formula = "RD = p\u0302\u2081 \u2212 p\u0302\u2082    RR = p\u0302\u2081 / p\u0302\u2082    OR = ad / bc",
    notes = "The gold standard for establishing temporal sequence between exposure and disease. Expensive and time-consuming, especially for rare diseases.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      p1 <- a/(a+b); p2 <- cc/(cc+d); RD <- p1-p2; RR <- p1/p2; OR <- (a*d)/(b*cc)
      se_or <- sqrt(1/a+1/b+1/cc+1/d); CI <- exp(log(OR)+c(-1,1)*qnorm(.975)*se_or)
      test_out <- capture.output(print(prop.test(c(a,cc), c(a+b,cc+d), correct=TRUE)))
      paste0("Risk (exposed):   ", round(p1,4), "\nRisk (unexposed): ", round(p2,4),
             "\n\nRisk Difference = ", round(RD,4), "\nRisk Ratio      = ", round(RR,4),
             "\nOdds Ratio      = ", round(OR,4), "\n  95% CI for OR: (", round(CI[1],4), ", ", round(CI[2],4), ")",
             "\n\n", paste(test_out, collapse="\n"))
    },
    r_template = 'mat <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)
p1 <- a/(a+b); p2 <- c/(c+d)
cat("RD =", round(p1-p2, 4), "\\nRR =", round(p1/p2, 4),
    "\\nOR =", round((a*d)/(b*c), 4), "\\n")
prop.test(c(a,c), c(a+b,c+d), correct=TRUE)', r_packages = "base R"),

  case_control = list(
    title = "Case-Control Study", icon = "\U0001F50D",
    category = "Study Design", cat_color = "#0077b6",
    description = "Subjects selected by disease status; exposure assessed retrospectively. Only the OR is directly estimable. If disease is rare (< 10%), OR \u2248 RR.",
    formula = "OR = ad / bc \u2248 RR (rare disease)",
    notes = "Cannot estimate RD or RR directly \u2014 sampling fractions are unknown. The exposure-odds ratio equals the disease-odds ratio regardless of sampling strategy.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      OR <- (a*d)/(b*cc); se <- sqrt(1/a+1/b+1/cc+1/d)
      CI <- exp(log(OR)+c(-1,1)*qnorm(.975)*se)
      paste0("Odds Ratio = ", round(OR,4), "\n  95% CI: (", round(CI[1],4), ", ", round(CI[2],4), ")",
             "\n  ln(OR) = ", round(log(OR),4), "  se(ln OR) = ", round(se,4),
             "\n\nIf disease incidence < 10%, this OR approximates the RR.")
    },
    r_template = 'OR <- (a*d)/(b*c); se <- sqrt(1/a+1/b+1/c+1/d)
CI <- exp(log(OR)+c(-1,1)*qnorm(0.975)*se)
cat("OR =", round(OR,2), "95% CI:", round(CI,2), "\\n")', r_packages = "base R"),

  cross_sectional = list(
    title = "Cross-Sectional Study", icon = "\U0001F4CA",
    category = "Study Design", cat_color = "#0077b6",
    description = "Exposure and disease measured simultaneously. Can estimate prevalence and OR. Cannot establish temporal sequence.",
    formula = "Prevalence ratio = prev\u2081 / prev\u2082    OR = ad / bc",
    notes = "Useful for estimating disease burden and generating hypotheses. Number of cases not fixed in advance.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      p1 <- a/(a+b); p2 <- cc/(cc+d); OR <- (a*d)/(b*cc)
      paste0("Prevalence (exposed):   ", round(p1,4), "\nPrevalence (unexposed): ", round(p2,4),
             "\nPrevalence ratio:       ", round(p1/p2,4), "\nOR = ", round(OR,4))
    },
    r_template = 'prev1 <- a/(a+b); prev2 <- c/(c+d)
cat("PR:", round(prev1/prev2,4), "OR:", round((a*d)/(b*c),4), "\\n")', r_packages = "base R"),

  risk_difference = list(
    title = "Risk Difference (RD)", icon = "\u2796",
    category = "Measures of Effect", cat_color = "#00B894",
    description = "Absolute difference in disease probability. Only from cohort or cross-sectional studies.",
    formula = "RD = p\u0302\u2081 \u2212 p\u0302\u2082\nse(RD) = \u221A(p\u0302\u2081q\u0302\u2081/n\u2081 + p\u0302\u2082q\u0302\u2082/n\u2082)",
    notes = "CI requires: n\u2081p\u0302\u2081(1\u2212p\u0302\u2081) \u2265 5 and n\u2082p\u0302\u2082(1\u2212p\u0302\u2082) \u2265 5.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      p1 <- a/(a+b); p2 <- cc/(cc+d); RD <- p1-p2
      se <- sqrt(p1*(1-p1)/(a+b)+p2*(1-p2)/(cc+d)); CI <- RD+c(-1,1)*qnorm(.975)*se
      paste0("RD = ", round(RD,4), "\nse = ", round(se,4), "\n95% CI: (", round(CI[1],4), ", ", round(CI[2],4), ")")
    },
    r_template = 'p1 <- a/(a+b); p2 <- c/(c+d); RD <- p1-p2
se <- sqrt(p1*(1-p1)/(a+b)+p2*(1-p2)/(c+d))
cat("RD =", round(RD,4), "95% CI:", round(RD+c(-1,1)*qnorm(.975)*se, 4), "\\n")', r_packages = "base R"),

  risk_ratio = list(
    title = "Risk Ratio (RR)", icon = "\u00F7",
    category = "Measures of Effect", cat_color = "#00B894",
    description = "Ratio of disease probabilities. CI on ln(RR) scale via delta method.",
    formula = "RR = p\u0302\u2081/p\u0302\u2082    se(ln RR) = \u221A(q\u0302\u2081/(n\u2081p\u0302\u2081) + q\u0302\u2082/(n\u2082p\u0302\u2082))",
    notes = "RR bounded above by 1/p\u2082. OR preferred for stratified analysis.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      p1 <- a/(a+b); p2 <- cc/(cc+d); RR <- p1/p2
      se <- sqrt((1-p1)/((a+b)*p1)+(1-p2)/((cc+d)*p2))
      CI <- exp(log(RR)+c(-1,1)*qnorm(.975)*se)
      paste0("RR = ", round(RR,4), "\n95% CI: (", round(CI[1],4), ", ", round(CI[2],4), ")")
    },
    r_template = 'p1 <- a/(a+b); p2 <- c/(c+d); RR <- p1/p2
se <- sqrt((1-p1)/((a+b)*p1)+(1-p2)/((c+d)*p2))
CI <- exp(log(RR)+c(-1,1)*qnorm(.975)*se)
cat("RR =", round(RR,4), "95% CI:", round(CI,4), "\\n")', r_packages = "base R"),

  odds_ratio = list(
    title = "Odds Ratio (OR)", icon = "\u2694",
    category = "Measures of Effect", cat_color = "#00B894",
    description = "Ratio of disease odds. Estimable from ANY study design. No upper-bound restriction.",
    formula = "OR = ad/bc    se(ln OR) = \u221A(1/a + 1/b + 1/c + 1/d)",
    notes = "When disease is rare (p < 0.10): OR \u2248 RR. The universal effect measure across study designs.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      OR <- (a*d)/(b*cc); se <- sqrt(1/a+1/b+1/cc+1/d)
      CI <- exp(log(OR)+c(-1,1)*qnorm(.975)*se)
      paste0("OR = ", round(OR,4), "\nln(OR) = ", round(log(OR),4), "  se = ", round(se,4),
             "\n95% CI: (", round(CI[1],4), ", ", round(CI[2],4), ")")
    },
    r_template = 'OR <- (a*d)/(b*c); se <- sqrt(1/a+1/b+1/c+1/d)
CI <- exp(log(OR)+c(-1,1)*qnorm(.975)*se)
cat("OR =", round(OR,4), "95% CI:", round(CI,4), "\\n")', r_packages = "base R"),

  attributable_risk = list(
    title = "Attributable Risk (AR%)", icon = "\u26A0",
    category = "Measures of Effect", cat_color = "#00B894",
    description = "Percentage of disease attributable to the exposure.",
    formula = "AR = 100% \u00D7 p(RR \u2212 1) / [p(RR \u2212 1) + 1]",
    notes = "p = prevalence of risk factor in population. Common exposure + modest RR can produce large AR.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'RR <- p1/p2; p <- prevalence_of_exposure
AR <- 100*p*(RR-1)/(p*(RR-1)+1)
cat("AR% =", round(AR,2), "\\n")', r_packages = "base R"),

  chi_square_2x2 = list(
    title = "Chi-Square Test (2\u00D72)", icon = "\u03C7\u00B2",
    category = "Hypothesis Tests", cat_color = "#FDCB6E",
    description = "Large-sample test of H\u2080: p\u2081 = p\u2082. Requires all expected cell counts \u2265 5.",
    formula = "X\u00B2 = \u03A3(|O \u2212 E| \u2212 0.5)\u00B2 / E  ~  \u03C7\u00B2\u2081 (Yates')",
    notes = "ALWAYS check expected counts first. If any < 5 \u2192 Fisher's exact test.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      mat <- matrix(c(a,b,cc,d), nrow=2, byrow=TRUE)
      ex <- chisq.test(mat, correct=FALSE)$expected; ok <- all(ex>=5)
      res <- if(ok) chisq.test(mat, correct=TRUE) else fisher.test(mat)
      paste0("Expected cell counts:\n", paste(capture.output(print(round(ex,1))), collapse="\n"),
             "\nAll \u2265 5? ", ok, if(!ok) "\n\u26A0 Using Fisher's exact test instead.\n" else "\n",
             "\n", paste(capture.output(print(res)), collapse="\n"))
    },
    r_template = 'mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
ex <- chisq.test(mat, correct=FALSE)$expected
cat("Expected:\\n"); print(round(ex,1))
if(all(ex>=5)) chisq.test(mat) else fisher.test(mat)', r_packages = "base R"),

  fisher_exact = list(
    title = "Fisher's Exact Test", icon = "\U0001F3AF",
    category = "Hypothesis Tests", cat_color = "#FDCB6E",
    description = "Exact test for 2\u00D72 tables when expected cells < 5. Based on the hypergeometric distribution.",
    formula = "P(X = a) = C(a+b,a) \u00B7 C(c+d,c) / C(n,a+c)",
    notes = "Preferred when any expected cell count < 5. Exact p-value without chi-square approximation.",
    has_input = TRUE, input_type = "2x2",
    r_code = function(a, b, cc, d) {
      mat <- matrix(c(a,b,cc,d), nrow=2, byrow=TRUE)
      paste(capture.output(print(fisher.test(mat))), collapse="\n")
    },
    r_template = 'mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
fisher.test(mat)', r_packages = "base R"),

  chi_square_trend = list(
    title = "Chi-Square Test for Trend", icon = "\u2197",
    category = "Hypothesis Tests", cat_color = "#FDCB6E",
    description = "For 2\u00D7k tables with ordered exposure categories. Tests linear dose-response.",
    formula = "X\u00B2_trend ~ \u03C7\u00B2\u2081", notes = "Scores assigned to ordered categories. 1 df \u2014 more powerful than heterogeneity test.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'events <- c(a1,a2,a3); totals <- c(n1,n2,n3)
prop.trend.test(events, totals, score=c(1,2,3))', r_packages = "base R"),

  chi_square_heterogeneity = list(
    title = "Chi-Square for Heterogeneity", icon = "\u2262",
    category = "Hypothesis Tests", cat_color = "#FDCB6E",
    description = "For unordered categories. Tests independence without assuming ordering.",
    formula = "X\u00B2 ~ \u03C7\u00B2_(R\u22121)(C\u22121)",
    notes = "Validity: \u2264 1/5 of expected counts < 5, none < 1.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'mat <- matrix(c(...), nrow=R, byrow=TRUE)\nchisq.test(mat)', r_packages = "base R"),

  mantel_haenszel = list(
    title = "Mantel-Haenszel Test", icon = "\u2696",
    category = "Confounding & Mantel-Haenszel", cat_color = "#E17055",
    description = "Tests E\u2192D association controlling for categorical confounder C across stratified 2\u00D72 tables.",
    formula = "X\u00B2_MH = [|\u03A3O\u1D62 \u2212 \u03A3E\u1D62| \u2212 0.5]\u00B2 / \u03A3Var(O\u1D62) ~ \u03C7\u00B2\u2081\nOR_MH = \u03A3(a\u1D62d\u1D62/n\u1D62) / \u03A3(b\u1D62c\u1D62/n\u1D62)",
    notes = "Requires: D binary, E binary, C single categorical. If E/C continuous or multiple confounders \u2192 logistic regression.",
    has_input = TRUE, input_type = "stratified",
    r_code = function(a1,b1,c1,d1,a2,b2,c2,d2) {
      tab <- array(c(a1,c1,b1,d1,a2,c2,b2,d2), dim=c(2,2,2),
                   dimnames=list(Exposure=c("Exp","Unexp"), Disease=c("D+","D-"), Stratum=c("S1","S2")))
      or1 <- (a1*d1)/(b1*c1); or2 <- (a2*d2)/(b2*c2)
      res <- capture.output(print(mantelhaen.test(tab, correct=TRUE)))
      paste0("Stratum-specific ORs:\n  Stratum 1: OR = ", round(or1,4),
             "\n  Stratum 2: OR = ", round(or2,4), "\n\n", paste(res, collapse="\n"))
    },
    r_template = 'tab <- array(c(a1,c1,b1,d1, a2,c2,b2,d2), dim=c(2,2,2),
  dimnames=list(Exposure=c("Exp","Unexp"), Disease=c("D+","D-"),
                Stratum=c("S1","S2")))
mantelhaen.test(tab, correct=TRUE)', r_packages = "base R"),

  homogeneity_test = list(
    title = "Test for Homogeneity of ORs", icon = "\u2263",
    category = "Confounding & Mantel-Haenszel", cat_color = "#E17055",
    description = "Tests H\u2080: OR\u2081 = OR\u2082 = \u2026 = OR\u2096. THE critical fork: determines common OR vs stratum-specific ORs.",
    formula = "X\u00B2_homog ~ \u03C7\u00B2_(K\u22121)", notes = "If rejected \u2192 effect modification (report separate ORs). If not rejected \u2192 estimate common OR.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'library(DescTools)\nBreslowDayTest(tab)  # tab = 2x2xK array', r_packages = "DescTools"),

  effect_modification = list(
    title = "Effect Modification", icon = "\u26A1",
    category = "Confounding & Mantel-Haenszel", cat_color = "#E17055",
    description = "Homogeneity rejected: C genuinely modifies E\u2192D. Report stratum-specific ORs \u2014 a pooled OR would be misleading.",
    formula = "Report OR\u2081, OR\u2082, \u2026 separately with individual CIs",
    notes = "This is a substantive biological finding, NOT a nuisance to be adjusted away.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'for(k in 1:K) {\n  t <- tab[,,k]; or <- (t[1,1]*t[2,2])/(t[1,2]*t[2,1])\n  se <- sqrt(sum(1/t))\n  ci <- exp(log(or)+c(-1,1)*qnorm(.975)*se)\n  cat("Stratum",k,": OR =",round(or,2),"CI:",round(ci,2),"\\n")\n}',
    r_packages = "base R"),

  logistic_regression = list(
    title = "Multiple Logistic Regression", icon = "\U0001F4C8",
    category = "Logistic Regression", cat_color = "#E84393",
    description = "Controls for many confounders simultaneously. Handles continuous and categorical predictors. e\u1D5D = adjusted OR.",
    formula = "log(p/(1\u2212p)) = \u03B1 + \u03B2\u2081x\u2081 + \u2026 + \u03B2\u2096x\u2096\nBinary/dummy x: OR = e^\u03B2    Continuous x: OR = e^(\u0394\u00B7\u03B2)    Categorical: k\u22121 dummies",
    notes = "RECOGNITION CUES:\n\u2022 Problem says 'logistic regression was run' or gives coefficients + SEs\n\u2022 Intercept + predictor coefficients (not a 2\u00D72 table)\n\u2022 Dummy-coded groups \u2192 \u03B2 = log(OR) vs. reference group\n\u2022 Continuous predictor \u2192 \u03B2 = log(OR) per 1-unit increase; use \u0394 for other increments\n\nADJUSTED vs. CRUDE ORs:\n\u2022 Adding covariates (e.g., BMI) to the model produces ADJUSTED ORs\n\u2022 If adjusted OR \u2260 crude OR, the covariate is a confounder\n\u2022 Each \u03B2 is interpreted holding all other predictors constant",
    has_input = TRUE, input_type = "logistic_coef",
    r_code = function(beta, se, delta=1) {
      lnOR <- beta * delta; OR <- exp(lnOR)
      se_delta <- se * abs(delta)
      CI <- exp(lnOR + c(-1,1) * qnorm(0.975) * se_delta)
      z <- beta / se; pval <- 2 * (1 - pnorm(abs(z)))
      header <- if (delta == 1) {
        paste0("Regression coefficient (\u03B2): ", round(beta, 4),
               "\nStandard error:             ", round(se, 4),
               "\n\n\u2192 OR = e^\u03B2 = e^(", round(beta,4), ") = ", round(OR, 4))
      } else {
        paste0("Regression coefficient (\u03B2): ", round(beta, 4),
               "\nStandard error:             ", round(se, 4),
               "\nUnit change (\u0394):            ", delta,
               "\n\n\u2192 OR per ", delta, "-unit change = e^(\u0394\u00B7\u03B2) = e^(", delta, " \u00D7 ", round(beta,4), ") = e^(", round(lnOR,4), ") = ", round(OR, 4))
      }
      paste0(header,
             "\n\n\u2500\u2500\u2500 95% CI CALCULATION (work on ln scale) \u2500\u2500\u2500",
             "\n  ln(OR) \u00B1 1.96 \u00D7 SE = ", round(lnOR,4), " \u00B1 ", round(qnorm(0.975)*se_delta,4),
             "\n  ln-scale CI = (", round(lnOR - qnorm(0.975)*se_delta, 4), ", ", round(lnOR + qnorm(0.975)*se_delta, 4), ")",
             if (delta != 1) paste0("\n  (* SE for \u0394-unit OR = |\u0394| \u00D7 SE = ", delta, " \u00D7 ", round(se,4), " = ", round(se_delta,4), ")") else "",
             "\n\n\u2500\u2500\u2500 ANSWER: 95% CI FOR OR \u2500\u2500\u2500",
             "\n  Exponentiate both bounds:",
             "\n  \u2192\u2192\u2192  95% CI = (", round(CI[1], 4), ", ", round(CI[2], 4), ")  \u2190\u2190\u2190",
             "\n\nWald test: z = \u03B2/SE = ", round(z, 4),
             "\n  p-value = ", format.pval(pval, digits=4),
             "\n\nInterpretation: This is the ", if(delta==1) "OR" else paste0("OR per ", delta, "-unit change"),
             " ADJUSTED for all other predictors in the model.")
    },
    r_template = 'model <- glm(disease ~ group + bmi, data=df,
             family=binomial(link="logit"))
summary(model)
exp(coef(model))           # ORs (per 1-unit for continuous)
exp(confint.default(model)) # 95% CIs

# --- From coefficients: dummy predictor ---
beta <- 0.186; se <- 0.060
OR <- exp(beta)  # vs reference group
CI <- exp(beta + c(-1,1)*qnorm(0.975)*se)

# --- From coefficients: continuous predictor (e.g., per 5-unit change) ---
beta_bmi <- 0.107; se_bmi <- 0.004; delta <- 5
OR_delta <- exp(delta * beta_bmi)
CI_delta <- exp(delta*beta_bmi + c(-1,1)*qnorm(.975)*abs(delta)*se_bmi)
cat("OR per", delta, "units =", round(OR_delta,4),
    "  95% CI:", round(CI_delta,4), "\\n")', r_packages = "base R (glm)"),

  logistic_predict = list(
    title = "Logistic Regression: Predicted Probability", icon = "\U0001F3AF",
    category = "Logistic Regression", cat_color = "#E84393",
    description = "Plug covariate values into a fitted logistic model to predict the probability of disease for a specific individual profile.",
    formula = "logit = \u03B1 + \u03B2\u2081x\u2081 + \u03B2\u2082x\u2082 + \u2026\np = 1 / (1 + e^(\u2212logit))",
    notes = "RECOGNITION CUES:\n\u2022 'Predict the probability of [disease] for a [specific person]'\n\u2022 You are given ALL coefficients (intercept + predictors) and specific covariate values\n\u2022 For dummy predictors: x = 1 if the person is in that group, x = 0 otherwise\n\u2022 For continuous predictors: x = the actual value (e.g., BMI = 23)\n\nSTEPS:\n1. Multiply each \u03B2 by its x value\n2. Sum everything including the intercept = logit\n3. Convert: p = 1/(1 + e^(\u2212logit))",
    has_input = TRUE, input_type = "logistic_predict",
    r_code = function(intercept, b1, x1, b2, x2, b3, x3, b4, x4, b5, x5) {
      terms <- list()
      logit <- intercept
      terms[[1]] <- paste0("  Intercept:     ", sprintf("%+.4f", intercept))
      pairs <- list(list(b1,x1), list(b2,x2), list(b3,x3), list(b4,x4), list(b5,x5))
      for (i in seq_along(pairs)) {
        bi <- pairs[[i]][[1]]; xi <- pairs[[i]][[2]]
        if (!is.null(bi) && !is.na(bi) && !is.null(xi) && !is.na(xi)) {
          prod_i <- bi * xi; logit <- logit + prod_i
          terms[[length(terms)+1]] <- paste0("  \u03B2", i, " \u00D7 x", i, ":    ",
            sprintf("%+.4f", bi), " \u00D7 ", sprintf("%.4g", xi), " = ", sprintf("%+.4f", prod_i))
        }
      }
      prob <- 1 / (1 + exp(-logit)); odds <- exp(logit)
      paste0("LOGIT CALCULATION:\n",
             paste(terms, collapse="\n"),
             "\n  ", paste(rep("\u2500", 34), collapse=""),
             "\n  logit = ", round(logit, 4),
             "\n\nCONVERSION TO PROBABILITY:",
             "\n  p = 1 / (1 + e^(\u2212logit))",
             "\n  p = 1 / (1 + e^(", sprintf("%+.4f", -logit), "))",
             "\n  p = 1 / (1 + ", round(exp(-logit), 4), ")",
             "\n  p = ", round(prob, 4),
             "\n\n  Odds = e^(logit) = ", round(odds, 4),
             "\n\nPredicted probability = ", round(prob, 4), " (", round(prob*100, 2), "%)")
    },
    r_template = '# From fitted model:
    # predict(model, newdata=data.frame(AA=0, Hispanic=1, BMI=23), type="response")

# From coefficients by hand:
logit <- -4.277 + 0.009*0 + 0.186*1 + 0.107*23
p <- 1 / (1 + exp(-logit))
cat("logit =", round(logit,4), "\\np =", round(p,4), "\\n")', r_packages = "base R"),

  confounding_assessment = list(
    title = "Confounding Assessment (Crude vs. Adjusted)", icon = "\U0001F50E",
    category = "Logistic Regression", cat_color = "#E84393",
    description = "Compare the crude OR (without covariate) to the adjusted OR (with covariate) to determine if confounding is present and its direction.",
    formula = "Positive confounder: crude OR farther from 1 than adjusted OR\n   (crude overstates the association)\nNegative confounder: crude OR closer to 1 than adjusted OR\n   (crude understates the association)\nNo confounding: crude OR \u2248 adjusted OR",
    notes = "RECOGNITION CUES:\n\u2022 Two models compared: one without a covariate, one with it\n\u2022 Question asks 'is X a confounder?' or 'what type of confounder?'\n\u2022 Must assume covariate is NOT in the causal pathway (otherwise it's a mediator, not confounder)\n\nDECISION RULE:\n\u2022 |ln(OR_crude)| > |ln(OR_adj)| \u2192 Positive confounder (inflates)\n\u2022 |ln(OR_crude)| < |ln(OR_adj)| \u2192 Negative confounder (masks)\n\u2022 |ln(OR_crude)| \u2248 |ln(OR_adj)| \u2192 Not a confounder\n\nCommon threshold: > 10% change in OR suggests meaningful confounding.",
    has_input = TRUE, input_type = "confounding_assess",
    r_code = function(or_crude, or_adj) {
      ln_crude <- log(or_crude); ln_adj <- log(or_adj)
      pct_change <- 100 * (or_crude - or_adj) / or_adj
      abs_change <- abs(pct_change)
      dist_crude <- abs(ln_crude); dist_adj <- abs(ln_adj)
      direction <- if (abs_change < 10) {
        "NOT A CONFOUNDER"
      } else if (dist_crude > dist_adj) {
        "POSITIVE CONFOUNDER"
      } else {
        "NEGATIVE CONFOUNDER"
      }
      explain <- if (abs_change < 10) {
        "The crude and adjusted ORs are similar (< 10% change), suggesting the covariate does not confound this association."
      } else if (dist_crude > dist_adj) {
        "The crude OR is farther from 1.0 than the adjusted OR. The covariate was INFLATING the apparent association \u2014 making it look stronger than it truly is. This is positive confounding."
      } else {
        "The crude OR is closer to 1.0 than the adjusted OR. The covariate was MASKING the true association \u2014 making it look weaker than it truly is. This is negative confounding."
      }
      paste0("CRUDE vs. ADJUSTED COMPARISON:",
             "\n  Crude OR:    ", round(or_crude, 4), "  (ln = ", round(ln_crude, 4), ", distance from null = ", round(dist_crude, 4), ")",
             "\n  Adjusted OR: ", round(or_adj, 4), "  (ln = ", round(ln_adj, 4), ", distance from null = ", round(dist_adj, 4), ")",
             "\n\n  % change: ", sprintf("%+.1f%%", pct_change),
             "\n  |% change| = ", round(abs_change, 1), "%",
             if (abs_change >= 10) " > 10% threshold \u2192 meaningful" else " < 10% threshold \u2192 minimal",
             "\n\n\u2500\u2500\u2500 VERDICT \u2500\u2500\u2500",
             "\n  \u2192 ", direction,
             "\n\n", explain)
    },
    r_template = '# Crude model (no covariate)
crude_model <- glm(disease ~ exposure, family=binomial, data=df)
OR_crude <- exp(coef(crude_model)["exposure"])

# Adjusted model (with covariate)
adj_model <- glm(disease ~ exposure + bmi, family=binomial, data=df)
OR_adj <- exp(coef(adj_model)["exposure"])

# Compare
cat("Crude OR:", round(OR_crude,4), "\\nAdjusted OR:", round(OR_adj,4),
    "\\n% change:", round(100*(OR_crude-OR_adj)/OR_adj, 1), "\\n")', r_packages = "base R"),

  kaplan_meier = list(
    title = "Kaplan-Meier Estimator", icon = "\U0001F4C9",
    category = "Survival Analysis", cat_color = "#00B894",
    description = "Nonparametric survival curve S(t). Handles censored data.",
    formula = "\u015C(t) = \u220F[1 \u2212 d\u2C7C/n\u2C7C]  for t\u2C7C \u2264 t", notes = "d\u2C7C = events at time t\u2C7C, n\u2C7C = at-risk just before.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'library(survival)\nkm <- survfit(Surv(time, status) ~ group, data=df)\nplot(km, col=c("blue","red"), xlab="Time", ylab="S(t)")', r_packages = "survival"),

  log_rank = list(
    title = "Log-Rank Test", icon = "\u2694",
    category = "Survival Analysis", cat_color = "#00B894",
    description = "Compares two+ survival curves. More powerful than comparing at specific time points.",
    formula = "X\u00B2_LR = (O \u2212 E)\u00B2 / Var ~ \u03C7\u00B2\u2081", notes = "Similar in spirit to MH test.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'library(survival)\nsurvdiff(Surv(time, status) ~ group, data=df)', r_packages = "survival"),

  cox_ph = list(
    title = "Cox Proportional Hazards", icon = "\u2693",
    category = "Survival Analysis", cat_color = "#00B894",
    description = "Semi-parametric regression for survival data. Hazard ratio e^\u03B2 interpreted like an OR.",
    formula = "h(t) = h\u2080(t) \u00B7 exp(\u03B2\u2081x\u2081 + \u2026 + \u03B2\u2096x\u2096)",
    notes = "Must verify PH assumption with cox.zph(). If violated \u2192 stratified Cox or time-varying coefficients.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'library(survival)\ncox <- coxph(Surv(time, status) ~ x1 + x2, data=df)\nsummary(cox)\ncox.zph(cox)', r_packages = "survival"),

  meta_analysis = list(
    title = "Meta-Analysis", icon = "\U0001F4DA",
    category = "Advanced Methods", cat_color = "#636E72",
    description = "Combines log-ORs across studies. Fixed vs random effects.",
    formula = "\u03BC\u0302 = \u03A3w\u1D62*y\u1D62 / \u03A3w\u1D62*", notes = "Test homogeneity with Q. If rejected \u2192 random effects.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'y <- c(ln_OR_1, ln_OR_2, ...); w <- 1/c(var1, var2, ...)
mu <- sum(w*y)/sum(w); Q <- sum(w*(y-mu)^2)
cat("Combined OR:", round(exp(mu),2), "Q =", round(Q,2), "\\n")', r_packages = "base R"),

  equivalence_study = list(
    title = "Equivalence Studies", icon = "\u2248",
    category = "Advanced Methods", cat_color = "#636E72",
    description = "Establish therapeutic equivalence within margin \u03B4.",
    formula = "Reject non-equivalence if CI \u2282 (\u2212\u03B4, +\u03B4)", notes = "\u03B4 chosen before the study on clinical grounds.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 'delta <- 0.10; diff <- p1-p2\nse <- sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)\nCI <- diff+c(-1,1)*qnorm(.975)*se\ncat("Equivalent?", CI[1]>-delta & CI[2]<delta, "\\n")', r_packages = "base R"),

  crossover = list(
    title = "Cross-Over Design", icon = "\u21C4",
    category = "Advanced Methods", cat_color = "#636E72",
    description = "Each subject receives both treatments with washout. Subject = own control.",
    formula = "", notes = "Watch for carry-over and period effects.",
    has_input = FALSE, input_type = "none", r_code = NULL,
    r_template = 't.test(outcome ~ period, data=df, paired=TRUE)', r_packages = "base R")
)

# =============================================================================
# CSS
# =============================================================================
css <- '
@import url("https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap");
body { background: linear-gradient(160deg, #0a0f1a, #0d1b2a 40%, #1b2838); color: #dde; min-height: 100vh; font-family: "Plus Jakarta Sans", sans-serif !important; font-size: 15px; }

/* ---- Navbar ---- */
.navbar { background: rgba(10,15,26,.97) !important; border-bottom: 2px solid rgba(0,180,216,.2); backdrop-filter: blur(12px); padding: 6px 12px !important; }
.navbar-brand { font-weight: 800 !important; font-size: 1.2rem !important; background: linear-gradient(135deg,#00b4d8,#55efc4); -webkit-background-clip: text; -webkit-text-fill-color: transparent; letter-spacing: -.3px; }
.navbar .nav-link { color: #8899aa !important; font-weight: 600; font-size: .88rem; border-radius: 12px; padding: 8px 16px !important; transition: all .25s ease; margin: 0 2px; font-family: "Plus Jakarta Sans", sans-serif !important; }
.navbar .nav-link:hover { color: #55efc4 !important; background: rgba(0,180,216,.1); }
.navbar .nav-link.active { color: #fff !important; background: rgba(0,180,216,.18); border: 1px solid rgba(0,180,216,.3); }

/* ---- Well / Navigator Panel ---- */
.well { background: rgba(13,27,42,.9) !important; border: 1px solid rgba(0,180,216,.18) !important; border-radius: 18px !important; padding: 22px 24px !important; }
.well h5 { color: #00b4d8; font-weight: 700; font-size: 1.05rem; letter-spacing: -.2px; }
.well .control-label { color: #00b4d8 !important; font-weight: 700; font-size: .88rem; }
.radio label { color: #bcc8d4 !important; font-size: .9rem; font-weight: 500; }
.radio input[type="radio"] { accent-color: #00b4d8; }
.well hr { border-color: rgba(0,180,216,.12); margin: 14px 0; }
.gate-title { color: #fdcb6e !important; font-weight: 700; font-size: .88rem; }

/* ---- Method Card ---- */
.method-card { background: rgba(13,27,42,.75); border: 1px solid rgba(0,180,216,.14); border-radius: 20px; padding: 28px 32px; backdrop-filter: blur(8px); animation: fadeIn .3s ease; }
@keyframes fadeIn { from {opacity:0; transform:translateY(6px)} to {opacity:1; transform:translateY(0)} }
.method-card h2 { font-size: 1.55rem; font-weight: 800; color: #fff; margin-bottom: 4px; letter-spacing: -.3px; }
.method-icon { font-size: 1.75rem; margin-right: 10px; }
.cat-badge { display:inline-block; padding:4px 13px; border-radius:20px; font-size:.72rem; font-weight:700; letter-spacing:.5px; text-transform:uppercase; margin-bottom:10px; }
.method-desc { color:#b0bec5; font-size:.95rem; line-height:1.65; margin:14px 0; }

/* ---- Formula & Notes ---- */
.formula-box { background:rgba(0,180,216,.06); border:1px solid rgba(0,180,216,.22); border-radius:12px; padding:14px 18px; font-family:"JetBrains Mono",monospace; font-size:.85rem; color:#dde; white-space:pre-wrap; margin:14px 0; }
.notes-box { background:rgba(253,203,110,.06); border-left:3px solid #fdcb6e; border-radius:0 10px 10px 0; padding:12px 16px; font-size:.84rem; color:#e0d4a8; margin:14px 0; white-space:pre-wrap; }
.notes-box strong { color:#fdcb6e; }

/* ---- Input Section ---- */
.input-section { background:rgba(0,184,148,.05); border:1px solid rgba(0,184,148,.18); border-radius:14px; padding:20px 22px; margin:18px 0; }
.input-section h5 { color:#55efc4; font-weight:700; font-size:.95rem; margin-bottom:12px; }
.input-section .form-control { background:rgba(10,15,26,.85) !important; border:1px solid rgba(0,184,148,.25) !important; color:#eee !important; border-radius:10px; font-family:"JetBrains Mono",monospace; font-size:.88rem; }
.input-section label { color:#8fd8c0 !important; font-size:.8rem; font-weight:600; }

/* ---- Results ---- */
.results-box { background:rgba(10,15,26,.92); border:1px solid rgba(0,180,216,.2); border-radius:12px; padding:16px 20px; font-family:"JetBrains Mono",monospace; font-size:.84rem; color:#55efc4; white-space:pre-wrap; overflow-x:auto; margin-top:12px; }

/* ---- Code Section ---- */
.code-section { margin-top:18px; }
.code-section h5 { color:#00b4d8; font-weight:700; font-size:.92rem; }
.code-box { background:#080e18; border:1px solid rgba(0,180,216,.15); border-radius:12px; padding:14px 18px; font-family:"JetBrains Mono",monospace; font-size:.8rem; color:#d0d6dc; white-space:pre-wrap; overflow-x:auto; margin-top:8px; }
.pkg-badge { display:inline-block; background:rgba(0,180,216,.1); color:#00b4d8; padding:2px 10px; border-radius:12px; font-size:.7rem; font-weight:700; margin-left:8px; }

/* ---- Buttons ---- */
.btn-run { background:linear-gradient(135deg,#00b894,#00cec9) !important; color:#0a0f1a !important; border:none !important; font-weight:700 !important; border-radius:10px !important; padding:9px 24px !important; font-size:.9rem !important; font-family:"Plus Jakarta Sans",sans-serif !important; letter-spacing:.2px; }
.btn-run:hover { box-shadow:0 4px 16px rgba(0,184,148,.35); transform:translateY(-1px); }
.btn-outline-secondary { border-color:rgba(0,180,216,.25) !important; color:#00b4d8 !important; background:transparent !important; border-radius:10px !important; font-family:"Plus Jakarta Sans",sans-serif !important; font-weight:600 !important; }
.btn-outline-secondary:hover { background:rgba(0,180,216,.08) !important; }

/* ---- Placeholders ---- */
.placeholder-panel { text-align:center; padding:80px 20px; }
.placeholder-panel h3 { color:#00b4d8; font-weight:800; font-size:1.4rem; }
.placeholder-panel p { color:#6a7b8b; max-width:480px; margin:10px auto; font-size:.95rem; }

/* ---- Empty State ---- */
.empty-state { text-align:center; padding:60px 20px; }
.empty-state h4 { color:#00b4d8; font-weight:700; font-size:1.15rem; }
.empty-state p { color:#6a7b8b; max-width:400px; margin:8px auto; font-size:.92rem; }

/* ---- Help Tab ---- */
.help-hero { text-align:center; padding:40px 20px 28px; }
.help-hero h2 { font-size:1.85rem; font-weight:800; background:linear-gradient(135deg,#00b4d8,#55efc4); -webkit-background-clip:text; -webkit-text-fill-color:transparent; }
.help-hero p { color:#6a7b8b; margin-top:6px; font-size:.95rem; }
.help-card { background:rgba(13,27,42,.75); border:1px solid rgba(0,180,216,.12); border-radius:16px; padding:22px 24px; margin-bottom:14px; transition:all .25s ease; }
.help-card:hover { border-color:rgba(0,180,216,.35); transform:translateY(-2px); box-shadow:0 4px 14px rgba(0,180,216,.08); }
.help-card h5 { color:#00b4d8; font-weight:700; font-size:.95rem; }
.help-card p { color:#8899aa; font-size:.88rem; margin-top:4px; }

/* ---- Dropdowns ---- */
.selectize-input { background:rgba(10,15,26,.9) !important; border-color:rgba(0,180,216,.25) !important; color:#ddd !important; border-radius:10px !important; }
.selectize-dropdown { background:rgba(13,27,42,.97) !important; border-color:rgba(0,180,216,.2) !important; border-radius:10px !important; }
.selectize-dropdown-content .option { color:#bcc8d4 !important; }
.selectize-dropdown-content .option.active { background:rgba(0,180,216,.2) !important; color:#fff !important; }

/* ---- Hint Boxes ---- */
.hint-box { border-radius:12px; padding:12px 14px; margin-bottom:12px; }
.hint-box p { margin:0; font-weight:600; font-size:.84rem; }
.hint-box ul { margin:6px 0 0 0; padding-left:18px; font-size:.8rem; }
.hint-box li { margin-bottom:2px; }
.hint-cohort { background:rgba(0,180,216,.08); border:1px solid rgba(0,180,216,.25); }
.hint-cohort p, .hint-cohort li { color:#5cc8e0; }
.hint-cc { background:rgba(0,184,148,.08); border:1px solid rgba(0,184,148,.25); }
.hint-cc p, .hint-cc li { color:#55efc4; }
.hint-cs { background:rgba(253,203,110,.08); border:1px solid rgba(253,203,110,.25); }
.hint-cs p, .hint-cs li { color:#fad390; }
.hint-reg { background:rgba(232,67,147,.08); border:1px solid rgba(232,67,147,.25); }
.hint-reg p, .hint-reg li { color:#e84393; }
'

# =============================================================================
# UI
# =============================================================================
ui <- navbarPage(
  title = "Biostatistical Methods",
  header = tags$head(
    tags$link(href = "https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap", rel = "stylesheet"),
    tags$style(HTML(css))
  ),

  # ---- Epidemiologic Studies ----
  tabPanel("Epidemiologic Studies",
    fluidRow(style = "padding: 20px 12px;",
      column(4,
        wellPanel(
          h5("\U0001F9ED Scenario navigator"),
          radioButtons("entry", "What are you starting with?",
            choices = c("Categorical outcome (disease yes/no)" = "categorical",
                        "Regression output (coefficients + SEs)" = "regression_output",
                        "Person-time data" = "persontime",
                        "Browse all methods" = "browse"),
            selected = character(0)),
          conditionalPanel("input.entry=='regression_output'", hr(),
            div(class="hint-box hint-reg",
              tags$p("\U0001F4A1 You likely have logistic regression output if you see:"),
              tags$ul(
                tags$li("An intercept + predictor coefficients"),
                tags$li("Standard errors for each coefficient"),
                tags$li("The problem mentions 'regression model was run'"),
                tags$li("Dummy-coded groups (reference category noted)"))),
            radioButtons("reg_type", "What type of regression?",
              choices = c("Logistic regression (binary outcome)" = "logistic",
                          "Linear regression (continuous outcome)" = "linear_regression"),
              selected = character(0))),
          conditionalPanel("input.entry=='regression_output' && input.reg_type=='logistic'", hr(),
            radioButtons("logistic_task", "What do you need to compute?",
              choices = c("OR from a single coefficient" = "logistic_regression",
                          "Predicted probability for a specific person" = "logistic_predict",
                          "Assess confounding (crude vs. adjusted OR)" = "confounding_assessment"),
              selected = character(0))),
          conditionalPanel("input.entry=='categorical'", hr(),
            radioButtons("study_type", "Study design?",
              choices = c("Cohort (prospective)"="cohort","Case-control"="case_control","Cross-sectional"="cross_sectional"),
              selected = character(0)),
            conditionalPanel("input.study_type=='cohort'",
              div(class="hint-box hint-cohort",
                tags$p("\U0001F4A1 You likely have a cohort study if you see:"),
                tags$ul(
                  tags$li("Subjects classified by EXPOSURE, then followed for disease"),
                  tags$li("Words like 'followed', 'prospective', 'incidence', 'risk'"),
                  tags$li("You can count disease cases in both exposed and unexposed"),
                  tags$li("All 3 measures available: RD, RR, and OR")))),
            conditionalPanel("input.study_type=='case_control'",
              div(class="hint-box hint-cc",
                tags$p("\U0001F4A1 You likely have a case-control study if you see:"),
                tags$ul(
                  tags$li("Subjects selected by DISEASE status (cases vs. controls)"),
                  tags$li("Exposure assessed retrospectively (looking back)"),
                  tags$li("Words like 'cases', 'controls', 'retrospective', 'matched'"),
                  tags$li("Only the OR is directly estimable (\u2248 RR if disease rare)")))),
            conditionalPanel("input.study_type=='cross_sectional'",
              div(class="hint-box hint-cs",
                tags$p("\U0001F4A1 You likely have a cross-sectional study if you see:"),
                tags$ul(
                  tags$li("Exposure and disease measured at the SAME TIME"),
                  tags$li("Words like 'prevalence', 'survey', 'at a single point in time'"),
                  tags$li("Cannot establish which came first (no temporal sequence)"),
                  tags$li("Prevalence ratio and OR estimable, but NOT true RR"))))),
          conditionalPanel("input.entry=='categorical' && input.study_type!=''", hr(),
            radioButtons("table_type", "Data structure?",
              choices = c("2 \u00D7 2 table"="2x2","2 \u00D7 k ordered"="2xk_ord","2 \u00D7 k unordered"="2xk_unord","R \u00D7 C table"="rxc"),
              selected = character(0))),
          conditionalPanel("input.entry=='categorical' && input.table_type!=''", hr(),
            radioButtons("confounding", "Confounding suspected?",
              choices = c("No \u2014 crude analysis"="no","Yes \u2014 control for confounder(s)"="yes"),
              selected = character(0))),
          conditionalPanel("input.confounding=='yes'", hr(),
            div(class="gate-title", "\u26A0 Mantel-Haenszel eligibility"),
            radioButtons("mh_eligible", NULL,
              choices = c("E binary, C single categorical \u2192 MH"="mh_yes",
                          "E/C continuous or multiple confounders \u2192 Logistic reg."="mh_no"),
              selected = character(0))),
          conditionalPanel("input.mh_eligible=='mh_yes'", hr(),
            radioButtons("homog_result", "Homogeneity of ORs?",
              choices = c("Fail to reject \u2192 Common OR"="homogeneous",
                          "Reject \u2192 Effect modification"="heterogeneous",
                          "Show me the test"="show_test"),
              selected = character(0))),
          conditionalPanel("input.entry=='persontime'", hr(),
            radioButtons("rates_constant", "Rates constant over time?",
              choices = c("Yes \u2014 rate methods"="yes","No \u2014 survival analysis"="no"),
              selected = character(0))),
          conditionalPanel("input.entry=='persontime' && input.rates_constant=='no'", hr(),
            radioButtons("surv_goal", "Goal?",
              choices = c("Estimate survival curve"="km","Compare two curves"="logrank","Multiple risk factors"="coxph"),
              selected = character(0))),
          conditionalPanel("input.entry=='browse'", hr(),
            selectInput("browse_method", "Select method:",
              choices = c("-- Select method --"="",
                "Cohort study"="cohort","Case-control"="case_control","Cross-sectional"="cross_sectional",
                "Risk Difference"="risk_difference","Risk Ratio"="risk_ratio","Odds Ratio"="odds_ratio","Attributable Risk"="attributable_risk",
                "Chi-square (2\u00D72)"="chi_square_2x2","Fisher's exact"="fisher_exact","Chi-square trend"="chi_square_trend","Chi-square heterogeneity"="chi_square_heterogeneity",
                "Mantel-Haenszel"="mantel_haenszel","Homogeneity of ORs"="homogeneity_test","Effect modification"="effect_modification",
                "Multiple logistic regression"="logistic_regression",
                "Logistic: Predicted probability"="logistic_predict",
                "Logistic: Confounding assessment"="confounding_assessment",
                "Kaplan-Meier"="kaplan_meier","Log-rank test"="log_rank","Cox PH"="cox_ph",
                "Meta-analysis"="meta_analysis","Equivalence studies"="equivalence_study","Cross-over design"="crossover"))),
          hr(),
          actionButton("reset_btn", "\u21BA  Reset selections", class="btn-sm btn-outline-secondary", width="100%")
        )
      ),
      column(8, uiOutput("method_display"))
    )
  ),

  # ---- Placeholders ----
  tabPanel("Person-Time & Survival", div(class="placeholder-panel",
    h3("\U0001F4C5 Person-Time & Survival Analysis"),
    p("Incidence rates, Kaplan-Meier, log-rank, Cox PH \u2014 next section to build."))),
  tabPanel("Categorical Data", div(class="placeholder-panel",
    h3("\U0001F4CA Categorical Data Analysis"),
    p("Binomial tests, contingency tables, McNemar's, Kappa \u2014 semester 1."))),
  tabPanel("Multisample Inference", div(class="placeholder-panel",
    h3("\U0001F4CB Multisample Inference"),
    p("ANOVA, Kruskal-Wallis, multiple comparisons \u2014 semester 1."))),
  tabPanel("Regression & Correlation", div(class="placeholder-panel",
    h3("\U0001F4C8 Regression & Correlation"),
    p("Linear regression, partial correlation, Fisher's z \u2014 semester 1."))),

  # ---- Help ----
  tabPanel("\u2753 Help",
    div(style = "max-width:620px; margin:0 auto; padding:20px;",
      div(class="help-hero", h2("Not sure where to start?"),
        p("Answer a few questions about your problem and we'll point you to the right method.")),
      div(class="help-card",
        h5("\U0001F3AF What kind of outcome are you measuring?"),
        p("Binary (disease yes/no), categorical (multiple groups), continuous (blood pressure), or time-to-event (survival)?")),
      div(class="help-card",
        h5("\U0001F522 How many groups are you comparing?"),
        p("One sample vs known value, two independent groups, two paired/matched groups, or more than two?")),
      div(class="help-card",
        h5("\u23F1 Is the timing of events important?"),
        p("Variable follow-up with censoring? You likely need survival analysis or person-time methods.")),
      div(class="help-card",
        h5("\U0001F500 Are there confounders to control for?"),
        p("Categorical or continuous? One or many? This determines MH vs logistic regression.")),
      div(class="help-card",
        h5("\U0001F50D What's the study design?"),
        p("Cohort, case-control, cross-sectional, or randomized trial? This constrains which effect measures are estimable.")),
      tags$p(style="text-align:center; color:#555; margin-top:20px; font-size:.82rem;",
        em("Interactive guided wizard coming in a future update."))
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  current_method <- reactive({
    if (!is.null(input$entry) && input$entry == "browse") {
      m <- input$browse_method; if (!is.null(m) && m != "") return(m); return(NULL)
    }
    if (!is.null(input$entry) && input$entry == "regression_output") {
      if (is.null(input$reg_type) || input$reg_type == "") return(NULL)
      if (input$reg_type == "linear_regression") return("linear_regression")
      if (input$reg_type == "logistic") {
        if (!is.null(input$logistic_task) && input$logistic_task != "") return(input$logistic_task)
        return(NULL)
      }
      return(NULL)
    }
    if (!is.null(input$entry) && input$entry == "persontime") {
      if (!is.null(input$rates_constant) && input$rates_constant == "no" && !is.null(input$surv_goal))
        return(switch(input$surv_goal, km="kaplan_meier", logrank="log_rank", coxph="cox_ph", NULL))
      return(NULL)
    }
    if (!is.null(input$entry) && input$entry == "categorical") {
      if (is.null(input$study_type) || input$study_type == "") return(NULL)
      if (is.null(input$table_type) || input$table_type == "") return(input$study_type)
      tm <- switch(input$table_type, "2x2"="chi_square_2x2","2xk_ord"="chi_square_trend",
                   "2xk_unord"="chi_square_heterogeneity","rxc"="chi_square_heterogeneity",NULL)
      if (is.null(input$confounding) || input$confounding == "") return(tm)
      if (input$confounding == "no") return(tm)
      if (input$confounding == "yes") {
        if (is.null(input$mh_eligible) || input$mh_eligible == "") return(NULL)
        if (input$mh_eligible == "mh_no") return("logistic_regression")
        if (input$mh_eligible == "mh_yes") {
          if (is.null(input$homog_result) || input$homog_result == "") return("mantel_haenszel")
          return(switch(input$homog_result, homogeneous="mantel_haenszel",
                        heterogeneous="effect_modification", show_test="homogeneity_test", NULL))
        }
      }
    }
    NULL
  })

  analysis_result <- reactiveVal(NULL)
  observeEvent(current_method(), { analysis_result(NULL) })

  output$method_display <- renderUI({
    mid <- current_method()
    if (is.null(mid)) return(div(class="empty-state",
      h4("\U0001F9ED Select a scenario to begin"),
      p("Use the navigator to walk the decision tree, or browse methods directly.")))
    m <- methods_db[[mid]]
    if (is.null(m)) return(p("Method not found."))
    tagList(div(class="method-card",
      div(class="cat-badge", style=paste0("background:",m$cat_color,"22;color:",m$cat_color,";border:1px solid ",m$cat_color,"44;"), m$category),
      h2(span(class="method-icon", m$icon), m$title),
      div(class="method-desc", m$description),
      if (!is.null(m$formula) && m$formula != "") div(class="formula-box", m$formula),
      if (!is.null(m$notes) && m$notes != "") div(class="notes-box", tags$strong("Key notes: "), m$notes),
      if (!is.null(m$has_input) && m$has_input && m$input_type == "2x2")
        div(class="input-section",
          h5("\U0001F4E5 Enter your 2 \u00D7 2 table"),
          fluidRow(
            column(6, fluidRow(
              column(6, numericInput("cell_a", "a (Exp, D+)", value=NULL, min=0, width="100%")),
              column(6, numericInput("cell_b", "b (Exp, D\u2212)", value=NULL, min=0, width="100%")))),
            column(6, fluidRow(
              column(6, numericInput("cell_c", "c (Unexp, D+)", value=NULL, min=0, width="100%")),
              column(6, numericInput("cell_d", "d (Unexp, D\u2212)", value=NULL, min=0, width="100%"))))),
          actionButton("run_2x2", "\u25B6  Run analysis", class="btn-run")),
      if (!is.null(m$has_input) && m$has_input && m$input_type == "stratified")
        div(class="input-section",
          h5("\U0001F4E5 Enter stratified 2 \u00D7 2 tables"),
          tags$p(style="color:#55efc4;font-size:.82rem;font-weight:700;", "Stratum 1"),
          fluidRow(
            column(3, numericInput("s1_a","a\u2081",value=NULL,min=0,width="100%")),
            column(3, numericInput("s1_b","b\u2081",value=NULL,min=0,width="100%")),
            column(3, numericInput("s1_c","c\u2081",value=NULL,min=0,width="100%")),
            column(3, numericInput("s1_d","d\u2081",value=NULL,min=0,width="100%"))),
          tags$p(style="color:#55efc4;font-size:.82rem;font-weight:700;margin-top:6px;", "Stratum 2"),
          fluidRow(
            column(3, numericInput("s2_a","a\u2082",value=NULL,min=0,width="100%")),
            column(3, numericInput("s2_b","b\u2082",value=NULL,min=0,width="100%")),
            column(3, numericInput("s2_c","c\u2082",value=NULL,min=0,width="100%")),
            column(3, numericInput("s2_d","d\u2082",value=NULL,min=0,width="100%"))),
          actionButton("run_strat", "\u25B6  Run MH analysis", class="btn-run")),
      if (!is.null(m$has_input) && m$has_input && m$input_type == "logistic_coef")
        div(class="input-section",
          h5("\U0001F4E5 Enter regression coefficient"),
          div(class="hint-box hint-reg",
            tags$p("\u03B2 = regression coefficient for the predictor of interest"),
            tags$ul(
              tags$li("For dummy-coded predictors: \u03B2 gives log(OR) vs. reference group, adjusted for all other predictors."),
              tags$li("For continuous predictors: \u03B2 gives log(OR) per 1-unit increase. Use \u0394 for other increments."))),
          radioButtons("pred_type", "Predictor type:",
            choices = c("Binary / dummy-coded (e.g., ethnic group)" = "dummy",
                        "Continuous (e.g., BMI, age)" = "continuous"),
            selected = "dummy", inline = FALSE),
          fluidRow(
            column(4, numericInput("coef_beta", "\u03B2 (coefficient)", value=NULL, step=0.001, width="100%")),
            column(4, numericInput("coef_se", "SE (standard error)", value=NULL, min=0, step=0.001, width="100%")),
            column(4, conditionalPanel("input.pred_type=='continuous'",
              numericInput("coef_delta", "\u0394 (unit change)", value=1, min=0.001, step=1, width="100%")))),
          actionButton("run_logistic", "\u25B6  Compute OR", class="btn-run")),
      if (!is.null(m$has_input) && m$has_input && m$input_type == "logistic_predict")
        div(class="input-section",
          h5("\U0001F4E5 Enter full model coefficients + covariate values"),
          div(class="hint-box hint-reg",
            tags$p("Enter the intercept, then each predictor\u2019s \u03B2 and the person\u2019s x value."),
            tags$ul(
              tags$li("For dummy predictors: x = 1 if person is in that group, x = 0 otherwise."),
              tags$li("For continuous: x = actual value (e.g., BMI = 23). Leave unused rows blank."))),
          fluidRow(column(6, numericInput("pred_intercept", "\u03B1 (intercept)", value=NULL, step=0.001, width="100%")),
                   column(6)),
          tags$p(style="color:#55efc4;font-size:.82rem;font-weight:700;margin-top:4px;", "Predictors"),
          fluidRow(column(2, tags$label(style="color:#8899aa;font-size:.78rem;", "")),
                   column(5, tags$label(style="color:#8899aa;font-size:.78rem;", "\u03B2 (coefficient)")),
                   column(5, tags$label(style="color:#8899aa;font-size:.78rem;", "x (value)"))),
          fluidRow(column(2, tags$label(style="color:#55efc4;font-size:.82rem;margin-top:8px;", "x\u2081")),
                   column(5, numericInput("pred_b1", NULL, value=NULL, step=0.001, width="100%")),
                   column(5, numericInput("pred_x1", NULL, value=NULL, step=0.001, width="100%"))),
          fluidRow(column(2, tags$label(style="color:#55efc4;font-size:.82rem;margin-top:8px;", "x\u2082")),
                   column(5, numericInput("pred_b2", NULL, value=NULL, step=0.001, width="100%")),
                   column(5, numericInput("pred_x2", NULL, value=NULL, step=0.001, width="100%"))),
          fluidRow(column(2, tags$label(style="color:#55efc4;font-size:.82rem;margin-top:8px;", "x\u2083")),
                   column(5, numericInput("pred_b3", NULL, value=NULL, step=0.001, width="100%")),
                   column(5, numericInput("pred_x3", NULL, value=NULL, step=0.001, width="100%"))),
          fluidRow(column(2, tags$label(style="color:#55efc4;font-size:.82rem;margin-top:8px;", "x\u2084")),
                   column(5, numericInput("pred_b4", NULL, value=NULL, step=0.001, width="100%")),
                   column(5, numericInput("pred_x4", NULL, value=NULL, step=0.001, width="100%"))),
          fluidRow(column(2, tags$label(style="color:#55efc4;font-size:.82rem;margin-top:8px;", "x\u2085")),
                   column(5, numericInput("pred_b5", NULL, value=NULL, step=0.001, width="100%")),
                   column(5, numericInput("pred_x5", NULL, value=NULL, step=0.001, width="100%"))),
          actionButton("run_predict", "\u25B6  Predict probability", class="btn-run")),
      if (!is.null(m$has_input) && m$has_input && m$input_type == "confounding_assess")
        div(class="input-section",
          h5("\U0001F50E Enter crude and adjusted ORs"),
          div(class="hint-box hint-reg",
            tags$p("Compare two models: one without the suspected confounder, one with it."),
            tags$ul(
              tags$li("Crude OR = from model WITHOUT the covariate."),
              tags$li("Adjusted OR = from model WITH the covariate."),
              tags$li("Compute each OR first using the \u2018OR from a single coefficient\u2019 tool, then enter both here."))),
          fluidRow(
            column(6, numericInput("or_crude", "Crude OR (no covariate)", value=NULL, min=0.001, step=0.01, width="100%")),
            column(6, numericInput("or_adj", "Adjusted OR (with covariate)", value=NULL, min=0.001, step=0.01, width="100%"))),
          actionButton("run_confound", "\u25B6  Assess confounding", class="btn-run")),
      uiOutput("results_panel"),
      div(class="code-section",
        h5("\U0001F4DD R code template", span(class="pkg-badge", m$r_packages %||% "base R")),
        div(class="code-box", m$r_template))
    ))
  })

  observeEvent(input$run_2x2, {
    mid <- current_method(); m <- methods_db[[mid]]
    req(m, m$r_code, input$cell_a, input$cell_b, input$cell_c, input$cell_d)
    tryCatch({
      analysis_result(m$r_code(input$cell_a, input$cell_b, input$cell_c, input$cell_d))
    }, error = function(e) analysis_result(paste("Error:", e$message)))
  })

  observeEvent(input$run_strat, {
    mid <- current_method(); m <- methods_db[[mid]]
    req(m, m$r_code, input$s1_a, input$s1_b, input$s1_c, input$s1_d,
        input$s2_a, input$s2_b, input$s2_c, input$s2_d)
    tryCatch({
      analysis_result(m$r_code(input$s1_a, input$s1_b, input$s1_c, input$s1_d,
                               input$s2_a, input$s2_b, input$s2_c, input$s2_d))
    }, error = function(e) analysis_result(paste("Error:", e$message)))
  })

  observeEvent(input$run_logistic, {
    mid <- current_method(); m <- methods_db[[mid]]
    req(m, m$r_code, input$coef_beta, input$coef_se)
    delta <- if (!is.null(input$pred_type) && input$pred_type == "continuous" && !is.null(input$coef_delta)) input$coef_delta else 1
    tryCatch({
      analysis_result(m$r_code(input$coef_beta, input$coef_se, delta))
    }, error = function(e) analysis_result(paste("Error:", e$message)))
  })

  observeEvent(input$run_predict, {
    mid <- current_method(); m <- methods_db[[mid]]
    req(m, m$r_code, input$pred_intercept)
    tryCatch({
      analysis_result(m$r_code(input$pred_intercept,
        input$pred_b1, input$pred_x1, input$pred_b2, input$pred_x2,
        input$pred_b3, input$pred_x3, input$pred_b4, input$pred_x4,
        input$pred_b5, input$pred_x5))
    }, error = function(e) analysis_result(paste("Error:", e$message)))
  })

  observeEvent(input$run_confound, {
    mid <- current_method(); m <- methods_db[[mid]]
    req(m, m$r_code, input$or_crude, input$or_adj)
    tryCatch({
      analysis_result(m$r_code(input$or_crude, input$or_adj))
    }, error = function(e) analysis_result(paste("Error:", e$message)))
  })

  output$results_panel <- renderUI({
    res <- analysis_result()
    if (is.null(res)) return(NULL)
    div(class="results-box", res)
  })

  observeEvent(input$reset_btn, {
    analysis_result(NULL)
    for (id in c("entry","study_type","table_type","confounding","mh_eligible","homog_result","rates_constant","surv_goal","reg_type","pred_type","logistic_task"))
      updateRadioButtons(session, id, selected = character(0))
    updateSelectInput(session, "browse_method", selected = "")
  })
}

shinyApp(ui, server)
