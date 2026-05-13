## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, 
  fig.height = 5
)


## ----setup_single-------------------------------------------------------------
library(ConfMatrix)

# 1. Setup audit data (The position in the vector marks the priority)
# Pos 1: Accepted, Pos 2: Below Std, Pos 3: Discrepant, Pos 4: Incompatible
observed_counts <- list(c(36, 8, 5, 1))
target_probs    <- list(c(0.75, 0.15, 0.08, 0.02))
compliance_labels <- c("Accepted", "Below_Standard", "Discrepant", "Incompatible")

# 2. Create the audit object
audit_single <- QCCS$new(
  Vectors = observed_counts, 
  Prob = target_probs,
  ClassNames = compliance_labels,
  Source = "Structural Model Audit",
  ID = "STR-001"
)

# 3. Run the Exact Test
exact_result <- audit_single$Exact.test(a = 0.05)
print(exact_result)


## ----multi_audit--------------------------------------------------------------
# Data for three different disciplines
obs_arch   <- c(40, 6, 3, 1) 
obs_struct <- c(42, 5, 2, 1) 
obs_mep    <- c(30, 12, 6, 2) 

# Shared quality targets for all disciplines
project_targets <- c(0.80, 0.12, 0.06, 0.02)

# Create a multi-audit object
multi_audit <- QCCS$new(
  Vectors = list(obs_arch, obs_struct, obs_mep),
  Prob = list(project_targets, project_targets, project_targets),
  ClassNames = compliance_labels,
  ID = "BIM_GLOBAL_PROJECT"
)

# Individual Chi-squared tests (with Bonferroni correction)
individual_results <- multi_audit$Ji.test()
print(individual_results) # Results for MEP

# Global Stability Test for the entire project
global_result <- multi_audit$JiGlobal.test()
print(global_result)


## ----comparison_plot----------------------------------------------------------
# Extract proportions for the chart
proportions <- do.call(rbind, lapply(multi_audit$Vectors, function(x) x/sum(x)))
rownames(proportions) <- c("Architecture", "Structural", "MEP")
colnames(proportions) <- multi_audit$ClassNames

# Grouped barplot (Order is preserved by vector position)
barplot(t(proportions), beside = TRUE, 
        col = c("#2c3e50", "#18bc9c", "#f39c12", "#e74c3c"),
        main = "Technical Compliance by Discipline",
        ylab = "Proportion",
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", cex = 0.8))

# Red dashed line shows the 80% target for 'Accepted' elements
abline(h = 0.80, col = "red", lty = 2, lwd = 1.5)


