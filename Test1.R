library("readr")
library("dplyr")
library("lubridate")
library("lmtest")
library("sandwich")
library("rdrobust")
library("modelsummary")
library("ggplot2")
library("showtext")
# =========================================================
# RDiT analysis: comments per user-day around downvote removal
# =========================================================
# This script implements a Regression Discontinuity in Time (RDiT) design.
# We (1) build a user-day panel of comment counts, (2) define a running variable
# (days relative to the cutoff) and a post indicator, (3) estimate a local-linear
# RDiT with day-clustered SEs, (4) run a one-sided test for a positive jump, and
# (5) estimate heterogeneous effects for new users and top users, including plots.

font_add(
  family = "Times New Roman",
  regular = "/Library/Fonts/Times New Roman.ttf",
  bold = "/Library/Fonts/Times New Roman Bold.ttf",
  italic = "/Library/Fonts/Times New Roman Italic.ttf",
  bolditalic = "/Library/Fonts/Times New Roman Bold Italic.ttf"
)

showtext_auto()
font_families()
theme_set(
  theme_minimal(base_family = "Times New Roman") +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text  = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14)
    )
)
# ---------- Load data ----------
df <- read_csv("data/Data_final.csv")
cutoff <- ymd_hm("2025-11-23 19:30", tz = "Europe/Zurich")
bw <- 7

# ---------- Clean + core time variables ----------
# Keep only the columns needed for the analysis, rename them for convenience,
# parse comment timestamps into Zurich time, and create a date-level variable
# ("day") used for aggregation and clustering.
df_clean <- df |>
  dplyr::select(`comment id`, userId, createdAt, `before change`, `user created`) |>
  rename(
    comment_id    = `comment id`,
    user_id       = userId,
    comment_date  = createdAt,
    before_change = `before change`,
    membership    = `user created`
  ) |>
  mutate(
    comment_date = ymd_hms(comment_date, tz = "Europe/Zurich"),
    day = as_date(comment_date)
  )

# Diagnostic: how many unique users are observed in the cleaned data?
length(unique(df_clean$user_id))

# ---------- User-day outcome ----------
# Collapse to a user-day panel where the outcome is the number of comments per
# user per day ("comments"). Define the running variable t_days as the distance
# (in days) from each day's midpoint (12:00) to the cutoff timestamp. Create
# post=1 for observations on/after the cutoff. Restrict to a symmetric window
# of +/- bw days around the cutoff to implement a local comparison.
df_userday <- df_clean |>
  group_by(user_id, day) |>
  summarise(comments = n(), .groups = "drop") |>
  mutate(
    day_mid = ymd_hms(paste(day, "12:00:00"), tz = "Europe/Zurich"),
    t_days  = as.numeric(difftime(day_mid, cutoff, units = "days")),
    post    = as.integer(t_days >= 0)
  ) |>
  filter(t_days >= -bw, t_days <= bw)

# ---------- Diagnostic RD plot ----------
# Visual check for a discontinuity at the cutoff. rdplot bins the running
# variable and overlays a polynomial fit (here p=1 for local linear), helping
# assess whether a "jump" at t_days=0 looks plausible.
rdplot(
  y = df_userday$comments,
  x = df_userday$t_days,
  c = 0,
  p = 1,
  title = "RDiT diagnostic: comments per user-day around downvote removal",
  x.label = "Days relative to cutoff (2025-11-23 19:30)",
  y.label = "Comments per user per day"
)


# ---------- Main RDiT model + clustered SEs by day ----------
# Estimate the canonical local-linear RDiT:
#   comments ~ t_days + post + t_days:post
# where:
# - post captures the level discontinuity ("jump") at the cutoff
# - t_days captures the pre-cutoff time trend
# - t_days:post allows the post-cutoff time trend to differ (slope change)
# Standard errors are clustered by day to allow correlation within each calendar day.
m_rdit <- lm(comments ~ t_days + post + t_days:post, data = df_userday)
V_day  <- vcovCL(m_rdit, cluster = ~ day, data = df_userday)
ct_main <- coeftest(m_rdit, vcov = V_day)
print(ct_main)

# LaTeX regression table printed to console (copy/paste into Overleaf).
# NOTE: By default, modelsummary may use tabularray environments; if you want
# classic LaTeX tabular output, set output="latex_tabular".
modelsummary(
  m_rdit,
  vcov = V_day,
  stars = TRUE,
  statistic = "({std.error})",
  output = "latex"
)

# ---------- One-sided test for H1: post > 0 ----------
# Hypothesis: the discontinuity at the cutoff is positive (commenting increases).
# We compute a one-sided p-value using df = (# unique days - 1), consistent with
# day-level clustering.
tau_hat <- ct_main["post", "Estimate"]
t_stat  <- ct_main["post", "t value"]
G <- length(unique(df_userday$day))
p_one_sided <- pt(t_stat, df = G - 1, lower.tail = FALSE)
cat("\nH1 (one-sided) tau_hat =", tau_hat, " p_one_sided =", p_one_sided, "\n")

# ---------- Raw + fitted RDiT plot (ggplot version) ----------

# Prediction grid
grid_main <- data.frame(
  t_days = seq(-bw, bw, by = 0.1)
)
grid_main$post <- as.integer(grid_main$t_days >= 0)
grid_main$pred <- predict(m_rdit, newdata = grid_main)

p_rdit_main <- ggplot(df_userday, aes(x = t_days, y = comments)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_line(
    data = grid_main,
    aes(y = pred),
    linewidth = 1.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Days relative to cutoff",
    y = "Comments per user-day",
    caption = "Data: Republik magazine, n = 884"
  ) +
  theme_minimal()

p_rdit_main +
  theme(text = element_text(family = "Times New Roman"))
ggsave(
  filename = "rdit_fitted_lines.pdf",
  plot = p_rdit_main + theme(text = element_text(family = "Times New Roman", size = 18)),
  width = 12,
  height = 8
)

# =========================================================
# HETEROGENEITY SETUP (FIXED: 1 membership per user)
# =========================================================
# We examine whether the discontinuity differs for:
# (1) new users (accounts created shortly before the cutoff)
# (2) top users (high pre-cutoff activity)
# The key idea is to create user-level indicators and merge them onto the user-day
# panel, then interact them with post (and optionally post-slope).

# ---------- New users (created within X months pre-cutoff) ----------
# Define "new user" as membership created within the last new_user_months before cutoff.
# Parse membership timestamps and collapse to one row per user using the earliest
# observed membership date (min) to avoid duplicates/inconsistencies.
new_user_months <- 3

df_clean <- df_clean |>
  mutate(membership = ymd_hms(membership, tz = "Europe/Zurich"))

user_info <- df_clean |>
  group_by(user_id) |>
  summarise(membership = min(membership, na.rm = TRUE), .groups = "drop") |>
  mutate(
    is_new_user = as.integer(
      membership >= (cutoff %m-% months(new_user_months)) &
        membership < cutoff
    )
  ) |>
  select(user_id, is_new_user)

# Merge new-user indicator onto the user-day panel; treat missing as 0 (not new).
df_userday <- df_userday |>
  left_join(user_info, by = "user_id") |>
  mutate(is_new_user = ifelse(is.na(is_new_user), 0L, is_new_user))

# ---------- Top users (top X% by pre-cutoff commenting volume) ----------
# Compute each user's total number of comments before the cutoff within the analysis
# window and classify the top top_share share as "top users". This captures whether
# highly active users respond differently at the cutoff.
top_share <- 0.05

user_activity <- df_userday |>
  filter(t_days < 0) |>
  group_by(user_id) |>
  summarise(pre_comments = sum(comments), .groups = "drop")

thr <- quantile(user_activity$pre_comments, probs = 1 - top_share, na.rm = TRUE)

user_top <- user_activity |>
  mutate(is_top_user = as.integer(pre_comments >= thr)) |>
  select(user_id, is_top_user)

# Merge top-user indicator; treat missing as 0 (not top).
df_userday <- df_userday |>
  left_join(user_top, by = "user_id") |>
  mutate(is_top_user = ifelse(is.na(is_top_user), 0L, is_top_user))

# Quick checks: counts by subgroup indicator
print(df_userday |> count(is_new_user))
print(df_userday |> count(is_top_user))

# =========================================================
# HETEROGENEOUS-EFFECT RDiT MODELS (clustered by day)
# =========================================================
# We extend the baseline RDiT by interacting subgroup indicators with post and the
# post-slope term. In each model:
# - post is the jump for the reference group (e.g., established or non-top users)
# - post:group is the *difference* in the jump for the subgroup
# - t_days:post:group is the *difference* in the post-cutoff slope change

# --- New users heterogeneity ---
m_het_new <- lm(
  comments ~ t_days + post + t_days:post +
    post:is_new_user + t_days:post:is_new_user,
  data = df_userday
)
V_new <- vcovCL(m_het_new, cluster = ~ day, data = df_userday)
ct_new <- coeftest(m_het_new, vcov = V_new)
cat("\n--- Heterogeneity: New users ---\n")
print(ct_new)
modelsummary(
  m_het_new,
  vcov = V_new,
  stars = TRUE,
  statistic = "({std.error})",
  output = "latex"
)

# Implied jumps:
# - jump for established users = coefficient on post
# - jump for new users        = post + post:is_new_user
jump_old <- ct_new["post", "Estimate"]
jump_new <- jump_old + ct_new["post:is_new_user", "Estimate"]
cat("Jump (non-new users):", jump_old, "\n")
cat("Jump (new users):    ", jump_new, "\n")

# --- Top users heterogeneity ---
m_het_top <- lm(
  comments ~ t_days + post + t_days:post +
    post:is_top_user + t_days:post:is_top_user,
  data = df_userday
)
V_top <- vcovCL(m_het_top, cluster = ~ day, data = df_userday)
ct_top <- coeftest(m_het_top, vcov = V_top)
cat("\n--- Heterogeneity: Top users ---\n")
print(ct_top)

modelsummary(
  m_het_top,
  vcov = V_top,
  stars = TRUE,
  statistic = "({std.error})",
  output = "latex"
)

# Implied jumps:
# - jump for non-top users = post
# - jump for top users     = post + post:is_top_user
jump_non_top <- ct_top["post", "Estimate"]
jump_top     <- jump_non_top + ct_top["post:is_top_user", "Estimate"]
cat("Jump (non-top users):", jump_non_top, "\n")
cat("Jump (top users):    ", jump_top, "\n")

# --- Optional combined model (both at once) ---
# Includes both tenure and activity interactions simultaneously, allowing each
# dimension to shift the jump and post-trend.
m_het_both <- lm(
  comments ~ t_days + post + t_days:post +
    post:is_new_user + t_days:post:is_new_user +
    post:is_top_user + t_days:post:is_top_user,
  data = df_userday
)
V_both <- vcovCL(m_het_both, cluster = ~ day, data = df_userday)
cat("\n--- Heterogeneity: Combined model ---\n")
print(coeftest(m_het_both, vcov = V_both))

# ---------- Plot: heterogeneity by tenure ----------
# Create a prediction grid for established vs. new users across the running
# variable, compute model predictions, and plot subgroup mean outcomes plus
# the fitted lines. This visualizes whether the discontinuity differs by tenure.
grid_new <- expand.grid(
  t_days = seq(-bw, bw, by = 0.1),
  is_new_user = c(0, 1)
)
grid_new$post <- as.integer(grid_new$t_days >= 0)
grid_new$pred <- predict(m_het_new, newdata = grid_new)

plot_new <- ggplot(df_userday, aes(x = t_days, y = comments, color = factor(is_new_user))) +
  stat_summary(fun = mean, geom = "point", size = 2, alpha = 0.7) +
  geom_line(data = grid_new, aes(y = pred), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "User type",
    values = c("0" = "#1b9e77", "1" = "#d95f02"),
    labels = c("Established users", "New users")
  ) +
  labs(
    x = "Days relative to cutoff",
    y = "Comments per user-day",
    caption = "Data: Republik magazine, n = 884"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  filename = "rdit_heterogeneity_new_users.pdf",
  plot = plot_new + theme(text = element_text(family = "Times New Roman", size = 18)),
  width = 12,
  height = 8
)

# ---------- Plot: heterogeneity by activity ----------
# Same approach as above, but comparing non-top vs. top users based on pre-cutoff
# comment volume.
grid_top <- expand.grid(
  t_days = seq(-bw, bw, by = 0.1),
  is_top_user = c(0, 1)
)
grid_top$post <- as.integer(grid_top$t_days >= 0)
grid_top$pred <- predict(m_het_top, newdata = grid_top)

plot_top <- ggplot(df_userday, aes(x = t_days, y = comments, color = factor(is_top_user))) +
  stat_summary(fun = mean, geom = "point", size = 2, alpha = 0.7) +
  geom_line(data = grid_top, aes(y = pred), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    name = "User type",
    values = c("0" = "#7570b3", "1" = "#e7298a"),
    labels = c("Non-top users", "Top users")
  ) +
  labs(
    x = "Days relative to cutoff",
    y = "Comments per user-day",
    caption = "Data: Republik magazine, n = 884"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  filename = "rdit_heterogeneity_top_users.pdf",
  plot = plot_top + theme(text = element_text(family = "Times New Roman", size = 18)),
  width = 12,
  height = 8
)


models <- list(
  "Base model"      = m_rdit,
  "New users (H3)"  = m_het_new,
  "Top users (H4)"  = m_het_top
)

vcovs <- list(
  "Base model"      = V_day,
  "New users (H3)"  = V_new,
  "Top users (H4)"  = V_top
)

modelsummary(
  models,
  vcov = vcovs,
  stars = TRUE,
  statistic = "({std.error})",
  output = "latex",          # or "latex_tabular" for classic tabular
  coef_map = c(
    "(Intercept)"            = "Intercept",
    "t_days"                 = "t (days)",
    "post"                   = "Post",
    "t_days:post"            = "t × Post",
    "post:is_new_user"       = "Post × New user",
    "t_days:post:is_new_user"= "t × Post × New user",
    "post:is_top_user"       = "Post × Top user",
    "t_days:post:is_top_user"= "t × Post × Top user"
  ),
  gof_map = c("nobs", "r.squared", "adj.r.squared")
)



# -----------------------------
# Helper: run RDiT at any cutoff
# -----------------------------
fit_rdit_cutoff <- function(df_clean, cutoff_dt, bw = 7, tz = "Europe/Zurich") {
  
  # Rebuild user-day panel relative to the placebo cutoff
  df_ud <- df_clean |>
    group_by(user_id, day) |>
    summarise(comments = n(), .groups = "drop") |>
    mutate(
      day_mid = ymd_hms(paste(day, "12:00:00"), tz = tz),
      t_days  = as.numeric(difftime(day_mid, cutoff_dt, units = "days")),
      post    = as.integer(t_days >= 0)
    ) |>
    filter(t_days >= -bw, t_days <= bw)
  
  # Fit your model
  m <- lm(comments ~ t_days + post + t_days:post, data = df_ud)
  
  # Day-clustered SEs (same as you do)
  V <- vcovCL(m, cluster = ~ day, data = df_ud)
  ct <- coeftest(m, vcov = V)
  
  # Return the discontinuity estimate ("post") and a few diagnostics
  out <- data.frame(
    cutoff = cutoff_dt,
    estimate = ct["post", "Estimate"],
    se = ct["post", "Std. Error"],
    t = ct["post", "t value"],
    p = ct["post", "Pr(>|t|)"],
    nobs = nrow(df_ud),
    ndays = length(unique(df_ud$day))
  )
  
  return(out)
}

# -----------------------------
# 1) TRUE cutoff estimate
# -----------------------------
true_cutoff <- ymd_hm("2025-11-23 19:30", tz = "Europe/Zurich")
true_res <- fit_rdit_cutoff(df_clean, true_cutoff, bw = bw)

true_res
# true_res$estimate is your baseline jump at the real cutoff

# -----------------------------
# 2) Placebo cutoffs (choose a window where nothing happened)
#    Common choice: ONLY pre-period cutoffs to avoid contamination.
# -----------------------------
placebo_days <- seq.Date(
  from = as.Date(true_cutoff) - 60,
  to   = as.Date(true_cutoff) - 15,
  by   = "day"
)

# mimic the real change time (19:30)
placebo_cutoffs <- ymd_hm(paste0(placebo_days, " 19:30"), tz = "Europe/Zurich")

placebo_res <- bind_rows(lapply(placebo_cutoffs, function(cd) {
  fit_rdit_cutoff(df_clean, cd, bw = bw)
}))

# -----------------------------
# 3) Plot placebo distribution + mark true estimate
# -----------------------------
plot_placebo <- ggplot(placebo_res, aes(x = estimate)) +
  geom_histogram(bins = 25) +
  geom_vline(xintercept = true_res$estimate, linewidth = 1) +
  labs(
    title = "Placebo RDiT discontinuity estimates (post coefficient)",
    subtitle = "Vertical line = true cutoff estimate",
    x = "Estimated discontinuity (jump) at placebo cutoff",
    y = "Number of placebo cutoffs",
    caption = "Data: Republik magazine, n = 884"
  )
ggsave(
  filename = "rdit_placebo_distribution.pdf",
  plot = plot_placebo + theme(text = element_text(family = "Times New Roman", size = 18)),
  width = 12,
  height = 8
)
# -----------------------------
# 4) Simple placebo-based "pseudo p-value"
#    (how often a placebo estimate is as large as the true estimate)
# -----------------------------
pseudo_p <- mean(abs(placebo_res$estimate) >= abs(true_res$estimate), na.rm = TRUE)
pseudo_p

