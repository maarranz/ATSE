# ==============================================================================
# adl_theory.R
#
# Teaching toolkit for exploring the theoretical dynamic properties of
# Autoregressive Distributed Lag (ADL) models. R port of the Python library
# adl_theory.py (https://github.com/maarranz/adl-theory).
#
# The model is
#
#   y_t = phi_1 y_{t-1} + ... + phi_p y_{t-p}
#         + beta_0 x_t + ... + beta_q x_{t-q} + u_t
#
# The functions do NOT estimate ADL models from data. They take the model
# coefficients as given and show what those coefficients imply: stability,
# dynamic and cumulative effects, long-run multiplier, mean and median lags.
#
# Usage (self-contained; only needs the plotly package):
#
#   source("adl_theory.R")
#
#   phi  <- c(0.6)              # phi_1, ..., phi_p
#   beta <- c(0.5, 0.3)         # beta_0, ..., beta_q
#
#   stab <- adl_stability(phi)
#   plot_inverse_roots(stab)
#
#   res <- adl_dynamic_effects(phi, beta)   # automatic 99% horizon
#   res                                      # printed summary + table
#   adl_summary(res)                         # summary as a data frame
#   plot_dynamic_effects(res)
#   plot_cumulative_effects(res)
#
# In a Quarto document, knitr::kable(adl_summary(res)) gives a formatted table.
#
# Author : Miguel A. Arranz
# License: MIT
# Version: 0.1.0
# ==============================================================================

if (!requireNamespace("plotly", quietly = TRUE)) {
  stop("adl_theory.R needs the 'plotly' package: install.packages(\"plotly\")",
       call. = FALSE)
}

# ------------------------------------------------------------------------------
# Internal helpers (not meant to be called by students)
# ------------------------------------------------------------------------------

# Equivalent of numpy's isclose(x, 0) with default tolerances
.adl_is_zero <- function(x) abs(x) <= 1e-8

# Format a number for annotations; NA / NaN / Inf become "Not defined"
.adl_fmt <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || !is.finite(x)) return("Not defined")
  formatC(x, format = "f", digits = digits)
}

.adl_check_results <- function(results) {
  if (!inherits(results, "adl_effects")) {
    stop("`results` must be the object returned by adl_dynamic_effects().",
         call. = FALSE)
  }
}

# Integer x-axis ticks: one per lag for short horizons, automatic otherwise
.adl_xaxis <- function(horizon) {
  ax <- list(title = "Lag", showgrid = FALSE, zeroline = FALSE)
  if (horizon <= 40) {
    ax$tickmode <- "linear"
    ax$tick0 <- 0
    ax$dtick <- 1
  }
  ax
}

# Horizontal line spanning the whole plot width
.adl_hline <- function(y, color = "black", width = 1, dash = "solid") {
  list(type = "line", xref = "paper", x0 = 0, x1 = 1, yref = "y", y0 = y, y1 = y,
       line = list(color = color, width = width, dash = dash), layer = "below")
}

# Vertical line spanning the whole plot height
.adl_vline <- function(x, color = "black", width = 1, dash = "solid") {
  list(type = "line", yref = "paper", y0 = 0, y1 = 1, xref = "x", x0 = x, x1 = x,
       line = list(color = color, width = width, dash = dash), layer = "below")
}

# Boxed text annotation in paper coordinates
.adl_box <- function(text, x, y, xanchor, yanchor, bordercolor = "gray") {
  list(text = text, x = x, y = y, xref = "paper", yref = "paper",
       xanchor = xanchor, yanchor = yanchor, align = "left", showarrow = FALSE,
       bordercolor = bordercolor, borderwidth = 1, borderpad = 8,
       bgcolor = "white")
}

.adl_summary_text <- function(s, lr_label = "Long-run effect") {
  paste0(
    "<b>Model summary</b><br>",
    "Stability: ", if (isTRUE(s$stable)) "Stable" else "Unstable", "<br>",
    lr_label, ": ", .adl_fmt(s$long_run_effect), "<br>",
    "Mean lag: ", .adl_fmt(s$mean_lag), "<br>",
    "Median lag: ", .adl_fmt(s$median_lag, digits = 0)
  )
}

# ==============================================================================
# Stability
# ==============================================================================

#' Stability of an ADL model
#'
#' The characteristic polynomial is 1 - phi_1 z - ... - phi_p z^p. The model
#' is stable if all its roots lie outside the unit circle, i.e. all inverse
#' roots lie inside it. Stability depends only on phi, not on beta.
#'
#' @param phi Autoregressive coefficients c(phi_1, ..., phi_p). Use
#'   numeric(0) for a finite distributed lag model (no autoregressive part).
#' @return An object of class "adl_stability": a list with phi, order (p),
#'   roots, inverse_roots, stable (TRUE/FALSE) and max_inverse_root.
#' @examples
#' adl_stability(c(0.8))$stable                 # TRUE
#' adl_stability(c(1.2, -0.5))                  # complex inverse roots
adl_stability <- function(phi) {
  phi <- as.numeric(phi)
  if (any(!is.finite(phi))) stop("All values in `phi` must be finite.", call. = FALSE)

  # Trailing zero coefficients do not add roots (same as numpy.roots)
  last <- max(c(0, which(phi != 0)))
  roots <- if (last == 0) complex(0) else polyroot(c(1, -phi[seq_len(last)]))
  inverse_roots <- 1 / roots

  structure(
    list(
      phi = phi,
      order = length(phi),
      roots = roots,
      inverse_roots = inverse_roots,
      stable = all(Mod(roots) > 1),
      max_inverse_root = if (length(roots) > 0) max(Mod(inverse_roots)) else 0
    ),
    class = "adl_stability"
  )
}

#' @export
print.adl_stability <- function(x, digits = 4, ...) {
  cat("ADL stability (AR order p =", x$order, ")\n")
  cat("Stable:", if (x$stable) "Yes" else "No",
      "  |  Max |inverse root| =", formatC(x$max_inverse_root, format = "f", digits = digits),
      "\n")
  if (length(x$inverse_roots) > 0) {
    tab <- data.frame(
      Real      = round(Re(x$inverse_roots), digits),
      Imaginary = round(Im(x$inverse_roots), digits),
      Modulus   = round(Mod(x$inverse_roots), digits),
      row.names = paste0("r", seq_along(x$inverse_roots))
    )
    cat("\nInverse roots:\n")
    print(tab)
  } else {
    cat("No autoregressive part: the model is a finite distributed lag.\n")
  }
  invisible(x)
}

# ==============================================================================
# Graphics: inverse roots
# ==============================================================================

#' Plot the inverse roots of the characteristic polynomial
#'
#' Interactive plot of the inverse roots against the unit circle. Points are
#' blue for a stable model and red for an unstable one; hovering shows the
#' real part, imaginary part and modulus of each root.
#'
#' @param stability Object returned by adl_stability().
#' @param title Figure title.
#' @param width,height Figure size in pixels.
#' @return A plotly htmlwidget.
plot_inverse_roots <- function(stability,
                               title = "Inverse Roots of the Characteristic Polynomial",
                               width = 700, height = 700) {
  if (!inherits(stability, "adl_stability")) {
    stop("`stability` must be the object returned by adl_stability().", call. = FALSE)
  }

  inv  <- stability$inverse_roots
  col  <- if (stability$stable) "royalblue" else "firebrick"
  theta <- seq(0, 2 * pi, length.out = 400)

  fig <- plotly::plot_ly(width = width, height = height)

  # Unit circle
  fig <- plotly::add_trace(
    fig, x = cos(theta), y = sin(theta),
    type = "scatter", mode = "lines", name = "Unit circle",
    line = list(color = "gray", dash = "dash"), hoverinfo = "skip"
  )

  # Inverse roots
  if (length(inv) > 0) {
    hover <- sprintf(
      "<b>Inverse root %d</b><br>Real: %.4f<br>Imaginary: %.4f<br>Modulus: %.4f",
      seq_along(inv), Re(inv), Im(inv), Mod(inv)
    )
    fig <- plotly::add_trace(
      fig, x = Re(inv), y = Im(inv),
      type = "scatter", mode = "markers+text", name = "Inverse roots",
      text = paste0("r", seq_along(inv)), textposition = "top center",
      marker = list(size = 10, color = col),
      customdata = hover, hovertemplate = "%{customdata}<extra></extra>"
    )
  }

  annotation <- paste0(
    "<b>", if (stability$stable) "Stable" else "Unstable", " model</b><br>",
    "ADL(", stability$order, ")<br>",
    "Max |inverse root| = ", formatC(stability$max_inverse_root, format = "f", digits = 3)
  )

  plotly::layout(
    fig,
    title = list(text = title),
    xaxis = list(title = "Real", zeroline = FALSE,
                 scaleanchor = "y", scaleratio = 1),
    yaxis = list(title = "Imaginary", zeroline = FALSE),
    shapes = list(.adl_hline(0), .adl_vline(0)),
    annotations = list(.adl_box(annotation, x = 0.02, y = 0.98,
                                xanchor = "left", yanchor = "top",
                                bordercolor = "black")),
    legend = list(x = 0.02, y = 0.02),
    margin = list(l = 70, r = 40, t = 80, b = 70),
    plot_bgcolor = "white", paper_bgcolor = "white"
  )
}

# ==============================================================================
# Dynamic effects
# ==============================================================================

#' Dynamic effects of an ADL(p, q) model
#'
#' Dynamic effects (dynamic multipliers) are computed recursively as
#'
#'   psi_0 = beta_0
#'   psi_h = beta_h + phi_1 psi_{h-1} + ... + phi_p psi_{h-p},
#'
#' with beta_h = 0 for h > q. The long-run effect is
#' sum(beta) / (1 - sum(phi)), and the mean lag is
#' sum(j beta_j) / sum(beta_j) + sum(i phi_i) / (1 - sum(phi_i)).
#' The median lag is the first lag at which the cumulative effect reaches
#' 50% of the long-run effect.
#'
#' For an unstable model the effects are still computed for a finite
#' horizon, but the long-run effect, percentages, mean and median lags are
#' not defined (NA), and a warning is issued.
#'
#' @param phi Autoregressive coefficients c(phi_1, ..., phi_p);
#'   numeric(0) for a finite distributed lag model.
#' @param beta Distributed-lag coefficients c(beta_0, ..., beta_q).
#' @param horizon Last lag to compute. If NULL (default), the first lag at
#'   which the cumulative effect reaches `tolerance` times the long-run
#'   effect is selected, searching up to `max_horizon`.
#' @param tolerance Share of the long-run effect that defines the automatic
#'   horizon (default 0.99).
#' @param max_horizon Upper limit for the automatic horizon (default 500).
#' @return An object of class "adl_effects": a list with phi, beta,
#'   stability (an adl_stability object), table (data frame with columns
#'   lag, effect, cumulative, pct_total, pct_cumulative) and summary (list
#'   with stable, horizon, long_run_effect, mean_lag, median_lag, tolerance,
#'   automatic_horizon).
#' @examples
#' res <- adl_dynamic_effects(phi = 0.6, beta = c(0.5, 0.3))
#' res$summary$long_run_effect      # 2
#' res$table
adl_dynamic_effects <- function(phi, beta, horizon = NULL,
                                tolerance = 0.99, max_horizon = 500) {
  phi  <- as.numeric(phi)
  beta <- as.numeric(beta)

  # ---- Input checks ---------------------------------------------------------
  if (length(beta) == 0) stop("`beta` must contain at least beta_0.", call. = FALSE)
  if (any(!is.finite(phi)))  stop("All values in `phi` must be finite.",  call. = FALSE)
  if (any(!is.finite(beta))) stop("All values in `beta` must be finite.", call. = FALSE)

  is_whole <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x)

  if (!is.null(horizon)) {
    if (!is_whole(horizon)) stop("`horizon` must be a whole number or NULL.", call. = FALSE)
    if (horizon < 0) stop("`horizon` must be non-negative.", call. = FALSE)
  }
  if (!(is.numeric(tolerance) && length(tolerance) == 1 &&
        tolerance > 0 && tolerance < 1)) {
    stop("`tolerance` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (!is_whole(max_horizon)) stop("`max_horizon` must be a whole number.", call. = FALSE)
  if (max_horizon < 0) stop("`max_horizon` must be non-negative.", call. = FALSE)

  # ---- Stability and long-run effect ---------------------------------------
  stability <- adl_stability(phi)
  stable <- stability$stable

  if (!stable) {
    warning("The ADL model is dynamically unstable. Dynamic effects can still ",
            "be calculated for a finite horizon, but the theoretical long-run ",
            "effect and lag-distribution summaries are not defined.",
            call. = FALSE)
  }

  denominator <- 1 - sum(phi)
  long_run_effect <- if (stable && !.adl_is_zero(denominator)) {
    sum(beta) / denominator
  } else {
    NA_real_
  }
  lr_usable <- stable && is.finite(long_run_effect) && !.adl_is_zero(long_run_effect)

  p <- length(phi)
  q <- length(beta) - 1

  # psi_h given psi_0, ..., psi_{h-1}  (psi is stored with index h + 1)
  effect_at <- function(psi, h) {
    beta_h <- if (h <= q) beta[h + 1] else 0
    ar <- 0
    for (i in seq_len(min(p, h))) ar <- ar + phi[i] * psi[h - i + 1]
    beta_h + ar
  }

  # ---- Recursion ------------------------------------------------------------
  automatic <- is.null(horizon)
  if (automatic) {
    psi <- numeric(0)
    for (h in 0:max_horizon) {
      psi[h + 1] <- effect_at(psi, h)
      if (lr_usable && sum(psi) / long_run_effect >= tolerance) break
    }
    selected_horizon <- length(psi) - 1
    if (lr_usable && selected_horizon == max_horizon &&
        sum(psi) / long_run_effect < tolerance) {
      warning(sprintf(
        "The cumulative effect did not reach %.1f%% of the theoretical long-run effect by lag %d.",
        100 * tolerance, as.integer(max_horizon)), call. = FALSE)
    }
  } else {
    selected_horizon <- horizon
    psi <- numeric(0)
    for (h in 0:selected_horizon) psi[h + 1] <- effect_at(psi, h)
  }

  cumulative <- cumsum(psi)

  if (lr_usable) {
    pct_total      <- 100 * psi / long_run_effect
    pct_cumulative <- 100 * cumulative / long_run_effect
  } else {
    pct_total <- pct_cumulative <- rep(NA_real_, length(psi))
  }

  # ---- Mean and median lag --------------------------------------------------
  mean_lag <- if (stable && !.adl_is_zero(sum(beta)) && !.adl_is_zero(denominator)) {
    sum((seq_along(beta) - 1) * beta) / sum(beta) + sum(seq_len(p) * phi) / denominator
  } else {
    NA_real_
  }

  median_lag <- NA_integer_
  if (lr_usable) {
    hit <- which(cumulative / long_run_effect >= 0.5)
    if (length(hit) > 0) median_lag <- as.integer(hit[1] - 1)
  }

  table <- data.frame(
    lag            = 0:selected_horizon,
    effect         = psi,
    cumulative     = cumulative,
    pct_total      = pct_total,
    pct_cumulative = pct_cumulative
  )

  summary <- list(
    stable            = stable,
    horizon           = as.integer(selected_horizon),
    long_run_effect   = long_run_effect,
    mean_lag          = mean_lag,
    median_lag        = median_lag,
    tolerance         = if (automatic) tolerance else NA_real_,
    automatic_horizon = automatic
  )

  structure(
    list(phi = phi, beta = beta, stability = stability,
         table = table, summary = summary),
    class = "adl_effects"
  )
}

#' @export
print.adl_effects <- function(x, digits = 4, n = 30, ...) {
  p <- length(x$phi)
  q <- length(x$beta) - 1
  cat(sprintf("ADL(%d, %d) dynamic effects\n", p, q))
  cat("  phi  =", paste(format(x$phi),  collapse = ", "), "\n")
  cat("  beta =", paste(format(x$beta), collapse = ", "), "\n\n")

  s <- adl_summary(x, digits = digits)
  for (i in seq_len(nrow(s))) cat(sprintf("  %-17s %s\n", s$Statistic[i], s$Value[i]))
  cat("\n")

  tab <- x$table
  names(tab) <- c("Lag", "Dynamic Effect", "Cumulative Effect",
                  "% Total Effect", "% Cumulative Effect")
  tab[-1] <- lapply(tab[-1], round, digits = digits)
  if (nrow(tab) > n) {
    print(utils::head(tab, n), row.names = FALSE)
    cat(sprintf("  ... %d more lags (see $table, or print(res, n = Inf))\n",
                nrow(tab) - n))
  } else {
    print(tab, row.names = FALSE)
  }
  invisible(x)
}

# ==============================================================================
# Summary table
# ==============================================================================

#' Summary table of an ADL model's dynamic properties
#'
#' @param results Object returned by adl_dynamic_effects().
#' @param digits Decimal places for numeric values.
#' @return A data frame with columns Statistic and Value (character), ready
#'   for knitr::kable() in a Quarto document.
adl_summary <- function(results, digits = 4) {
  .adl_check_results(results)
  s <- results$summary

  horizon_method <- if (isTRUE(s$automatic_horizon)) {
    sprintf("Automatic (%.0f%% threshold)", 100 * s$tolerance)
  } else {
    "User specified"
  }

  data.frame(
    Statistic = c("Stability", "Long-run effect", "Mean lag", "Median lag",
                  "Selected horizon", "Horizon method"),
    Value = c(
      if (isTRUE(s$stable)) "Stable" else "Unstable",
      .adl_fmt(s$long_run_effect, digits),
      .adl_fmt(s$mean_lag, digits),
      .adl_fmt(s$median_lag, 0),
      as.character(s$horizon),
      horizon_method
    ),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# Graphics: dynamic and cumulative effects
# ==============================================================================

#' Plot the dynamic effects (dynamic multipliers) of an ADL model
#'
#' Bar chart of the dynamic effect at each lag. Hovering shows the effect and
#' its share of the long-run effect; a box summarises stability, long-run
#' effect, mean lag and median lag.
#'
#' @param results Object returned by adl_dynamic_effects().
#' @param width,height Figure size in pixels.
#' @return A plotly htmlwidget.
plot_dynamic_effects <- function(results, width = 900, height = 520) {
  .adl_check_results(results)
  tab <- results$table
  s   <- results$summary

  share <- ifelse(is.finite(tab$pct_total),
                  sprintf("%.2f%%", tab$pct_total), "Not defined")
  hover <- sprintf("<b>Lag %d</b><br>Dynamic effect: %.4f<br>Share of LR effect: %s",
                   tab$lag, tab$effect, share)

  fig <- plotly::plot_ly(width = width, height = height)
  fig <- plotly::add_trace(
    fig, x = tab$lag, y = tab$effect, type = "bar", name = "Dynamic effect",
    marker = list(color = "#636EFA"),
    customdata = hover, hovertemplate = "%{customdata}<extra></extra>"
  )

  plotly::layout(
    fig,
    title = list(text = "Dynamic Effects", x = 0.5, xanchor = "center"),
    xaxis = .adl_xaxis(s$horizon),
    yaxis = list(title = "Dynamic Effect", zeroline = FALSE, gridcolor = "lightgray"),
    shapes = list(.adl_hline(0)),
    annotations = list(.adl_box(.adl_summary_text(s), x = 0.99, y = 0.98,
                                xanchor = "right", yanchor = "top")),
    showlegend = FALSE,
    hovermode = "x",
    margin = list(l = 70, r = 40, t = 80, b = 70),
    plot_bgcolor = "white", paper_bgcolor = "white"
  )
}

#' Plot the cumulative dynamic effects of an ADL model
#'
#' Line chart of the cumulative effect at each lag, with the long-run effect
#' (dashed line), the median lag (highlighted point) and, when the horizon
#' was chosen automatically, the tolerance horizon (dotted line).
#'
#' @param results Object returned by adl_dynamic_effects().
#' @param width,height Figure size in pixels.
#' @return A plotly htmlwidget.
plot_cumulative_effects <- function(results, width = 900, height = 520) {
  .adl_check_results(results)
  tab <- results$table
  s   <- results$summary

  share <- ifelse(is.finite(tab$pct_cumulative),
                  sprintf("%.2f%%", tab$pct_cumulative), "Not defined")
  hover <- sprintf("<b>Lag %d</b><br>Cumulative effect: %.4f<br>Share of LR effect: %s",
                   tab$lag, tab$cumulative, share)

  shapes <- list()
  annotations <- list()

  fig <- plotly::plot_ly(width = width, height = height)
  fig <- plotly::add_trace(
    fig, x = tab$lag, y = tab$cumulative,
    type = "scatter", mode = "lines+markers", name = "Cumulative effect",
    line = list(color = "#636EFA"), marker = list(color = "#636EFA"),
    customdata = hover, hovertemplate = "%{customdata}<extra></extra>"
  )

  # Long-run effect
  if (is.finite(s$long_run_effect)) {
    shapes[[length(shapes) + 1]] <- .adl_hline(s$long_run_effect, width = 1.5, dash = "dash")
    annotations[[length(annotations) + 1]] <- list(
      text = "Long-run effect", x = 0, xref = "paper", xanchor = "left",
      y = s$long_run_effect, yref = "y", yanchor = "bottom", showarrow = FALSE
    )
  }

  # Median lag
  if (!is.na(s$median_lag)) {
    m_eff <- tab$cumulative[tab$lag == s$median_lag]
    if (length(m_eff) == 1) {
      fig <- plotly::add_trace(
        fig, x = s$median_lag, y = m_eff,
        type = "scatter", mode = "markers", name = "Median lag",
        marker = list(size = 10, color = "#EF553B"),
        hovertemplate = sprintf("<b>Median lag: %d</b><br>Cumulative effect: %.4f<extra></extra>",
                                s$median_lag, m_eff)
      )
    }
  }

  # Automatic horizon
  if (isTRUE(s$automatic_horizon) && is.finite(s$tolerance)) {
    shapes[[length(shapes) + 1]] <- .adl_vline(s$horizon, color = "gray", dash = "dot")
    annotations[[length(annotations) + 1]] <- list(
      text = sprintf("%.0f%% horizon", 100 * s$tolerance),
      x = s$horizon, xref = "x", xanchor = "left",
      y = 1, yref = "paper", yanchor = "top", showarrow = FALSE
    )
  }

  annotations[[length(annotations) + 1]] <- .adl_box(
    .adl_summary_text(s, lr_label = "LR effect"),
    x = 0.99, y = 0.05, xanchor = "right", yanchor = "bottom"
  )

  plotly::layout(
    fig,
    title = list(text = "Cumulative Dynamic Effects", x = 0.5, xanchor = "center"),
    xaxis = .adl_xaxis(s$horizon),
    yaxis = list(title = "Cumulative Effect", zeroline = TRUE,
                 zerolinecolor = "black", zerolinewidth = 1, gridcolor = "lightgray"),
    shapes = shapes,
    annotations = annotations,
    hovermode = "x",
    margin = list(l = 70, r = 40, t = 80, b = 70),
    plot_bgcolor = "white", paper_bgcolor = "white"
  )
}
