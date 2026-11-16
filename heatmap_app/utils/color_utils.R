annotation_palettes <- list(
  c("#5B9BD5","#E07B54","#6DBF8A","#F4C842","#A87CC7","#E05C7A"),
  c("#A8A8A8","#6DBF8A","#F4C842","#E07B54","#5B9BD5","#A87CC7"),
  c("#7B9E87","#C97D6E","#7D9CC0","#D4A85A","#8E7BB5","#C9A87C"),
  c("#E8C47A","#7ABDE8","#E87AB5","#7AE8A0","#E8907A","#907AE8")
)

get_cat_colors <- function(lvls, palette_index = 1) {
  pal <- annotation_palettes[[((palette_index - 1) %% length(annotation_palettes)) + 1]]
  setNames(pal[seq_along(lvls)], lvls)
}

# ---- The three anchor points of the heatmap colour ramp -----------------------
# colorRamp2(breaks, colors) needs `breaks` STRICTLY INCREASING. It does not
# check: given c(0.5, 0, 2) it silently builds a scale that runs
# white -> dark blue -> brown -> yellow, i.e. the mid colour lands BELOW the low
# colour. The figure looks deliberate and nothing warns.
#
# That was reachable in one drag. The "Z-score scale limits" slider spans
# -10..10 with both handles free, so a user zooming into up-regulated genes
# (0.5 .. 2) got exactly those breaks, because the midpoint was hardcoded to 0.
#
# The rule, in one place because two callers need it and must agree -- the live
# draw in app.r and the `rerun_heatmap_NN.R` script written by repro_module.R.
# If they disagreed, the bundle would stop reproducing the figure it ships with.
#
# zlim  : length-2 numeric from the slider (or anything at all -- it is validated)
# fixed : the anchors the ANALYSIS chose (2%/98% quantiles for none/log2, the
#         data-driven arms for zscore_col). When present they win: they are
#         already ordered by construction and are not the user's to override.
#
# Returns list(breaks = numeric(3), note = NULL | one sentence for the user).
# `note` is non-NULL only when the request had to be reinterpreted, so a caller
# can show it without spamming the normal path.
heatmap_color_breaks <- function(zlim, fixed = NULL) {
  # 1. Analysis-chosen anchors win, provided they are usable.
  if (!is.null(fixed) && length(fixed) == 3 && all(is.finite(fixed)) &&
      fixed[1] < fixed[2] && fixed[2] < fixed[3]) {
    return(list(breaks = as.numeric(fixed), note = NULL))
  }

  note <- NULL
  zl   <- suppressWarnings(as.numeric(zlim))

  # 2. A cleared slider sends NA, and %||% does not catch that (global.R:42):
  #    it only replaces NULL / length-0 / "". Fall back to the widget default.
  if (length(zl) != 2 || !all(is.finite(zl))) {
    zl <- c(-1, 1)

  } else {
    # 3. Reversed handles: honour the range the user drew rather than refusing.
    if (zl[1] > zl[2]) zl <- rev(zl)

    # 4. Both handles on the same value: there is no range to draw. Three equal
    #    breaks make colorRamp2() stop outright ("at least two distinct break
    #    values"), which the caller could only report as a generic render error.
    if (zl[1] == zl[2]) {
      note <- sprintf(paste0(
        "Both z-score scale limits are set to %s, which is not a range: the ",
        "colour scale needs a low and a high value. Falling back to the default ",
        "(-1, 1) — drag the two handles apart to choose your own."),
        format(zl[1]))
      zl <- c(-1, 1)
    }
  }

  # 5. THE FIX. The mid colour sits at 0 only when 0 is actually inside the
  #    range. For a window that lies entirely above or below zero, 0 is not a
  #    midpoint at all and using it would put the mid colour outside the scale,
  #    producing the non-monotonic ramp described above. Use the centre of the
  #    chosen window instead, and say so -- the mid colour normally reads as
  #    "no change", and here it no longer does.
  mid <- if (zl[1] < 0 && zl[2] > 0) {
    0
  } else {
    m <- (zl[1] + zl[2]) / 2
    note <- c(note, sprintf(paste0(
      "The z-score scale limits (%s to %s) do not include 0, so the middle ",
      "colour is placed at %s — the centre of the window you chose — ",
      "instead of at z = 0. The middle colour therefore no longer means ",
      "'no change'."), format(zl[1]), format(zl[2]), format(round(m, 3))))
    m
  }

  list(breaks = c(zl[1], mid, zl[2]),
       note   = if (length(note)) paste(note, collapse = " ") else NULL)
}
