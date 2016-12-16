plotmbc <-
function (x, parameters, z, groups = NULL, p = NULL)
{
  p <- if (is.null(p))
    ncol(x)
  else p
  groups <- if (is.null(groups))
    ncol(parameters$mean)
  else groups
  close.screen(all.screens = TRUE)
  scr <- split.screen(c(p + 2, p + 2))
  scr2 <- as.vector(matrix(scr, ncol = p + 2)[c(-1, -(p + 2)), c(-1, -(p + 2))])
  usr.matrix <- mar.matrix <- fig.matrix <- matrix(ncol = 4L, nrow = length(
    scr2))
  rows <- rep(1:p, each = p)
  cols <- rep(1:p, p)
  groupcolour <- rep(rainbow(min(max(3L, groups), 10L)), length.out = groups)

  dev.hold()
  for (i in seq_along(scr2)){
    screen(scr2[i])
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    par(mgp = c(3, 0.25, 0.15))
    plot(x[, cols[i]], x[, rows[i]], xlab = "", ylab = "", xaxt = "n", yaxt =
      "n", col = if (identical(rows[i], cols[i])) NULL else "gray", cex = 0.5)
    cn <- colnames(x)[c(cols[i], rows[i])]
    if (all(cn %in% rownames(parameters$mean))){
      for (k in 1:groups){
        mclust::mvn2plot(mu = parameters$mean[cn, k], sigma = parameters$sigma[
          cn, cn, k], col = groupcolour[k], lwd = 2)
      }
    }
    if (identical(rows[i], 1L) & (2 * (round(cols[i] / 2)) == cols[i]))
      axis(3, cex.axis = 0.7, tcl = -0.2)
    if (identical(rows[i], p) & !(2 * (round(cols[i] / 2)) == cols[i]))
      axis(1, cex.axis = 0.7, tcl = -0.2)
    if (identical(cols[i], 1L) & (2 * (round(rows[i] / 2)) == rows[i]))
      axis(2, cex.axis = 0.7, tcl = -0.2)
    if (identical(cols[i], p) & !(2 * (round(rows[i] / 2)) == rows[i]))
      axis(4, cex.axis = 0.7, tcl = -0.2)
    if (identical(rows[i], cols[i]))
      text(x = mean(range(x[, rows[i]])), y = mean(range(x[, cols[i]])),
        labels = colnames(x)[rows[i]])
    mar.matrix[i, ] <- par("mar")
    usr.matrix[i, ] <- par("usr")
    fig.matrix[i, ] <- par("fig")
  }
  coords <- data.frame(fig.matrix)
  names(coords) <- c("xleft", "xright", "ybottom", "ytop")
  coords$xcplots.index <- scr2
  dev.flush()
  invisible(structure(list(x = x, p = p, groups = groups, groupcolour =
    groupcolour, rows = rows, cols = cols, scr2 = scr2, coords = coords, device
    = dev.cur(), mar.matrix = mar.matrix, usr.matrix = usr.matrix), class =
    "mbcplot"))
}

update.mbcplot <-
function (object, parameters)
{
  dev.set(object$device)
  dev.hold()
  for (i in seq_along(object$scr2)){
    cn <- colnames(object$x)[c(object$cols[i], object$rows[i])]
    if (all(cn %in% rownames(parameters$mean))){
      screen(object$scr2[i])
      par(mar = object$mar.matrix[i, ])
      par(usr = object$usr.matrix[i, ])
      par(fig = object$fig.matrix[i, ])
      rect(object$usr.matrix[i, 1], object$usr.matrix[i, 3], object$usr.matrix[
        i, 2], object$usr.matrix[i, 4], col = "white")
      points(object$x[, cn, drop = FALSE], col = "gray", cex = 0.5)
      for (k in 1:object$groups){
        mclust::mvn2plot(mu = parameters$mean[cn, k], sigma = parameters$sigma[
          cn, cn, k], col = object$groupcolour[k], lwd = 2)
      }
    }
  }
  dev.flush()
  invisible(object)
}
