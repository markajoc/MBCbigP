plotmbc <-
function (x, parameters, z, groups = NULL, p = NULL)
{
  dev.flush()
  p <- if (is.null(p))
    ncol(x)
  else p
  groups <- if (is.null(groups))
    ncol(parameters$mean)
  else groups
  close.screen(all.screens = TRUE)
  scr <- split.screen(c(p + 2, p + 2))
  scr2 <- (matrix(scr, ncol = p + 2)[c(-1, -(p + 2)), c(-1, -(p + 2))])
  scr2 <- scr2[lower.tri(scr2, diag = TRUE)]
  usr.matrix <- mar.matrix <- fig.matrix <- matrix(ncol = 4L, nrow = length(
    scr2))
  rows <- matrix(rep(1:p, each = p), ncol = p)
  cols <- matrix(rep(1:p, p), ncol = p)
  rows <- as.vector(rows[lower.tri(rows, diag = TRUE)])
  cols <- as.vector(cols[lower.tri(cols, diag = TRUE)])
  groupcolour <- rep(c("dodgerblue2", "red3", "green3", "slateblue",
    "darkorange", "skyblue1", "violetred4", "forestgreen", "steelblue4",
    "slategrey", "brown", "black", "darkseagreen", "darkgoldenrod3",
    "olivedrab", "royalblue", "tomato4", "cyan2", "springgreen2"), length.out
    = groups)
  clustering <- mclust::map(z)
  dev.hold()
  for (i in seq_along(scr2)){
    screen(scr2[i])
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    par(mgp = c(3, 0.25, 0.15))
    plot(x[, cols[i]], x[, rows[i]], xlab = "", ylab = "", xaxt = "n", yaxt =
      "n", col = if (identical(rows[i], cols[i])) NULL else groupcolour[clustering], cex = 0.5)
    if (!identical(rows[i], cols[i])){
      cn <- colnames(x)[c(cols[i], rows[i])]
      if (all(cn %in% rownames(parameters$mean))){
        for (k in 1:groups){
          mclust::mvn2plot(mu = parameters$mean[cn, k], sigma =
            parameters$sigma[cn, cn, k], col = "gray", lwd = 2)
        }
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
function (object, parameters, z, ...)
{
  clustering <- mclust::map(z)
  dev.set(object$device)
  dev.hold()
  for (i in seq_along(object$scr2)){
    cn <- colnames(object$x)[c(object$cols[i], object$rows[i])]
    if (all(cn %in% rownames(parameters$mean))){
      if (!identical(cn[1L], cn[2L])){
        screen(object$scr2[i])
        par(mar = object$mar.matrix[i, ])
        par(usr = object$usr.matrix[i, ])
        par(fig = object$fig.matrix[i, ])
        rect(object$usr.matrix[i, 1], object$usr.matrix[i, 3],
          object$usr.matrix[i, 2], object$usr.matrix[i, 4], col = "white")
        points(object$x[, cn, drop = FALSE], col = object$groupcolour[
          clustering], cex = 0.5)
        for (k in 1:object$groups){
          mclust::mvn2plot(mu = parameters$mean[cn, k], sigma =
            parameters$sigma[cn, cn, k], col = "gray", lwd = 2)
        }
      }
    }
  }
  dev.flush()
  invisible(object)
}

plotmbcbigp <-
function (x_A, x_B, mean_A, mean_B, sigma_AA, sigma_AB, sigma_BB, z,
  groups = NULL, p = NULL)
{
  x <- cbind(x_A, x_B)
  parameters <- list()
  parameters$mean <- reform_mean(mean_A = mean_A, mean_B = mean_B, groups =
    groups)
  parameters$sigma <- reform_sigma(sigma_AA = sigma_AA, sigma_AB = sigma_AB,
    sigma_BB = sigma_BB, groups = groups)
  plotmbc(x = x, parameters = parameters, z = z, groups = groups, p = p)
}

update.mbcbigpplot <-
function (object, mean_A, mean_B, sigma_AA, sigma_AB, sigma_BB, z, groups, ...)
{
  parameters <- list()
  parameters$mean <- reform_mean(mean_A = mean_A, mean_B = mean_B, groups =
    groups)
  parameters$sigma <- reform_sigma(sigma_AA = sigma_AA, sigma_AB = sigma_AB,
    sigma_BB = sigma_BB, groups = groups)
  update.mbcplot(object, parameters, z)
}
