Break = "\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   Welcome to 'EDP 480C.6 Statistical Analysis for Experimental Data'.
   Programs developed by Reza Norouzian, Copyright (C) 2020-present"

message(Break, notice, Break)

#================================================================================================================================

d2f <- function(d, n1, n2) sqrt(d2peta(d, n1, n2) / (1 - d2peta(d, n1, n2) ))
f2d <- function(f, n1, n2) peta2d(f2peta(f), n1, n2)                  
peta2f <- function(peta) sqrt(peta / (1 - peta))
f2peta <- function(f) (f^2) / (1 + f^2)
peta2F <- function(peta, df1, df2) (peta / df1) / ((1 - peta)/df2)
F2peta <- function(F.value, df1, df2) (F.value*df1) / ((F.value*df1) + df2)
d2r <- function(d, n1 = 300, n2 = 300) sqrt((d^2) / ((d^2) + (((n1 + n2)^2) - (2*(n1 + n2))) / (n1 * n2)))
r2d <- function(r, n1 = 300, n2 = 300) sqrt((r^2)*(((n1 + n2)^2)-(2*(n1 + n2)))/(n1 * n2)/(1-(r^2)))
d2t <- function(d, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  d*sqrt(N)
}

t2d <- function(t, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  t/sqrt(N)
}

ncp2peta <- function(ncp, N) { ncp / (ncp + N) }


peta2ncp <- function(peta, N) { (peta*N) / (1 - peta) }

peta2N <- function(peta, ncp) { (ncp - (peta * ncp)) / peta }

d2peta <- function(d, n1 = 300, n2 = 300) (d^2) / ((d^2) + (((n1 + n2)^2) - (2*(n1 + n2))) / (n1 * n2))

peta2d <- function(peta, n1 = 300, n2 = 300) sqrt((peta)*(((n1 + n2)^2)-(2*(n1 + n2)))/(n1 * n2)/(1-(peta)))

F2pomega <- function(F.value, df1, N){
  (df1 * (F.value - 1)) / ((df1 * (F.value - 1)) + N)
}

pomega2F <- function(pomega, df1, N) {
  1 - ( (N * pomega )/(df1 * (pomega - 1)) )
}


peta2pomega <- function(peta, df1, df2, N){
  f <- peta2F(peta, df1, df2)  
  F2pomega(f, df1, N)
}


pomega2peta <- function(pomega, df1, df2, N){
  f <- pomega2F(pomega, df1, N)  
  F2peta(f, df1, df2)
}


exp.pov <- exp.peta <- Vectorize(function(pbase = 0, df1, df2, N){
  
  integrate(function(x, df1, df2, pbase, N){
    
    x * dpeta(x = x, df1 = df1, df2 = df2, pbase = pbase, N = N)
    
  }, 0, 1, df1 = df1, df2 = df2, pbase = pbase, N = N)[[1]]
  
})


exp.d <- Vectorize(function(dbase = 0, n1, n2 = NA){
  
  integrate(function(x, n1, n2, dbase){
    
    x * dcohen(x = x, dbase = dbase, n1 = n1, n2 = n2)
    
  }, -Inf, Inf, n1 = n1, n2 = n2, dbase = dbase)[[1]]
  
})


exp2peta <- Vectorize(function(exp.val, df1, df2, N){
  
  optimize(function(x){
    
    abs(exp.val - exp.peta(pbase = x, df1 = df1, df2 = df2, N = N))
    
  }, 0:1, tol = 1e-9)[[1]]
  
})


exp2d <- Vectorize(function(exp.val, n1, n2 = NA){
  
  uniroot(function(x){
    
    exp.val - exp.d(dbase = x, n1 = n1, n2 = n2)
    
  }, c(-4, 4), extendInt = "yes")[[1]]
})

#==================================================================================================================================


f.balance <- function(F.unbalance, cell.makeup, df1, df2, N, conf.level = .9)
{
  
  fbalance <- F.unbalance * (mean(cell.makeup) / harmonic(cell.makeup))
  
  ci <- peta.ci(f = c(fbalance, F.unbalance), df1 = df1, df2 = df2, N = N, conf.level = conf.level)
  l <- length(F.unbalance)
  rownames(ci) <- paste((2*l):1, c(rep("balanced", l), rep("Unbalace", l)))
  ci[nrow(ci):1,]
}


#====================================================================================================================================


is.whole <- function(x)  abs(x - round(x)) < .Machine$double.eps^.5


#====================================================================================================================================




gpower.peta <- function(spss, df2, N, design){
  
  (spss * df2) / (N - (spss * design))
  
}

gpower.peta.bw <- function(peta, rho = .5, N, m, n.group){
  
  ((1 - rho)*peta*(N - n.group)*(m -1)) / ((1 - peta)*(m*N) + (1-rho)*peta*(N - n.group)*(m -1))
  
}



gpower.peta.b <- function(peta, rho = .5, N, m, n.group){
  
  (peta*(N - n.group)*(1 + (m-1)*rho)) / ((1 - peta)*(m*N) + peta*(N - n.group)*(1 + (m-1)*rho))
  
}




#===========================================================================================================================


dens.plot <- function(x, adjust = 1, na.rm = TRUE, n = 1e3, from = min(x), to = max(x), add = FALSE, hdi = FALSE, ci = FALSE, level = .95, xlab = deparse(substitute(x)), main = NA, lwd = 2, lty = 1, ...){
  
  UseMethod("dens.plot")
}


dens.plot.default <- function(x, adjust = 1, na.rm = TRUE, n = 1e3, from = min(x), to = max(x), add = FALSE, hdi = FALSE, ci = FALSE, level = .95, xlab = deparse(substitute(x)), main = NA, lwd = 2, lty = 1, ...){
  
  d <- density(x, adjust = adjust, na.rm = na.rm, n = n, from = from, to = to)
  
  if(!add){
    
    graphics.off()                            
    
    plot(d, zero.line = FALSE, xlab = xlab, main = main, lwd = lwd, lty = lty, ...)
    
  } else {
    
    lines(d, lwd = lwd, lty = lty, ...)
    
  }
  
  
  alpha <- (1 - level)/2
  q <- if(ci) quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE) else c(NA, NA)
  i <- if(hdi) hdir(x, level = level) else c(NA, NA)
  
  
  mode <- d$x[which.max(d$y)]
  mean <- mean(x)
  median <- median(x)
  sd <- sd(x)
  mad <- mad(x)
  
  if(hdi){
    h <- min(d$y)
    lines(i, rep(h, 2), lend = 1, lwd = 6, lty = 1, xpd = NA, ...)
    text(i, h, round(i, 3), pos = 3, cex = .8, font = 2, xpd = NA)
    points(mode, h, pch = 21, bg = "cyan", col = "magenta", cex = 1.7, xpd = NA)
  }
  
  invisible(list(lower = i[1], upper = i[2], level = level, mean = mean, mode = mode, median = median, 
                 mad = mad, sd = sd, q1 = q[[1]], q2 = q[[2]], x = d$x, y = d$y))
}


#===================================================================================================================




set.margin <- function() 
{
  par.mf <- par("mfrow", "mfcol")
  if (all(unlist(par.mf) == 1)) {
    par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
        tck = -0.02)
  }
}



#==================================================================================================================

HDI <- function(fun, lower = 0, upper = 1, level = .95, eps = 1e-3)
{
  UseMethod("HDI")
}

HDI.default <- function(fun, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  x <- formals(fun)
  FUN <- function(x) fun(x)
  
  posterior <- function(x) FUN(x)/integrate(FUN, lower, upper)[[1]]
  mode <- optimize(posterior, c(lower, upper), maximum = TRUE)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                 right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("You may change prior parameters or 'lower' & 'upper'.", call. = FALSE)
    return(ur[[1]])
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- 1
  find.lims <- function(a) {
    ur <- uniroot(function(h) areafun(h) / post.area - a,
                  c(eps, posterior(mode) - eps))
    return(ur[[1]])
  }
  f <- find.lims(level)
  return(c(inverse.posterior(f, "left"),
           inverse.posterior(f, "right")))
}


#==================================================================================================================

hdi <- function(x, y, level = .95)
{
  UseMethod("hdi")
}

hdi.default <- function(x, y, level = .95){
  
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  areas <- diff(x) * .5 * (head(y, -1) + tail(y, -1))
  peak <- which.max(areas)
  range <- c(peak, peak)
  found <- areas[peak]
  while(found < level) {
    if(areas[range[1]-1] > areas[range[2]+1]) {
      range[1] <- range[1]-1
      found <- found + areas[range[1]-1]
    } else {
      range[2] <- range[2]+1
      found <- found + areas[range[2]+1]
    }
  }
  val <- x[range]
  attr(val, "indexes") <- range
  attr(val, "area") <- found
  return(val)
}

#==================================================================================================================

hdir <- function(sample, level = .95)
{
  UseMethod("hdir")
}

hdir.default <- function(sample, level = .95){
  
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  #if(length(sample) < 1e3) message("Warning: \n\tInsufficient sample to produce reliable 'interval' estimates.")  
  sorted <- sort(sample)
  index <- ceiling(level*length(sorted))
  n <- length(sorted)- index
  width <- numeric(n)
  for(i in 1:n){
    width[i] <- sorted[i+ index]- sorted[i]
  }
  lower <- sorted[which.min(width)]
  upper <- sorted[which.min(width)+ index]
  return(c(lower, upper))
}

#==================================================================================================================

cip <- function(fun, lower = 0, upper = 1, level = .95){
  
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  
  p <- (1 - level) / 2
  
  inv.cdf(c(p, 1-p), fun, lower, upper)
}

#====================================================================================================================


hdiq <- function(qdist, level = .95, ...)
{
  UseMethod("hdiq")
}


hdiq.default <- function(qdist, level = .95, ...)
{
  
  alpha <-  1L - level
  width <- function(lower, qdist, level, ...){
    qdist(level + lower, ...) - qdist(lower, ...)
  }
  
  low <- optimize(width, c(0, alpha), qdist = qdist, level = level, ...)[[1]]
  
  return(c(qdist(low, ...), qdist(level + low, ...)))
}

#==================================================================================================================


cdf <- Vectorize(function(q, fun, lower, upper){
  
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  
  x <- formals(fun)
  f <- function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  
  integrate(f, lower, q)[[1]]
  
}, c("q", "lower", "upper"))


#==================================================================================================================


inv.cdf <- Vectorize(function(p, fun, lower, upper){
  
  if(!is.function(fun)) stop("Error: 'fun' must be a function.")
  if(length(formals(fun)) > 1) stop("Error: 'fun' must be a 'single-argument' function.")
  
  uniroot(function(q) cdf(q, fun, lower, upper) - p, c(lower, upper), extendInt = "yes")[[1]]
  
}, c("p", "lower", "upper"))

#==========================================================================================================================


typem.anova <- function(peta.h1, df1, df2, N, alpha = .05, peta.h0 = 0, peta.obs = .1, xlab = bquote(eta[p]^2), from = 0, to = .2){
  
  graphics.off()  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), xpd = TRUE)
  
  h0 = curve(dpeta(x, df1, df2, peta.h0, N), from = from, to = to, lwd = 2, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i")
  a <- qpeta(alpha, df1, df2, peta.h0, N, lower.tail = FALSE)
  x = seq(a, 1, l = 1e3) ; y = dpeta(x, df1, df2, peta.h0, N)
  polygon(c(a, x, 1), c(0, y, 0), col = 2, border = NA)
  lines(h0, lwd = 2)
  abline(v = peta.h0, col = 2, xpd = FALSE) 
  
  h1 = curve(dpeta(x, df1, df2, peta.h1, N), from = from, to = to, lwd = 2, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i")
  x = seq(a, 1, l = 1e3) ; y = dpeta(x, df1, df2, peta.h1, N)
  polygon(c(a, x, 1), c(0, y, 0), col = 2, border = NA)
  lines(h1, lwd = 2)
  abline(v = peta.h1, col = 4, xpd = FALSE)
  
  points(peta.obs, 0, pch = 23, bg = "cyan", col = "magenta", cex = 1.5)
  
  abline(v = a, col = 2, lty = 2, xpd = NA)
  
  power <- ppeta(a, df1, df2, peta.h1, N, lower.tail = FALSE)
  p.value <- ppeta(peta.obs, df1, df2, peta.h0, N, lower.tail = FALSE)
  
  random.p <- rpeta(1e6, df1, df2, peta.h1, N)
  sig <- random.p > a
  exaggeration <- mean(random.p[sig]) / peta.h1
  
  data.frame(power = power, p.value = p.value, exaggeration = exaggeration, row.names = "Result:")
}

#=====================================================================================================================================


typem.anova.fun <- function(df1, df2, N, peta.h0 = 0, peta.min = 0, peta.max = .5, alpha = .05){
  
  peta <- function(df1, df2, peta.h1, peta.h0, N, alpha){
    
    a <- qpeta(alpha, df1, df2, peta.h0, N, lower.tail = FALSE)
    power <- ppeta(a, df1, df2, peta.h1, N, lower.tail = FALSE)
    random.p <- rpeta(1e4, df1, df2, peta.h1, N)
    sig <- random.p > a
    exaggration <- mean(random.p[sig]) / peta.h1
    
    list(power = power, exaggration = exaggration)
  }
  
  peta_range = seq(peta.min, peta.max, by = .001)
  n = length(peta_range)
  power = numeric(n)
  exaggration = numeric(n)
  
  for(i in 1L:n){
    g = peta(peta.h1 = peta_range[i], df1 = df1, df2 = df2, N = N, alpha = alpha, peta.h0 = peta.h0)
    power[i] = g$power
    exaggration[i] = g$exaggration
  }
  
  plot(power, exaggration, type = "l", ylim = c(1, 10), xaxt = "n", yaxt = "n", lwd = 2, font.lab = 2, col = 4)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  axis(2, at = seq(1, 10, by = 2), las = 1)
  abline(h = 1, v = alpha, col = 8)
}


#========================================================================================================================


type.sm <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA, digits = 6)
{
  UseMethod("type.sm")
}

type.sm.default <- function(d = .1, obs.d = .6, n1 = 20, n2 = NA, digits = 6){
  
  graphics.off()  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), xpd = TRUE)  
  
  N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
  
  min.d = qt(1e-4, df)*d.SE  ;  max.d = qt(0.9999, df)*d.SE
  
  `d|H0` = curve( dt(x/d.SE, df)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, type = "n", yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i")
  
  CI = qt(c(.025, .975), df)*d.SE
  
  x = seq(min.d, CI[1], l = 1e4) ;  y = dt(x /d.SE, df)/d.SE
  xx = seq(max.d, CI[2], l = 1e4) ; yy = dt(xx/d.SE, df)/d.SE
  
  polygon(c(min.d,  x, CI[1]), c( 0,  y, 0), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(0, yy, 0), col = 2, border = NA)   
  
  lines(`d|H0`, lwd = 2)
  
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  
  legend("topright", "Observed \nS.S. Effect", pch = 23, pt.bg = 3, pt.cex = 1.2, bty = "n", text.font = 2, adj = c(0, .3))   
  abline(v = 0, col = 2, xpd = FALSE) 
  
  par(mar = c(5, 4, 1, 2))
  
  `d|H1` = curve( dt(x/d.SE, df, ncp)/d.SE, min.d, max.d, n = 1e4, xlab = "Effect Size", 
                  ylab = NA, font = 2, font.lab = 2, yaxt = "n", bty = "n",
                  cex.axis = 1, cex.lab = 1, yaxs = "i", ty = "n")
  
  x = seq(min.d, CI[1], l = 1e4)   ;  y = dt(x /d.SE, df, ncp)/d.SE
  xx = seq(max.d, CI[2], l = 1e4)   ; yy = dt(xx/d.SE, df, ncp)/d.SE 
  
  polygon(c(min.d,  x, CI[1]), c( 0,  y, 0), col = 2, border = NA)  
  polygon(c(max.d, xx, CI[2]), c(0, yy, 0), col = 2, border = NA)  
  
  lines(`d|H1`, lwd = 2)
  
  axis(1, at = d, col = 4, col.axis = 4, font = 2)
  points(obs.d, 0, pch = 23, bg = 3, cex = 1.4, xpd = TRUE)
  abline(v = d, col = 4, xpd = FALSE)
  
  segments(c(CI[1], CI[2]), 0, c(CI[1], CI[2]), 20, lty = 2, col = 2, xpd = NA)
  
  type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
  power = type.s.area + pt(ifelse(d > 0, CI[2]/d.SE, CI[1]/d.SE), df, ncp, lower.tail = ifelse(d > 0, FALSE, TRUE))
  type.s = type.s.area / power
  p.value = 2*pt(abs(obs.d)/d.SE, df, lower.tail = FALSE)
  random.d = rt(n = 1e6, df, ncp)*d.SE
  sig = if(d > 0) abs(random.d) > CI[2] else -abs(random.d) < CI[1]
  exaggration = if(d > 0) mean(abs(random.d)[sig])/ d else mean(-abs(random.d)[sig])/ d
  
  round(data.frame(exaggration = exaggration, type.s = type.s, power = power, Crit.d = CI[2], p.value = p.value, row.names = "Results:"), digits = digits)
}


#=======================================================================


type.sm.fun <- function(n1, n2 = NA, d.min = 0, d.max = 1.4, alpha = .05)
{
  UseMethod("type.sm.fun")
}

type.sm.fun.default <- function(n1, n2 = NA, d.min = 0, d.max = 1.4, alpha = .05){
  
  type.sm <- function(n1, n2, d, alpha){
    
    N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    d.SE = 1/sqrt(N) ; ncp = d*sqrt(N)
    CI = qt(c(alpha/2, 1-(alpha/2)), df)*d.SE
    
    type.s.area = pt(ifelse(d > 0, CI[1]/d.SE, CI[2]/d.SE), df, ncp, lower.tail = ifelse(d > 0, TRUE, FALSE))
    power = type.s.area + pt(ifelse(d > 0, CI[2]/d.SE, CI[1]/d.SE), df, ncp, lower.tail = ifelse(d > 0, FALSE, TRUE))
    type.s = type.s.area / power
    random.d = rt(1e4, df, ncp)*d.SE
    sig = if(d > 0) abs(random.d) > CI[2] else -abs(random.d) < CI[1]
    exaggration = if(d > 0) mean(abs(random.d)[sig])/ d else mean(-abs(random.d)[sig])/ d
    
    list(exaggration = exaggration, type.s = type.s, power = power)
  }
  
  graphics.off()  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mfrow = c(2, 1), mgp = c(2, .5, 0), mar = c(4, 4, 3, 2), las = 1)  
  
  d_range = seq(d.min, d.max, by = 5e-3)
  n = length(d_range)
  power = numeric(n)
  type.s = numeric(n)
  exaggration = numeric(n)
  
  for(i in 1L:n){
    a = type.sm(d = d_range[i], n1 = n1, n2 = n2, alpha = alpha)
    power[i] = a$power
    type.s[i] = a$type.s
    exaggration[i] = a$exaggration
  }
  
  plot(power, type.s, type = "l", xaxt = "n", lwd = 2, font.lab = 2, col = 2)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  abline(v = alpha, col = 8)
  plot(power, exaggration, type = "l", ylim = c(1, 10), xaxt = "n", yaxt = "n", lwd = 2, font.lab = 2, col = 4)
  axis(1, at = c(alpha, seq(.2, 1, by = .2)))
  axis(2, at = seq(1, 10, by = 2))
  abline(h = 1, v = alpha, col = 8)
  
}


#=================================================================================================================


model.ci <- function(fit, level = .95){
  
  est <- coef(fit)
  se <- summary(fit)$coef[, 2]
  
  n <- length(est)
  mat <- matrix(c(rep(-1, n), rep(1, n)), nrow = n)
  
  p <- (1 - level)/2
  z <- -qnorm(p)
  ci <- est + mat * (se * z)
  rownames(ci) <- names(est)
  col1 <- paste0(format(p * 1e2, nsmall = 1), "%")
  col2 <- paste0(format((1 - p) * 1e2, nsmall = 1), "%")
  colnames(ci) <- c(col1, col2)
  ci
}     



#===========================================================================================================================


pov.curve <- povcurve <- function(pov = seq(0, .5, .1), df1 = 3, df2 = 73, N = 78, biased = TRUE, labels = TRUE){
  
  if(!biased) pov <- exp2peta(pov, df1, df2, N)
  
  options(warn = -1) ; p = sort(pov)
  max.p = qpeta(.999999, df1, df2, max(p), N)  
  
  for(i in 1:length(p)){      
    H = curve(dpeta(x, df1, df2, p[i], N), 0, max.p, n = 1e3, xlab = "Effect Size (POV)",
              ylab = NA, type = "n", add = i!= 1, bty = "n", axes = FALSE, font.lab = 2, yaxs = "i")
    
    polygon(H, col = adjustcolor(i, .7), border = NA, xpd = NA)
    if(labels) text(p[i], max(H$y), bquote(bolditalic(H[.(i-1)])), pos = 3, xpd = NA)
    axis(1, at = round(p[i], 3), col = i, col.axis = i, font = 2)
    segments(p[i], 0, p[i], dpeta(p[i], df1, df2, p[i], N), lty = 3)
  }
}               


#===========================================================================================================================


d.curve <- dcurve <- function(d = seq(0,2,.5), n1 = 30, n2 = NA, biased = TRUE, labels = TRUE){
  
  
  if(!biased) d <- exp2d(d, n1, n2)
  options(warn = -1) ; d = sort(d)
  min.d = qcohen(1e-5, min(d), n1, n2)  ;  max.d = qcohen(.99999, max(d), n1, n2)  
  
  for(i in 1:length(d)){      
    H = curve(dcohen(x, d[i], n1, n2), min.d, max.d, n = 1e3, xlab = "Effect Size (d)", 
              ylab = NA, type = "n", add = i!= 1, bty = "n", axes = FALSE, font.lab = 2, yaxs = "i")
    
    polygon(H, col = adjustcolor(i, .7), border = NA, xpd = NA)
    if(labels) text(d[i], max(H$y), bquote(bolditalic(H[.(i-1)])), pos = 3, xpd = NA)
    axis(1, at = round(d[i], 3), col = i, col.axis = i, font = 2)
    segments(d[i], 0, d[i], dcohen(d[i], d[i], n1, n2), lty = 3)
  }
}

#===========================================================================================================================


compare.model <- function(..., digits = 1e2){
  
  L <- list(...)
  if(length(L) < 2) stop("You need to have a least '2' fitted models for a comparison.", call. = FALSE)
  names(L) <- substitute(...())
  combs <- t(combn(x = names(L), m = 2))
  
  p.value <- round(apply(combs, 1, function(i) pchisq(2 * abs(logLik(L[[i[2]]]) - logLik(L[[i[1]]])), df = abs(L[[i[1]]]$df.residual - L[[i[2]]]$df.residual), lower.tail = FALSE)), digits)
  result <- data.frame(combs, p.value)
  names(result) <- c("model.1", "model.2", "p.value")
  
  Sig. <- symnum(result$p.value, cut = c(0, .001, .01, .05, .1, 1), na = FALSE, symbols = c("***", "**", "*", ":-(", ":-(("), corr = FALSE)
  output <- cbind(result, Sig.)
  names(output)[4] <- " "
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 ':-(' 0.1 ':-((' 1\n-------------\n") 
  return(output)
}




#===========================================================================================================================      

CI.d <- function(d, n1, n2 = NA, conf.level = .95, CI){
  
  min.d <- min(qcohen(1e-4, CI[1], n1, n2), qcohen(1e-4, CI[2], n1, n2 ))
  max.d <- max(qcohen(.9999, CI[1], n1, n2), qcohen(.9999, CI[2], n1, n2))
  
  ylim <- c(0, max(dcohen(seq(min.d, max.d, l = 5e2), CI[1], n1, n2), dcohen(seq(min.d, max.d, l = 5e2), CI[2], n1, n2), na.rm = TRUE))
  
  L <- curve( dcohen(x, CI[1], n1, n2), min.d, max.d, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = "Cohen's d", font.lab = 2, mgp = c(1.5, .5, 0), ylim = ylim)
  U <- curve( dcohen(x, CI[2], n1, n2), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4) 
  abline(v = c(CI[1], CI[2], d), col = c(4, 2, 1), lty = 2 ) ; points(d, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
  
}


#===========================================================================================================================


CI.peta <- function(peta, df1, df2, N, conf.level = .95, CI){
  
  min.p <- min(qpeta(1e-5, df1, df2, CI[1], N), qpeta(1e-5, df1, df2, CI[2], N))
  max.p <- max(qpeta(.99999, df1, df2, CI[1], N), qpeta(.99999, df1, df2, CI[2], N))
  
  ylim <- c(0, max(dpeta(seq(0, 1, l = 5e2), df1, df2, CI[1], N), dpeta(seq(0, 1, l = 5e2), df1, df2, CI[2], N), na.rm = TRUE))
  
  L <- curve( dpeta(x, df1, df2, CI[1], N), min.p, max.p, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = bquote(eta[p]^2), font.lab = 2, mgp = c(1.8, .5, 0), ylim = ylim)
  U <- curve( dpeta(x, df1, df2, CI[2], N), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4) 
  abline(v = c(CI[1], CI[2], peta), col = c(4, 2, 1), lty = 2 ); points(peta, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
  
}


#===========================================================================================================================


CI.R2 <- function(R2, df1, df2, N, conf.level = .95, CI){
  
  min.r <- min(qpeta(1e-5, df1, df2, CI[1], N), qpeta(1e-5, df1, df2, CI[2], N))
  max.r <- max(qpeta(.99999, df1, df2, CI[1], N), qpeta(.99999, df1, df2, CI[2], N))
  
  ylim <- c(0, max(dpeta(seq(0, 1, l = 5e2), df1, df2, CI[1], N), dpeta(seq(0, 1, l = 5e2), df1, df2, CI[2], N), na.rm = TRUE))
  
  L <- curve( dpeta(x, df1, df2, CI[1], N), min.r, max.r, n = 5e2, col = 4, lwd = 2, xpd = TRUE, ylab = "Density", xlab = bquote(R^2), font.lab = 2, mgp = c(1.75, .5, 0), ylim = ylim)
  U <- curve( dpeta(x, df1, df2, CI[2], N), n = 5e2, col = 2, add = TRUE, lwd = 2, xpd = TRUE)
  lines(CI, c(0, 0), lend = 1, lwd = 4)
  abline(v = c(CI[1], CI[2], R2), col = c(4, 2, 1), lty = 2 ) ; points(R2, 0, pch = 21, bg = "cyan", col = "magenta", cex = 2)
  text(CI, c(max(L$y)/2, max(U$y)/2), round(CI, 3) , srt = 90, pos = 3, col = c(4, 2), font = 2)
}                       


#===========================================================================================================================


R2.ci <- function(R2, n.pred, N, f = NA, df1 = NA, df2 = NA, conf.level = .95, digits = 20, show = FALSE){ 
  
  if(is.na(df1)) df1 <- n.pred 
  if(missing(n.pred) & df1) n.pred <- df1
  if(is.na(df2)) df2 <- N - n.pred - 1
  if(missing(N)) N <- df1 + df2 + 1  
  if(missing(df2) & N) df2 <- N - df1 - 1
  
  a <- if(!missing(R2)){ peta.ci(peta = R2, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits)
  } else { peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits) }
  
  names(a)[1] <- "R2"
  
  if(show){
    
    r <- nrow(a)
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
    
    I <- eq(a$R2, df1, df2, N, conf.level)
    
    R2 <- I[[1]] ; df1 <- I[[2]] ; df2 <- I[[3]]; N <- I[[4]] ; conf.level <- I[[5]]
    
    for(i in 1:r) CI.R2(R2 = R2[i], df1 = df1[i], df2 = df2[i], N = N[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
    
  }
  
  return(a)
}

#===========================================================================================================================


peta.ci <- function(peta, f = NA, df1, df2, N, conf.level = .95, digits = 1e2, show = FALSE)
{
  UseMethod("peta.ci")
} 

peta.ci.default <- function(peta, f = NA, df1, df2, N, conf.level = .95, digits = 1e2, show = FALSE){
  
  ci <- Vectorize(function(peta, f, N, df1, df2, conf.level){
    
    q <- ifelse(is.na(f), peta2F(peta, df1, df2), f) 
    
    alpha <- (1 - conf.level)/2
    
    u <- function (ncp, alpha, q, df1, df2) {
      suppressWarnings(pf(q = q, df1 = df1, df2 = df2, ncp, lower.tail = FALSE)) - alpha
    }
    
    g <- try(uniroot(u, c(0, 1e7), alpha = alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(g, "try-error")) g <- 0
    h <- try(uniroot(u, c(0, 1e7), alpha = 1-alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(h, "try-error")) h <- 0
    I <- c(g, h)
    
    I <- I / (I + N)
    
    P.eta.sq <- if(is.na(f)) peta else F2peta(f, df1, df2)
    
    return(c(P.eta.sq = P.eta.sq, lower = I[1], upper = I[2], conf.level = conf.level, ncp = peta2ncp(P.eta.sq, N), F.value = q))
  })
  
  peta <- if(missing(peta)) NA else peta
  
  a <- round(data.frame(t(ci(peta = peta, f = f, N = N, df1 = df1, df2 = df2, conf.level = conf.level))), digits = digits)
  
  if(show){
    
    r <- nrow(a)
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
    
    I <- eq(a$P.eta.sq, df1, df2, N, conf.level)
    
    peta <- I[[1]] ; df1 <- I[[2]] ; df2 <- I[[3]]; N <- I[[4]] ; conf.level <- I[[5]]
    
    for(i in 1:r) CI.peta(peta = peta[i], df1 = df1[i], df2 = df2[i], N = N[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
    
  }
  
  return(a)
  
}

#===========================================================================================================================


d.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95, digits = 1e2, show = FALSE)
{
  UseMethod("d.ci")
}

d.ci.default <- function(d, t = NA, n1, n2 = NA, conf.level = .95, digits = 1e2, show = FALSE){
  
  ci <- Vectorize(function(d, t, n1, n2, conf.level){
    
    options(warn = -1)  
    alpha = (1 - conf.level)/2
    N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    d.SE = 1/sqrt(N)
    q = ifelse(is.na(t), d/d.SE, t)
    
    f <- function(ncp, alpha, q, df){
      alpha - suppressWarnings(pt(q, df, ncp, lower.tail = FALSE))
    }
    
    CI <- sapply(c(alpha, 1-alpha),
                 function(x) uniroot(f, interval = c(-1e7, 1e7), alpha = x, q = q, df = df, extendInt = "yes")[[1]]*d.SE)
    
    Cohen.d = ifelse(is.na(t), d, t*d.SE)
    
    return(c(Cohen.d = Cohen.d, lower = CI[1], upper = CI[2], conf.level = conf.level, ncp = q))
  })
  
  d <- if(missing(d)) NA else d
  
  a <- round(data.frame(t(ci(d = d, t = t, n1 = n1, n2 = n2, conf.level = conf.level))), digits = digits)
  
  if(show){
    
    r <- nrow(a)
    graphics.off()
    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    if(r > 1) { par(mfrow = n2mfrow(r)) ; set.margin2() }
    
    I <- eq(a$Cohen.d, n1, n2, conf.level)
    
    d <- I[[1]] ; n1 <- I[[2]] ; n2 <- I[[3]]; conf.level <- I[[4]]
    
    for(i in 1:r) CI.d(d = d[i], n1 = n1[i], n2 = n2[i], conf.level = conf.level[i], CI = c(a$lower[i], a$upper[i]))
    
  }
  
  return(a)
  
}                               


#===========================================================================================================================


dcohen <- function(x, dbase = 0, n1, n2 = NA, log = FALSE){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  dt(x*sqrt(N), df, ncp, log = log)*sqrt(N)
}

#=======================================================================================================================================


qcohen <- function(p, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, dbase, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- dbase*sqrt(N)
    
    qt(p, df, ncp, lower.tail = lower.tail, log.p = log.p)/sqrt(N)
  })
  q(p = p, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}


#=======================================================================================================================================


pcohen <- function(q, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, dbase, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- dbase*sqrt(N)
    
    pt(q*sqrt(N), df, ncp, lower.tail = lower.tail, log.p = log.p)
  })
  p(q = q, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}


#=======================================================================================================================================


rcohen <- function(n, dbase = 0, n1, n2 = NA){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  rt(n, df, ncp)/sqrt(N)
}


#=========================================================================================================================================

dpeta <- function(x, df1, df2, pbase = 0, N, log = FALSE){
  x[x > .9999999] <- .9999999
  x[x < 0] <- 0
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- x / (1 - x) * d
  stats::df(f, df1, df2, ncp, log = log) * d * ( 1 / (1 - x) + x / (1 - x)^2 )
}


#=========================================================================================================================================


ppeta <- function(q, df1, df2, pbase = 0, N, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, df1, df2, pbase, N, lower.tail, log.p){
    
    q[q > .9999999] <- .9999999
    q[q < 0] <- 0
    pbase[pbase > .9999999] <- .9999999
    pbase[pbase < 0] <- 0
    ncp <- (pbase * N) / (1 - pbase)
    d <- df2 / df1
    f <- q / (1 - q) * d
    stats::pf(f, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
  })
  p(q = q, df1 = df1, df2 = df2, pbase = pbase, N = N, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


qpeta <- function(p, df1, df2, pbase = 0, N, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, df1, df2, pbase, N, lower.tail, log.p){
    
    p[p > 1] <- 1
    p[p < 0] <- 0
    pbase[pbase > .9999999] <- .9999999
    pbase[pbase < 0] <- 0
    ncp <- (pbase * N) / (1 - pbase)
    d <- df2 / df1
    f <- stats::qf(p, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
    f / (f + d)
  })
  q(p = p, df1 = df1, df2 = df2, pbase = pbase, N = N, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


rpeta <- function(n, df1, df2, pbase = 0, N){
  pbase[pbase > .9999999] <- .9999999
  pbase[pbase < 0] <- 0
  ncp <- (pbase * N) / (1 - pbase)
  d <- df2 / df1
  f <- stats::rf(n, df1, df2, ncp)
  f / (f + d)
}


#==================================================================================================================

dpetab <- function(x, df1, df2, ncp = 0, log = FALSE){
  x[x > .9999999] <- .9999999
  x[x < 0] <- 0
  d <- df2 / df1
  f <- x / (1 - x) * d
  stats::df(f, df1, df2, ncp, log = log) * d * ( 1 / (1 - x) + x / (1 - x)^2 )
}


#=========================================================================================================================================


ppetab <- function(q, df1, df2, ncp = 0, lower.tail = TRUE, log.p = FALSE){
  q[q > .9999999] <- .9999999
  q[q < 0] <- 0
  d <- df2 / df1
  f <- q / (1 - q) * d
  stats::pf(f, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
}


#=========================================================================================================================================


qpetab <- function(p, df1, df2, ncp = 0, lower.tail = TRUE, log.p = FALSE){
  p[p > 1] <- 1
  p[p < 0] <- 0
  d <- df2 / df1
  f <- stats::qf(p, df1, df2, ncp, lower.tail = lower.tail, log.p = log.p)
  f / (f + d)
}


#=========================================================================================================================================


rpetab <- function(n, df1, df2, ncp = 0){
  d <- df2 / df1
  f <- stats::rf(n, df1, df2, ncp)
  f / (f + d)
}



#=========================================================================================================================================

anova.es <- function(fit = NULL, f, df1, df2, N, conf.level = .9, digits = 6)
{
  UseMethod("anova.es")
}


anova.es.default <- function(fit = NULL, f, df1, df2, N, conf.level = .9, digits = 6){
  
  if(!is.null(fit)){
    
    if(class(fit)[1] != "lm" && class(fit)[1] != "aov") { stop("Error: 'fit' must be a fitted model from base R's 'aov()' or 'lm()' commands.") }        
    N <- nobs(fit)
    sit <- anova(fit)
    #if(!("F value" %in% names(sit))) { stop("Error: Fitted model does not include any 'F value'.") } 
    f <- head(sit[,4], -1)
    df1 <- head(sit$Df, -1)
    df2 <- tail(sit$Df, 1)
  }
  
  if(length(f) != length(df1)){message("Warning: The length of 'f' and 'df1' must be equal. Check your inputted values.\n")}
  I <- eq(f, df1) 
  f = I[[1]] ; df1 = I[[2]]
  
  omega <- (df1 * (f - 1)) / as.numeric(crossprod(df1, f) + df2 + 1)
  eta <- (df1 * f) / as.numeric(crossprod(df1, f) + df2)
  pomega <- (df1 * (f - 1)) / ((df1 * (f - 1)) + N)
  peta <- peta.ci(f = f, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = digits)
  
  result <- round(data.frame(F.value = f, eta.sq = eta, P.eta.sq = peta[,1], lower.P.eta.sq = peta[,2], 
                             upper.P.eta.sq = peta[,3], conf.level = conf.level, omega.sq = omega, 
                             P.omega.sq = pomega, row.names = paste0("effect ", 1:length(f), ":")), digits = digits)
  
  message("Note: If analysis includes random-effects, manually input the right 'df2' to obtain correct 'P.eta- or P.omega-sq.'")
  
  if(is.null(fit)){  
    
    return(result)
    
  }else{
    
    rownames(result) <- head(rownames(sit), -1)
    
    return(result)
  } 
}


#=========================================================================================================================================
                  

plan.t.test <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, d.range = seq(.1, .5, .05),
                        two.tailed = TRUE, xlab = "Cohen's d", xlim = c(NULL, NULL), ylim = NULL)
{
  
  UseMethod("plan.t.test")
}


plan.t.test.default <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, d.range = seq(.1, .5, .05),
                                two.tailed = TRUE, xlab = "Cohen's d", xlim = c(NULL, NULL), ylim = NULL){
  
  graphics.off()  
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  d[d == 0] <- 1e-4
  if(d < -6) d <- -6  
  if(d > 6) d <- 6
  if(power == 0) power <- sig.level
  
  from <- xlim[1]
  to <- xlim[2]
  
  d2 <- d
  d3 <- d.range
  d <- abs(d)
  sig.level <- if(two.tailed) sig.level/2 else sig.level
  k <- base.rate / (1 + base.rate)
  options(warn = -1)
  
  f <- if(two.tailed){ function(x){
    
    power - (pt(qt(sig.level, df = x), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
  }
    
  } else {
    
    function(x){
      
      power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
      
    }
  }
  
  df <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "yes")[[1]])
  
  n1 <- df + 1
  n2 <- if(paired) NA else round(base.rate*n1)
  
  loop <- length(d3)
  
  dfb <- numeric(loop)
  n1b <- numeric(loop)
  n2b <- numeric(loop)
  
  for(i in 1:loop){
    
    f <- if(two.tailed){ function(x){
      
      power - (pt(qt(sig.level, df = x), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
    }
      
    } else {
      
      function(x){
        
        power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d3[i]*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
        
      }
    }
    
    dfb[i] <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
    
    n1b[i] <- dfb[i] + 1
    n2b[i] <- if(paired) NA else round(base.rate*n1b[i])
  }
  
  base.rate <- if(paired) NA else base.rate
  method <- paste(if(paired) "One- or Paired sample" else "Two-sample", "t test power analysis")
  note <- paste("Use 'base.rate' to specify how many times one group might be larger than the other e.g.,'base.rate = 1.1' ")
  
  a <- qcohen(sig.level, 0, n1, n2)
  b <- -a
  
  est.power <- if(two.tailed) pcohen(a, d, n1, n2) + pcohen(b, d, n1, n2, lower.tail = FALSE) else pcohen(b, d, n1, n2, lower.tail = FALSE)
  
  from <- if(is.null(from)) min(qcohen(1e-5, 0, n1, n2), qcohen(1e-5, d2, n1, n2), na.rm = TRUE) else from
  to <- if(is.null(to)) max(qcohen(.99999, 0, n1, n2), qcohen(.99999, d2, n1, n2), na.rm = TRUE) else to
  
  x <- seq(from, to, 1e-4)
  ylimb <- c(0, max(c(dcohen(x, 0, n1, n2), dcohen(x, d2, n1, n2)), na.rm = TRUE) )
  
  ylim <- if(is.infinite(ylimb[2]) & is.null(ylim)) NULL else if(is.null(ylim)) ylimb else ylim
  
  par(mfrow = c(2, 1), mgp = c(2.5, .5, 0), mar = c(4, 4, 2, 2), tck = -.02)
  
  h0 <- curve(dcohen(x, 0, n1, n2), from, to, n = 1e3, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", ylim = ylim, font.lab = 2)
  
  x1 <- seq(from, a, length.out = 1e3) ; y1 <- dcohen(x1, 0, n1, n2) 
  x2 <- seq(b, to, length.out = 1e3) ; y2 <- dcohen(x2, 0, n1, n2)
  
  if(d2 < 0 & !two.tailed || two.tailed) polygon(c(from, x1, a), c(0, y1, 0), col = adjustcolor(2, .25), border = NA)
  if(d2 > 0 & !two.tailed || two.tailed) polygon(c(b, x2, to), c(0, y2, 0), col = adjustcolor(2, .25), border = NA) 
  lines(h0, lwd = 2, col = 2, xpd = TRUE)
  
  g <- if(d2 < 0 & !two.tailed) a else if(d2 > 0 & !two.tailed) b else c(a, b)
  p <- if(two.tailed) rep(par('usr')[4], 2) else par('usr')[4]
  abline(v = g, col = 2, lty = 2)
  
  crit <- round(g, 4) 
  
  points(g, p, pch = 19, col = 2, xpd = NA)  
  
  text(g, par('usr')[4], paste("critical d =", crit), pos = 3, cex = .7, font = 2, xpd = TRUE)
  
  h1 <- curve(dcohen(x, d2, n1, n2), from, to, n = 1e3, add = TRUE) 
  x1 <- seq(from, a, length.out = 1e3) ; y1 <- dcohen(x1, d2, n1, n2) 
  x2 <- seq(b, to, length.out = 1e3) ; y2 <- dcohen(x2, d2, n1, n2)
  if(d2 < 0 & !two.tailed || two.tailed) polygon(c(from, x1, a), c(0, y1, 0), border = NA, density = 15, col = 4, xpd = TRUE)
  if(d2 > 0 & !two.tailed || two.tailed) polygon(c(b, x2, to), c(0, y2, 0), border = NA, density = 15, col = 4, xpd = TRUE) 
  lines(h1, lwd = 2, col = 4, xpd = TRUE)
  
  legend("topleft", legend = c("Sig. Area(s)", "Power"), inset = c(-.15, 0), density = c(NA, 35), x.intersp = c(.3, .3),
         bty = "n", xpd = NA, cex = .7, text.font = 2, angle = c(NA, 45), fill = c(adjustcolor(2, .4), 4), border = c(2, 4), adj = c(0, .4))
  
  plot(d3, n1b, type = "b", pch = 19, lwd = 2, xlab = xlab, las = 1, col = 4, font.lab = 2, ylab = "Group Sample Size", ylim = c(0, max(n1b, n2b, na.rm = TRUE)), xaxt = "n")
  axis(1, at = d3)
  if(!paired)lines(d3, n2b, col = 2, lty = 3, lwd = 2, type = "b")
  
  points(if(!paired)rep(d, 2) else d, if(!paired) c(n1, n2) else n1, col = "magenta", bg = "cyan", pch = 21, cex = 1.5)
  
  text(d3, n1b, n1b, pos = 1, font = 2, col = 4, cex = .7, xpd = TRUE)
  if(!paired) text(d3, n2b, n2b, pos = 3, font = 2, col = 2, cex = .7, xpd = TRUE)
  
  if(paired){
    
    legend("topright", legend = "Group 1", col = 4, pch = 19, cex = .7, text.font = 2, lwd = 1,
           pt.cex = 1, bty = "n")
  } else {
    
    legend("topright", paste("Group", 1:2), col = c(4, 2), pch = c(19, 1), cex = .7, text.font = 2, x.intersp = c(.6, .6),
           lwd = 1, adj = c(0, .4), pt.cex = 1, pt.lwd = c(1, 2), lty = c(1, 3), bty = "n")
  }
  
  box()
  sig.level <- if(two.tailed) sig.level*2 else sig.level
  two.tailed <- if(two.tailed) "Yes" else "No"
  
  structure(list(n1 = n1, n2 = n2, base.rate = base.rate, d = d, est.power = est.power, sig.level = sig.level, 
                 two.tailed = two.tailed, method = method, note = note), class = "power.htest")
}



#=====================================================================================================================================



plan.f.test <- function(pov, n.level, design, sig.level = .05, n.covar = 0, n.pred = NULL, power = .8, peta.range = seq(1e-1, .9, 1e-1),
                        xlab = NULL, ylim = NULL, to = NULL, d = NA)
{
  
  UseMethod("plan.f.test")
}


plan.f.test.default <- function(pov, n.level, design, sig.level = .05, n.pred = NULL, n.covar = 0, power = .8, peta.range = seq(1e-1, .9, 1e-1),
                                xlab = NULL, ylim = NULL, to = NULL, d = NA){
  
  graphics.off()  
  original.par <- par(no.readonly = TRUE)
  on.exit(par(original.par))
  options(warn = -1)
  regress <- if(!is.null(n.pred)) TRUE else FALSE
  if(regress) n.level <- n.pred  
  
  peta2 <- peta.range
  
  peta2[peta2 == 0] <- 1e-2
  peta2[peta2 == 1] <- .99
  
  if(!is.na(d) & !regress) { message("\nNote: You are doing reseach planning for 'pairwise' comparisons.") ;  n.level <- design <- 2 }
  if(!is.na(d)) pov <- d2peta(d = d, n1 = 300, n2 = 300) 
  peta <- pov
  
  if(any(n.level < 1)) stop("Error: You must have at least '2 levels' or '1 predictor' for regression.", call. = FALSE)
  if(any(n.level <= 1) & !regress) stop("Error: You must have at least '2 levels'.", call. = FALSE)
  xlab <- if(is.null(xlab) && !regress) bquote(eta[p]^2) else if(is.null(xlab) && regress) bquote(bold(R^2)) else xlab
  if(!regress && missing(design)) stop("'design' must be numerically specified e.g., 'design = 2*4'.", call. = FALSE)
  if(regress){ n.level <- n.level + 1 ; design <- n.level }
  df1 <- n.level - 1
  n.covar[n.covar < 0] <- 0
  x <- sapply(list(n.level, design, n.covar), round)
  n.level <- x[1] ; design <- x[2] ; n.covar <- x[3]
  
  f <- function(x){
    
    power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta * (x + design + n.covar) ) /(1 - peta), lower.tail = FALSE))
  }
  
  df2 <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
  
  N <- df2 + design + n.covar
  
  loop <- length(peta2)
  
  Nb <- numeric(loop)
  df2b <- numeric(loop)
  
  for(i in 1:loop){
    
    f <- function(x){
      
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta2[i] * (x + design + n.covar) ) /(1 - peta2[i]), lower.tail = FALSE))
    }
    
    df2b[i] <- ceiling(uniroot(f, c(1e-8, 1e6), extendInt = "downX")[[1]])
    
    
    Nb[i] <- df2b[i] + design + n.covar
    
  }
  
  a <- qpeta(sig.level, df1, df2, 0, N, lower.tail = FALSE)
  
  to <- if(is.null(to)) max(qpeta(.999999, df1, df2, 0, N), qpeta(.999999, df1, df2, peta, N), na.rm = TRUE) else to
  x <- seq(0, to, 1e-4)
  ylimb <- c(0, max(dpeta(x, df1, df2, 0, N), dpeta(x, df1, df2, peta, N), na.rm = TRUE))
  
  ylim <- if(is.infinite(ylimb[2]) & is.null(ylim)) NULL else if(is.null(ylim)) ylimb else ylim
  
  est.power <- ppeta(a, df1, df2, peta, N, lower.tail = FALSE)
  
  par(mfrow = c(2, 1), mgp = c(2.5, .5, 0), mar = c(4, 4, 2, 2), tck = -.02)
  
  h0 <- curve(dpeta(x, df1, df2, 0, N), from = 0, to = to, n = 1e4, xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", ylim = ylim, font.lab = 2)
  
  x = seq(a, to, length.out = 1e3) ; y = dpeta(x, df1, df2, 0, N)
  polygon(c(a, x, to), c(0, y, 0), col = adjustcolor(2, .25), border = NA)
  lines(h0, lwd = 2, col = 2, xpd = TRUE)
  abline(v = a, col = 2, lty = 2) ; crit <- round(a, 4) ; points(a, par('usr')[4], pch = 19, col = 2, xpd = NA)
  
  es <- if(regress) bquote(R^2) else bquote(eta[p]^2)
  
  text(a, par('usr')[4], bquote(bold("critical"~ .(es) == .(crit)~"or"~.(crit*1e2)*"%")), pos = 3, cex = .7, font = 4, xpd = TRUE)
  
  h1 <- curve(dpeta(x, df1, df2, peta, N), from = 0, to = to, n = 1e4, add = TRUE) # xlab = xlab, ylab = NA, yaxt = "n", bty = "n", yaxs = "i", main = bquote(bolditalic(H[1]))
  x <- seq(a, to, length.out = 1e3) ; y <- dpeta(x, df1, df2, peta, N)
  polygon(c(a, x, to), c(0, y, 0), border = NA, density = 15, col = 4, xpd = TRUE)
  lines(h1, lwd = 2, col = 4, xpd = TRUE)
  
  legend("topleft", legend = c("Sig. Area", "Power"), inset = c(-.15, 0), density = c(NA, 35), x.intersp = c(.3, .3),
         bty = "n", xpd = NA, cex = .7, text.font = 2, angle = c(NA, 45), fill = c(adjustcolor(2, .4), 4), border = c(2, 4), adj = c(0, .4))
  
  plot(peta2, Nb, type = "b", lwd = 2, xlab = xlab, las = 1, col = "green4", font.lab = 2, xaxt = "n", ylab = "Total Sample Size")
  axis(1, at = peta2)
  
  legend("topright", legend = bquote(bold("Current required"~ bolditalic("\"N\""))), col = "magenta", pt.bg = "cyan", pch = 21, cex = .7,
         pt.cex = 1.2, bty = "n")
  box()
  points(peta, N, col = "magenta", bg = "cyan", pch = 21, cex = 1.5)
  
  text(peta2, Nb, Nb, pos = 3, font = 2, col = "gray40", cex = .8, xpd = TRUE)
  
  method <- paste("fixed-effects", if(regress) "Regression" else if(n.covar == 0) "ANOVA" else "ANCOVA", "power analysis") 
  
  balannced.N <- if(!regress) ceiling(N/design) * design else NA
  
  n.level <- if(regress) n.level-1 else n.level
  design <- if(regress) n.level else design
  
  r  <- structure(list(method, pov, est.power, a, sig.level, n.covar, design, n.level, df1, df2, N, balannced.N), class = "power.htest")
  
  setNames(r, c("method", ifelse(regress, "R-squared", "peta squared"), "est.power", ifelse(regress, "crit.Rsq", "crit.peta"), 
                "sig.level", "n.covar", "design", ifelse(regress, "n.pred", "n.level"), "df1", "df2", "total.N", "balanced.N"))
}


#=====================================================================================================================================


plan.mrm <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                     rho = .5, d = NA)
{
  
  UseMethod("plan.mrm")
}


plan.mrm.default <- function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level = .05, n.covar = 0, power = .8, eps = .9,
                             rho = .5, d = NA){
  
  if(!anyNA(d)) peta <- d2peta(d = d, n1 = 300, n2 = 300)
  
  options(warn = -1)
  
  G <- Vectorize(function(peta, n.rep, n.group, factor.type = c("between", "within", "bw"), sig.level, n.covar, power, eps, rho, d){
    
    factor.type <- match.arg(factor.type)
    
    m <- n.rep
    if(rho <= 0) rho <- 1e-7 else if(rho >= 1) rho <-.9999999
    if(eps < .5) eps <- .5 else if(eps > 1) eps <- 1
    if(any(n.group < 1)) stop("You must have at least '1 group' in your design.", call. = FALSE)
    if(any(m < 1)) stop("Incorrect # of measurements, change 'n.rep'.", call. = FALSE)
    if(any(factor.type != "between" & m < 2)) stop("You must have at least '2 repeated measurements' in your design.", call. = FALSE)
    if(any(factor.type == "between" & n.group < 2)) stop("You must have at least '2 groups' in your design.", call. = FALSE)
    if(any(missing(n.group))) stop("'n.group' must be numerically specified.", call. = FALSE)
    peta <- if(missing(peta)) NA else peta
    if(n.covar < 0) n.covar <- 0
    g <- sapply(list(n.group, n.covar, m), round)
    n.group <- g[1] ; n.covar <- g[2] ; m <- g[3]
    
    
    df1 <- switch(factor.type, between = n.group - 1, within = (m - 1)*eps, bw = (n.group - 1)*(m - 1)*eps)
    
    u <- if(factor.type == "between") m / (1 + (m - 1)*rho) else m / (1 - rho)
    
    f <- if(factor.type == "between"){ function(x){
      
      power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( x + n.group + n.covar) ) /(1 - peta))*u, lower.tail = FALSE))
    } 
      
    } else {
      
      function(x){ 
        power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = ((peta * ( ((x)/(m-1)) + n.group + n.covar) ) /(1 - peta))*eps*u, lower.tail = FALSE))
      }
    }
    
    df2 <- uniroot(f, c(1e-8, 1e3), extendInt = "yes")[[1]]
    
    df2 <- if(factor.type == "between") ceiling(df2) else df2
    
    N <- if(factor.type == "between") ceiling(df2 + n.group + n.covar)  else ceiling((df2 / ((m - 1)*eps)) + n.group + n.covar) 
    
    balanced.N <- if(factor.type == "between") ceiling(N/n.group) * n.group else NA
    
    a <- qpetab(sig.level, df1, df2, 0, lower.tail = FALSE)
    
    ncp <- if(factor.type == "between") (peta2f(peta)^2)*N*u else (peta2f(peta)^2)*N*u*eps
    
    est.power <- ppetab(a, df1, df2, ncp, lower.tail = FALSE)
    
    ro <- sapply(list(a, est.power), round, 4)
    
    list(peta = peta, total.N = N, balanced.N = balanced.N, factor.type = factor.type, n.group = n.group, n.rep = n.rep, n.covar = n.covar, sig.level = sig.level, crit.peta = ro[1], est.power = ro[2])
  })
  
  if(!anyNA(d) & any(n.group == 2)) message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.\n")    
  if(!anyNA(d) & any(n.group == 1)) message("\nNote: For 'pairwise' comparisons, 'total.N' is for '1' group.\n") 
  
  a <- data.frame(t(G(peta = peta, n.rep = n.rep, n.group = n.group, factor.type = factor.type, sig.level = sig.level, 
                      n.covar = n.covar, power = power, eps = eps, rho = rho, d = d)))
  
  names(a)[1] <- if(!is.na(d)) "d" else "peta"
  a[, 1] <- if(is.na(d)) peta else d
  a
}

#===========================================================================================================================


power.t <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, two.tailed = TRUE)
{
  
  pwr <- Vectorize(function(d, sig.level, power, base.rate, paired, two.tailed)
  {
    
    d <- abs(d)
    sig.level <- if(two.tailed) sig.level/2 else sig.level
    k <- base.rate / (1 + base.rate)
    options(warn = -1)
    
    f <- if(two.tailed){ function(x){
      
      power - (pt(qt(sig.level, df = x), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
    }
      
    } else {
      
      function(x){
        
        power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
        
      }
    }
    
    df <- ceiling(uniroot(f, c(1e-8, 1e7), extendInt = "downX")[[1]])
    
    n1 <- df + 1
    n2 <- if(paired) NA else round(base.rate*n1)
    
    return(c(n1 = n1, n2 = n2))
  })
  
  data.frame(t(pwr(d = d, sig.level = sig.level, power = power, base.rate = base.rate, paired = paired, two.tailed = two.tailed)))
}



#=================================================================================================================

plan.t.ci <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%")
{
  UseMethod("plan.t.ci")
}


plan.t.ci.default <- function(d, t = NA, n1, n2 = NA, conf.level = .95, width = NA, base.rate = 1, paired = FALSE, assure = .99, expect = FALSE, reduce.by = "0%", increase.by = "0%"){
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(is.na(width) & missing(n1) || is.na(width) & is.na(t) & missing(d) || is.na(width) & !paired & missing(n2)) stop("Either provide 'width' or provide 't or d', 'n1' and/or 'n2' from prior study.", call. = FALSE)  
  if(!is.na(t)) d <- t2d(t = t, n1 = n1, n2 = n2)
  if(is.na(width)) width <- d.width(d = d, t = t, n1 = n1, n2 = n2, conf.level = conf.level)
  if(expect) assure <- .5
  
  inc <- if(is.character(increase.by)) as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 else increase.by
  red <- if(is.character(reduce.by)) as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2 else reduce.by
  
  fac <- if(inc != 0 & red == 0) { 1 + inc
  } else if(red != 0 & inc == 0) { 1 - red 
  } else { 1 }
  
  
  if(fac <= 0 || inc == 0 & fac > 1) fac <- 1
  
  width <- width * fac
  
  G <- Vectorize(function(d, conf.level, width, base.rate, paired, assure, expect){
    
    n.d <- function(d, conf.level, width, base.rate, paired, assure){
      
      alpha <- (1 - conf.level)/2
      k <- base.rate 
      
      f <- function(ncp, alpha, d, df){
        alpha - suppressWarnings(pt(d*sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)), df, ncp, lower.tail = FALSE))
      }
      
      dbase <- function(df){
        sapply(c(alpha, 1 - alpha),
               function(x) uniroot(f, c(-5e1, 5e1), alpha = x, d = d, df = df, extendInt = "yes")[[1]]/sqrt(if(paired) df + 1 else ((k/(1 + k))^2)*(df + 2)))
      }
      
      m <- function(df, width){
        abs(abs(diff(dbase(df))) - width)
      }
      
      df <- optimize(m, c(1, 1e7), width = width)
      
      if(round(df$objective, 4) != 0) stop("Impossible planning: change input values.", call. = FALSE)
      
      n1 <- ceiling(if(paired) df[[1]] + 1 else (df[[1]] + 2)/(1 + k))
      n2 <- if(paired) NA else round(k * n1)
      
      list(d = d, n1 = n1, n2 = n2, base.rate = base.rate, width = width, conf.level = conf.level, assure = assure, paired = paired)
    }
    
    n <- n.d(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
    
    a <- d.ci(d = d, n1 = n$n1, n2 = n$n2, conf.level = c(assure, 2*assure - 1))$upper
    
    dnew <- function(dnew = dnew, n1 = n$n1, n2 = n$n2, d = d, assure = assure){
      total <- sum(pcohen(c(-dnew, dnew), d, n1 = n1, n2 = n2, lower.tail = c(TRUE, FALSE)))
      return(abs(total - (1 - assure)))
    }
    
    dnew <- optimize(dnew, a, d = d, assure = assure)[[1]]
    n.d(d = dnew, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure)
  })
  
  if(paired) base.rate <- NA
  a <- data.frame(t(G(d = d, conf.level = conf.level, width = width, paired = paired, base.rate = base.rate, assure = assure, expect = expect)), row.names = NULL)
  a[,1] <- d
  a                                                                                                          
}


#===================================================================================================================================================


plan.f.ci <- function(pov, design = 2 * 2, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95, width = NA, assure = .99, expect = FALSE, reduce.by = "0%", d = NA, lower, upper, increase.by = "0%", tol = 1e3)
{
  
  UseMethod("plan.f.ci")
  
}

plan.f.ci.default <- function(pov, design = 2 * 2, f = NA, n.level = 2, n.pred = NULL, n.covar = 0, conf.level = .95, width = NA, assure = .99, expect = FALSE, reduce.by = "0%", d = NA, lower, upper, increase.by = "0%", tol = 1e3){
  
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  if(expect) assure <- .5
  regress <- if(!is.null(n.pred)) TRUE else FALSE
  if(regress) n.level <- n.pred
  # if(!is.na(d)) { pov <- d2peta(d = d, n1 = 500, n2 = 500) ; n.level <- 2 ;
  # message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.") }
  #if(!is.na(d)) { pov <- d2peta(d = d, n1 = 500, n2 = 500) ; n.level <- 2 }
  # if(!is.na(d) & is.na(width)) width <- d.width.meta(lower = lower, upper = upper)
  # if(!is.na(d) & width >= .3) width <- .3
  # if(!is.na(d) & pov <= .15) pov <- .15
  
  #------------- 
  
  if(!anyNA(d)) { pov <- d2peta(d = d, n1 = 500, n2 = 500) ; n.level <- 2 ;
  message("\nNote: For 'pairwise' comparisons, 'total.N' is for '2' groups.") }
  
  width <- ifelse(!is.na(d) & is.na(width), d.width.meta(lower = lower, upper = upper), width)
  width <- ifelse(!is.na(d) & width >= .3,  .3, width)
  pov <- ifelse(!is.na(d) & pov <= .15,  .15, pov)   
  
  #--------------
  
  
  peta <- pov
  
  inc <- if(is.character(increase.by)) as.numeric(substr(increase.by, 1, nchar(increase.by)-1))/ 1e2 else increase.by
  red <- if(is.character(reduce.by)) as.numeric(substr(reduce.by, 1, nchar(reduce.by)-1))/ 1e2 else reduce.by
  
  fac <- if(inc != 0 & red == 0) { 1 + inc
  } else if(red != 0 & inc == 0) { 1 - red 
  } else { 1 }
  
  
  if(fac <= 0 || inc == 0 & fac > 1) fac <- 1
  
  width <- width * fac
  
  
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, expect){
    
    
    n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress){
      
      alpha <- (1 - conf.level)/2
      if(regress){ n.level <- n.level + 1 ; design <- n.level }
      df1 <- n.level - 1
      if(n.covar < 0) n.covar <- 0
      options(warn = -1)
      
      f <- function(alpha, q, df1, df2, ncp){
        alpha - suppressWarnings(pf(peta2F(peta, df1, df2), df1, df2, ncp, lower.tail = FALSE))
      }
      
      pbase <- function(df2){      
        
        b <- sapply(c(alpha, 1 - alpha), function(x) 
          tryCatch(uniroot(f, c(0, 1e7), alpha = x, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], error = function(e) NA))
        if(any(is.na(b))) b <- c(1, tol)     
        ncp2peta(b, df2 + design + n.covar)
      }
      
      m <- function(df2, width){
        abs(diff(pbase(df2))) - width
      }
      
      df2 <- uniroot(m, c(0, 1e3), width = width, extendInt = "yes")
      
      
      if(round(df2$f.root, 3) != 0) return( c(1, message("\nimpossible planning found (denoted 'NA'): You may change your 'width'."))[1])
      
      
      df2 <- ceiling(df2[[1]])
      
      N <- ceiling(df2 + design + n.covar)
      n.covar <- if(n.covar == 0) NA else n.covar
      n.level <- if(regress) n.level-1 else n.level
      design <- if(regress) n.level else design
      df1 <- if(regress) n.level else df1
      
      list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2, design = design)
    }
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)  
    
    peta <- exp.peta(pbase = n$peta, df1 = n$df1, df2 = n$df2, N = n$total.N)
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)
    
    peta.max <- root(pov = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = conf.level)$m
    
    a <- peta.ci(peta = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 2*assure - 1)
    
    nLU <- try(sapply(c(a$lower, a$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N), silent = TRUE)
    
    if(inherits(nLU, "try-error")) return(c(peta = peta, total.N = NA, width = width, n.level = n.level, design = design, conf.level = conf.level, max.width = NA))
    
    NN1 <- max(nLU, na.rm = TRUE)
    
    b <- peta.ci(peta = peta.max, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 1 - assure)
    
    nLU <- sapply(c(b$lower, b$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar)$total.N)
    
    NN2 <- max(nLU, na.rm = TRUE)
    
    NN3 <- if(!(peta.max %inn% c(a$lower, a$upper))) NN1 else max(NN1, NN2)
    
    max.w <- round(root(pov = peta, df1 = n$df1, N = NN3, df2 = NN3 - design - n.covar, conf.level = conf.level)$w, 3)
    
    return(c(peta = peta, total.N = NN3, width = width, n.level = n.level, design = design, conf.level = conf.level, max.width = max.w))
    
  })
  
  a <- data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, assure = assure, regress = regress, expect = expect)), regress = regress, assure = assure, row.names = NULL)
  names(a)[4] <- if(regress) "n.pred" else "n.level"
  names(a)[1] <- if(regress) "R2" else if(!is.na(d)) "d" else "peta"
  a[, 1] <- if(is.na(d)) pov else d
  a
}                                                 


#=====================================================================================================================================================

d.unbias <- function(d, n1, n2 = NA, t = NA){
  
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  d <- if(!is.na(t)) t/sqrt(N) else d
  d * exp(lgamma(df/2)-log(sqrt(df/2)) - lgamma((df-1)/2))
}


#=======================================================================================================================================


graph <- function(x, y = NULL, type = "p", xlim = NULL, ylim = NULL, 
                  log = "", main = NULL, sub = NULL, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), 
                  ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, 
                  panel.last = NULL, asp = NA, add = FALSE, show = TRUE, ...)
{
  
  if(!add){ 
    
    if(!show){type <- "n"; axes <- FALSE; ann <- FALSE} 
    
    graphics.off()  
    
    plot(x = x, y = y, type = type, xlim = xlim, ylim = ylim, 
         log = log, main = main, sub = sub, xlab = xlab, ylab = ylab, 
         ann = ann, axes = axes, frame.plot = frame.plot, panel.first = panel.first, 
         panel.last = panel.last, asp = asp, ...)
    
  }else{
    
    lines(x = x, y = y, type = type, ...)
  }
}


#========================================================================================================================================================


peta2d.fun <- function(peta = seq(.1, .5, .1), n = seq(30, 300, 10), base.rate = 1, xlab = "Group Sample Size", ylab = "Cohen's d", ylim = NA, ...)
{
  
  n <- sort(n)
  peta <- sort(peta)  
  
  d <- lapply(peta, function(x) peta2d(x, n1 = n, n2 = base.rate*n))
  ylim <- if(is.na(ylim)) range(d) else ylim
  for(i in 1:length(peta)){
    graph(n, d[[i]], type = "l", add = i!= 1, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    text(mean(n), mean(d[[i]]), bquote(eta[p]^2 == .(round(peta[i], 3))), col = 2, pos = 3, xpd = NA, cex = .8)
  }
}



#========================================================================================================================================================



"%inn%" <- function(x = 3.5, interval = c(3, 5)){
  
  r <- range(interval, na.rm = TRUE)
  
  x >= r[1] & x <= r[2] 
}

#=============================================================================================================================================================================================================================


root <- function(pov = .6, df1 = 3, df2 = 108, N = 100, conf.level = .95, show = FALSE, ...){
  
  f <- function(x){ 
    
    ci <- peta.ci(peta = x, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = 1e2)
    
    abs(ci$upper - ci$lower)
  }
  
  m <- optimize(f, 0:1, maximum = TRUE)[[1]]
  
  est <- uniroot(function(x) f(pov) - f(x), if(pov >= m) c(0, m) else c(m, 1))[[1]]
  
  if(show) curve(f, panel.f = abline(v = c(pov, est), h = f(pov), col = 2, lty = c(2, 1, 1)), ...) 
  
  list(m = m, est = est, w = f(m))
}


#=============================================================================================================================================================================================================================


power.f <- Vectorize(function(peta = .2, n.level = 2, design = 2 * 3, sig.level = .05, n.covar = 0, power = .8, regress = FALSE, pair.design = NULL){
  
  if(n.level <= 1) stop("Error: You must have at least '2 levels' or '2 predictors'.")
  #if(!regress & missing(design)) stop("Error: 'design' must be numerically specified e.g., 'design = 2 * 4'.")
  if(regress){ n.level <- n.level + 1 ; design <- n.level }
  df1 <- n.level - 1
  if(n.covar < 0) n.covar <- 0
  x <- sapply(list(n.level, design, n.covar), round)
  n.level <- x[1] ; design <- x[2] ; n.covar <- x[3]
  
  f <- function(x){
    
    power - suppressWarnings(pf(qf(sig.level, df1 = df1, df2 = x, lower.tail = FALSE), df1 = df1, df2 = x, ncp = (peta * (x + design) ) /(1 - peta), lower.tail = FALSE))
  }
  
  df2 <- ceiling(uniroot(f, c(1, 1e6), extendInt = "yes")[[1]]) - n.covar
  
  N <- df2 + design + n.covar
  
  bal <- ceiling(N/design) * design
  
  N <- if(!is.null(pair.design)) pair.design * (bal/2) else N
  
  return(N)
  
})


#====================================================================================================================================


Anova2 <- function(eta.sq = .25, n = 5, min.score = 0, max.score = 25, coef = 1.2, sig.level = .05, ...){
  
  beta = qnorm(c(1e-16, .9999999999999999))
  q = c(min.score, max.score)
  musd = solve(cbind(1L, beta), q)  
  m1 = musd[[1]]  
  sdw = musd[[2]] 
  
  x = ((sdw^2) * eta.sq) / (1 - eta.sq)
  m2 = coef * m1
  A = m1 + m2
  
  a = function(m3) abs((m1 - (A + m3)/3)^2 + (m2 - (A + m3)/3)^2 + (m3 - (A + m3)/3)^2 - 3*x)
  
  m3 = optimize(a, c(-3*max.score, 3*max.score))[[1]]
  
  mus = c(m1, m2, m3)
  k = length(mus)
  sigb = var(mus)*(k-1)/k
  eta.p = sigb / (sigb + sdw^2)
  group = gl(k, n, labels = paste0("GROUP ", 1:k))
  y = as.vector(mapply(rnorm, n = rep(n, k), mean = mus, sd = rep(sdw, k)))
  sim = anova(aov(y ~ group))
  eta = sim[1, 2] / sum(sim[, 2])
  msb = sim[1, 3]
  mse = sim[2, 3]
  omega = (sim[1, 2] - sim[1, 1]*mse) / ((sim[1, 2] + sim[2, 2]) + mse)
  eps = (sim[1, 2] - sim[1, 1]*mse) / (sim[1, 2] + sim[2, 2])
  lab = c(paste0("subj #", rev(n)[1]), paste0(rep(".", n - 2)), paste0("subj #", 1L))
  
  
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  
  set.margin2()
  
  par(font.lab = 2, font = 2, mgp = c(1.4, .2, 0), ...)
  dotchart(y, groups = group, gcol = 2:(k+1), pch = 19, color = (2:(k+1))[group],
           xlab = "Participants' Scores", labels = rep(lab, k)) 
  
  ms = unique(ave(y, group))
  sd = round(tapply(y, group, sd), 3)
  g = rev(cumsum(rev(tapply(group, group, length)) + 2) - 1)
  u = par("usr")
  
  segments(ms, c(u[4], g[2:3]), ms, c(g[2:3], u[3]), col = 2:(k+1), lty = 2)
  arrows(ms[1:2], g[2:3], ms[2:3], g[2:3], code = 3, length = .12, angle = 40, lwd = 2)
  
  ms = round(ms, 3)
  legend("topright", paste0("Mean = ", ms[1], "\n", "sd = ", sd[1]), bty = "n", text.col = "red4")
  legend("left", paste0("Mean = ", ms[2], "\n", "sd = ", sd[2]), bty = "n", text.col = "darkgreen")
  legend("bottomleft", paste0("Mean = ", ms[3], "\n", "sd = ", sd[3]), bty = "n", text.col = 4)
  legend("right",  paste0(" Eta Sq. = ", signif(eta, 3)), bty = "n")
  
  uneq = function(a, b, sig = 4) round(a, sig) != round(b, sig)
  if(uneq(eta.sq, eta.p)) message("\nNote: You may need to change 'coef' to get requested 'eta.sq'.")
  p = power.anova.test(groups = k, n = n, between.var = msb, within.var = mse, sig.level = sig.level)[[6]]
  
  cbind(sim, Pop.Eta.Sq = c(eta.p, NA), Eta.Sq = c(eta, NA), Omega.Sq = c(omega, NA), Epsilon.Sq = c(eps, NA), Power = c(p, NA))
}

#======================================================================================================================================


random <- function(N, n.groups, keep = FALSE){
  
  r <- round(c(N, n.groups))
  N <- r[1]
  n.groups <- r[2]
  
  n <- N
  N <- ceiling(N/n.groups) * n.groups
  n.set <- N/n.groups
  
  if(keep) set.seed(9036)
  a <- sample(rep(seq_len(n.groups), n.set))
  b <- table(if(n != N) head(a, -(N - n)) else a, dnn = NULL)
  if(n != N) a[(n+1):N] <- NA
  a <- matrix(a, n.groups, n.set)
  
  a[] <- sprintf("%s%d:~%d", "p", seq_len(N), a)
  rownames(a) <- rep(" ", n.groups)
  colnames(a) <- rep(" ", n.set)
  
  list(Assignments = noquote(a), Groups = b)
}


#===========================================================================================================================


set.margin2 <- function() 
{
  par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + .1, 
      tck = -0.02)
}

#===========================================================================================================================

denscurve <- function(..., adjust = 1, na.rm = TRUE, n = 1e3, hdi = FALSE, level = .95, xlab = "x", ylim = NA, xlim = NA, labels = NA, bottom = 1, top = 1, scale = 1){
  
  L <- if(all(sapply(list(...), inherits, "data.frame"))) as.list(...) else list(...)
  lab <- if(all(sapply(list(...), inherits, "data.frame"))) names(L) else substitute(...())
  
  loop <- length(L)
  soop <- seq_len(loop)
  
  a <- lapply(L, function(x) density(x, adjust = adjust, na.rm = na.rm, n = n))
  
  from <- numeric(loop)
  to <- numeric(loop)
  hi <- numeric(loop)
  if(hdi) CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
  
  for(i in soop){
    from[i] <- min(a[[i]]$x)
    to[i] <- max(a[[i]]$x)
    hi[i] <- max(a[[i]]$y)
    if(hdi) CI[i,] <- hdir(L[[i]], level = level)
    mode[i] <- a[[i]]$x[which.max(a[[i]]$y)]
  }
  
  f = hi + soop
  m = scale*hi + soop
  
  plot(rep(soop, 2), rep(soop, 2), type = "n", xlim = if(is.na(xlim)) c(min(from), max(to)) else xlim, ylim = if(is.na(ylim)) c(bottom*1, top*max(f)) else ylim, ylab = NA, yaxt = "n", xlab = xlab, mgp = c(2, .3, 0))
  axis(2, at = soop, labels = if(is.na(labels)) lab else labels, font = 2, las = 1, cex.axis = .8, tck = -.012, mgp = c(2, .3, 0), padj = rep(.35, loop))
  abline(h = soop, col = 8, lty = 3)
  
  for(i in soop){
    polygon(x = a[[i]]$x, y = scale*a[[i]]$y +i, col = adjustcolor(i, .4), border = NA, xpd = NA)
  }
  
  if(hdi){   
    segments(CI[, 1], soop, CI[, 2], soop, lend = 1, lwd = 4, col = soop, xpd = NA)                            
    segments(mode, soop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, soop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = decimal(CI, 2); o = decimal(mode, 2)
    text(c(CI[,1], o, CI[,2]), soop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
  }  
  return(invisible(a))
}        


#=====================================================================================

isFALSE <- function(x) identical(FALSE, x)
    
       
#=====================================================================================       
 
       
samp.dist <- function(n, pop.dist = c('nor','exp','uni','poi','bin','gam','chi','tds', 'bet'), param1 = NULL, param2 = NULL, times = 1e4, others = FALSE, xlab = NA, ylab = "Density", obs.mean = NULL, rev.page = FALSE, seed.pop = 0, digits = 3, xaxt = "s", reset = FALSE, unknown = FALSE, obs.sd = NULL, ...){
  
  pop.dist <- match.arg(pop.dist)
  
  if(unknown & is.null(obs.sd)) stop("Provide 'obs.sd' to see the place of your observed test statistic.", call. = FALSE)
  
  samples <- switch(pop.dist,
                    "exp" = replicate(times, rexp(n, param1)),
                    "nor" = replicate(times, rnorm(n, param1, param2)),
                    "uni" = replicate(times, runif(n, param1, param2)),
                    "poi" = replicate(times, rpois(n, param1)),
                    "bin" = replicate(times, rbinom(n, param1, param2)),
                    "gam" = replicate(times, rgamma(n, param1, param2)),
                    "chi" = replicate(times, rchisq(n, param1)),
                    "tds" = replicate(times, rt(n, param1)),
                    "bet" = replicate(times, rbetab(n, param1, param2)))
  
  set.seed(seed.pop)
  pop <- switch(pop.dist,
                "exp" = rexp(1e4, param1),
                "nor" = rnorm(1e4, param1, param2),
                "uni" = runif(1e4, param1, param2),
                "poi" = rpois(1e4, param1),
                "bin" = rbinom(1e4, param1, param2),
                "gam" = rgamma(1e4, param1, param2),
                "chi" = rchisq(1e4, param1),
                "tds" = rt(1e4, param1),
                "bet" = rbetab(1e4, param1, param2))
  
  
  all.sample.means <- colMeans(samples, na.rm = TRUE)   
  all.sample.sd <- apply(samples, 2, sd, na.rm = TRUE)
  
  if(reset){
  graphics.off()
  org.par <- par(no.readonly = TRUE)
  on.exit(par(org.par))
  }
  
  h <- if(others) 5 else 3
  dev <- if(!rev.page) n2mfrow(h) else rev(n2mfrow(h))
  par(mfcol = dev, mar = c(2.5, 2.6, 1.8, .5), mgp = c(1.65, .4, 0))
  
  
  pt.curve(pop, main = "Population", cex.main = 1, col = 4, ylab = ylab, xlab = xlab, pch = '.', xaxt = xaxt, ...)
  sd.pop <- sd(pop, na.rm = TRUE)
  m <- mean(pop, na.rm = TRUE)
  abline(v = c(m, m-sd.pop, m+sd.pop), col = 3, lty = c(2, 1,1))
  at1 <- axTicks(1)
  
  pt.curve(all.sample.means,main = "Sampling Distribution of Means (Normal)", cex.main = 1, xlim = range(pop, na.rm = TRUE), xlab = xlab, ylab = ylab, pch = '.', xaxt = "n", ...)
  at2 <- axTicks(1)
  if(xaxt == "s") axis(1, at = union(at1, at2), ...)
  
  se <- sd(all.sample.means, na.rm = TRUE)
  m <- mean(all.sample.means, na.rm = TRUE)
  if(!is.null(obs.mean)) {  
    points(obs.mean, 0, pch = 23, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
    obs.z <-  if(!unknown) (obs.mean-m) / (sd.pop/sqrt(n)) else (obs.mean-m) / (obs.sd/sqrt(n))
  }
  
  abline(v = c(m, m-se, m+se), col = 3, lty = c(2, 1,1))
  
  z <- if(!unknown)(all.sample.means - m) / (sd.pop/sqrt(n)) else (all.sample.means - m) / (all.sample.sd/sqrt(n))
  
  pt.curve(z, main = paste("Sampling Distribution of",if(unknown) "'t'" else "'z'","(Std. Normal)"), cex.main = 1, xlab = xlab, ylab = ylab, col = 'purple', pch = '.', xaxt = xaxt, xlim = if(!is.null(obs.mean)) range(z, obs.z, na.rm = TRUE) else range(z, na.rm = TRUE), ...)
  se.z <- sd(z, na.rm = TRUE)
  m.z <- mean(z, na.rm = TRUE)
  if(!is.null(obs.mean)) points(obs.z, 0, pch = 23, bg = "cyan", col = 'magenta', cex = 2, xpd = NA)
  qs <- quantile(z, names = F, probs = c(.025, .975), na.rm = TRUE)
  arrows(qs[1], .1, -3, .1, code = 2, length = .12, angle = 20)
  arrows(qs[2], .1, 3, .1, code = 2, length = .12, angle = 20)
  abline(v = c(m.z, qs), col = 3, lty = c(2, 1,1))
  
  if(others){
    
    all.sample.sums <- colSums(samples, na.rm = TRUE)
    #all.sample.vars <- apply(samples,2,var, na.rm = TRUE) 
    pt.curve(all.sample.sums, main = "Sampling Distribution of Sum", cex.main = 1, xlab = xlab, pch = ".", ylab = ylab, xaxt = xaxt, ...)
    pt.curve(all.sample.sd, main = "Sampling Distribution of Sd", cex.main = 1, xlab = xlab, pch = ".", ylab = ylab, xaxt = xaxt, ...) 
    #m.sd <- mean(all.sample.sd) ; med.sd <- median(all.sample.sd)
    abline(v = mean(all.sample.sd, na.rm = TRUE), col = 3)
  }
  
  obj <- round(c(sd.pop = sd.pop, se = se, clt.se = sd.pop / sqrt(n), obs.z = if(!is.null(obs.mean)) obs.z else NA), digits)
  names(obj)[4] <- if(unknown) "obs.t" else "obs.z"
  setNames(obj, names(obj))
}  
                 
#=====================================================================================
       
need22 <- c('car','psych','tidyverse','effects','ez', 'haven', 'rstatix', 'ggpubr')

not.have11 <- need22[!(need22 %in% installed.packages()[,"Package"])]
if(length(not.have11)) install.packages(not.have11)


suppressWarnings(
suppressMessages({ 
  
  for(i in need22){
    library(i, character.only = TRUE)
  }
}))       
