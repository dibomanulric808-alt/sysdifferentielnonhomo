# =============================================================================
#  NON-HOMOGENEOUS DIFFERENTIAL SYSTEM:  X'(t) = A X(t) + B(t)
#  Presentation – Differential Equations
#
#  Contents
#    0. Utility functions (matrix exponential, Runge-Kutta 4, Duhamel)
#    1. Example 1: constant forcing term   (x' = y, y' = -2x - 3y + 2)
#    2. Example 2: forced damped oscillator (x'' + 0.4 x' + x = cos(1.5 t))
#    3. Application: two-compartment pharmacokinetics (infusion)
#
#  No external package required (base R only).
#  Usage: open in RStudio, then Ctrl+Shift+S (Source) or Ctrl+Enter.
#  Every computation is checked with stopifnot(): if the script runs to the
#  end, all consistency tests have passed.
# =============================================================================

rm(list = ls())
SAVE_FIGS <- TRUE      # TRUE: also saves the figures in ./figures/

# ------------------------------------------------------------------ 0. TOOLS
# Spectral decomposition A = V diag(l) V^-1 (A assumed diagonalizable)
spectre <- function(A) {
  e <- eigen(A)
  list(l = e$values, V = e$vectors, Vinv = solve(e$vectors))
}

# Matrix exponential: e^(tA) = V diag(e^(l t)) V^-1
expA <- function(t, sp) {
  Re(sp$V %*% (exp(sp$l * t) * sp$Vinv))
}

# Solution of the homogeneous system X' = AX with X(0) = X0
sol_homogene <- function(t, sp, X0) as.vector(expA(t, sp) %*% X0)

# 4th-order Runge-Kutta integrator: f(t, X) -> vector
# Returns a matrix: one row per time point of the grid t.
rk4 <- function(f, X0, t) {
  n <- length(t); X <- matrix(NA_real_, n, length(X0)); X[1, ] <- X0
  for (i in 1:(n - 1)) {
    h  <- t[i + 1] - t[i]
    k1 <- f(t[i],         X[i, ])
    k2 <- f(t[i] + h / 2, X[i, ] + h / 2 * k1)
    k3 <- f(t[i] + h / 2, X[i, ] + h / 2 * k2)
    k4 <- f(t[i] + h,     X[i, ] + h * k3)
    X[i + 1, ] <- X[i, ] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  X
}

# Duhamel's formula: X(t) = e^(tA) X0 + int_0^t e^((t-s)A) B(s) ds
# (integral computed with the trapezoidal rule, n subintervals)
duhamel <- function(t, sp, X0, Bfun, n = 4000) {
  s  <- seq(0, t, length.out = n + 1)
  g  <- sapply(s, function(si) as.vector(expA(t - si, sp) %*% Bfun(si)))
  g  <- matrix(g, ncol = n + 1)                      # dim x (n+1)
  h  <- s[2] - s[1]
  I  <- h * (rowSums(g) - 0.5 * (g[, 1] + g[, n + 1]))
  as.vector(expA(t, sp) %*% X0) + I
}

# Plots on screen and (optionally) to a PNG file
figure <- function(nom, expr) {
  e <- substitute(expr); env <- parent.frame()
  eval(e, env)
  if (SAVE_FIGS) {
    dir.create("figures", showWarnings = FALSE)
    png(file.path("figures", paste0(nom, ".png")), 900, 600, res = 120)
    eval(e, env); invisible(dev.off())
  }
}

cat("=================== EXAMPLE 1 ===================\n")
# ----------------------------------------------------- 1. EXAMPLE 1
#   x' = y
#   y' = -2x - 3y + 2        X(0) = (0, 0)
A  <- matrix(c(0, 1,
              -2, -3), nrow = 2, byrow = TRUE)
Bc <- c(0, 2)                      # constant forcing term
X0 <- c(0, 0)

cat("Matrix A:\n"); print(A)
sp <- spectre(A)
cat("\nEigenvalues:", sp$l, "\n")
cat("Eigenvectors (columns):\n"); print(sp$V)

# Step 1: constant particular solution  A Xp + B = 0  =>  Xp = -A^-1 B
Xp <- as.vector(-solve(A, Bc))
cat("\nParticular solution Xp = -A^-1 B =", Xp, "\n")
stopifnot(all(abs(A %*% Xp + Bc) < 1e-12))          # Xp is indeed an equilibrium

# Step 2: constants C such that X(0) = Xh(0) + Xp = X0
Cst <- as.vector(sp$Vinv %*% (X0 - Xp))
cat("Constants C (normalized eigenbasis):", Cst, "\n")

# General solution: X(t) = Xp + V diag(e^(l t)) C
sol_exacte <- function(t, X0 = c(0, 0)) {
  C <- as.vector(sp$Vinv %*% (X0 - Xp))
  Re(Xp + sp$V %*% (exp(sp$l * t) * C))
}

# Analytical formula derived in class:
x_an <- function(t) 1 - 2 * exp(-t) + exp(-2 * t)
y_an <- function(t) 2 * exp(-t) - 2 * exp(-2 * t)

tps <- seq(0, 5, by = 0.01)
Xspec <- t(sapply(tps, sol_exacte))
cat("\nMax error (spectral method / analytical formula):",
    max(abs(Xspec[, 1] - x_an(tps)), abs(Xspec[, 2] - y_an(tps))), "\n")
stopifnot(max(abs(Xspec[, 1] - x_an(tps))) < 1e-10,
          max(abs(Xspec[, 2] - y_an(tps))) < 1e-10)

# Check 2: numerical integration (RK4) of the full system
Xrk <- rk4(function(t, X) as.vector(A %*% X) + Bc, X0, tps)
err_rk <- max(abs(Xrk[, 1] - x_an(tps)), abs(Xrk[, 2] - y_an(tps)))
cat("Max error (RK4 / analytical formula)                :", err_rk, "\n")
stopifnot(err_rk < 1e-6)

# Check 3: Duhamel's formula at a few time points
for (tk in c(1, 2.5, 5)) {
  Xd <- duhamel(tk, sp, X0, function(s) Bc)
  cat(sprintf("Duhamel t = %.1f : (%.6f, %.6f)  | analytical: (%.6f, %.6f)\n",
              tk, Xd[1], Xd[2], x_an(tk), y_an(tk)))
  stopifnot(max(abs(Xd - c(x_an(tk), y_an(tk)))) < 1e-4)
}

# Check 4: the function satisfies  x' = y  (numerical derivative, midpoints)
dx <- diff(x_an(tps)) / diff(tps)
stopifnot(max(abs(dx - y_an(head(tps, -1) + diff(tps) / 2))) < 1e-4)

# --- Plot 1: x(t) and y(t)
figure("ex1_courbes", {
  plot(tps, x_an(tps), type = "l", lwd = 3, col = "#0F8B8D", ylim = c(0, 1.1),
       xlab = "t", ylab = "value",
       main = "Example 1: solutions x(t) and y(t)")
  lines(tps, y_an(tps), lwd = 3, col = "#E07A2F")
  abline(h = Xp[1], lty = 2, col = "grey50")
  grid()
  legend("right", c("x(t)", "y(t)", "equilibrium x = 1"),
         col = c("#0F8B8D", "#E07A2F", "grey50"), lwd = c(3, 3, 1), lty = c(1, 1, 2), bty = "n")
})

# --- Plot 2: decomposition Xh + Xp
figure("ex1_decomposition", {
  xh <- x_an(tps) - Xp[1]
  plot(tps, x_an(tps), type = "l", lwd = 3, col = "#0F8B8D", ylim = c(-1.1, 1.1),
       xlab = "t", ylab = "x", main = "Example 1: x(t) = xh(t) + xp")
  lines(tps, xh, lwd = 3, col = "#E07A2F")
  abline(h = Xp[1], lty = 2, col = "grey50")
  grid()
  legend("right", c("x(t) = xh + xp", "xh(t) (transient)", "xp = 1"),
         col = c("#0F8B8D", "#E07A2F", "grey50"), lwd = c(3, 3, 1), lty = c(1, 1, 2), bty = "n")
})

# --- Plot 3: phase portrait (several initial conditions)
figure("ex1_portrait_phase", {
  CI <- list(c(0, 0), c(0, 1), c(2, 1), c(-0.5, -1), c(2.2, -0.8), c(-1, 1.5))
  cols <- c("#0F8B8D", "#E07A2F", "#7A5CC8", "#3C9D5D", "#C2455A", "#1F2A44")
  tp <- seq(0, 6, by = 0.02)
  plot(NA, xlim = c(-1.2, 2.5), ylim = c(-1.2, 1.8), xlab = "x", ylab = "y",
       main = "Example 1: phase portrait")
  grid()
  for (k in seq_along(CI)) {
    tr <- t(sapply(tp, sol_exacte, X0 = CI[[k]]))
    lines(tr, col = cols[k], lwd = 2)
    points(CI[[k]][1], CI[[k]][2], pch = 16, col = cols[k])
  }
  points(Xp[1], Xp[2], pch = 4, cex = 2, lwd = 3)
  text(Xp[1], Xp[2], "equilibrium (1, 0)", pos = 3)
})

cat("\n=================== EXAMPLE 2 ===================\n")
# ----------------------------------------------------- 2. EXAMPLE 2
#   x'' + 0.4 x' + x = cos(1.5 t)      x(0) = x'(0) = 0
#   Let y = x'  =>  X' = A X + B(t),  A = [[0,1],[-1,-0.4]],  B(t) = (0, cos(w t))
w  <- 1.5
A2 <- matrix(c(0, 1,
              -1, -0.4), nrow = 2, byrow = TRUE)
b  <- c(0, 1)
B2 <- function(t) b * cos(w * t)
X02 <- c(0, 0)

sp2 <- spectre(A2)
cat("Eigenvalues:", format(sp2$l, digits = 4), "\n")
stopifnot(all(Re(sp2$l) < 0))                       # stable system

# Particular solution  Xp(t) = W1 cos(wt) + W2 sin(wt)
#   cos : w W2 = A W1 + b        sin : -w W1 = A W2
#   => (w I + A^2 / w) W2 = b   then   W1 = -A W2 / w
W2 <- as.vector(solve(w * diag(2) + A2 %*% A2 / w, b))
W1 <- as.vector(-A2 %*% W2 / w)
cat("W1 =", W1, "\nW2 =", W2, "\n")
Xp2 <- function(t) W1 * cos(w * t) + W2 * sin(w * t)

# Check: Xp2 satisfies the full equation (analytical derivative)
tt0 <- 0.7
dXp <- -w * W1 * sin(w * tt0) + w * W2 * cos(w * tt0)
stopifnot(max(abs(dXp - (A2 %*% Xp2(tt0) + B2(tt0)))) < 1e-12)

# Complete solution: X(t) = Xp(t) + e^(tA) (X0 - Xp(0))
sol2 <- function(t) as.vector(Xp2(t) + expA(t, sp2) %*% (X02 - Xp2(0)))

t2 <- seq(0, 40, by = 0.005)
Xan2 <- t(sapply(t2, sol2))
Xrk2 <- rk4(function(t, X) as.vector(A2 %*% X) + B2(t), X02, t2)
err2 <- max(abs(Xan2 - Xrk2))
cat("Max error (exact solution / RK4):", err2, "\n")
stopifnot(err2 < 1e-6)

# Theoretical amplitude of the steady state
amp_th <- 1 / sqrt((1 - w^2)^2 + (0.4 * w)^2)
amp_Xp <- sqrt(W1[1]^2 + W2[1]^2)
amp_num <- max(abs(Xan2[t2 > 30, 1]))
cat(sprintf("Steady-state amplitude: theoretical %.4f | via Xp %.4f | simulated %.4f\n",
            amp_th, amp_Xp, amp_num))
stopifnot(abs(amp_th - amp_Xp) < 1e-10, abs(amp_th - amp_num) < 1e-2)

figure("ex2_reponse", {
  plot(t2, Xan2[, 1], type = "l", lwd = 2.5, col = "#0F8B8D", ylim = c(-1.2, 1.2),
       xlab = "t", ylab = "x(t)", main = "Example 2: response of the forced oscillator")
  lines(t2, cos(w * t2), col = "grey70", lwd = 1.5)
  abline(h = c(-amp_th, amp_th), lty = 2, col = "#E07A2F")
  grid()
  legend("topright", c("x(t)", "forcing cos(1.5 t)", "steady-state amplitude"),
         col = c("#0F8B8D", "grey70", "#E07A2F"), lwd = c(2.5, 1.5, 1), lty = c(1, 1, 2), bty = "n")
})

figure("ex2_portrait_phase", {
  plot(Xan2[, 1], Xan2[, 2], type = "l", lwd = 1.5, col = "#7A5CC8",
       xlab = "x", ylab = "x'", main = "Example 2: phase portrait (limit cycle)")
  tp <- seq(0, 2 * pi / w, length.out = 200)               # steady-state cycle
  lines(t(sapply(tp, Xp2)), col = "#C2455A", lwd = 3)
  points(0, 0, pch = 16)
  grid()
  legend("topright", c("trajectory from (0,0)", "limit cycle Xp"),
         col = c("#7A5CC8", "#C2455A"), lwd = c(1.5, 3), bty = "n")
})

cat("\n=================== PHARMACOKINETICS ===================\n")
# ----------------------------------------------------- 3. PHARMACOKINETICS
#   x1: amount in the central compartment, x2: peripheral compartment
#   x1' = -(k10 + k12) x1 + k21 x2 + R      (R: infusion rate)
#   x2' =  k12 x1 - k21 x2
k10 <- 0.3; k12 <- 0.2; k21 <- 0.1; R <- 1
A3 <- matrix(c(-(k10 + k12),  k21,
                k12,         -k21), nrow = 2, byrow = TRUE)
B3 <- c(R, 0)
sp3 <- spectre(A3)
Xeq <- as.vector(-solve(A3, B3))
cat("Eigenvalues:", sp3$l, "\n")
cat("Equilibrium (plateau):", Xeq, "\n")
stopifnot(all(sp3$l < 0), max(abs(Xeq - c(10 / 3, 20 / 3))) < 1e-12)

t3 <- seq(0, 60, by = 0.1)
sol3 <- function(t) Re(Xeq + sp3$V %*% (exp(sp3$l * t) * as.vector(sp3$Vinv %*% (c(0, 0) - Xeq))))
X3 <- t(sapply(t3, sol3))
X3rk <- rk4(function(t, X) as.vector(A3 %*% X) + B3, c(0, 0), t3)
stopifnot(max(abs(X3 - X3rk)) < 1e-5)
cat("Max error (exact / RK4):", max(abs(X3 - X3rk)), "\n")

figure("pharmaco", {
  matplot(t3, X3, type = "l", lty = 1, lwd = 3, col = c("#0F8B8D", "#E07A2F"),
          xlab = "time", ylab = "amount of drug",
          main = "Constant infusion: two-compartment model")
  abline(h = Xeq, lty = 2, col = c("#0F8B8D", "#E07A2F"))
  grid()
  legend("right", c("central x1", "peripheral x2"),
         col = c("#0F8B8D", "#E07A2F"), lwd = 3, bty = "n")
})

cat("\n>>> All tests passed: the script works correctly.\n")

