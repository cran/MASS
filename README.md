1. Fixed glm.nb() behavior for data that are effectively Poisson: when the dispersion parameter (theta) diverges, the routine now triggers a Poisson fallback implemented via theta.ml2.

2. Updated theta.ml() so that when estimation indicates a Poisson-like degenerate case (theta → ∞), it returns NA. The glm.nb2() wrapper consumes this NA and applies the Poisson fallback, preventing numeric errors and misleading non-convergence messages.
