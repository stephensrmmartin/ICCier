ICCier
================
Stephen R. Martin
2019-08-22

# Introduction

`ICCier` is a package for estimating and predicting group-specific or
observation-specific intraclass correlation coefficients (ICCs). In this
vignette, we describe the rationale for specific ICCs, and walk the
reader through the package with various examples. First, we provide some
background on ICCs and their variety of forms. In this section, we
introduce the reader to the notation used throughout the rest of the
vignette. Second, we describe the notion and rationale for varying ICCs
(e.g., group-specific ICCs). Third, we describe the Bayesian mixed
effects location scale model (MELSM) that underlies the estimation of
these ICCs. The MELSM is the “engine” that makes varying ICCs possible,
and we encourage all readers to read this section to best understand the
model output. Fourth, we walk through the `ICCier` functions, followed
by numerous examples for a variety of scenarios. Finally, we end with
some best practices and tips for using the `ICCier` package.

# Background

## Intraclass Correlation Coefficients

The intraclass correlation coefficient (ICC) is the proportion of total
variance that is due to between-group variance. It can also be thought
of as a reliability estimate, or the expected correlation between any
two scores within the same group. In its most basic form, it can be
expressed as: \[
\text{ICC(1)} = \frac{\sigma^2_b}{\sigma^2_b + \sigma^2_{\epsilon}},
\] where \(\sigma^2_b\) is the variance between groups, and
\(\sigma^2_\epsilon\) is the residual variance.

The ICC ranges from 0 to 1, with 0 meaning none of the variance is due
to between-group variation, and 1 meaning all is due to between-group
variation. The *groups* of data may be many things, such as trials
within *persons*, persons within *schools*, regions within *countries*,
and so on. For example, if one has repeated trials within persons, and
observes and ICC(1) of .45, then 45% of the variance in the observations
is due to differences between persons.

There are several ICCs, and `ICCier` supports many of them. The ICC(1)
is specified above, and represents the proportion of variability in the
*observations* due to the grouping factor. The ICC(2) represents the
proportion of variability in the *mean score* estimates that is due to
true group mean variability. It is specified as: \[
\text{ICC(2)} = \frac{\sigma^2_b}{\sigma^2_b + \frac{\sigma^2_\epsilon}{n}},
\] where \(n\) is the number of observations within each group. The
ICC(2) is the same as the ICC(1), except the second variance component
is the squared standard error of the mean (i.e., the sampling variance
of the mean). When \(n\) varies across group (i.e., unbalanced groups),
it is common to estimate the expected sampling variance by using the
average \(n_i\) (E.g., See Raudenbush & Bryk, 2002).

Thus far, we have only described the *unconditional* ICC. Unconditional
ICCs describe the proportion of variance due to group-varying *means*.
In unconditional ICCs, \(\sigma^2_b\) describes the variance of group
means, and \(\sigma^2_\epsilon\) is the within-group variance.
Conditional ICCs include a model on the expected values of the outcome.
That is, \(\sigma^2_\epsilon\) is no longer the within-group variance,
but the *residual variance* after accounting for the linear predictors.

Moreover, the \(\sigma^2_b\) is ill-defined for conditional ICCs,
because there is no longer a simple mean estimate for each group. There
are two broad approaches for the \(\sigma^2_b\) in conditional ICCs.
One, which we call *unadjusted* (and is default in `ICCier`), replaces
\(\sigma^2_b\) with \(\sigma^2_{u_0i}\), the variance of the
*intercept*. The corresponding ICC, labeled uICC(1), then represents the
proportion of variance (after controlling for the predictors) due to
groups, assuming the predictors are all at zero. If one is interested in
the ICC or reliability of conditional group means (i.e., the group means
when the predictors are zero), then this is the conditional ICC of
interest.

The second, which we call *adjusted* (and is used when `adjusted =
TRUE`), replaces \(\sigma^2_b\) with \(\sigma^2_{\hat u_{ij}}\), the
variance due to *all random coefficients*. That is,
\(\hat u_{ij} = u_{0i} + x_{ij}u_{1i} + \ldots\), and
\(\sigma^2_b = \text{Var}(\hat u_{ij})\). By comparison, the
*unadjusted* conditional ICC merely uses
\(\sigma^2_b = \text{Var}(u_{0i})\). This means that the aICC(1)
represents the proportion of total residual variance that is due to
group-specific variation, as a whole. Alternatively stated, it
represents the proportion of variance that is *predicted* by random
effects. If one is interested in individual differences (i.e., in
intercepts and slopes), then this is the conditional ICC of interest.

When no random slopes are present, adjusted and unadjusted ICCs are
equivalent. The ICC(1) and ICC(2) can be estimated for unconditional
models, and for *unadjusted* conditional models.

In sum, `ICCier` can handle unconditional and conditional ICCs with
single grouping factors. These include the ICC(1), ICC(2), uICC(1),
uICC(2), and aICC(1). This alone is not novel — various packages permit
the estimation of these ICCs. `ICCier`, however, can estimate these
quantities *per* group, and predict them *per* observation.
Consequently, `ICCier` can be used to probe sources of
reliability.

## Group and observation-specific ICCs

# The Underlying Model

# Function Overview

## `ICCier`

## `summary`

## `predict` and `fitted`

## `ranef` and `coef`

## `loo`

# Examples

## Unconditional, Intercept-only

## Unconditional, Scale Model

## Conditional, Intercept-only

## Conditional, Scale Model

## Level 2 Variables

# Model Comparison

# Best Practices

## Diagnostic Failures and You

## Scaling and Centering

## Pre-Filtering (Self: as in, using complete cases only, so fitted can be combined easily)

# References

<div id="refs" class="references">

<div id="ref-Raudenbush2002">

Raudenbush, S. W., & Bryk, A. S. (2002). *Hierarchical linear models:
applications and data analysis methods* (2nd ed.). Thousand Oaks.

</div>

</div>
