ICCier
================
Stephen R. Martin
2019-10-15

Introduction
============

Researchers often employ grouped, hierarchical, or repeated measures designs when measuring the construct of interest. Knowing that groups and persons vary, researchers may wish to express the proportion of variance in the measured construct that is due to the grouping variables. One such measure is the intraclass correlation coefficient (ICC).

Although the ICC can take many forms, it broadly represents the proportion of some variance due to underlying group variation. Consequently, the ICC can help answer many pertinent statistical questions. In multilevel models, it can indicate how much variance could be explained by group-level variables in hierarchical models (Raudenbush & Bryk, 2002). In psychometrics, it can represent reliability in repeatedly measured constructs (Shieh, 2016). In individual-difference research, it can represent the proportion of predictive variance due to between-person model variation (Johnson, 2014).

To compute the ICC, the variance components of interest must be estimated. For example, the ICC(1) is composed of the between-group variance and within-group variance:
$$
\\text{ICC(1)} = \\frac{\\sigma^2\_b}{\\sigma^2\_b + \\sigma^2\_\\epsilon}
$$
 A one-way mixed effects model (Raudenbush & Bryk, 2002) can estimate the corresponding variance components. These estimates are simply plugged into the equation to produce the estimated ICC(1).

However, this approach makes the strong assumption that the variance components are homogenous across groups and possible covariates. Namely, it assumes the residual variance is equal across all groups. This assumption is needlessly restrictive, and threatens the generalizability of estimated ICCs.

To relax the homogeneity assumption, the mixed effects location scale model (MELSM) can be used(Hedeker, Mermelstein, & Demirtas, 2008, 2012; Li & Hedeker, 2012; Rast, Hofer, & Sparks, 2012). The MELSM extends the mixed effects model to include random effects models on variance components. For example, a one-way mixed effects model includes randomly varying intercepts, but a single fixed residual variance. By contrast, the MELSM one-way mixed effects model includes both randomly varying intercepts, and randomly varying residual variances.

The MELSM does not only loosen the homogeneity assumption, but can explicitly model each variance component. In the context of the ICC, the MELSM allows one to predict both *σ*<sub>*b*</sub><sup>2</sup> and *σ*<sub>*ϵ*</sub><sup>2</sup>. Therefore, the ICC(1) can be modeled and predicted from covariates:
$$
\\text{ICC(1)}\_{ij} = \\frac{\\hat\\sigma^2\_{b,i}}{\\hat\\sigma^2\_{b,i} + \\hat\\sigma^2\_{\\epsilon,ij}} = \\frac{\\exp(\\boldsymbol{X}\_i\\boldsymbol{\\eta})^2}{\\exp(\\boldsymbol{X}\_i\\boldsymbol{\\eta})^2 + \\exp(\\boldsymbol{z}\_{ij}\\boldsymbol{\\gamma}\_{i})^2}
$$

This approach therefore allows researchers to model and predict the ICC, and with it access to new questions and design considerations. For example, researchers can examine which conditions, sampling paradigm, demographics, or number of trials produce the greatest or least reliability. Therefore, the researcher can plan future studies more efficiently. One may examine the extent to which ICCs vary across groups, and therefore whether reliability is equivalent across multiple groups. As this approach is relatively recent (Brunton-Smith, Sturgis, & Leckie, 2017; Hedeker et al., 2008; Williams, Martin, & Rast, 2019), we are excited to see how researchers employ it.

To this end, we developed `ICCier`. `ICCier` is a user-friendly package for estimating and predicting ICCs using the MELSM. <!-- `ICCier` is a package for estimating and predicting group-specific or observation-specific intraclass correlation coefficients (ICCs). --> <!-- In this vignette, we describe the rationale for specific ICCs, and walk the reader through the package with various examples. --> In this vignette, we walk the reader through the package using various examples, and provide further descriptions about the supported ICCs and the underlying model. First, we provide background on the supported ICCs. Second, we walk through the `ICCier` functions through numerous examples. Third, we describe the Bayesian mixed effects location scale model (MELSM) that underlies the estimation of these ICCs. Finally, we end with some best practices and tips for using the `ICCier` package.

Intraclass Correlation Coefficients
===================================

The intraclass correlation coefficient (ICC) represents the proportion of variance that is due to between-group variance. The ICC ranges from 0 to 1, with 0 meaning none of the variance is due to between-group variation, and 1 meaning all is due to between-group variation. The *groups* of data may be many things, such as trials within *persons*, persons within *schools*, regions within *countries*, and so on. It can also be thought of as a reliability estimate, or the expected correlation between any two scores within the same group. There are several ICC variants, and `ICCier` supports many of them. The ICC(1) represents the proportion of variability in the *observations* due to the grouping factor. The ICC(2) represents the proportion of variability in the *mean score* estimates that is due to true group mean variability. The uICC(1) and aICC(1) are conditional ICCs, which are ICCs after controlling for some covariates. These are all discussed in turn below.

ICC(1)
------

The ICC(1) can be expressed as:
$$
\\text{ICC(1)} = \\frac{\\sigma^2\_b}{\\sigma^2\_b + \\sigma^2\_{\\epsilon}},
$$
 where *σ*<sub>*b*</sub><sup>2</sup> is the between-group variance, and *σ*<sub>*ϵ*</sub><sup>2</sup> is the within-group variance.

The ICC(1) is the proportion of total variance in the observations that is due to between-group variance. For example, if one has repeated trials within persons, and observes and ICC(1) of .45, then 45% of the variance in the observations is due to differences between persons. The ICC(1) can be interpreted as the proportion of total variance that can be attributed to, and explained by, group-level differences and covariates (E.g., See Raudenbush & Bryk, 2002). Therefore, a low ICC(1) indicates that most explainable variance is due to within-group, and not between-group, variables. Additionally, the ICC(1) can be interpreted as a reliability measure for individual observations. If ICC(1) is high, then most variance in scores is due to between-person differences, rather than within-person fluctuations or measurement error.

`ICCier` extends this so that the between-group variance (*σ*<sub>*b*</sub><sup>2</sup>) and the within-group variance (*σ*<sub>*ϵ*</sub><sup>2</sup>) can be modeled. Specifically, `ICCier` defines the ICC(1) more generally as:

$$
\\text{ICC(1)}\_{ij} = \\frac{\\hat\\sigma^2\_{b,i}}{\\hat\\sigma^2\_{b,i} + \\hat\\sigma^2\_{\\epsilon,ij}} = \\frac{\\exp(\\boldsymbol{X\_i\\eta})^2}{\\exp(\\boldsymbol{X\_i\\eta})^2 + \\exp(\\boldsymbol{z\_{ij}\\gamma\_i})^2}
$$
 Therefore, an ICC(1) can be predicted for each group (i) and each observation therein (j). More details about the model implementation are available below (\[\#model\]).

ICC(2)
------

The ICC(2) can be expressed as:
$$
\\text{ICC(2)} = \\frac{\\sigma^2\_b}{\\sigma^2\_b + \\frac{\\sigma^2\_\\epsilon}{n}},
$$
 where *n* is the number of observations within each group. The ICC(2) is the same as the ICC(1), except the second variance component is the approximate sampling variance of the mean (i.e., the squared standard error of the mean). When *n* varies across groups (i.e., unbalanced groups), it is common to estimate the expected sampling variance by using the average *n*<sub>*i*</sub> (E.g., See Raudenbush & Bryk, 2002).

The ICC(2) is the proportion of variance in group mean estimates that is due to true between-group variance. Therefore, the ICC(2) is a reliability measure for group mean estimates. If ICC(2) is high, then most variance in group mean estimates is due to true group differences, rather than sampling error.

`ICCier` extends this similarly to the ICC(1) above.
$$
\\text{ICC(1)}\_{ij} = \\frac{\\hat\\sigma^2\_{b,i}}{\\hat\\sigma^2\_{b,i} + \\hat\\sigma^2\_{\\epsilon,ij}} = \\frac{\\exp(\\boldsymbol{X\_i\\eta})^2}{\\exp(\\boldsymbol{X\_i\\eta})^2 + \\frac{\\exp(\\boldsymbol{z\_{ij}\\gamma\_i})^2}{n\_{ij}}}
$$
 The *n*<sub>*i**j*</sub> can either be the number of observations within the group (*n*<sub>*i**j*</sub> = *n*<sub>*i*</sub>), or the number of observations up until the *j*-th observation. The former provides the predicted ICC(2) at the given covariates assuming the same number of observations are made. However, the latter allows researchers to examine how mean reliability changes across the number of observations. For example, the trial number can predict variance components, and also be used as *n*<sub>*i**j*</sub>. By doing so, the experimenter can estimate how mean score reliability changes over the course of the experiment.

uICC(1)
-------

Thus far, we have only described the *unconditional* ICC. Unconditional ICCs are called "unconditional" because the means are not conditioned on any covariates. Therefore, the between-group variance is the variance in unconditional means, and the within-group variance is the variance around each group mean. In contrast, *conditional* ICCs estimate the variance in *conditional* means, and the residual variance is the variance around the *conditional* mean. That is, conditional ICCs include a model on the expected values of the outcome. Therefore, *σ*<sub>*ϵ*</sub><sup>2</sup> is no longer the within-group variance, but the *residual variance* after accounting for the linear predictors. Moreover, the *σ*<sub>*b*</sub><sup>2</sup> is ill-defined for conditional ICCs, because there is no longer a simple mean estimate for each group.

There are two broad approaches for the *σ*<sub>*b*</sub><sup>2</sup> in conditional ICCs. One, which we call *unadjusted* (default in `ICCier`) replaces the between-group variance of the means (*σ*<sub>*b*</sub><sup>2</sup>) with the between-group variance in the intercept (*σ*<sub>*u*<sub>0*i*</sub></sub><sup>2</sup>). The corresponding ICC, labeled uICC(1), can be expressed as:
$$
\\text{uICC(1)}\_{ij} = \\frac{\\hat\\sigma^2\_{u\_0,i}}{\\hat\\sigma^2\_{u\_0,i} + \\hat\\sigma^2\_{\\epsilon,ij}} = \\frac{\\exp(\\boldsymbol{X\_i\\eta\_1})^2}{\\exp(\\boldsymbol{X\_i\\eta\_1})^2 + \\exp(\\boldsymbol{z\_{ij}\\gamma\_i})^2}
$$

The uICC(1) represents the proportion of variance (after controlling for the predictors) due to groups, assuming the predictors are all at zero. If one is interested in the ICC when predictor variables are at zero, then this is the conditional ICC of interest.

aICC(1)
-------

The second conditional ICC which we call *adjusted* (`adjusted = TRUE`) replaces the between-group variance of the means (*σ*<sub>*b*</sub><sup>2</sup>) with the between-group variance of the conditional means (*σ*<sub>*û*<sub>*i**j*</sub></sub><sup>2</sup>). Whereas the unadjusted ICC only includes the variance between intercepts, the adjusted ICC includes all variance due to between-group differences in the intercepts and slopes. That is, if *û*<sub>*i**j*</sub> = *u*<sub>0*i*</sub> + *x*<sub>*i**j*</sub>*u*<sub>1*i*</sub> + …, then *σ*<sub>*b*, *i*</sub><sup>2</sup> = Var(*û*<sub>*i**j*</sub>). The aICC(1) can be expressed as:
$$
\\text{aICC(1)}\_{ij} = \\frac{\\hat\\sigma^2\_{\\hat u,i}}{\\hat\\sigma^2\_{\\hat u,i} + \\hat\\sigma^2\_{\\epsilon,ij}} = \\frac{\\boldsymbol{\\mathrm{Tr}(m\_{i}\\Sigma\_i^{(L)}m\_{i}')}}{\\mathrm{Tr}(\\boldsymbol{m\_{i}\\Sigma\_i^{(L)}m\_{i}'}) + n\_i\\exp(\\boldsymbol{z\_{ij}\\gamma\_i)^2}}
$$
 , where **m**<sub>**i**</sub> is a matrix of conditional mean predictors for group *i*, and *Σ*<sub>*i*</sub><sup>(*L*)</sup> is the group-specific covariance matrix of mean coefficient random effects. Further details are available in the [model description](#model) section (See also Johnson, 2014). When no random slopes are present, adjusted and unadjusted ICCs are equivalent.

The aICC(1) represents the proportion of total residual variance that is due to group-specific variation, as a whole. Alternatively stated, it represents the proportion of variance that is *predicted* by all random effects. If one is interested in individual differences (i.e., in intercepts and slopes), then this is the conditional ICC of interest.

Examples
========

`ICCier`
--------

`ICCier` is the model-fitting function. It has two interfaces for specifying the model --- a multiple formula and a single-formula interface. This vignette focuses only on the multiple-formula interface, because it is more intuitive. With the multiple-formula interface, only three arguments are mandatory:

1.  *x*
    The name of the variable for which the ICC is to be computed.
2.  *group*
    The name of the grouping variable.
3.  *data*
    A data.frame containing the measured and grouping variables.

> Example
> -------
>
> `ICCier(x = MathAch, group = School, data = MathAchieve)`

When only these three arguments are provided, `ICCier` estimates one between-group variance, and per-group residual variances.

In order to model between-group variance, another formula is required. The formula should include group-level covariates on the right-hand side, with 'between' on the left-hand side, and is specified anywhere after the data argument.

> Example
> -------
>
> `ICCier(MathAch, School, MathAchieve, between ~ MEANSES)`

In this example, the between-group variance can vary across the `MEANSES` variable. Consequently, the ICC can be larger or smaller depending on the mean SES of the school.

To model the within-group variance, another formula is required. The within-group variance formula has 'within' on the LHS, and covariates on the RHS. The RHS formula is partitioned into two formulas: within-group variables and between-group formulas, separated by a pipe (|). Within-group variables are variables whose values vary within groups (i.e., level 1 variables), whereas between-group variables only vary between groups (i.e., level 2 variables).

> Example
> -------
>
> `ICCier(MathAch, School, MathAchieve, within ~ Sex | MEANSES)`

In this example, the residual variance for each observation is modeled from Sex (level 1), MEANSES (level 2), and their interaction (Sex:MEANSES). All level 1 variable coefficients are assumed to randomly vary and covary, and are predicted from level 2 variables. Consequently, each group is expected to have a different residual variance. Moreover, when predicting from within-group variables, an ICC is estimated for each observation. To exclude level 1 or level 2 predictors, simply include an intercept term instead. E.g., to model within-group variance only from MEANSES:

> Example
> -------
>
> `ICCier(MathAch, School, MathAchieve, within ~ 1 | MEANSES)`

Both between-group and within-group variances can be modeled simultaneously by including both formulas. E.g.,

> Example
> -------
>
> `ICCier(MathAch, School, MathAchieve, within ~ Sex + Minority | MEANSES, between ~ MEANSES + Sector)`

In order to estimate conditional ICCs, a mean model must be included. The mean model is specified using 'mean' on the LHS, and covariates on the RHS. Similar to the within-group variance model, the RHS is split into within-group and between-group predictors.

> Example
> -------
>
> `ICCier(MathAch, School, MathAchieve, mean ~ Sex | Sector)`

In this example, `MathAch` is predicted by Sex (level 1), Sector (Level 2), and their interaction (Sex:Sector). For simplicity, no between or within-group variance models were included. The within-group variance component is now the residual error variance from the conditional mean model. Because `adjusted` was not specified, the between-group variance component is the between-group variance of the intercepts (i.e., variance between groups when Sex and Sector are zero).

Other notable arguments to `ICCier` include: \* adjusted
When estimating conditional ICCs, to use the adjusted ICC (aICC, `TRUE`) or not (uICC, `FALSE`, default). \* iter
Number of MCMC iterations to run (including warmup). \* cores
Number of cores to use. By default, it will use the number of cores available. \* chains
Number of MCMC chains to run. By default, four MCMC chains are used.

Unconditional, Intercept-only
-----------------------------

Unconditional, Scale Model
--------------------------

Conditional, Intercept-only
---------------------------

Conditional, Scale Model
------------------------

Level 2 Variables
-----------------

Model Comparison
----------------

Best Practices
==============

Diagnostic Failures and You
---------------------------

Scaling and Centering
---------------------

Pre-Filtering (Self: as in, using complete cases only, so fitted can be combined easily)
----------------------------------------------------------------------------------------

The Underlying Model
====================

This section details the model implemented by the `ICCier` package. The ICC requires a between-group variance and a within-group variance. The mixed-effect location scale model (MELSM) is used to model these components. The MELSM employed by `ICCier` has three submodels.

The first submodel is the mean model. The mean model is a mixed-effects model for the mean given the within-group and between-group predictors. It is "maximal" because all within-group predictors have random coefficients that covary. When an unconditional ICC is estimated, the mean submodel is an intercept-only model, with group-varying means (i.e., a one-way random effects model).

The second submodel is the within-group (residual) variance model. Unlike traditional multilevel models, the MELSM does not assume a homogenous residual standard deviation (SD). Instead, the expected SD is modeled using a log-linear mixed effects model from within-group and between-group predictors. At minimum, `ICCier` uses an intercept-only model, with group-varying residual SDs.

Finally, the third submodel is the between-group variance model. The between-group random effect SDs are log-linearly modeled from between-group predictors. At minimum, `ICCier` uses an intercept-only model, such that each random effect SD is constant.

Model details
-------------

Let *y*<sub>*i**j*</sub> be the *j*th observation in the *i*th group response vector **y**<sub>*i*</sub>. Let **m**<sub>**i**</sub> represent the observation-level design matrix for the response variable, and **M**<sub>*i*</sub> represent a group-level row vector. The response variable is modeled as follows:
$$
\\begin{align\*}
    \\boldsymbol{y\_i} &= \\boldsymbol{m}\_i\\boldsymbol{\\beta}\_i + \\boldsymbol{\\epsilon}\_i \\\\
    \\boldsymbol\\beta\_i' &= \\boldsymbol{M}\_i\\boldsymbol{B} + \\boldsymbol{u}\_i^{(\\mu)'}
\\end{align\*}
$$
 This defines a maximal random effects model, such that all within-group variables have covarying random effects. All random coefficients are multivariate-linearly modeled from level-2 variables in the coefficient matrix **B** If *q* level-1 predictors are present, and *Q* level-2 predictors are present, then **B** is a *Q* × *q* matrix of coefficients, and **u**<sub>*i*</sub><sup>(*μ*)′</sup> is a *q*-length row vector of residual random effects.

The errors are assumed to be normally distributed with a modeled residual standard deviation:
$$
\\begin{align\*}
    \\epsilon\_{ij} &\\sim \\mathcal{N}(0,\\hat\\sigma\_{\\epsilon,ij}) \\\\
\\end{align\*}
$$
 The standard deviation, *σ̂*<sub>*ϵ*, *i**j*</sub>, is modeled similarly to the mean. Let **z**<sub>*i*</sub> represent the observation-level design matrix for the predicted residual SD, and **Z**<sub>*i*</sub> represent a group-level row vector. The residual SD is modeled as follows:
$$
\\begin{align\*}
    \\log\\boldsymbol{\\hat\\sigma}\_i &= \\boldsymbol{z}\_i\\boldsymbol{\\gamma}\_i \\\\
    \\boldsymbol{\\gamma}\_i' &= \\boldsymbol{Z}\_i\\boldsymbol{\\Gamma} + \\boldsymbol{u}\_i^{(\\sigma)}
\\end{align\*}
$$
 All level-1 variables have covarying random effects. All random coefficients are multivariate-linearly modeled from level-2 variables in the coefficient matrix **Γ**. If there are *p* level-1 predictors for the residual SD, and *P* group-level predictors, then **Γ** is a *P* × *p* matrix of coefficients, and **u**<sub>*i*</sub><sup>(*σ*)′</sup> is a *p*-length row vector of residual random effects.

The random effects of the mean model and residual SD model (**u**<sub>*i*</sub><sup>(*μ*)</sup>, **u**<sub>*i*</sub><sup>(*σ*)</sup>) are distributed multivariate normal with mean zero and a dense covariance matrix. However, because between-group SDs are modeled, the covariance matrix depends on the group. The between-group random-effect SDs are multivariate log-linearly modeled from a level-2 design matrix (**X**<sub>**i**</sub>) and a matrix of coefficients (**η**). The estimated between-group SDs are then included into a diagonal scale matrix, and combined with the correlation matrix (**Ω**) to produce the group-specific covariance matrix from which random effects are distributed.
$$
\\begin{align\*}
    \\begin{bmatrix}
        \\boldsymbol{u}\_i^{(\\mu)} \\\\
        \\boldsymbol{u}\_i^{(\\sigma)}
    \\end{bmatrix} &\\sim 
    \\text{MVN}(\\boldsymbol{0}, \\boldsymbol{\\Sigma}\_i) \\\\
    \\boldsymbol{\\Sigma\_i} &= \\boldsymbol{D}\_i\\boldsymbol{\\Omega}\\boldsymbol{D}\_i \\\\
    \\text{Diag}(\\boldsymbol{D\_i})' &= \\exp{\\boldsymbol{X}\_i\\boldsymbol{\\eta}}
\\end{align\*}
$$

Altogether, the model is specified as follows:
$$
\\begin{align\*}
    \\boldsymbol{y}\_i &\\sim \\mathcal{N}(\\boldsymbol{\\hat\\mu\_i}, \\boldsymbol{\\hat\\sigma}\_i) \\\\
    \\boldsymbol{\\hat\\mu}\_i &= \\boldsymbol{m}\_i\\boldsymbol{\\beta}\_i \\tag{Mean Model}\\\\ 
    \\boldsymbol{\\beta}\_i' &= \\boldsymbol{M}\_i\\boldsymbol{B} + \\boldsymbol{u}\_i^{(\\mu)} \\\\
    \\log\\boldsymbol{\\hat\\sigma}\_i &= \\boldsymbol{z}\_i\\boldsymbol{\\gamma}\_i \\tag{Residual SD Model}\\\\
    \\boldsymbol{\\gamma}\_i' &= \\boldsymbol{Z}\_i\\boldsymbol{\\Gamma} + \\boldsymbol{u}\_i^{(\\sigma)} \\\\
    \\begin{bmatrix}
        \\boldsymbol{u}\_i^{(\\mu)} \\\\
        \\boldsymbol{u}\_i^{(\\sigma)} 
    \\end{bmatrix} &= \\text{MVN}(\\boldsymbol{0}, \\boldsymbol{\\Sigma}\_i) \\\\
    \\boldsymbol{\\Sigma}\_i &= \\boldsymbol{D}\_i\\boldsymbol{\\Omega}\\boldsymbol{D}\_i \\\\
    \\log\\text{Diag}(\\boldsymbol{D}\_i)' &= \\boldsymbol{X}\_i\\boldsymbol{\\eta} \\tag{Between-group SD Model}
\\end{align\*}
$$

As previously noted, the implemented model employs a maximal structure, such that every group-level variable predicts each within-group coefficient. This allows researchers to explore the impact of group-level variables on the within-group variance model. In practice, this does not affect the estimated ICCs. The ICC calculation only requires teh predicted between-group and within-group variance. If one only had a level 1 within-group variance model, with no level 2 predictors, then the random effect variances would be larger than with level 2 covariates. Consequently, the random effects would change, but the predicted within-group variance would remain approximately the same. Therefore, although the model is maximal, the dense **Γ** matrix provides additional insight into varying ICCs without affecting the estimated ICCs themselves.

Relationship to ICC
-------------------

### Unconditional ICC

### Conditional ICC

Priors
------

<!-- This section details the model as implemented in the \texttt{ICCier} package. -->
<!-- Let $y_{ij}$ be the $j$th observation in the $i$th group response vector $\Vector{y}_i$. -->
<!-- Because we are using an unconditional one-way model, the expected response for any $y_{ij}$ is $\mu_i = \beta_0 + u_{0i}$. -->
<!-- The expected residual standard deviation for $y_{ij}$ is log-linearly modeled from a Level 1 design matrix ($\Vector{x_i}$) and coefficient vector ($\Vector{\gamma}_i$). -->
<!-- All Level 1 scale coefficients are modeled using a multivariate linear model from a Level 2 design matrix ($\Vector{X}$), coefficient \emph{matrix} ($\Vector{\Gamma}$), and random effect row vectors ($\Vector{u}^{\sigma\intercal}_{i}$). -->
<!-- Here, $\Vector{X}_i$ is the row vector containing group-level predictors for group $i$. -->
<!-- All random effects ($u_{0i},\Vector{u}^\sigma_i$) are assumed to vary and covary. -->
<!-- The standard deviations of the random effects are similarly modeled using a multivariate linear model from a Level 2 design matrix, and coefficient matrix $\Vector{\eta}$. -->
<!-- Altogether, the model is defined as follows. -->
<!-- \begin{align*} -->
<!--     \Vector{y}_{i} &\sim \mathcal{N}(\mu_{i}, \Vector{\hat\sigma}_{\epsilon,i}) \\ -->
<!--     \mu_{i} &= \beta_0 + u_{0i} \tag{Location Model} \\ -->
<!--     %\Vector{y}_{i} &= \beta_0 + u_{0i} + \Vector{e}_{i} \tag{Location Model}\\ -->
<!--     \log\Vector{\hat\sigma}_{\epsilon,i} &= \Vector{x}_{i}\Vector{\gamma}_{i} \tag{Level 1 Scale Model}\\ -->
<!--     %\Vector{e}_{i} &\sim \mathcal{N}(0,\Vector{\hat\sigma}_{\epsilon,i}) \\ -->
<!--     %\log\Vector{\hat\sigma}_{\epsilon,i} &= \Vector{x}_{i}\Vector{\gamma}_{i} \tag{Level 1 Scale Model}\\ -->
<!--     \Vector{\gamma}_{i}^\intercal &= \Vector{X}_{i}\Vector{\Gamma} + \Vector{u}^{\sigma\intercal}_i \tag{Level 2 Scale Model}\\ -->
<!--     \begin{bmatrix} -->
<!--         u_{0i} \\ -->
<!--         \Vector{u}^{\sigma}_i -->
<!--     \end{bmatrix} &\sim \mathcal{N}(\Vector{0},\Vector{\Sigma}_i = \Vector{\sigma}_i\Vector{\Omega}\Vector{\sigma}_i) \\ -->
<!--     \diag(\log\Vector{\sigma}_i)^\intercal &= \Vector{X}_i\Vector{\eta} \tag{Between-group Scale Model}  -->
<!-- %    \text{ICC(1)}_{ij} &= \frac{\exp(\Vector{X}_i\Vector{\eta}_{\bullet 1})^2}{\exp(\Vector{X}_i\Vector{\eta}_{\bullet 1})^2 + \exp(\Vector{x}_{ij}\Vector{\gamma}_i)^2} -->
<!-- \end{align*} -->
<!-- The expected unconditional ICC(1)$_{ij}$ is therefore -->
<!-- \begin{align*} -->
<!--     \text{ICC(1)}_{ij} &= \frac{\sigma_{b,i}^2}{\sigma_{b,i}^2 + \hat\sigma^2_{\epsilon,ij}} \\ -->
<!--     &= \frac{\exp(\Vector{X}_i\Vector{\eta}_{\bullet 1})^2}{\exp(\Vector{X}_i\Vector{\eta}_{\bullet 1})^2 + \exp(\Vector{x}_{ij}\Vector{\gamma}_i)^2}, -->
<!-- \end{align*} -->
<!-- where $\Vector{\eta}_{\bullet 1}$ refers to the first column of $\Vector{\eta}$, which specifically contains the predictive coefficients for the group mean standard deviation (i.e., for $\sigma_{b,i}$). -->
<!-- The model implemented by \texttt{ICCier} employs a maximal structure, such that every level 2 variable predicts each level 1 coefficient. -->
<!-- This allows researchers to explore the impact of level 2 variables on the level 1 scale model. -->
<!-- In practice, this does not affect the estimated ICCs. -->
<!-- The ICC calculation only requires the predicted between-group and within-group variance. -->
<!-- If one only had a level 1 scale model, with no level 2 predictors, then the random effect variances would be larger than with level 2 covariates. -->
<!-- Consequently, the random effects would change, but the predicted within-person variance would be the same. -->
<!-- Therefore, although the model is maximal, the dense $\Vector{\Gamma}$ matrix provides additional insight into varying ICCs without affecting the estimated ICCs themselves. -->
<!-- \subsubsection{Priors} -->
Broadly, there are four sets of parameters for which priors must be defined: \* Fixed mean model parameters (**B**) \* Fixed residual variance model parameters (**Γ**) \* Between-group variance model parameters (**η**) \* Random effect correlations (**Ω**)

In practice, one should choose priors based on the a-priori known characteristics of the data, such that priors imply reasonable expected values and variances for the given scenario (Gelman, Simpson, & Betancourt, 2017). Nevertheless, the following priors were specified to accommodate a wide range of possible data.

$$
\\begin{align\*}
  B\_{pq} &\\sim \\mathcal{N}(0,10) \\\\
  \\Gamma\_{pq} &\\sim \\mathcal{ST}(\\nu = 3, \\mu = 0, \\sigma = 5) \\\\
  \\eta\_{pq} &\\sim \\mathcal{ST}(\\nu = 3, \\mu = 0, \\sigma = 5) \\\\
  \\boldsymbol{\\Omega} &\\sim \\text{LKJ}(1)
\\end{align\*}
$$

Note that the individual elements in each coefficient matrix (,,) are assigned the corresponding priors. The two scale models' coefficients are assigned a heavy-tailed Student T prior. The heavy-tailed prior accommodates variances approaching zero better than a normal distribution permits. In our testing of this model, between-group variances in particular may be small or effectively zero, and the wide-tailed priors effectively attenuate any sampling problems. The LKJ prior is a spherical prior over correlation matrices that takes a single parameter. We set the parameter to one, implying a uniform prior over correlation matrices (<span class="citeproc-not-found" data-reference-id="lewandowski2009b">**???**</span>)

<!-- Broadly, there are four sets of parameters for which priors must be defined: Fixed location intercept ($\beta_0$), fixed scale coefficients ($\Vector{\Gamma}$), between-person scale coefficients ($\Vector{\eta}$), and the random effect correlations ($\Vector{\Omega}$). -->
<!-- In practice, one should choose priors based on the a-priori known characteristics of the data, such that the priors imply reasonable expected values and variances for the given scenario \citep{Gelman2017_prior_likelihood}. -->
<!-- Nevertheless, the following priors were specified to accommodate a wide range of possible data. -->
<!-- \begin{align*} -->
<!--     \beta_0 &\sim \mathcal{N}(0,10) \\ -->
<!--     \Gamma_{pq} &\sim \mathcal{ST}(\nu = 3, \mu = 0, \sigma = 5), \forall p,q \\ -->
<!--     \eta_{pq} &\sim \mathcal{ST}(\nu = 3, \mu = 0, \sigma = 5), \forall p,q \\ -->
<!--     \Vector{\Omega} &\sim LKJ(1) -->
<!-- \end{align*} -->
<!-- Note that $\Gamma_{pq}$ and $\eta_{pq}$ represent elements in the respective $\Vector{\Gamma}$ and $\Vector{\eta}$ coefficient matrices. -->
<!-- Each element in these coefficient matrices is given a heavy-tailed Student T prior. -->
<!-- %We recommend using a heavy-tailed prior for scale model parameters. -->
<!-- The heavy-tailed prior accommodates variances approaching zero better than a normal distribution permits. -->
<!-- In our testing of this model, between-group variances in particular may be small or effectively zero, and the wide-tailed priors effectively attenuated any sampling problems. -->
<!-- The LKJ prior is a spherical prior over correlation matrices that takes a single parameter. -->
<!-- We set the parameter to one, implying a uniform prior over correlation matrices [@lewandowski2009b]. -->
References
==========

Brunton-Smith, I., Sturgis, P., & Leckie, G. (2017). Detecting and understanding interviewer effects on survey data by using a cross-classified mixed effects location-scale model. *Journal of the Royal Statistical Society: Series A (Statistics in Society)*, *180*(2), 551–568. <https://doi.org/10.1111/rssa.12205>

Gelman, A., Simpson, D., & Betancourt, M. (2017). The prior can often only be understood in the context of the likelihood. *Entropy*. <https://doi.org/10.3390/e19100555>

Hedeker, D., Mermelstein, R. J., & Demirtas, H. (2008). An application of a mixed-effects location scale model for analysis of ecological momentary assessment (EMA) data. *Biometrics*, *64*(2), 627–634. <https://doi.org/10.1111/j.1541-0420.2007.00924.x>

Hedeker, D., Mermelstein, R. J., & Demirtas, H. (2012). Modeling between-subject and within-subject variances in ecological momentary assessment data using mixed-effects location scale models. *Statistics in Medicine*, *31*(27), 3328–3336. <https://doi.org/10.1002/sim.5338>

Johnson, P. C. (2014). Extension of nakagawa & schielzeth’s r2glmm to random slopes models. *Methods in Ecology and Evolution*, *5*(9), 944–946. <https://doi.org/10.1111/2041-210X.12225>

Li, X., & Hedeker, D. (2012). A three-level mixed-effects location scale model with an application to ecological momentary assessment data. *Statistics in Medicine*, *31*(26), 3192–3210. <https://doi.org/10.1002/sim.5393>

Rast, P., Hofer, S. M., & Sparks, C. (2012). Modeling Individual Differences in Within-Person Variation of Negative and Positive Affect in a Mixed Effects Location Scale Model Using BUGS/JAGS. *Multivariate Behavioral Research*, *47*(2), 177–200. <https://doi.org/10.1080/00273171.2012.658328>

Raudenbush, S. W., & Bryk, A. S. (2002). *Hierarchical linear models: applications and data analysis methods* (2nd ed.). Thousand Oaks.

Shieh, G. (2016). Choosing the best index for the average score intraclass correlation coefficient. *Behavioral Research Methods*, *48*(3), 994–1003. <https://doi.org/10.3758/s13428-015-0623-y>

Williams, D. R., Martin, S. R., & Rast, P. (2019). Putting the Individual into Reliability: Bayesian Testing of Homogeneous Within-Person Variance in Hierarchical Models. *PsyArXiv*. <https://doi.org/10.31234/OSF.IO/HPQ7W>
