---
title: "Bayesian Model Selection with SLGF"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{slgf_vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: "REFERENCES.bib"
# csl: journal-of-the-american-statistical-association.csl
---





# 1 Introduction

The package ``slgf`` implements the suspected latent grouping factor (SLGF) methodology of @metzger2019. In this vignette, we first provide a brief theoretical and mathematical review of the method. Next, through the use of real observational data and a designed experiment, we illustrate the usage of ``slgf`` in the context of two common linear model layouts: one-way analysis of variance (ANOVA) and a two-way unreplicated layout. These examples will illustrate the user specification of the models, group-based regression and variance structure, SLGF selection, and regression effect prior specification. 

# 2 Functions

* `ms.slgf()` Computes posterior model probabilities for each user-specified model based on the SLGF grouping methodology 

This package requires a data frame containing a response variable and at least one categorical factor with three or more levels. In many cases this data frame may also contain additional covariates, either continuous or categorical.

* `extract.hats()` Returns the concentrated maximum likelihood estimates for regression coefficients, variance(s), and $g$ from a user-specified model. 

* `groupings()` Computes the unique groupings available for a given $K$. 

The function ``ms.slgf`` requires several inputs to compute and output posterior model probabilities for all models and model classes of interest. The user begins with a data frame containing a continuous response, at least one categorical predictor, and any other covariates of interest. This data frame must not contain column names with the character string ``group``, as this specific string is used by `ms.slgf` to create the group-based regression effect structure. The user must first identify a suspected latent grouping factor, usually by plotting the data and noting a latent structure within the levels of a categorical predictor as illustrated in Section 3.1. The user indicates, via the arguments ``response`` and ``lgf``, character strings corresponding to the response and the suspected latent grouping factor variable names, respectively. 

Next the user determines the model classes they wish to evaluate. We note the distinction between these model classes and the *R* class of a variable. The argument ``usermodels`` is a list where each element contains a string of *R* class ``formula`` or ``character``. The user also specifies which classes should also be considered in a heteroscedastic context via the argument ``het``, which provides an indicator ``1`` or ``0`` corresponding to each model class specified in ``usermodels``. Together the arguments ``usermodels`` and ``het`` create the full set of model classes considered. 

Next the user chooses a prior to place on the regression effects. As described in Section 3.2, ``prior="flat"`` (the default) implements a constant prior and ``prior="zs"`` imposes the Zellner-Siow mixture $g$-prior. 

Finally the user must specify the minimum number of levels of the slgf that can comprise a group, via the argument ``min.levels``, which defaults to ``1``. Because the number of possible grouping schemes increases exponentially with $K$, the number of levels of the SLGF, the user can reduce the number of candidate models by increasing ``min.levels``, limiting the computational burden. Additionally, when considering data with limited degrees of freedom, increasing ``min.levels`` can also ensure estimability of the parametrized ``usermodels``; see Section 4.2 for more detail. 

# 3 Analytic Details

## 3.1 Model Specification

For a thorough review of the model specification see @metzger2019. We must define a model flexible enough to account for all possible model classes. We begin with the model

\equation
$\boldsymbol{Y}=X\boldsymbol{\beta}+\boldsymbol{\varepsilon}$,
\equation

where 

* $\boldsymbol{Y}$ is an $N\times 1$ vector of continuous observations;
* $X$ is an $N\times P$ design matrix;
* $\boldsymbol{\beta}$ is a $P\times 1$ vector of regression effects; and
* $\boldsymbol{\varepsilon}$ is an $N\times 1$ residual vector where $\boldsymbol{\varepsilon}\sim N(\boldsymbol{0},\, \Sigma)$.

We must augment our notation to account for several additional components: the full slgf with $K$ degrees of freedom or a 2-degree of freedom group effect; interactions with the slgf or group effect; other effects of interest unrelated to the slgf; and potential group-based heteroscedasticity.Thus we let $X=(\boldsymbol{1}^T|W|V|U)$ and $\boldsymbol{\beta}=(\alpha,\boldsymbol{\nu},\boldsymbol{\tau},\boldsymbol{\rho})$ to obtain

\begin{equation}
\boldsymbol{Y} = \boldsymbol{1}^T \alpha + W\boldsymbol{\nu} + V\boldsymbol{\tau} + U\boldsymbol{\rho} + \boldsymbol{\varepsilon}
\end{equation}

where 

* $\boldsymbol{1}^T$ is an $N\times 1$ vector of 1s;
* $\alpha$ is a scalar intercept common to all models; 
* $\boldsymbol{\nu}$ is the $K$-dimensional vector of slgf effects;
* $W$ is the corresponding $K\times N$ design matrix;
* $\boldsymbol{\tau}$ is the $J$-dimensional vector of additional effects; 
* $V$ is the corresponding $J\times N$ design matrix;
* $\boldsymbol{\rho}$ is the $L$-dimensional vector of slgf-interaction effects; and
* $U$ is the corresponding $L\times N$ design matrix.

Not all linear model classes will incorporate each term; for example, the one-way ANOVA example of Section 4.1 contains only a single categorical predictor so $\bf{\tau}:=\boldsymbol{0}$ and $\bf{\rho}:=\boldsymbol{0}$. When the slgf is modeled as a group effect, denote this effect as $\tilde{\boldsymbol{\nu}}$; when there is an interaction involving the group effect, denote it as $\tilde{\boldsymbol{\rho}}$.  

In heteroscedastic contexts, let $\sigma^2_1$ and $\sigma^2_2$ represent the error variances of groups 1 and 2, respectively. Let $\tilde{\Sigma}$ denote the covariance matrix where the $i$th diagonal element is $\sigma^2_1$ if $y_i$ belongs to group 1, or $\sigma^2_2$ if $y_i$ belongs to group 2. 

## 3.2 Parameter Priors

In the interest of objectivity, we implement noninformative priors on the regression effects and error variance(s). For homoscedastic models,

\equation
$P(\boldsymbol{\beta},\sigma^2|m)\propto \frac{1}{\sigma^2}$
\equation

and for a heteroscedastic models, 
\equation
$P(\boldsymbol{\beta},\sigma^2_1, \sigma^2_2|m)\propto \frac{1}{\sigma^2_1\cdot\sigma^2_2}$
\equation

The user implements these priors with the argument `prior="flat"`. 

However, in contexts with limited data, such as the two-way unreplicated layout of Section 4.2, we recommend the Zellner-Siow mixture $g$-prior [@Liangetal], which reduces the minimal training sample size necessary for the computation of the fractional Bayes factor. For homoscedastic models,

\equation
$P(\alpha, \sigma^2|m)\propto \frac{1}{\sigma^2}$ and $\boldsymbol{\beta}_{-\alpha}|\Sigma,g,m \sim N(\boldsymbol{0},\,g(X^T\Sigma X)^{-1})$;
\equation

for heteroscedastic models,

\equation
$P(\alpha, \sigma^2_1, \sigma^2_2|m)\propto \frac{1}{\sigma^2_1\cdot\sigma^2_2}$ and $\boldsymbol{\beta}_{-\alpha}|\tilde{\Sigma},g,m \sim N(\boldsymbol{0},\,g(X^T\tilde{\Sigma} X)^{-1})$;
\equation

and in both cases, 

\equation
$g\sim \texttt{IG}\big(\frac{1}{2},\, \frac{N}{2}\big)$.
\equation

The user implements these priors with the argument `prior="zs"`. 

## 3.3 Fractional Bayes Factors and Posterior Model Probabilities

We invoke a fractional Bayes factor approach to compute well-defined posterior model probabilities for each model; for a thorough review and justification see @OHaganfbfs and @metzger2019. 

Let $\mathcal{M}$ represent the full set of models under consideration, representing all classes and grouping schemes of interest. Denote $\boldsymbol{\theta}$ as the full set of unknown parameters associated with a model $m\in \mathcal{M}$ and $\pi(\boldsymbol{\theta})$ as the joint prior on these parameters. A fractional Bayes factor is a ratio of two fractional marginal model probabilites, where a fractional marginal likelihood is defined as

\equation
$q^b(\boldsymbol{Y}|\boldsymbol{\theta})=\frac{\int P(\boldsymbol{Y}|\boldsymbol{\theta},m)\pi(\boldsymbol{\theta})d\boldsymbol{\theta}}{\int P(\boldsymbol{Y}|\boldsymbol{\theta},m)^b\pi(\boldsymbol{\theta})d\boldsymbol{\theta}}$
\equation

for a some fractional exponent $b$. We exclusively implement $b=\frac{m_0}{N}$ where $m_0$ is the minimal training sample size required for $P(\boldsymbol{Y}|m)$ to be a proper distribution. 

Thus `slgf` must approximate the integrals $\int P(\boldsymbol{Y}|\boldsymbol{\theta},m)\pi(\boldsymbol{\theta})d\boldsymbol{\theta}$ and $\int P(\boldsymbol{Y}|\boldsymbol{\theta},m)^b\pi(\boldsymbol{\theta})d\boldsymbol{\theta}$ for all $m\in \mathcal{M}$. In the case of noninformative regression priors, $\boldsymbol{\beta}$ is integrated analytically, and $\sigma^2$ or $\sigma^2_1, \sigma^2_2$ are integrated using a Laplace approximation after a log-variance transformation. In the Zellner-Siow mixture $g$-prior case, $\alpha$ and $\boldsymbol{\beta}_{-\alpha}$ are integrated analytically, and $\sigma^2$ or $\sigma^2_1, \sigma^2_2$ and $g$ are again integrated using a Laplace approximation with a log-variance transformation. Let $\boldsymbol{\tilde{\theta}}$ represent the set of unknown parameters after the regression effects $\boldsymbol{\beta}$ have been integrated out. Then for dimensionality $d=2$ in the noninformative prior case and $d=3$ in the Zellner-Siow mixture $g$-prior case,

\equation
$\log\big(\int P(\boldsymbol{Y}|\boldsymbol{\tilde{\theta}},m)\pi(\boldsymbol{\tilde{\theta}})d\boldsymbol{\tilde{\theta}}\big)\approx \frac{d}{2}\log(2\pi)-\frac{1}{2}\log|-H^{\star}|+\log(P(\boldsymbol{Y}|\boldsymbol{\tilde{\theta}}^{\star},m))$
\equation

where $\boldsymbol{\tilde{\theta}}^{\star}$ is the mode of $P(\boldsymbol{Y}|\boldsymbol{\tilde{\theta}},m)\pi(\boldsymbol{\tilde{\theta}})$ and  $H^{\star}$ is the Hessian matrix evaluated at $\boldsymbol{\tilde{\theta}}^{\star}$. These values are computed with the functions ``optim`` and ``numDeriv::hessian``, respectively. We make a similar computation for $\int P(\boldsymbol{Y}|\boldsymbol{\tilde{\theta}},m)^b\pi(\boldsymbol{\tilde{\theta}})d\boldsymbol{\tilde{\theta}}$ to compute the fractional marginal model probability $q^b(\boldsymbol{Y}|\boldsymbol{\theta})$ for all $m\in \mathcal{M}$, well defined for both homoscedastic and heteroscedastic models. Once log-fractional marginal likelihoods have been computed for all models, we subtract the maximum from this set so that the set of log-fractional marginal likelihoods has been rescaled to have a maximum of 0. Each value is exponentiated to obtain a set of fractional marginal likelihoods with maximum 1 to avoid numerical underflow when computing posterior model probabilities. 

Model priors are imposed uniformly prior by model class, and for classes containing multiple models, the prior on each class is uniformly divided among the models it contains.  

We finally compute posterior model probabilities for each model:

\equation
$P(m^{\prime}|\boldsymbol{Y})=\frac{P(\boldsymbol{Y}|m^{\prime})P(m^{\prime})}{\underset{\mathcal{M}}{\sum}P(\boldsymbol{Y}|m)P(m)}$
\equation


# 4 Examples

## 4.1 One-way Analysis of Variance (ANOVA)

First consider the data of @ObrienHeft, who studied olfactory function by age ($y$-axis), where age was divided into five categories ($x$-axis). Load the dataset:


```r
data(smell)
```

A simple boxplot suggests potential heteroscedasticity, with latent grouping structure 1,2,3:4,5. We propose the alternative notation {1,2,3}{4,5} to avoid ambiguity with the "`:`" symbol used with strings of class `formula` in R.  


```r
boxplot(smell$olf ~ smell$agecat, pch = 16, xlab = "Age Category", 
    ylab = "Olfactory Score", xaxt = "n", col = c("gray30", "gray30", 
        "gray30", "white", "white"), border = c("black", "black", 
        "black", "red", "red"), main = "O'Brien and Heft (1995) Smell Data")
axis(1, labels = c("1\n(n=38)", "2\n(n=35)", "3\n(n=21)", "4\n(n=43)", 
    "5\n(n=42)"), at = 1:5, lwd.tick = 0)
```

<img src="figure/unnamed-chunk-3-1.png" title="O'Brien and Heft (1995) studied olfactory function by age ($y$-axis), where age was divided into five categories ($x$-axis). The data suggests potential heteroscedasticity, with latent grouping structure {1,2,3}{4,5}." alt="O'Brien and Heft (1995) studied olfactory function by age ($y$-axis), where age was divided into five categories ($x$-axis). The data suggests potential heteroscedasticity, with latent grouping structure {1,2,3}{4,5}." style="display: block; margin: auto;" />

We thus consider ``agecat`` as the suspected latent grouping factor. The apparent latent grouping scheme we observe is denoted {1,2,3}{4,5} (or equivalently, {4,5}{1,2,3}), but all possible grouping schemes are considered. The means may also differ by level, but this is more difficult to distinguish by the plot. Thus we first consider the following model and heteroscedasticity structures:


```r
smell.models <- list("olf~1", "olf~agecat", "olf~group")
smell.het <- c(0, 0, 1)
```

By specifying `smell.models` and `smell.het` as done above, we consider the third model specified with group-based heteroscedasticity. This elicits four model classes: a homoscedastic global mean, homoscedastic age level means, homoscedastic group-based means, and group-based means with group-based heteroscedasticity. Finally we note that with a relatively large sample size, we prefer the use of noninformative priors via ``prior="flat"``, and we specify ``min.levels=1``, as we have no prior information on the number of levels of ``agecat`` that may be grouped together and we wish to consider a comprehensive set of candidate models. With the model specification list argument `usermodels` and the heteroscedasticity vector argument `het`, the number of unique classes can be obtained as ``length(usermodels)+sum(het)``. 


```r
smell.out <- ms.slgf(dataf=smell, response="olf", lgf="agecat", usermodels=smell.models, 
                     prior="flat", het=smell.het, min.levels=1)
```

Note by specifying the argument `min.levels=1`, we consider all possible grouping schemes containing at least one level of the SLGF. We could specify `min.levels=2` to consider a smaller model space with at least two levels per group, but `min.levels=3` or greater is not possible with five levels of the SLGF. The output `smell.out` is a list of class ``slgf``, with six elements when ``prior="flat"`` and seven when ``prior="zs"``; see the help file for full details. 

We summarize the five most probable models of 32 considered:


```r
smell.out$results[1:5, c(1:3,6,11)]
#>        Model       Scheme Variance   Fmodprob               Class
#> 1  olf~group {4,5}{1,2,3} Heterosk 0.99987407 olf~group, Heterosk
#> 2  olf~group {1,2}{3,4,5} Heterosk 0.00012589 olf~group, Heterosk
#> 3 olf~agecat         None   Homosk 0.00000003  olf~agecat, Homosk
#> 4  olf~group {5}{1,2,3,4} Heterosk 0.00000000 olf~group, Heterosk
#> 5  olf~group {4,5}{1,2,3}   Homosk 0.00000000   olf~group, Homosk
```

We strongly favor the model with group-based mean effects and variances via scheme {4,5}{1,2,3}. 

The function ``extract.hats`` provides the estimates for the coefficient(s) and variance(s) for a given ``model.index``, as well as $g$ if ``prior="zs"``. 


```r
smell.hats <- extract.hats(smell.out, model.index=1)
print(smell.hats$`sigsq.{4,5}`)
#> [1] 0.01211023
print(smell.hats$`sigsq.{1,2,3}`)
#> [1] 0.05868753
```

That is, we compute $\hat{\sigma}^2_{\text{1,2,3}}=$ 0.05869 and $\hat{\sigma}^2_{\text{4,5}}=$ 0.01211. Let us also consider the case where ``het=c(1,1,1)``; that is, we include two additional classes: group-based variances and a single global mean, and group-based variances with means by ``agecat``. 


```r
smell.out <- ms.slgf(dataf=smell, response="olf", lgf="agecat", 
                     usermodels=smell.models, prior="flat", 
                     het=c(1,1,1), min.levels=1)

smell.out$results[1:5, c(1:3,6,11)]
#>        Model       Scheme Variance   Fmodprob                Class
#> 1  olf~group {4,5}{1,2,3} Heterosk 0.54204570  olf~group, Heterosk
#> 2 olf~agecat {4,5}{1,2,3} Heterosk 0.44946817 olf~agecat, Heterosk
#> 3 olf~agecat {1,2}{3,4,5} Heterosk 0.00840126 olf~agecat, Heterosk
#> 4  olf~group {1,2}{3,4,5} Heterosk 0.00007605  olf~group, Heterosk
#> 5 olf~agecat {1,3}{2,4,5} Heterosk 0.00000566 olf~agecat, Heterosk
```

Now the most probable models are a bit less conclusive, as the distinct category-means model with scheme {4,5}{1,2,3} group-based heteroscedasticity accounts for a meaningful amount of posterior model probability. We can easily summarize the scheme and class probabilities, which strongly favor scheme {4,5}{1,2,3} and moderately favor the group-based means and variances model class:


```r
smell.out$scheme.Probs
#>              Scheme.Prob
#> {4,5}{1,2,3}  0.99151391
#> {1,2}{3,4,5}  0.00847731
#> {1,3}{2,4,5}  0.00000566
#> {2,3}{1,4,5}  0.00000219
#> {1}{2,3,4,5}  0.00000047
#> {5}{1,2,3,4}  0.00000025
#> {2}{1,3,4,5}  0.00000018
#> None          0.00000002
#> {1,4}{2,3,5}  0.00000000
#> {1,5}{2,3,4}  0.00000000
#> {2,4}{1,3,5}  0.00000000
#> {2,5}{1,3,4}  0.00000000
#> {3,4}{1,2,5}  0.00000000
#> {3,5}{1,2,4}  0.00000000
#> {3}{1,2,4,5}  0.00000000
#> {4}{1,2,3,5}  0.00000000
smell.out$class.Probs
#>                      Class.Prob
#> olf~group, Heterosk  0.54212175
#> olf~agecat, Heterosk 0.45787818
#> olf~1, Heterosk      0.00000004
#> olf~agecat, Homosk   0.00000002
#> olf~1, Homosk        0.00000000
#> olf~group, Homosk    0.00000000
```

## 4.3 Unreplicated Two-way Layouts

Next we analyze a two-way unreplicated layout. Consider the data analyzed by @franckdogs, where six dogs with lymphoma were studied. Two individual samples were taken from healthy and tumor tissue within each dog, and the copy number variation was measured for each sample. We arrange dogs into rows and tissue types into columns of a two-way layout. We first plot the data to determine whether a latent grouping structure underlies the data:

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

We strongly suspect that dogs 1, 2, and 5 behave distinctly from dogs 3, 4, and 6. The telltale non-parallel lines suggest an underlying interaction, but with only a single observation in each cell, we lack the degrees of freedom to fit a standard row by column interaction term. Instead, we let ``dog`` represent the SLGF, exclude the column effect, and parametrize a group-by-column interaction term, isomorphic to fitting distinct column effects by group. Thus we consider four reasonable model classes: a dog and tissue effect, a dog and tissue-by-column interaction, and the heteroscedastic counterparts of each class. We accomplish this with the arguments ``usermodels=list("gene~dog+tissue", "gene~dog+group:tissue")`` and ``het=c(1,1)``. Additionally, we must specify ``min.levels=2`` or ``min.levels=3`` to ensure sufficient degrees of freedom to estimate the ``group:col`` effect. Here we consider ``min.levels=2`` in the interest of assessing a more complete set of candidate models. 

Because of the limited amount of data, our choice of prior is more impactful in this case. We impose ``prior="zs"`` to utilize the Zellner-Siow mixture $g$-prior, as the fractional Bayes factor exponent would require a prohibitively high proportion of the data for model training.  

We first put the two-way layout into a ``data.frame`` format compatible with the ``slgf`` function, and then implement this approach. Note that warnings are related to the optimization initialization, and are suppressed here. 


```r
lymphoma.tall <- maketall(lymphoma)
lymphoma.tall <- data.frame("gene"=lymphoma.tall[,1], "dog"=lymphoma.tall[,2], 
                            "tissue"=lymphoma.tall[,3])

lymphoma.models <- list("gene~dog+tissue", "gene~dog+group:tissue")
lymphoma.out <- ms.slgf(dataf=lymphoma.tall, response="gene", lgf="dog", 
                        usermodels=lymphoma.models, prior="zs", 
                        het=c(1,1), min.levels=2)
```

Note the ``:`` operator in the ``usermodels`` syntax, which does not automatically include the main effects ``group`` and ``tissue`` which are not both estimable. As expected, we conclude with high probability that scheme {1,2,5}{3,4,6} underlies the data. The five most probable models are given by:


```r
lymphoma.out$results[1:5, c(1:3,6)]
#>                   Model         Scheme Variance   Fmodprob
#> 1 gene~dog+group:tissue {1,2,5}{3,4,6}   Homosk 0.73643491
#> 2 gene~dog+group:tissue {1,2,5}{3,4,6} Heterosk 0.24833051
#> 3       gene~dog+tissue           None   Homosk 0.00311082
#> 4 gene~dog+group:tissue {1,2}{3,4,5,6} Heterosk 0.00211760
#> 5       gene~dog+tissue {1,2,5}{3,4,6} Heterosk 0.00109052
```

while the class and five highest scheme probailities are


```r
lymphoma.out$class.Probs
#>                                 Class.Prob
#> gene~dog+group:tissue, Homosk   0.74081898
#> gene~dog+group:tissue, Heterosk 0.25326343
#> gene~dog+tissue, Homosk         0.00311082
#> gene~dog+tissue, Heterosk       0.00280677
head(lymphoma.out$scheme.Probs, 5)
#>                Scheme.Prob
#> {1,2,5}{3,4,6}  0.98585594
#> None            0.00311082
#> {1,2}{3,4,5,6}  0.00288222
#> {1,5}{2,3,4,6}  0.00145899
#> {3,6}{1,2,4,5}  0.00081783
```

The ``group.datafs`` element of ``lymphoma.out`` contains ``data.frames`` associated with each model and grouping scheme. We first determine which element of ``lymphoma.out$group.datafs`` contains the ``data.frame`` of interest via the column ``dataf.Index``:


```r
lymphoma.out$results[1,c(1:3,6,8)]
#>                   Model         Scheme Variance  Fmodprob dataf.Index
#> 1 gene~dog+group:tissue {1,2,5}{3,4,6}   Homosk 0.7364349          18
```

This tells us that element 18 of ``lymphoma.out$group.datafs`` contains the ``data.frame`` with the {1,2,5}{3,4,6} group effect:


```r
lymphoma.out$group.datafs[[18]]
#>      gene dog tissue group
#> 1  9.3278   1      1 1,2,5
#> 2  9.2168   1      2 1,2,5
#> 3  9.5108   2      1 1,2,5
#> 4  9.3942   2      2 1,2,5
#> 5  8.7535   3      1 3,4,6
#> 6  9.4158   3      2 3,4,6
#> 7  8.6372   4      1 3,4,6
#> 8  9.2480   4      2 3,4,6
#> 9  9.4981   5      1 1,2,5
#> 10 9.4626   5      2 1,2,5
#> 11 8.7322   6      1 3,4,6
#> 12 9.3439   6      2 3,4,6
```

# References
