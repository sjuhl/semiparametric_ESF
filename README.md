# When Space is Just a Nuisance: A Simple Solution for Spatial Dependence

## Summary
Spatial autocorrelation is a ubiquitous phenomenon in cross-sectional data that poses notable challenges for statistical inference including biased and inconsistent parameter estimates. While political scientists so far predominantly rely on parametric spatial regression models (such as a spatial lag or spatial error model), semiparametric spatial filtering techniques constitute a valuable alternative. By using the eigenfunction decomposition of a transformed connectivity matrix, the filtering approach generates a synthetic proxy variable from a linear combination of judiciously selected eigenvectors. This synthetic variable acts as a surrogate for omitted spatial effects and removes spatial autocorrelation from the model residuals. This study introduces eigenvector-based spatial filtering to political science and discusses its strengths and limitations in comparison to parametric spatial models. Analytical results and Monte Carlo simulations show that spatial filtering outperforms non-spatial OLS and incorrectly specified spatial models, and performs comparably under common conditions to spatial models that correctly reflect the data-generating process. These techniques can easily be implemented with an easy-to-use R package. This study concludes that spatial filtering constitutes a very flexible and intuitive alternative to parametric spatial models for avoiding damaging misspecification problems.  Indeed, since we rarely know the true spatial process in applied research, spatial filtering may be the best option when space is just a nuisance.

## Contact
- [Sebastian Juhl](http://www.sebastianjuhl.com) (SAP SE & University of Missouri)
- [Laron K. Williams](https://williamslaro.github.io/) (University of Missouri)
