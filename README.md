# Odds-Ratio-AB-Analysis
Reliability of Odds Ratio Confidence Intervals: Simulation Study

Background: 
An odds ratio (OR) is a statistical measure that measures the association between an 
exposure and an outcome. Typically, it's used to compare two groups, one with some exposure 
and one without. It then measures the likeliness of an outcome between the different exposure 
groups. It is particularly useful in case-control studies, cohort studies, and clinical trials.

Goal: 
* Examine the reliability of three commonly used confidence intervals for 
  ORs – Woolf’s, Gart’s, and Agresti’s – by assessing their empirical coverage rates through statistical analysis. 
* Evaluate how well these confidence intervals maintain their nominal coverage levels across 
  different combinations of large, medium, and small sample sizes, varying ORs, and varying probabilities.

Tools/Methodology:
* Monte Carlo Simulation with 2x2 contingency tables simulated under predefined conditions
* R (programming software)
* ggplot2 (data visualization package in R)

Key Findings:
* Convergence of all tests towards a flat empirical coverage rate as p1 (probability of success in group 1) 
  values approach extremes. The tests converged as p1 values approached 1 when n1 (sample size of group 1)
  was significantly larger than n2 (sample size of group 2), and converged when p1 values approached 0 when n2 was
  significantly larger than n1.
    <img width="542" height="221" alt="image" src="https://github.com/user-attachments/assets/ccdd96c3-fbf5-4a36-b165-955220130124" />
* Downward trends showing high to low empirical coverage rates as p1 values increase. Concave and convex parabolic
  (symmetrical) results, which were even on both tails but with empirical coverage rates straying
  lower/higher, respectively, were common for odds ratios closer to 1 and similar sample sizes.
    <img width="238" height="221" alt="image" src="https://github.com/user-attachments/assets/bf3c6282-7767-4b43-bd46-366ae3ea9cce" />
* Extremely high empirical coverage rates, which then dropped drastically, when OR was high (over 10), n1 was small (3, 10, 20),
  and n2 was moderately sized (30, 35, 75)
    <img width="238" height="221" alt="image" src="https://github.com/user-attachments/assets/80dbe69a-44d1-42bc-94a5-9932f2ef0ffa" />
