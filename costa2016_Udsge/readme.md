# Archives of codes from Understanding DSGE (Costa, 2016)
Code description:
- `rbc_log.mod`: log-linearized basic RBC
- `nk_basic.mod`: log-linearized basic NK
- `nk_wage_sticky.mod`: log-linearized NK with wage stickiness (nominal rigidity)
- `nk_habit.mod`: log-linearized NK with habit formation and non-Ricardian agents (household rigidity)
- `nk_adjcost.mod`: log-linearized NK with capital adjustment and under-utilization costs (production rigidity)
- `nk_gov.mod`: log-linearized NK with fiscal policy and monetary policy
- `nk_gov_fiscal_adj_1.mod` to `nk_gov_fiscal_adj_4.mod`: NK with fiscal adjustments. 1 is no adjustment, 2 to 4 are small, medium, and big adjustments. Parameters governming these are $\phi_{\tau^c}, \phi_{\tau^l}, \phi_{\tau^k}$. After running these mod files, we can use the `rep_fig76_taxpol.m` to show the policy differences.
