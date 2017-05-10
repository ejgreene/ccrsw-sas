CCR_SW for SAS<sup>&reg;</sup>
===

This SAS macro for covariate-constrained randomization of stepped-wedge designs is released under the GNU Lesser General Public License.  It hasn't been peer-reviewed yet, but it's offered for use as-is with no implied warranty (not even for merchantability or fitness for a particular purpose).

The macro is very similar to [CCR](https://github.com/ejgreene/ccr-sas) but assigns to steps rather than arms.  The randomization aims to balance the covariates across both arms, weighted by the number of periods each cluster spends in each arm, within strata and overall.

This usage information is lifted from the macro comments; a more user-friendly article is in the works.

```SAS
/****************************/
/* ARGUMENTS AND PARAMETERS */
/****************************/

/* REQUIRED ARGUMENTS:
   dataset = the input dataset, which should have one observation per cluster
   stratid = the variable in the dataset identifying which stratum a cluster
             is in; can leave blank for an unstratified design
   clustid = the variable in the dataset identifying each cluster
   x       = list of covariate variables
   vtype   = 'd' -> dummy covariate
                    (display checks as frequency tables of steps & differences)
             'c' -> continuous/many-valued covariate
                    (display checks as mean/min/q1/median/q3/max of difference)
   slimit  = stratum balancing criteria; set to any for unstratified design
   olimit  = overall balancing criteria
             formatting for criteria: {type}{difference} -or- ANY
             ANY imposes no constraint (allows for constraining only at stratum
               or overall level)
             {type} = 'm' for mean or 's' for sum
             if {difference} starts with f, it's treated as an allowable
               fraction of the control value
             FOR EXAMPLE:
             s100 -> sums must differ by no more than 100
             mf.5 -> means must differ by no more than half the control mean
             NOTE: values in vtype, slimit, and olimit should be separated by
                   whitespace (CCR_V1.0 used asterisks)

   OPTIONAL PARAMETERS:
   steps    = the number of steps to allocate to [default: 3]
   stepsize = the variables in the dataset identifying how many clusters in a
                stratum to allocate to each step
              if not set, clusters are split as evenly as possible between steps
                (at both levels)
   ssample  = number of possible stratum allocations to sample and check in
	            each stratum [default: 100,000]
   osample  = number of combinations of acceptable stratum allocations to
                sample and check [default: 100,000]
              if not set, all combinations will be generated and checked
   seed     = seed for sampling allocations and choosing final randomization
              if not set, SAS will seed from the system clock
                (rand function and proc surveyselect default behavior)
   select   = # of final randomizations to display [default: 1]
   samestephi = percentage at which clusters appear in the same step "a lot"
                 [default: 75]
   samesteplo = percentage at which clusters appear in the same step "rarely"
                 [default: undefined -> 100/(2*steps), which is half chance]
   adjstephi = percentage at which clusters appear in adjacent steps "a lot"
                 [default: 75]
   adjsteplo = percentage at which clusters appear in adjacent steps "rarely"
                 [default: undefined -> 100*(steps-1)/(steps)^2,
                                                       which is half chance]
   verbose controls how much is sent to output and log
           (though everything is available afterwards in temporary datasets
            regardless of setting)
      0/false -> only print summaries and checks
      1/true -> print lists of allocations that satisfy the constraints,
                plus the final randomization
      2 -> print lists of all possible within-stratum allocations
      3 -> print lists of all possible (or all sampled) overall allocations
   setting binomsig prints all pairs of clusters whose rates of appearance in
     the same step differ significantly from chance (parameter sets alpha)
   setting debug prints a lot of macro variables to the log
   setting sizecheck to 0 disables checking whether a sample size
     is larger than the number of overall allocations

   Probability of two clusters being in same step: 1/s
   Probability of two clusters being x steps apart: 2*(s-x)/s^2
*/
```
