Version 1.2.0 (2023-07-21)
-------------------------

	* Include new function elicit_surv. This opens a Shiny application which interfaces with functions from the package and can be used to produce reports. 


   * Fixed issue with error_mod_normal and other functions; argument model was called when it should have been mode
   * For the penalized maximum likelihood model spline models resulted in an error ; DLSsurvspline was accessed in the wrong order with test.deriv
   * Error in the lik_lno function; resulted in incorrect likelihood for the expert opinion in the lognormal model
   


