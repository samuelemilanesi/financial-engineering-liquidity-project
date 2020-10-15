ReadMe Final project Group 6
Milano 14.05.19  
Financial Engineering, Politecnico
---------------------------------------------------------------------------

Script: runFinalProject_Group6_1_2
- readExcelData: reads data from excel
- clean_date: Take dates and check if they are business dates: if so, they do not
              change, otherwise they become business dates according to Modified
              Following day convention

---------------------------------------------------------------------------

Functions for point i):
- bootstrap_OIS: Perform bootstrap for EONIA curve; 
- bootstrap_pseudo: Perform bootstrap for Euribor-6m curve.

---------------------------------------------------------------------------

Functions for point ii):
- tenordates: Find relevant subsequent dates from a starting date according to a
              regular timestep between dates and adjusted according to a specified
              business day convention and a calendar of holidays.
- interp_discounts: Auxiliary function that computes the vector of B required at dates

- calibrate_model: Calibrates volatility parameters
- normal_vols_swaption_MHJM: Computes market price for a complete set of ATM Cash-Settlement 
                             diagonal reciver swaptions vs Euribor6m (fixed leg paid annually, 
                             floating semiannually), given a corresponding set of normal volatilities. 
                             The function uses the normal-Black formula (a.k.a. Bachelier formula) to 
                             convert normal volatilities into prices and works under the hypothesis of 
                             multicurve HJM framework.
- compute_model_swaption_MHJM:  Computes MHJM model price for a complete set of ATM Cash-Settlement 
                                diagonal reciver swaptions vs Euribor6m (fixed leg paid annually, 
                                floating semiannually).
---------------------------------------------------------------------------
Script: runFinalProject_Group6_3_4

---------------------------------------------------------------------------
Functions for point iii):
- readExcelData: reads data from excel
- buildStruct: Create the struct containing data of BNPP and Santander on 
               10 September 2015, according to paper [2]
- dirty_from_clean: Compute dirty prices from clean prices
- liquidity_spread_bounds: Compute upper and lower bounds for the 
                           illiquidity price
- bootstrap_Z_curve: Perform the bootstrap on Z-spread: Zeta-spread curve 
                     is assumed to be constant up to the maturity of the 
                     bond with the lowest maturity and it is linearly 
                     interpolated afterwards
- build_system_of_equation:  Auxiliary function that builds the system of 
                             equations that are solved in the 
                             bootstrap_Z_curve. The system is meant to be a
                             non-linear system related to a unique bond,
                             the i-th bond.
- do_plot_error: Plot the difference between the upper and lower bounds
                      for the illiquidity price Delta_tau for obligor bonds
---------------------------------------------------------------------------

Functions for point iv):
- bond_yield:    Compute bond yields for the liquid bond
- do_plot_yield: Plot obligor yields and yields obtained for illiquid bonds 
                 with ttl equal to two weeks and two months
---------------------------------------------------------------------------
