# surface_tension_isotherms
fit an isotherm to surface tension vs concentration data

Documentation for Surface Tension Isotherm Fitting V.1 (Jan 2023)
Alison Bain, 2023

function takes in surface tension vs concntration data 
(aranged from lowest concentration to highest concentration) 
and fits the selected isotherm (Gibs or Langmuir) to the low
 concnetration region and a straight line to the high 
concnetration region. 

User selects while points are included in the isotherm region,
the remaining points are assumed to be in the plateau region.

fit_isotherm(st_data,con_data,min_pt, max_pt,n,σ_sol,T,isotherm)

st_data: surface tenion data (N/m)

con_data: cooresponding concentration data (mM = mol/m^3)

min_pt: index of lowest concnetration point to include for isotherm fitting, 0 to include all points

max_pt: index of highest concetration point to include for isotherm fitting must be >2

n: van hofft factor for isotherm, generally 1 for nonionic surfactants and 2 for ionic surfactnats

σ_sol: surface tension of solvent (N/m) = 0.0728 for water

T:Temperature at which the measurements were taken (K)

isotherm: string, either 'gibs' or 'langmuir'

Function outputs two csv files. 
(successive runs in the same directory will save over old files)

...Isotherm_Parameter_Output.csv

fitting parameters for isotherms

CMC (mM)

error on CMC (mM)

Gamma_max (mol/m^2)

error on Gamma_max (mol/m^2)

a (m^3/mol) - Langumuir only 

error on a (m^3/mol) - Langmuir only

...Isotherm_Fit_Output.csv

plotting data for isotherms:

surfactnatn concntration (mM)

surface tension for isotherm region fit (N/m)

surface tension for plateau region fit (N/m)

Function will also plot the measurements and the fits of the two regions
for a visual check. Number of points included in the fit
can be changed as necessary. 

Note, errors are calculated using the covarrience matrix from 
fitting except in the case of the CMC for the Langmuir isotherm 
where the minimum difference between the two fitted curves is taken.
There error is then caclulated as half of the concnetration step. 
