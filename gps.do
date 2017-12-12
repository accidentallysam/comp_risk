
set more off

use "/Volumes/SVBROWN 1/Hazard Model/GPSmatchAfH.dta", clear

quietly bysort CO:  gen dup = cond(_N==1,0,_n)

drop if dup>1
drop dup

keep CO latitude longitude

outsheet using "/Volumes/SVBROWN 1/Hazard Model/MATLAB/data/gps.csv", c nol replace

