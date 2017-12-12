set scheme s1mono
graph set eps logo off
graph set eps preview off
graph set eps orientation portrait
graph set eps mag 100

set matsize 1000

set more off
pause on

local symptom = 99 // illness
// local symptom = 7 // diarrhea

local datadir = "/Volumes/SVBROWN 1/Mali-AfH/AfH-data-Stata/Stata"
// local datadir = "E:/Mali-AfH/AfH-data-Stata/Stata"
// local datadir = "F:/Mali-AfH/AfH-data-Stata/Stata"

cd "/Volumes/SVBROWN 1/Hazard Model"
// cd "E:/Hazard Model"
// cd "F:/Hazard Model"

/////////
// load in the spells data
use Spells-imp.dta, clear

local symptom_label : label symptom_label `symptom'
local symptom_label = proper("`symptom_label'")
local fname = strtoname("`symptom_label'")
local fname = substr("`fname'",1,4)

cap log close Hazard_`fname'
log using Hazard_`fname', replace text name(Hazard_`fname')

/////////
// load in covariates
// child age and sex
merge m:1 ID_conc ID_HH ID ID_enfant using "`datadir'/1-Kids-AfH-only.dta", keepus(child_age_Aug20 sex) force
drop if _merge==2
drop _merge
rename child_age_Aug20 age
replace age = 0 if age < 0
rename sex male
replace male = 2-male // recode to 1 = male
// wealth quartile
merge m:1 ID_conc ID_HH using "`datadir'/Assets.dta", keepus(value_rank) force
drop if _merge==2
drop _merge
gen quart = 2
replace quart =1 if value_rank<1077/4
replace quart =4 if value_rank>=(1077/4)*3
replace quart =3 if quart !=4 &   value_rank>= 1077/2
drop value_rank
// region (and true IDs)
preserve
use "`datadir'/FullListAfH.dta", clear
keep ID* True* strat_region
collapse (max) True* strat_region, by(ID*)
save temp-1.dta, replace
restore
merge m:1 ID_conc ID_HH ID ID_enfant using "temp-1.dta", force
drop if _merge==2
drop _merge
rm temp-1.dta
rename strat_region region
replace ID_conc = True_ID_conc if True_ID_conc ~= .
replace ID_HH = True_HH if True_HH ~= .
replace ID = True_mother if True_mother ~= .
drop True*

/////////
// drop unreasonable observations
keep if age <= 5

/////////
// format the data for the hazard model
// keep only observations of the symptom of interest
keep if symptom == `symptom'
// reshape so we can grab the first visit
reshape long visit vstart vlen vtype, ///
	i(intid inter ID_conc ID_HH ID ID_enfant symptom spstart) j(j)
// drop if vtype == 3 // drop informal visits
drop if visit == 9 | visit == 10 // drop "other" and "don't know"
drop j
// sort by visit date (and then type, so "better" types dominate)
sort intid inter ID_conc ID_HH ID ID_enfant symptom spstart vstart vtype
// generate a new "j" variable that corresponds to the new sort order
qui by intid inter ID_conc ID_HH ID ID_enfant symptom spstart:  gen dup = _n
// reshape back to the original format
reshape wide visit vstart vlen vtype, ///
	i(intid inter ID_conc ID_HH ID ID_enfant symptom spstart) j(dup)

// generate exit states
gen tvisit = vstart1-(spstart-1)
replace vtype1 = . if min(splen,tvisit) > 7
egen time = rowmin(splen tvisit)
replace time = 7 if time>7
gen visit = (vtype1<3)*1+(vtype1==3)*2+(vtype1 == .)*3

order intid interviewername ID_conc ID_HH ID ID_enfant symptom ///
	time visit spstart tvisit splen visit1 vstart1 vlen1 vtype1

egen id = group(intid inter ID_conc ID_HH ID ID_enfant symptom), m
sort id spstart
qui by id:  gen spell_num = _n

/////////
// re-format from spell data into discrete-time data
rename time spell_len
egen mtime = max(spell_len)
expand mtime
drop mtime
bysort intid inter ID_conc ID_HH ID ID_enfant spstart: gen time=_n
bysort intid inter ID_conc ID_HH ID ID_enfant spstart: gen N=_N
// qui gen y=0
// replace y=1 if visit==1 & time==N
// replace time = time-1 // day 1 --> t = 0

/////////
// now get the data into a format that matlab can use
rename visit exit_state
rename time t

keep  id spell_num spell_len exit_state t age male quart region
sort  id spell_num t

rename age    age1
rename male   male1
rename quart  quart1
rename region region1

gen age2    = age1
gen male2   = male1
gen quart2  = quart1
gen region2 = region1

gen age21 = age1*age1
gen age22 = age2*age2

order id spell_num spell_len exit_state t ///
	age1 age21 male1 quart1 region1 ///
	age2 age22 male2 quart2 region2

// create dummies
/*
foreach str in age { // start numbering at 0
	forvalues j = 1/2 {
		qui tab `str'`j', gen(temp)
		drop `str'`j'
		forvalues i = 1/10 {
			local k = `i'-1
			cap rename temp`i' `str'`k'`j'
		}
	}
}
*/
foreach str in quart region { // start numbering at 1
	forvalues j = 1/2 {
		qui tab `str'`j', gen(temp)
		drop `str'`j'
		forvalues i = 1/10 {
			cap rename temp`i' `str'`i'`j'
		}
	}
}
// drop age01 age02
drop quart11 quart12
drop region11 region12

egen empty = rowmiss(*)
drop if empty
drop empty
outsheet using "/Volumes/SVBROWN 1/Hazard Model/MATLAB/data/data.csv", c nol replace

cap log close Hazard_`fname'
