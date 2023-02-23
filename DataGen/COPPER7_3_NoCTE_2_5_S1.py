# COPPER model for Canada electricity system, national and provincial scale.
# Written by Reza Arjmand, Ph.D candidate at UVic (rezaarjmand@uvic.ca; rezaarjmand@yahoo.com). Linkedin: www.linkedin.com/in/rezaarjmand
# Please CITE this paper (Arjmand, R., & McPherson, M. (2021). Canada’s electricity system transition under alternative policy scenarios. Energy Policy, under review)
# Version 7.2 Nov 2021
# Converted windg solarg from equality to less than equal constraints
# Done: technology evolution, limit new thermal generation, thermal phase out, just small hydro, national carbon limit
# Done: new storage technology, fix pumped continues to include new build
# Done: limited new capacity of a specific generation type
# Done: Solved tranmission  cap problem
# Done: Transmission expansion limit
# New: supply by source
# New: added a provincial breakdown of new installed capacity
# New: added supply by source
# discount
# renewable min max limit per AP
# version 7.3
# New: added cost breakdown


from pyomo.environ import *
import pandas as pd
import numpy as np
import os
import csv
import gc
import time
import sys
from random import randint

start = time.time()
# uncomment when running on CC
#os.chdir('/home/arjmand/COPPER7')
# carbon Tax in each period
ctax={'2025':95,'2030':170,'2035':220,'2040':270,'2045':320,'2050':370}
### planning reserve data
reserve_margin=0.15 ####[CERI 174 p37(21), ReEDS model doc 2018,]
###percent of available hydro that can be upgraded to  store energy (retrofit to add pumped hydro storage)
pump_ret_limit=0.2
discount=0.05
inflation=0
#downsampling or hierarchical clustring
downsampling=False
hierarchical=True
# for test run
test=False
hydro_development=True
autrarky=False
storage_continous=True
provincial_emission_limit=False
emission_limit_ref_year=2017
national_emission_limit=True
nat_em_limit={'2025':100,'2030':11,'2035':0,'2040':0,'2045':0,'2050':0}
## local or national natural gas price
local_gas_price=True
## activate OBPS
OBPS_on=False
## SMR and CCS technologies
SMR_CCS=True
## GPS and CPO regulations
GPS=True
CPO=True
##constrainted transmission expansion, constraints the tranmission expansion capacity to a coefficient of the current capacity
CTE_extant=False
CTE_coef=0
CTE_custom=True
transexpansionlimits=pd.read_csv(r'trans_expansion_limits_NoCTE.csv',header=None)
trans_expansion_limits=dict(transexpansionlimits.values)
## Thermal unit phase out, type:year
thermal_phase_out=True
phase_out_type_year={'gasSC':'2040','gasCC':'2040','coal':'2030','diesel':'2040'}
ph_out_t=list(phase_out_type_year.keys())
## for linkage feedback
min_installed_LB_PHP=False
## new thermal unit restriction, types
new_thermal_limit=False
limited_tplants=['coal','gasCC','gasSC','nuclear','diesel','biomass']
## just consider small hydro projects, under 100 MW
just_small_hydro=False
DSF_tranmission=False
non_emitting_limit=False
## limit electrcity generation emission
nonemitting_limit={'2025':0,'2030':0.9,'2035':0,'2040':0,'2045':0,'2050':0}
## technology evolutiion level
technology_evolution_on={'base':True, 'evolving':False} #False #True
## limited capacity addition of a generation type
limited_new_thermal_gen={}
renewable_portfolio_standards={'wind.max':True,'wind.min':False,'solar.max':False,'solar.min':False}
wind_max_ap_limit = pd.read_csv(r'wind_max_ap_limit.csv',header=0,index_col=0)
wind_min_ap_limit = pd.read_csv(r'wind_min_ap_limit.csv',header=0,index_col=0)
solar_max_ap_limit = pd.read_csv(r'solar_max_ap_limit.csv',header=0,index_col=0)
solar_min_ap_limit = pd.read_csv(r'solar_min_ap_limit.csv',header=0,index_col=0)
DSF_CCS_limit=False
DSF_demand=False
### New parameters (Jake)
max_nuclear_aba = True
max_smr_aba = True
max_bio_aba = True

##### Reading the configuration excel sheet #####

configuration = pd.read_excel (r'COPPER_configuration.xlsx',header=0)
config=dict(zip(list(configuration.iloc[:]['Parameter']),list(configuration.iloc[:]['Value'])))

#Initializing all sets
### All regions
ap=["British Columbia", "Alberta",  "Saskatchewan", "Manitoba", "Ontario","Quebec", "New Brunswick", "Newfoundland and Labrador",  "Nova Scotia","Prince Edward Island"] #all provinces
aba=["British Columbia.a", "Alberta.a",  "Saskatchewan.a", "Manitoba.a", "Ontario.a","Ontario.b","Quebec.a","Quebec.b", "New Brunswick.a", "Newfoundland and Labrador.a","Newfoundland and Labrador.b",  "Nova Scotia.a","Prince Edward Island.a"] #all possible balancing areas
###auxiliary, just for creatin input files####
aba1=['a','b']
pds=['2030','2040','2050'] #pds=['2025','2030','2035','2040','2045','2050']
season=['winter','summer']
capacity_val=pd.read_csv(r'wind_solar_capacity_value.csv',header=None)
header_A=list(capacity_val.iloc[0,:])
header_B=list(capacity_val.iloc[1,:])
del(header_A[0])
del(header_B[0])
ind=list(capacity_val.iloc[:,0])
del(ind[0])
del(ind[0])
capacity_value = pd.DataFrame(np.array(capacity_val.loc[2:12,1:5]), columns = pd.MultiIndex.from_tuples(zip(header_A,header_B)), index=ind)


tech_evolution=pd.read_csv(r'technology_evolution.csv',header=0,index_col=0)
technology_evolution=dict()
for PD in pds:
    for GT in ['wind','solar','LB']:
        if technology_evolution_on['base']:
            technology_evolution[GT+'.'+PD]=tech_evolution[PD][GT+'.'+'base']   ### canada energy future 2020 base scenario
        elif technology_evolution_on['evolving']:
            technology_evolution[GT+'.'+PD]=tech_evolution[PD][GT+'.'+'evolving']
        else:
            technology_evolution[GT+'.'+PD]=1




#### gas PHP installed capacity limit
LB_PHP_installed_limit=dict()
for ABA in aba:
    LB_PHP_installed_limit[ABA]=0

##we can skip runnig some days using smaple_rat, for example if sample_rate=3 it will run days 1,4,7,10,... 
# number of days that we want to run
if downsampling:
   sample_rate=int(config['sample rate'])
   rundaynum=int(config['run day number'])
   rundays=list(range(1,rundaynum+1,sample_rate))
   cap_cost_alter=(365/len(rundays))
### uses hierarchical clustring results   
elif hierarchical:
     if test:
         run_days = pd.read_csv(r'run_days_test.csv',header=None)
     elif DSF_demand:
         run_days = pd.read_csv(r'run_days_dsf.csv',header=None)
     else:
         run_days = pd.read_csv(r'run_days-24.csv',header=None)
     rundays=list(run_days.values)
     rundays=[int(RD) for RD in rundays]
     cap_cost_alter=(365/len(rundays))

#####Storage dataol(x) -
st=['PHS', 'LB']
storage_hours={'PHS':8,'LB':4};
storage_capital={'LB':215664,'PHS':141023}
storage_cost=dict()
for PD in pds:
    for ST in st:
        if ST+'.'+PD in technology_evolution:
            storage_cost[ST+'.'+PD]=storage_capital[ST]*technology_evolution[ST+'.'+PD]/cap_cost_alter
        else:
            storage_cost[ST+'.'+PD]=storage_capital[ST]*1/cap_cost_alter

# storage_cost={'PHS.2030':141026/cap_cost_alter,'PHS.2040':141026/cap_cost_alter,'PHS.2050':141026/cap_cost_alter,'LB.2030':121097/cap_cost_alter,'LB.2040':105979/cap_cost_alter,'LB.2050':90860/cap_cost_alter,
#               'PHS.2025':141026/cap_cost_alter,'PHS.2035':141026/cap_cost_alter,'PHS.2045':141026/cap_cost_alter,'LB.2025':148815/cap_cost_alter,'LB.2035':113538/cap_cost_alter,'LB.2045':98419/cap_cost_alter}
#storage_cost=storage_cost*0.3
store_fix_o_m={'PHS.2030':9000/cap_cost_alter,'PHS.2040':9000/cap_cost_alter,'PHS.2050':9000/cap_cost_alter,'LB.2030':27580/cap_cost_alter,'LB.2040':24140/cap_cost_alter,'LB.2050':20680/cap_cost_alter,
               'PHS.2025':9000/cap_cost_alter,'PHS.2035':9000/cap_cost_alter,'PHS.2045':9000/cap_cost_alter,'LB.2025':33885/cap_cost_alter,'LB.2035':25852/cap_cost_alter,'LB.2045':22410/cap_cost_alter}
storage_efficiency={'PHS':0.8,'LB':0.85}#0.80

runhours=rundays[-1]*24
foryear=int(config['forecast year'])
refyear=int(config['refrence year'])

## converts GJ to MWh
GJtoMWh=config['GJtoMWh']# 3.6
#autonomy_pct=config['autonomy_pct']#0
autonomy_pct={"British Columbia" :0.85, "Alberta" :0.60,  "Saskatchewan" :0.60, "Manitoba" :0.80, "Ontario" :0.80, "Quebec" :0, "New Brunswick" :0.50, "Newfoundland and Labrador" :0,  "Nova Scotia" :0.50,"Prince Edward Island" :0}
#share carbon reduced
carbon_reduction=config['carbon_reduction']#0
#reference case carbon emissions in electricity sector in 2005 in Mt by province (source: https://www.cer-rec.gc.ca/en/data-analysis/energy-markets/provincial-territorial-energy-profiles/provincial-territorial-energy-profiles-explore.html)
carbon_2005_ref={"British Columbia" :1.04, "Alberta" :48.83,  "Saskatchewan" :14.82, "Manitoba" :.36, "Ontario" :33.9, "Quebec" :0.65, "New Brunswick" :7.8, "Newfoundland and Labrador" :0.82,  "Nova Scotia" :10.77,"Prince Edward Island" :0}

#reference case carbon emissions in electricity sector in 2017 in Mt by province (source: https://www.cer-rec.gc.ca/en/data-analysis/energy-markets/provincial-territorial-energy-profiles/provincial-territorial-energy-profiles-explore.html)
carbon_2017_ref={"British Columbia" :0.15, "Alberta" :44.33,  "Saskatchewan" :15.53, "Manitoba" :.07, "Ontario" :1.99, "Quebec" :0.26, "New Brunswick" :3.65, "Newfoundland and Labrador" :1.53,  "Nova Scotia" :6.5,"Prince Edward Island" :0.01}
#maximum carbon emissions in the target year in Mt
carbon_limit=dict()
if emission_limit_ref_year==2017:
    for AP in ap:
        carbon_limit[AP]=carbon_2017_ref[AP]*(1-carbon_reduction)
        
elif emission_limit_ref_year==2005:
    for AP in ap:
        carbon_limit[AP]=carbon_2005_ref[AP]*(1-carbon_reduction)

h=list(range(1,8761))

gltoba=pd.read_csv(r'map_gl_to_ba.csv',header=None)
map_gl_to_ba=dict(gltoba.values)        #map grid locations to balancing areas
del gltoba

gltopr=pd.read_csv(r'map_gl_to_pr.csv',header=None)
map_gl_to_pr=dict(gltopr.values)        #map grid locations to provinces
del gltopr

##### Generation fleets data
if SMR_CCS:
    gendata = pd.read_excel (r'Generation_type_data_SMR_CCS.xlsx',header=0 )
else:    
    gendata = pd.read_excel (r'Generation_type_data.xlsx',header=0 )
allplants=list(gendata.iloc[:]['Type'])

tplants=list()
isthermal=list(gendata.iloc[:]['Is thermal?'])
cc=0
for i in isthermal:
    if i:
        tplants.append(allplants[cc])
    cc+=1

## set of non-emitting generation technologies
non_emitting_tplants=list()
non_emitting=list(gendata.iloc[:]['non-emitting thermal?'])
cc=0
for i in non_emitting:
    if i:
        non_emitting_tplants.append(allplants[cc])
    cc+=1

## Reading thermal units specifications
max_cap_fact=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['max_cap_fact'])))#{'gasCC':0.8, 'gasSC': .3, 'nuclear': 0.95, 'coal': 0.9, 'diesel': 0.95, 'biomass': 0.9} #annual maximum capacity factor
min_cap_fact=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['min_cap_fact'])))#{'gasCC': 0.2, 'gasSC': .02, 'nuclear': 0.75,  'coal': 0.5, 'diesel': 0.05, 'biomass': 0.2} #annual minimum capacity factor
ramp_rate_percent=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['ramp_rate_percent'])))#{'gasCC': 0.1, 'gasSC':0.1 , 'nuclear': 0.05,  'coal': 0.05, 'diesel': 0.1, 'biomass': 0.05} #ramp rate in percent of capacity per hour
efficiency=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['efficiency'])))#{'gasCC': 0.509, 'gasSC': 0.28, 'nuclear': 0.327, 'coal': 0.39, 'diesel': 0.39, 'biomass': 0.39}

### Reading cost data
fixed_o_m=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['fixed_o_m'])))#(fixedom.values) 
for k in fixed_o_m:
    fixed_o_m[k]=fixed_o_m[k]/cap_cost_alter
  
variable_o_m=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['variable_o_m'])))#dict(variableom.values)
fuelprice=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['fuelprice'])))#dict(fuel_price.values) 


capitalcost1=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['capitalcost'])))#dict(capital_cost.values)
if technology_evolution_on['evolving'] and SMR_CCS:
    capitalcost1['SMR']=629574
capitalcost=dict()
for ABA in aba:
    for G in capitalcost1:
        capitalcost[ABA+'.'+G]=capitalcost1[G]


if GPS:
    for ABA in aba:
        capitalcost[ABA+'.gasSC']=1000000000 
        capitalcost[ABA+'.coal']=1000000000

for k in capitalcost:
    capitalcost[k]=capitalcost[k]/cap_cost_alter
    
    
#### Transmission sytem costs
trans_o_m=config['trans_o_m']/cap_cost_alter #10860
transcost=config['transcost']/cap_cost_alter #184
intra_ba_transcost=config['intra_ba_transcost']/cap_cost_alter #557
## renewables (wind and solar) data
merrawindsetall=pd.read_csv(r'merra_wind_set_all.csv',header=None)
wl=list(merrawindsetall.iloc[:,0])

allsolarlocations=pd.read_csv(r'all_solar_locations.csv',header=None)
sl=list(allsolarlocations.iloc[:,0])

gridlocations=pd.read_csv(r'grid_locations.csv',header=None)
gl=list(gridlocations.iloc[:,0])
gl=[str(GL) for GL in gl]

distancetogrid=pd.read_csv(r'distance_to_grid.csv',header=None)
distance_to_grid=dict(distancetogrid.values)
windcost=dict()
solarcost=dict()
for PD in pds:
    for GL in gl:
        windcost[PD+'.'+GL]=capitalcost[map_gl_to_ba[int(GL)]+'.wind']*technology_evolution['wind.'+PD]+distance_to_grid[int(GL)]*intra_ba_transcost
        solarcost[PD+'.'+GL] = capitalcost[map_gl_to_ba[int(GL)]+'.solar']*technology_evolution['solar.'+PD] + distance_to_grid[int(GL)]*intra_ba_transcost



#thermal plant co2 emissions in tonne per MWh of electricity generated
carbondioxide=dict(zip(list(gendata.iloc[:]['Type']),list(gendata.iloc[:]['fuel_co2'])))

#### Adding carbon price to the fuelcost, fuelcost=fuelprice + carbonprice ($/MWh)
## ref CERI study 168 table 3.11 from Statistics Canada (CANSIM Table 129-0003 Sales of natural gas, annual)
if local_gas_price:
    gasprice={"British Columbia.a":2.69, "Alberta.a":2.60,  "Saskatchewan.a":2.55, "Manitoba.a":2.73, "Ontario.a":6.77,"Ontario.b":6.77,"Quebec.a":6.73,"Quebec.b":6.73, "New Brunswick.a":6.21, "Newfoundland and Labrador.a":7.39,"Newfoundland and Labrador.b":7.39,  "Nova Scotia.a":7.39,"Prince Edward Island.a":7.39}
else:
    gasprice=dict()
    for ABA in aba:
        gasprice[ABA]=fuelprice['gasCC']

if OBPS_on:        
    if 'SMR'in tplants:
        OBPS=dict(zip(tplants,[0.37,0.37,0.52,0.37,0.37,0.37,0.37,0.37,0.37]))
    else:
        OBPS=dict(zip(tplants,[0.37,0.52,0.37,0.37,0.37,0.37]))
        
    fuelcost=dict()
    for PD in pds:
        for TP in tplants:
            
            for ABA in aba:
                if TP=='gasCC' or TP=='gasSC':
                    fuelcost[PD+'.'+TP+'.'+ABA] = (gasprice[ABA]/efficiency[TP])*GJtoMWh
                else:
                    fuelcost[PD+'.'+TP+'.'+ABA] = (fuelprice[TP]/efficiency[TP])*GJtoMWh
                if ABA=='British Columbia.a'or TP=='biomass' or TP=='SMR' or TP=='nuclear':
                    fuelcost[PD+'.'+TP+'.'+ABA]+=carbondioxide[TP]*ctax[PD]
                else:
                    
                    fuelcost[PD+'.'+TP+'.'+ABA]+=(carbondioxide[TP]-OBPS[TP])*ctax[PD]
else:
    fuelcost=dict()
    for PD in pds:
        for TP in tplants:
        
            for ABA in aba:
                if TP=='gasCC' or TP=='gasSC':
                    fuelcost[PD+'.'+TP+'.'+ABA] = (gasprice[ABA]/efficiency[TP])*GJtoMWh
                else:
                    fuelcost[PD+'.'+TP+'.'+ABA] = (fuelprice[TP]/efficiency[TP])*GJtoMWh
              
                fuelcost[PD+'.'+TP+'.'+ABA]+=carbondioxide[TP]*ctax[PD]

### Reading input file for nuclear development (By balancing area)
maxnuclear = pd.read_csv(r'max_nuclear.csv',header=0)
max_nuclear=dict(zip(list(maxnuclear.iloc[:]['ABA']),list(maxnuclear.iloc[:]['MAX_nuclear'])))

### Reading input file for SMR development (By balancing area)
maxsmr = pd.read_csv(r'max_smr.csv',header=0)
max_smr=dict(zip(list(maxsmr.iloc[:]['ABA']),list(maxsmr.iloc[:]['MAX_SMR'])))

### Reading input file for biomass development (By balancing area) - Added by Jake - 18 SEP 2022
maxbio = pd.read_csv(r'biomass_max_aba_limit.csv',header=0)
max_bio=dict(zip(list(maxbio.iloc[:]['ABA']),list(maxbio.iloc[:]['BIO_MAX'])))

### Reading transmission expansion routes            
if not DSF_tranmission:
    transmapba=pd.read_csv(r'transmission_map_ba.csv',header=None)
else:
    transmapba=pd.read_csv(r'transmission_map_ba_DSF.csv',header=None)
transmap=list(transmapba.iloc[:,0])        #map grid locations to provinces
del transmapba

transmapdis = pd.read_csv(r'transmission_map_distance.csv',header=None)
distance=dict(transmapdis.values)
del transmapdis
### Reading extant wind and solar data
extantwindsolar = pd.read_csv(r'extant_wind_solar.csv',header=0)
extant_wind_solar=list()
for PD in pds:
    extant_wind_solar.append(dict(zip(list(extantwindsolar.iloc[:]['location']),list(extantwindsolar.iloc[:][PD]))))

del extantwindsolar

### Reading extant generation capacity data
extantcapacity = pd.read_csv(r'extant_capacity.csv',header=0)
extant_capacity=dict()
for PD in pds:
    label_ABA=list(extantcapacity.iloc[:]['ABA'])
    label_ABA=[PD+'.'+LL for LL in label_ABA]
    excap=dict(zip(label_ABA,list(extantcapacity.iloc[:][PD])))
    extant_capacity.update(excap)

### Reading extant tranmission lines capacity data
extanttrans = pd.read_csv(r'extant_transmission.csv',header=0)
extant_transmission=list()
for PD in pds:
    extant_transmission.append(dict(zip(list(extanttrans.iloc[:]['ABA']),list(extanttrans.iloc[:][PD]))))
del extanttrans

##### Reading hydro units capacity factor data
hydrocf = pd.read_csv(r'hydro_cf.csv',header=None)
hydro_cf=dict(hydrocf.values)
del hydrocf

## ### Reading demand related data
demand_growth = pd.read_csv(r'annual_growth_2_5.csv',header=0,index_col=0)

demandall1 = pd.read_csv(r'demand.csv',header=None)
demandall=dict(demandall1.values)
del demandall1

##### Reading population data to disaggregate demand data
population1 = pd.read_csv(r'population.csv',header=None)
population=dict(population1.values)
del population1

demandus = pd.read_csv(r'us_demand.csv',header=None)
demand_us=dict(demandus.values)
# if DSF_demand:
#     for K in demand_us.keys():
#         demand_us[K]=0
del demandus

maphd = pd.read_csv(r'set_map_days_to_hours.csv',header=None)
map_hd=dict(maphd.values)
del maphd

maphm = pd.read_csv(r'set_map_months_to_hours.csv',header=None)
map_hm=dict(maphm.values)
del maphm

surfacearea = pd.read_csv(r'surface_area.csv',header=None)
surface_area=dict(surfacearea.values)
del surfacearea

demand13to18 = pd.read_csv(r'growth_2013_2018.csv',header=None)
demand_13_to_18=dict(demand13to18.values)
del demand13to18

## according to the sample_rate these lines remove the days that we don't want to run

demandall_2018=dict()
demand_all=dict()
if DSF_demand:
    demanddsf=maphm = pd.read_csv(r'demand_dsf_dict.csv',header=None)
    demand_all=dict(demanddsf.values)
## Calculating seasonal peak demand
national_demand=np.zeros((len(h),len(pds)))
peak_demand=dict()
peak_days=dict()
for PD in pds:
    
    for H in h:
        ND=0
        for AP in ap:
            if not DSF_demand:
                demandall_2018[AP+'.'+str(H)]=demandall[AP+'.'+str(H)]*(1+demand_13_to_18[AP])
                if pds.index(PD)==0:
                    demand_all[PD+'.'+AP+'.'+str(H)]=demandall_2018[AP+'.'+str(H)]*(1+demand_growth[PD][AP])**(int(PD)-refyear)
                else:
                    demand_all[PD+'.'+AP+'.'+str(H)]=demand_all[pds[pds.index(PD)-1]+'.'+AP+'.'+str(H)]*(1+demand_growth[PD][AP])**(int(PD)-int(pds[pds.index(PD)-1]))
                    
            ND+=demand_all[PD+'.'+AP+'.'+str(H)]
        national_demand[h.index(H),pds.index(PD)]=ND
    peak_demand[PD+'.'+'winter']= max(max(national_demand[:2160,pds.index(PD)]),max(national_demand[6480:8760,pds.index(PD)]))   
    peak_demand[PD+'.'+'summer']= max(national_demand[2160:6480,pds.index(PD)])
    peak_days[PD+'.'+'summer']=map_hd[int(np.where(national_demand == peak_demand[PD+'.'+'summer'])[0])]
    peak_days[PD+'.'+'winter']=map_hd[int(np.where(national_demand == peak_demand[PD+'.'+'winter'])[0])]


del h[runhours:]
hours=len(h)

nummonths=map_hm[rundays[-1]*24]
m=list(range(1,nummonths+1))
       
d=rundays.copy()
h3=h.copy()
for H in h3:
    if map_hd[H] not in rundays:
        h.remove(H)

del h3
hours=len(h)
h2=h.copy()
del h2[hours-1]


## Reading hourly wind and solar CF data in each grid cell 
with open('windcf.csv') as csv_file:
    reader = csv.reader(csv_file)
    windcf = dict(reader)

for k in windcf:
    windcf[k]=float(windcf[k])

with open('solarcf.csv') as csv_file:
    reader = csv.reader(csv_file)
    solarcf = dict(reader)

for k in solarcf:
    solarcf[k]=float(solarcf[k])


meanwincf=sum(list(windcf.values()))/len(list(windcf.values()))
meansolcf=sum(list(solarcf.values()))/len(list(solarcf.values()))
meanwincf_reduced=sum(windcf[str(H)+'.'+GL] for H in h for GL in gl)/(len(gl)*24*len(rundays))
meansolarcf_reduced=sum(solarcf[str(H)+'.'+GL] for H in h for GL in gl)/(len(gl)*24*len(rundays))


##calculte the diffrence between to hours in a row
time_diff=dict()
for I in list(range(len(h)-1)):
    time_diff[h[I]]=h[I+1]-h[I]



gl2=gl.copy()
for GL in gl2:
    if map_gl_to_pr[int(GL)] not in ap:
        gl.remove(GL)
        
###maximum wind and solar capacity that can be installed in each grid cell in MW per square km
maxwindperkmsq=config['maxwindperkmsq'] #2
maxsolarperkmsq=config['maxsolarperkmsq'] #31.28
maxwind=dict()
maxsolar=dict()
for GL in gl:
    maxwind[GL]=surface_area[int(GL)]*maxwindperkmsq
    maxsolar[GL]=surface_area[int(GL)]*maxsolarperkmsq
## Extant PHS data
ba_storage_capacity=dict()
for ABA in aba:
    for ST in st:
        ba_storage_capacity[ABA+'.'+ST]=0
ba_storage_capacity['Ontario.a'+'.'+'PHS']=174

#Calculating transmission loss coefficients
translossfixed=config['translossfixed'] #0.02
translossperkm=config['translossperkm'] #0.00003
transloss=dict()
for ABA in aba:
    for ABBA in aba:
        if ABA+'.'+ABBA in distance:
            transloss[ABA+'.'+ABBA]=distance[ABA+'.'+ABBA]*translossperkm+translossfixed
## disaggregating demand based on population
populationaba=dict()
for ABA in aba:
    populationaba[ABA]=0
for ABA in aba:
    for GL in population:
        if map_gl_to_ba[int(GL)]==ABA:
                populationaba[ABA]=populationaba[ABA]+population[GL]
populationap=dict()
demand=dict()

for PD in pds:
    for AP in ap:
        populationap[AP]=sum(populationaba[ABA] for ABA in aba if AP in ABA)
    for ABA in aba:
        pvba=ABA.replace('.a','')
        pvba=pvba.replace('.b','')
        for H in h:
            demand[PD+'.'+ABA+'.'+str(H)]=demand_all[PD+'.'+pvba+'.'+str(H)]*(populationaba[ABA]/populationap[pvba])
## CREating extant generation dic
extant_thermal=dict()
for AP in ap:
    for ABA in aba1:
        for TP in tplants:
            for PD in pds:
                extant_thermal[PD+'.'+AP+'.'+ABA+'.'+TP]=0
                if PD+'.'+AP+'.'+ABA+'.'+TP in extant_capacity:
                    extant_thermal[PD+'.'+AP+'.'+ABA+'.'+TP]=extant_capacity[PD+'.'+AP+'.'+ABA+'.'+TP]
                if CPO and TP=='coal' and int(PD)>=2030:
                    extant_thermal[PD+'.'+AP+'.'+ABA+'.'+TP]=0

hydro_capacity=dict()
extant_wind_gen=dict()
extant_solar_gen=dict()
for PD in pds:
    for AP in ap:
        for ABA in aba1:
    
            for H in h:
                extant_wind_gen[PD+'.'+AP+'.'+ABA+'.'+str(H)]=0
                extant_solar_gen[PD+'.'+AP+'.'+ABA+'.'+str(H)]=0

#CALCULATE OUTPUT POWER FOR EXTANT WIND FARM AND SOLAR POWER PLANTS  
for PD in pds:          
    for GL in gl:
        for EXG in extant_wind_solar[0]:
            #if int(PD)<=2030:
    #            if EXG==str(GL)+'.'+'hydro':
    #                hydro_capacity[map_gl_to_ba[int(GL)]+'.'+'hydro']+=extant_generation[EXG]
            if EXG==str(GL)+'.'+'wind':
                for H in h:
                    extant_wind_gen[PD+'.'+map_gl_to_ba[int(GL)]+'.'+str(H)]+=(extant_wind_solar[pds.index(PD)][EXG]*windcf[str(H)+'.'+str(GL)])
            if EXG==str(GL)+'.'+'solar':
                for H in h:
                    extant_solar_gen[PD+'.'+map_gl_to_ba[int(GL)]+'.'+str(H)]+=(extant_wind_solar[pds.index(PD)][EXG]*solarcf[str(H)+'.'+str(GL)])


## calculating hydro CF for different type of hydro power plant

hydro_minflow=0.1
ror_hydroout=dict()
day_hydroout=dict()
month_hydroout=dict()
day_hydro_historic=dict()
month_hydro_historic=dict()
ror_hydro_capacity=dict()
day_hydro_capacity=dict()
month_hydro_capacity=dict()
day_minflow=dict()
month_minflow=dict()
for PD in pds:
    for AP in ap:
       for ABA in aba1:
           ror_hydro_capacity[PD+'.'+AP+'.'+ABA] = 0
           day_hydro_capacity[PD+'.'+AP+'.'+ABA] = 0
           month_hydro_capacity[PD+'.'+AP+'.'+ABA] = 0
           if PD+'.'+AP+'.'+ABA+'.'+'hydro_run' in extant_capacity:
               ror_hydro_capacity[PD+'.'+AP+'.'+ABA] = extant_capacity[PD+'.'+AP+'.'+ABA+'.'+'hydro_run']
               day_hydro_capacity[PD+'.'+AP+'.'+ABA] = extant_capacity[PD+'.'+AP+'.'+ABA+'.'+'hydro_daily']
               month_hydro_capacity[PD+'.'+AP+'.'+ABA] = extant_capacity[PD+'.'+AP+'.'+ABA+'.'+'hydro_monthly']
           day_minflow[PD+'.'+AP+'.'+ABA]=day_hydro_capacity[PD+'.'+AP+'.'+ABA]*hydro_minflow
           month_minflow[PD+'.'+AP+'.'+ABA]=month_hydro_capacity[PD+'.'+AP+'.'+ABA]*hydro_minflow
           for D in d:
               day_hydro_historic[PD+'.'+str(D)+'.'+AP+'.'+ABA]=0
               
           for M in m:
               month_hydro_historic[PD+'.'+str(M)+'.'+AP+'.'+ABA]=0
           for H in h:

             if AP+'.'+str(H) in hydro_cf:
                  ror_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA]=ror_hydro_capacity[PD+'.'+AP+'.'+ABA]*hydro_cf[AP+'.'+str(H)]
                  day_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA] = day_hydro_capacity[PD+'.'+AP+'.'+ABA]*hydro_cf[AP+'.'+str(H)]
                  month_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA] = month_hydro_capacity[PD+'.'+AP+'.'+ABA]*hydro_cf[AP+'.'+str(H)]
             else:
                  ror_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA] =0
                  day_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA] =0
                  month_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA] =0
                
           
             day_hydro_historic[PD+'.'+str(map_hd[H])+'.'+AP+'.'+ABA]=day_hydro_historic[PD+'.'+str(map_hd[H])+'.'+AP+'.'+ABA]+day_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA]
             month_hydro_historic[PD+'.'+str(map_hm[H])+'.'+AP+'.'+ABA]=month_hydro_historic[PD+'.'+str(map_hm[H])+'.'+AP+'.'+ABA]+month_hydroout[PD+'.'+str(H)+'.'+AP+'.'+ABA]
            


## hydro renewal and greenfield development necessary inputs 
if hydro_development:
    if storage_continous and not just_small_hydro:
        #hydro_new = pd.read_excel (r'hydro_new_recon_nopump.xlsx',header=0)
        #Read from the input file with the correct number of splits
        hydro_new = pd.read_excel (r'hydro_new_recon_nopump.xlsx',header=0)
    elif not just_small_hydro:
        hydro_new = pd.read_excel (r'hydro_new_recon.xlsx',header=0)
    elif just_small_hydro:
        hydro_new = pd.read_excel (r'hydro_new_recon_nopump_justsmall.xlsx',header=0)
        
              
    
    hydro_renewal=list(hydro_new.iloc[:]['Short Name'])
    cost_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Annualized Capital Cost ($M/year)'])))
    capacity_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Additional Capacity (MW)'])))
    devperiod_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Development Time (years)'])))
    location_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Balancing Area'])))
    distance_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Distance to Grid (km)'])))
    type_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Type'])))
    fixed_o_m_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Fixed O&M ($/MW-year)'])))
    variable_o_m_renewal=dict(zip(hydro_renewal,list(hydro_new.iloc[:]['Variable O&M ($/MWh)'])))
    
    hr_ror=list()
    cost_ror_renewal=dict()
    capacity_ror_renewal=dict()
    hr_ror_location=dict()
    
    hr_day=list()
    cost_day_renewal=dict()
    capacity_day_renewal=dict()
    hr_day_location=dict()
    
    hr_mo=list()
    cost_month_renewal=dict()
    capacity_month_renewal=dict()
    hr_month_location=dict()
    
    hr_pump=list()
    cost_pump_renewal=dict()
    capacity_pump_renewal=dict()
    hr_pump_location=dict()
    for k in hydro_renewal:
        if foryear-2020>=devperiod_renewal[k]:
            if type_renewal[k]=='hydro_run':
                hr_ror.append(k)
                cost_ror_renewal[k]=cost_renewal[k]*1000000+distance_renewal[k]*intra_ba_transcost*capacity_renewal[k]
                capacity_ror_renewal[k]=capacity_renewal[k]
                hr_ror_location[k]=location_renewal[k]
    
            if type_renewal[k]=='hydro_daily':
                hr_day.append(k)
                cost_day_renewal[k]=cost_renewal[k]*1000000+distance_renewal[k]*intra_ba_transcost*capacity_renewal[k]
                capacity_day_renewal[k]=capacity_renewal[k]
                hr_day_location[k]=location_renewal[k]
    
            if type_renewal[k]=='hydro_monthly':
                hr_mo.append(k)
                cost_month_renewal[k]=cost_renewal[k]*1000000+distance_renewal[k]*intra_ba_transcost*capacity_renewal[k]
                capacity_month_renewal[k]=capacity_renewal[k]
                hr_month_location[k]=location_renewal[k]
            if type_renewal[k]=='hydro_pump':
                hr_pump.append(k)
                cost_pump_renewal[k]=cost_renewal[k]*1000000+distance_renewal[k]*intra_ba_transcost*capacity_renewal[k]
                capacity_pump_renewal[k]=capacity_renewal[k]
                hr_pump_location[k]=location_renewal[k]
    
    
                  
               
    ror_renewalout=dict()
    for HR_ROR in hr_ror:
        for H in h:
            province_loc=hr_ror_location[HR_ROR].replace('.a','')
            province_loc=province_loc.replace('.b','')
            ror_renewalout[str(H)+'.'+HR_ROR]=capacity_ror_renewal[HR_ROR]*hydro_cf[province_loc+'.'+str(H)]
    
    for k in cost_ror_renewal:
        cost_ror_renewal[k]=cost_ror_renewal[k]/cap_cost_alter        
            
    
    day_renewal_historic=dict()
    day_renewalout=dict()
    for HR_DAY in hr_day:
        for D in d:
            day_renewal_historic[str(D)+'.'+HR_DAY]=0 
        for H in h:
            province_loc=hr_day_location[HR_DAY].replace('.a','')
            province_loc=province_loc.replace('.b','')
            day_renewalout[str(H)+'.'+HR_DAY]=capacity_day_renewal[HR_DAY]*hydro_cf[province_loc+'.'+str(H)]
            day_renewal_historic[str(map_hd[H])+'.'+HR_DAY]=day_renewal_historic[str(map_hd[H])+'.'+HR_DAY]+day_renewalout[str(H)+'.'+HR_DAY]
    
    for k in cost_day_renewal:
        cost_day_renewal[k]=cost_day_renewal[k]/cap_cost_alter
    
    
    month_renewal_historic=dict()
    month_renewalout=dict()
    for HR_MO in hr_mo:
        for M in m:
            month_renewal_historic[str(M)+'.'+HR_MO]=0
        for H in h:
            province_loc=hr_month_location[HR_MO].replace('.a','')
            province_loc=province_loc.replace('.b','')
            month_renewalout[str(H)+'.'+HR_MO]=capacity_month_renewal[HR_MO]*hydro_cf[province_loc+'.'+str(H)]
            month_renewal_historic[str(map_hm[H])+'.'+HR_MO]=month_renewal_historic[str(map_hm[H])+'.'+HR_MO]+month_renewalout[str(H)+'.'+HR_MO]
    
    for k in cost_month_renewal:
        cost_month_renewal[k]=cost_month_renewal[k]/cap_cost_alter
    
    
    for k in cost_pump_renewal:
        cost_pump_renewal[k]=cost_pump_renewal[k]/cap_cost_alter
    
#####recontract input data ###########33
windcost_recon=dict()
solarcost_recon=dict()

for PD in pds:
    windcost_recon[PD]=config['windcost_recon']*technology_evolution['wind.'+PD]/cap_cost_alter  #118662.7485
    solarcost_recon[PD]=config['solarcost_recon']*technology_evolution['solar.'+PD]/cap_cost_alter #94204.76075


windsolarrecon = pd.read_csv(r'wind_solar_location_recon.csv',header=0)

wind_recon_capacity=dict()
solar_recon_capacity=dict()
for PD in pds:
    for GL in gl:
        wind_recon_capacity[PD+'.'+GL]=0
        solar_recon_capacity[PD+'.'+GL]=0
        
        if str(GL)+'.'+'wind' in list(windsolarrecon[:]['location']):
            wind_recon_capacity[PD+'.'+GL]=windsolarrecon[PD][list(windsolarrecon[:]['location']).index(str(GL)+'.'+'wind')]
        if str(GL)+'.'+'solar' in list(windsolarrecon[:]['location']):
            solar_recon_capacity[PD+'.'+GL]=windsolarrecon[PD][list(windsolarrecon[:]['location']).index(str(GL)+'.'+'solar')]
        

cleared_data=gc.collect()
######### Doublicating some sets ###############
ttplants=tplants
ggl=gl
hh=h
app=ap
abba=aba

end = time.time()
print(f'\n==================================================\n\
Initializing input data time (Sec): {round((end-start)/60)} Min and {round((end-start)%60)} Sec \
\n==================================================')
start=time.time()

########Creating the Optimization Model######################

model = ConcreteModel()

#### Defineing variables####

model.capacity_therm=Var( pds, aba,tplants, within=NonNegativeReals,initialize=1) #new thermal plant capacity in MW
model.retire_therm=Var( pds, aba,tplants, within=NonNegativeReals,initialize=0)  #retire extant thermal capacity in MW
model.capacity_wind=Var( pds, gl, within=NonNegativeReals,initialize=0)  #wind plant capacity in MW
model.capacity_solar=Var( pds, gl, within=NonNegativeReals,initialize=0)  #solar plant capacity in MW
if storage_continous:
    model.capacity_storage=Var(pds, st, aba, within=NonNegativeReals,initialize=0)  #storage plant capacity in MW
    
model.supply=Var(pds, h,aba,tplants, within=NonNegativeReals,initialize=0)  #fossil fuel supply in MW
model.windout=Var(pds, h,aba, within=NonNegativeReals,initialize=0)  #wind hourly power output
model.solarout=Var(pds, h,aba, within=NonNegativeReals,initialize=0)  #solar hourly power output
model.storageout=Var(pds,st,h,aba, within=NonNegativeReals,initialize=0)  #pumped hydro hourly output
model.storagein=Var(pds,st, h,aba, within=NonNegativeReals,initialize=0)  #pumped hydro hourly input
model.storageenergy=Var(pds,st, h,aba, within=NonNegativeReals,initialize=0)  #total stored pump hydro energy in MWh
model.daystoragehydroout=Var(pds, h,aba, within=NonNegativeReals,initialize=0)  #day storage hydro output in MW
model.monthstoragehydroout=Var(pds, h,aba, within=NonNegativeReals,initialize=0)  #month storage hydro output in MW
model.transmission=Var(pds, h,aba,abba, within=NonNegativeReals,initialize=0)  #hourly transmission in MW from ap,ba to apa,abba
model.capacity_transmission=Var(pds, aba,abba, within=NonNegativeReals,initialize=0)  #transmission capacity in MW from ap,aba to app,abba
#model.carbon=Var(ap,aba, within=NonNegativeReals,initialize=0)  #carbon emissions annual in Mt
if hydro_development:
    model.ror_renewal_binary=Var(pds, hr_ror, within=Binary,initialize=1)
    model.day_renewal_binary=Var(pds, hr_day, within=Binary,initialize=1)
    model.month_renewal_binary=Var(pds, hr_mo, within=Binary,initialize=1)
    model.dayrenewalout=Var(pds, h,hr_day, within=NonNegativeReals,initialize=0)
    model.monthrenewalout=Var(pds, h,hr_mo, within=NonNegativeReals,initialize=0)
    if not storage_continous:
        model.pumphydro=Var(pds, hr_pump, within=Binary,initialize=0)
model.capacity_wind_recon=Var( pds, gl, within=NonNegativeReals,initialize=0)  #wind recontract capacity in MW
model.capacity_solar_recon=Var( pds, gl, within=NonNegativeReals,initialize=0)  #solar recontract capacity in MW

## Intializing storages' energy
for ABA in aba:
    for PD in pds:
        for ST in st:
            model.storageenergy[PD,ST,h[0],ABA].fix(0)

m_counter=2
for H in h:
    if map_hm[H]==m_counter:
        for ABA in aba:
            for PD in pds:
                for ST in st:
                    model.storageenergy[PD,ST,H,ABA].fix(0)
        m_counter+=1
if DSF_CCS_limit:
    for PD in pds:
        for ABA in aba:
            if ABA!='Alberta.a' and ABA!='Saskatchewan.a' and ABA!='Ontario.a' and ABA!='Ontario.b':
                model.capacity_therm[PD,ABA,'gasccs'].fix(0)
         

###Objective function total cost minimization###
def obj_rule(model):
    capcost=dict()
    fcost=dict()
    fixedOM=dict()
    variableOM=dict()
    hydrorenewalccost=dict()
    hydrorenewalomcost=dict()
    newstorage_ccost=dict()
    newstorage_omcost=dict()
    tcapcost=0
    tfcost=0
    tfixedOM=0
    tvariableOM=0
    thydrorenewalccost=0
    thydrorenewalomcost=0
    tnewstorage_ccost=0
    tnewstorage_omcost=0
    for PD in pds:
        ind=pds.index(PD)
        capcost[PD]=sum(model.capacity_therm[PDD,ABA,TP]*capitalcost[ABA+'.'+TP]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1] for TP in tplants for ABA in aba)\
              +sum(model.capacity_wind[PDD,GL] * windcost[PDD+'.'+GL]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1]  for GL in gl)\
              +sum(model.capacity_solar[PDD,GL] * solarcost[PDD+'.'+GL]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1]  for GL in gl)\
              +sum(model.capacity_wind_recon[PDD,GL] * windcost_recon[PDD]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1]  for GL in gl)\
              +sum(model.capacity_solar_recon[PDD,GL] * solarcost_recon[PDD]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1]  for GL in gl)\
              +sum(model.capacity_transmission[PDD,ABA,ABBA]*transcost*distance[ABA+'.'+ABBA]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1]  for ABA in aba for ABBA in aba if ABA+'.'+ABBA in transmap)
           
        fcost[PD]=sum(model.supply[PD,H,ABA,TP] * fuelcost[PD+'.'+TP+'.'+ABA]  for H in h for ABA in aba for TP in tplants)
           
        fixedOM[PD]=sum((extant_thermal[pds[0]+'.'+ABA+'.'+TP]+sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1]))* fixed_o_m[TP]  for ABA in aba for TP in tplants)\
               +sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL])* fixed_o_m['wind']  for PDD in pds[:ind+1] for GL in gl)\
               +sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'wind'] * fixed_o_m['wind']  for GL in gl if str(GL)+'.'+'wind' in extant_wind_solar[0])\
               +sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL])* fixed_o_m['solar']  for PDD in pds[:ind+1] for GL in gl)\
               +sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'solar'] * fixed_o_m['solar']  for GL in gl if str(GL)+'.'+'solar' in extant_wind_solar[0])\
               +sum(model.capacity_transmission[PDD,ABA,ABBA]*trans_o_m  for PDD in pds[:ind+1] for ABA in aba for ABBA in aba if ABA+'.'+ABBA in transmap)\
               +sum(extant_transmission[pds.index(PD)][ABA+'.'+ABBA]*trans_o_m  for ABA in aba for ABBA in aba if ABA+'.'+ABBA in extant_transmission[pds.index(PD)])\
               +sum(ror_hydro_capacity[PD+'.'+ABA] * fixed_o_m['hydro']  for ABA in aba)\
               +sum(day_hydro_capacity[PD+'.'+ABA] * fixed_o_m['hydro']  for ABA in aba)\
               +sum(month_hydro_capacity[PD+'.'+ABA] * fixed_o_m['hydro']  for ABA in aba)
           
        variableOM[PD]=sum(model.supply[PD,H,ABA,TP] * variable_o_m[TP]  for H in h for ABA in aba for TP in tplants)\
              +sum(model.windout[PD,H,ABA]*variable_o_m['wind']  for H in h for ABA in aba)\
              +sum(model.solarout[PD,H,ABA]*variable_o_m['solar']  for H in h for ABA in aba)\
              +sum(ror_hydroout[PD+'.'+str(H)+'.'+ABA]*variable_o_m['hydro']  for H in h for ABA in aba)\
              +sum(model.daystoragehydroout[PD,H,ABA]*variable_o_m['hydro']  for H in h for ABA in aba)\
              +sum(model.monthstoragehydroout[PD,H,ABA]*variable_o_m['hydro']  for H in h for ABA in aba)
           
        hydrorenewalccost[PD]=0
        hydrorenewalomcost[PD]=0
        
        if hydro_development:
            hydrorenewalccost[PD]=sum(cost_ror_renewal[HR_ROR]*model.ror_renewal_binary[PDD,HR_ROR]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1] for HR_ROR in hr_ror)\
                  +sum(cost_day_renewal[HR_DAY]*model.day_renewal_binary[PDD,HR_DAY]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1] for HR_DAY in hr_day)\
                  +sum(cost_month_renewal[HR_MO]*model.month_renewal_binary[PDD,HR_MO]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1] for HR_MO in hr_mo)
            hydrorenewalomcost[PD]=sum(capacity_ror_renewal[HR_ROR]*model.ror_renewal_binary[PDD,HR_ROR]*fixed_o_m_renewal[HR_ROR] for PDD in pds[:ind+1] for HR_ROR in hr_ror)\
                  +sum(capacity_day_renewal[HR_DAY]*model.day_renewal_binary[PDD,HR_DAY]*fixed_o_m_renewal[HR_DAY] for PDD in pds[:ind+1] for HR_DAY in hr_day)\
                  +sum(capacity_month_renewal[HR_MO]*model.month_renewal_binary[PDD,HR_MO]*fixed_o_m_renewal[HR_MO] for PDD in pds[:ind+1] for HR_MO in hr_mo)\
                  +sum(ror_renewalout[str(H)+'.'+HR_ROR]*model.ror_renewal_binary[PDD,HR_ROR]*variable_o_m_renewal[HR_ROR] for PDD in pds[:ind+1] for H in h for HR_ROR in hr_ror)\
                  +sum(model.dayrenewalout[PD,H,HR_DAY]*variable_o_m_renewal[HR_DAY]  for H in h for HR_DAY in hr_day)\
                  +sum(model.monthrenewalout[PD,H,HR_MO]*variable_o_m_renewal[HR_MO]  for H in h for HR_MO in hr_mo)
                  
            if not storage_continous:
                
                hydrorenewalccost[PD]+=sum(cost_pump_renewal[HR_PUMP]*model.pumphydro[PDD,HR_PUMP] for PDD in pds[:ind+1] for HR_PUMP in hr_pump)
                hydrorenewalomcost[PD]+=sum(capacity_pump_renewal[HR_PUMP]*model.pumphydro[PDD,HR_PUMP]*fixed_o_m_renewal[HR_PUMP] for PDD in pds[:ind+1] for HR_PUMP in hr_pump)
        newstorage_ccost[PD]=0
        newstorage_omcost[PD]=0
        if storage_continous:
            newstorage_ccost[PD]=sum(model.capacity_storage[PDD,ST,ABA] * storage_cost[ST+'.'+PDD]*((1+inflation)**(int(PDD)-2018)) for PDD in pds[:ind+1] for ST in st for ABA in aba)
            newstorage_omcost[PD]=sum(store_fix_o_m[ST+'.'+PDD] * model.capacity_storage[PDD,ST,ABA] for PDD in pds[:ind+1] for ST in st for ABA in aba)
        disc_coef=(1/(1+discount)**(int(PD)-2018))
        inf_coef=((1+inflation)**(int(PD)-2018))
        tcapcost+=capcost[PD]*disc_coef
        tfcost+=fcost[PD]*disc_coef*inf_coef
        tfixedOM+=fixedOM[PD]*disc_coef*inf_coef
        tvariableOM+=variableOM[PD]*disc_coef*inf_coef
        thydrorenewalccost+=hydrorenewalccost[PD]*disc_coef
        thydrorenewalomcost+=hydrorenewalomcost[PD]*disc_coef*inf_coef
        tnewstorage_ccost+=newstorage_ccost[PD]*disc_coef
        tnewstorage_omcost+=newstorage_omcost[PD]*disc_coef*inf_coef
    return (tcapcost+tfcost+tvariableOM+thydrorenewalccost+thydrorenewalomcost+tnewstorage_ccost+tnewstorage_omcost+tfixedOM)

model.obj = Objective(rule=obj_rule,sense=minimize)

######Planning reserve requirment#####

def planning_reserve(model,PD,SEAS):
    ind=pds.index(PD)
    cap_val=sum((extant_thermal[pds[0]+'.'+ABA+'.'+TP]+sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1] )) for ABA in aba for TP in tplants)\
            +sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL])*float(capacity_value[SEAS]['wind'][map_gl_to_pr[int(GL)]]) for PDD in pds[:ind+1] for GL in gl)+sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'wind'] * float(capacity_value[SEAS]['wind'][map_gl_to_pr[int(GL)]])  for GL in gl if str(GL)+'.'+'wind' in extant_wind_solar[0] )\
            +sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL])*float(capacity_value[SEAS]['solar'][map_gl_to_pr[int(GL)]]) for PDD in pds[:ind+1] for GL in gl)+sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'solar'] * float(capacity_value[SEAS]['solar'][map_gl_to_pr[int(GL)]]) for GL in gl if str(GL)+'.'+'solar' in extant_wind_solar[0] )\
            +sum(ror_hydro_capacity[PD+'.'+ABA] for ABA in aba)\
            +sum(day_hydro_capacity[PD+'.'+ABA] for ABA in aba)\
            +sum(month_hydro_capacity[PD+'.'+ABA] for ABA in aba)\
    
    if hydro_development:
        cap_val+=sum(capacity_ror_renewal[HR_ROR]*model.ror_renewal_binary[PDD,HR_ROR] for PDD in pds[:ind+1] for HR_ROR in hr_ror)\
                +sum(capacity_day_renewal[HR_DAY]*model.day_renewal_binary[PDD,HR_DAY] for PDD in pds[:ind+1] for HR_DAY in hr_day)\
                +sum(capacity_month_renewal[HR_MO]*model.month_renewal_binary[PDD,HR_MO] for PDD in pds[:ind+1] for HR_MO in hr_mo)
    if storage_continous:
        cap_val+=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[:ind+1] for ST in st for ABA in aba)
            
    return cap_val >= peak_demand[PD+'.'+SEAS]*(reserve_margin+1)
model.planning_reserve=Constraint(pds,season, rule=planning_reserve)

###constrain retirements to extant plants
def retire(model,PD,ABA,TP):
    ind=pds.index(PD)
    return model.retire_therm[PD,ABA,TP] <= extant_thermal[pds[0]+'.'+ABA+'.'+TP]+sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind])
model.retire=Constraint(pds,aba,tplants, rule=retire)

#### forces the model to retire the plants that their lifes ended
def lifetime(model,PD,ABA,TP):
    ind=pds.index(PD)
    ex_thermal=extant_thermal[PD+'.'+ABA+'.'+TP]
    
    #######for under construction units, this will prevent infeasibility
    if ind>=1:
        if extant_thermal[PD+'.'+ABA+'.'+TP]-extant_thermal[pds[ind-1]+'.'+ABA+'.'+TP]>0:
            ex_thermal=100000

    return extant_thermal[pds[0]+'.'+ABA+'.'+TP]-sum(model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1]) <= ex_thermal
model.lifetime=Constraint(pds,aba,tplants, rule=lifetime)

###wind generation limit
def windg(model,PD,H,ABA):
    ind=pds.index(PD)
    return model.windout[PD,H,ABA]==sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL])*windcf[str(H)+'.'+str(GL)] for PDD in pds[:ind+1] for GL in gl if ABA==map_gl_to_ba[int(GL)])+ extant_wind_gen[PD+'.'+ABA+'.'+str(H)]
model.windg=Constraint(pds,h,aba, rule=windg)



####solar generation limit
def solarg(model,PD,H,ABA):
    ind=pds.index(PD)
    return model.solarout[PD,H,ABA]==sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL])*solarcf[str(H)+'.'+str(GL)] for PDD in pds[:ind+1] for GL in gl if ABA==map_gl_to_ba[int(GL)])+ extant_solar_gen[PD+'.'+ABA+'.'+str(H)]
model.solarg=Constraint(pds,h,aba, rule=solarg)

# if maxcurtailment!=1:
#     ###wind generation limit
#     def minwindg(model,PD,H,ABA):
#         ind=pds.index(PD)
#         return model.windout[PD,H,ABA]>=(1-maxcurtailment)*(sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL])*windcf[str(H)+'.'+str(GL)] for PDD in pds[:ind+1] for GL in gl if ABA==map_gl_to_ba[int(GL)])+ extant_wind_gen[PD+'.'+ABA+'.'+str(H)])
#     model.minwindg=Constraint(pds,h,aba, rule=minwindg)
    
#     def minsolarg(model,PD,H,ABA):
#         ind=pds.index(PD)
#         return model.solarout[PD,H,ABA]>=(1-maxcurtailment)*(sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL])*solarcf[str(H)+'.'+str(GL)] for PDD in pds[:ind+1] for GL in gl if ABA==map_gl_to_ba[int(GL)])+ extant_solar_gen[PD+'.'+ABA+'.'+str(H)])
#     model.minsolarg=Constraint(pds,h,aba, rule=minsolarg)

def wind_recon(model,PD,GL):
    ind=pds.index(PD)
    return model.capacity_wind_recon[PD,GL]<=wind_recon_capacity[PD+'.'+GL]+sum(wind_recon_capacity[PDD+'.'+GL]-model.capacity_wind_recon[PDD,GL] for PDD in pds[:ind])
model.wind_recon=Constraint(pds,gl, rule=wind_recon)

def solar_recon(model,PD,GL):
    ind=pds.index(PD)
    return model.capacity_solar_recon[PD,GL]<=solar_recon_capacity[PD+'.'+GL]+sum(solar_recon_capacity[PDD+'.'+GL]-model.capacity_solar_recon[PDD,GL] for PDD in pds[:ind])
model.solar_recon=Constraint(pds,gl, rule=solar_recon)

##provincial supply and demand balance
if autrarky:
    def autrarky(model,PD,AP):
        
        TP_supply=sum(model.supply[PD,H,ABA,TP] for H in h for ABA in aba for TP in tplants if AP in ABA)
        wind_solar_supply=sum(model.windout[PD,H,ABA] for H in h for ABA in aba if AP in ABA)\
            +sum(model.solarout[PD,H,ABA] for H in h for ABA in aba if AP in ABA)
        hydro_supply=sum(ror_hydroout[PD+'.'+str(H)+'.'+ABA] for H in h for ABA in aba if AP in ABA)\
            +sum(model.daystoragehydroout[PD,H,ABA] for H in h for ABA in aba if AP in ABA)\
            +sum(model.monthstoragehydroout[PD,H,ABA] for H in h for ABA in aba if AP in ABA)
        return TP_supply+wind_solar_supply+hydro_supply >=autonomy_pct[AP]*sum(demand[PD+'.'+ABA+'.'+str(H)] for ABA in aba for H in h if AP in ABA)
    model.autrarky=Constraint(pds,ap, rule=autrarky)

###supply and demand balance
aux=[1]
def demsup(model,PD,H,ABA):
    TP_supply=sum(model.supply[PD,H,ABA,TP] for TP in tplants)
    wind_solar_supply=model.windout[PD,H,ABA]+model.solarout[PD,H,ABA]
    hydro_supply=ror_hydroout[PD+'.'+str(H)+'.'+ABA]+model.daystoragehydroout[PD,H,ABA]+model.monthstoragehydroout[PD,H,ABA]
    storage_supply=sum(model.storageout[PD,ST,H,ABA]-model.storagein[PD,ST,H,ABA] for ST in st)
    renewal_supply=0
    if hydro_development:
        renewal_supply=sum(ror_renewalout[str(H)+'.'+HR_ROR]*model.ror_renewal_binary[PD,HR_ROR] for HR_ROR in hr_ror if ABA==hr_ror_location[HR_ROR])\
            +sum(model.dayrenewalout[PD,H,HR_DAY] for HR_DAY in hr_day if ABA==hr_day_location[HR_DAY])\
            +sum(model.monthrenewalout[PD,H,HR_MO] for HR_MO in hr_mo if ABA==hr_month_location[HR_MO])
    return TP_supply+wind_solar_supply+hydro_supply+storage_supply+renewal_supply>=demand[PD+'.'+ABA+'.'+str(H)]+sum(demand_us[ABA+'.'+str(H)] for i in aux if ABA+'.'+str(H) in demand_us)\
                                                                        +sum(model.transmission[PD,H,ABA,ABBA]-(1-transloss[ABBA+'.'+ABA])*model.transmission[PD,H,ABBA,ABA] for ABBA in aba if ABA+'.'+ABBA in transmap)
model.demsup=Constraint(pds,h,aba, rule=demsup)

###maximum annual capacity factor for thermal plants
def maxcapfactor(model,PD,ABA,TP):
    ind=pds.index(PD)
    return sum(model.supply[PD,H,ABA,TP] for H in h)<=(sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1])+extant_thermal[pds[0]+'.'+ABA+'.'+TP])*hours * max_cap_fact[TP]
model.maxcapfactor=Constraint(pds,aba,tplants, rule=maxcapfactor)

###minimum annual capacity factor for thermal plants

def mincapfactor(model,PD,ABA,TP):
    ind=pds.index(PD)
    return sum(model.supply[PD,H,ABA,TP] for H in h)>=(sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1])+extant_thermal[pds[0]+'.'+ABA+'.'+TP])*hours * min_cap_fact[TP]
model.mincapfactor=Constraint(pds,aba,tplants, rule=mincapfactor)

##transmission capacity constraint
def transcap(model,PD,H,ABA,ABBA):
    ind=pds.index(PD)
    return model.transmission[PD,H,ABA,ABBA]+model.transmission[PD,H,ABBA,ABA]<=sum(model.capacity_transmission[PDD,ABA,ABBA]+model.capacity_transmission[PDD,ABBA,ABA] for PDD in pds[:ind+1] for i in aux if ABA+'.'+ABBA in transmap)+sum(extant_transmission[pds.index(PD)][ABA+'.'+ABBA] for i in aux if ABA+'.'+ABBA in extant_transmission[pds.index(PD)])
model.transcap=Constraint(pds,h,aba,abba, rule=transcap)

if CTE_extant:
    def TEcap(model,ABA,ABBA):
        return sum(model.capacity_transmission[PD,ABA,ABBA]+model.capacity_transmission[PD,ABBA,ABA] for PD in pds)<=0+CTE_coef*sum(extant_transmission[0][ABA+'.'+ABBA] for i in aux if ABA+'.'+ABBA in extant_transmission[0])
    model.TEcap=Constraint(aba,abba, rule=TEcap)
elif CTE_custom:
    def TEcap(model,ABA,ABBA):
        
        if ABA+'.'+ABBA in trans_expansion_limits:
            TEL=trans_expansion_limits[ABA+'.'+ABBA]
        elif ABBA+'.'+ABA in trans_expansion_limits:
            TEL=trans_expansion_limits[ABBA+'.'+ABA]
        else:
            TEL=100000
                
        return sum(model.capacity_transmission[PD,ABA,ABBA]+model.capacity_transmission[PD,ABBA,ABA] for PD in pds)<=0+TEL
    
    model.TEcap=Constraint(aba,abba, rule=TEcap)
    


#####capacity constraints for thermal plants
def cap(model,PD,H,ABA,TP):
    ind=pds.index(PD)
    return model.supply[PD,H,ABA,TP]<=sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1])+extant_thermal[pds[0]+'.'+ABA+'.'+TP]
model.cap=Constraint(pds,h,aba,tplants, rule=cap)

### this constraint limits the PHS retrofit capacity to percentage of available hydro reservior facility in each BA
if storage_continous:
    def pumpretrofitlimit(model,PD,ABA):
        ind=pds.index(PD)
        return sum(model.capacity_storage[PDD,'PHS',ABA] for PDD in pds[:ind+1])<=(day_hydro_capacity[pds[0]+'.'+ABA]+month_hydro_capacity[pds[0]+'.'+ABA])*pump_ret_limit
    model.pumpretrofitlimit=Constraint(pds,aba, rule=pumpretrofitlimit)
    


#####pumped hydro energy storage
def pumpen(model,PD,ST,H,ABA):
    return model.storageenergy[PD,ST,H+time_diff[H],ABA]==model.storageenergy[PD,ST,H,ABA]-model.storageout[PD,ST,H,ABA]+model.storagein[PD,ST,H,ABA]*storage_efficiency[ST]
model.pumpen=Constraint(pds,st,h2,aba, rule=pumpen)

######pumped hydro energy storage
def pumcap(model,PD,ST,H,ABA):
    ind=pds.index(PD)
    pump_new_con=0
    pump_integer_cap=0
    
    if storage_continous:
        if ST=='LB' and int(PD)>=2040:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[ind-2:ind+1])
        else:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[:ind+1]) 
    if hydro_development and not storage_continous:
        pump_integer_cap=sum(model.pumphydro[PDD,HR_PUMP]*capacity_pump_renewal[HR_PUMP] for PDD in pds[:ind+1] for HR_PUMP in hr_pump if hr_pump_location[HR_PUMP]==ABA)
    return model.storageenergy[PD,ST,H,ABA]<=(ba_storage_capacity[ABA+'.'+ST]+pump_integer_cap+pump_new_con)*storage_hours[ST]
model.pumcap=Constraint(pds,st,h,aba, rule=pumcap)

######pump hydro power capacity
def storageoutmax(model,PD,ST,H,ABA):
    ind=pds.index(PD)
    pump_new_con=0
    pump_integer_cap=0
    
    if storage_continous:
        if ST=='LB' and int(PD)>=2040:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[ind-2:ind+1])
        else:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[:ind+1])            
    if hydro_development and not storage_continous:
        pump_integer_cap=sum(model.pumphydro[PDD,HR_PUMP]*capacity_pump_renewal[HR_PUMP] for PDD in pds[:ind+1] for HR_PUMP in hr_pump if hr_pump_location[HR_PUMP]==ABA)

    
    return model.storageout[PD,ST,H,ABA]<=ba_storage_capacity[ABA+'.'+ST]+pump_integer_cap+pump_new_con

model.storageoutmax=Constraint(pds,st,h,aba, rule=storageoutmax)

#####pump hydro storageing capacity
def storageinmax (model,PD,ST,H,ABA):
    ind=pds.index(PD)
    
    pump_new_con=0
    pump_integer_cap=0
    
    if storage_continous:
        if ST=='LB' and int(PD)>=2040:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[ind-2:ind+1])
        else:
            pump_new_con=sum(model.capacity_storage[PDD,ST,ABA] for PDD in pds[:ind+1]) 
    if hydro_development and not storage_continous:
        pump_integer_cap=sum(model.pumphydro[PDD,HR_PUMP]*capacity_pump_renewal[HR_PUMP] for PDD in pds[:ind+1] for HR_PUMP in hr_pump if hr_pump_location[HR_PUMP]==ABA)

    
    return model.storagein[PD,ST,H,ABA]*storage_efficiency[ST]<=ba_storage_capacity[ABA+'.'+ST]+pump_integer_cap+pump_new_con
model.storageinmax=Constraint(pds,st,h,aba, rule=storageinmax)

##### forces the first hour of months energy to zero
m_counter=2
fhom=list()
for H in h:
    if map_hm[H]==m_counter:
        fhom.append(H)
        m_counter+=1
        
def storage_rest_monthly (model,PD,ST,FHOM,ABA):
    return model.storageenergy[PD,ST,FHOM,ABA]==0
model.storage_rest_monthly=Constraint(pds,st,fhom,aba, rule=storage_rest_monthly)
#####hydro storage for systems with intra-day storage
def hydro_daystorage (model,PD,D,ABA):
    return sum(model.daystoragehydroout[PD,H,ABA] for H in h if map_hd[H]==D)<=day_hydro_historic[PD+'.'+str(D)+'.'+ABA]
model.hydro_daystorage=Constraint(pds,d,aba, rule=hydro_daystorage)

if hydro_development:
    def hydro_dayrenewal (model,PD,D,HR_DAY):
        return sum(model.dayrenewalout[PD,H,HR_DAY] for H in h if map_hd[H]==D)<=day_renewal_historic[str(D)+'.'+HR_DAY]
    model.hydro_dayrenewal=Constraint(pds,d,hr_day, rule=hydro_dayrenewal)

###hydro storage for systems with intra-month storage
def hydro_monthstorage (model,PD,M,ABA):
    return sum(model.monthstoragehydroout[PD,H,ABA] for H in h if map_hm[H]==M)<=month_hydro_historic[PD+'.'+str(M)+'.'+ABA]
model.hydro_monthstorage=Constraint(pds,m,aba, rule=hydro_monthstorage)

if hydro_development:
    def hydro_monthrenewal (model,PD,M,HR_MO):
        return sum(model.monthrenewalout[PD,H,HR_MO] for H in h if map_hm[H]==M)<=month_renewal_historic[str(M)+'.'+HR_MO]
    model.hydro_monthrenewal=Constraint(pds,m,hr_mo, rule=hydro_monthrenewal)

####hydro minimum flow constraints for systems with intra-day storage
def hydro_dayminflow (model,PD,H,ABA):
    return model.daystoragehydroout[PD,H,ABA]>=day_minflow[PD+'.'+ABA]
model.hydro_dayminflow=Constraint(pds,h,aba, rule=hydro_dayminflow)

if hydro_development:
    def renewal_dayminflow (model,PD,H,HR_DAY):
        return model.dayrenewalout[PD,H,HR_DAY]>=capacity_day_renewal[HR_DAY]*hydro_minflow*model.day_renewal_binary[PD,HR_DAY]
    model.renewal_dayminflow=Constraint(pds,h,hr_day, rule=renewal_dayminflow)

####hydro minimum flow constraints for systems with intra-month storage
def hydro_monthminflow (model,PD,H,ABA):
    return model.monthstoragehydroout[PD,H,ABA]>=month_minflow[PD+'.'+ABA]
model.hydro_monthminflow=Constraint(pds,h,aba, rule=hydro_monthminflow)

if hydro_development:
    def renewal_monthminflow (model,PD,H,HR_MO):
        return model.monthrenewalout[PD,H,HR_MO]>=capacity_month_renewal[HR_MO]*hydro_minflow*model.month_renewal_binary[PD,HR_MO]
    model.renewal_monthminflow=Constraint(pds,h,hr_mo, rule=renewal_monthminflow)

##hydro capacity constraints for systems with intra-day storage
def hydro_daycap (model,PD,H,ABA):
    return model.daystoragehydroout[PD,H,ABA]<=day_hydro_capacity[PD+'.'+ABA]
model.hydro_daycap=Constraint(pds,h,aba, rule=hydro_daycap)

if hydro_development:
    def renewal_daycap (model,PD,H,HR_DAY):
        return model.dayrenewalout[PD,H,HR_DAY]<=capacity_day_renewal[HR_DAY]*model.day_renewal_binary[PD,HR_DAY]
    model.renewal_daycap=Constraint(pds,h,hr_day, rule=renewal_daycap)

#####hydro capacity constraints for systems with intra-month storage
def hydro_monthcap (model,PD,H,ABA):
    return model.monthstoragehydroout[PD,H,ABA]<=month_hydro_capacity[PD+'.'+ABA]
model.hydro_monthcap=Constraint(pds,h,aba, rule=hydro_monthcap)

if hydro_development:
    def renewal_monthcap (model,PD,H,HR_MO):
        return model.monthrenewalout[PD,H,HR_MO]<=capacity_month_renewal[HR_MO]*model.month_renewal_binary[PD,HR_MO]
    model.renewal_monthcap=Constraint(pds,h,hr_mo, rule=renewal_monthcap)


##### The following constraints ensure that the model does not build a hydro renewal or greenfield project more that one time during all periods (pds)
if hydro_development:
    def ror_onetime (model,HR_ROR):
        return sum(model.ror_renewal_binary[PD,HR_ROR] for PD in pds)<=1
    model.ror_onetime=Constraint(hr_ror, rule=ror_onetime)
    
    def day_onetime (model,HR_DAY):
        return sum(model.day_renewal_binary[PD,HR_DAY] for PD in pds)<=1
    model.day_onetime=Constraint(hr_day, rule=day_onetime)
    
    def month_onetime (model,HR_MO):
        return sum(model.month_renewal_binary[PD,HR_MO] for PD in pds)<=1
    model.month_onetime=Constraint(hr_mo, rule=month_onetime)
#up ramp limit
def ramp_up (model,PD,H,ABA,TP):
    ind=pds.index(PD)
    return model.supply[PD,H+time_diff[H],ABA,TP]<=model.supply[PD,H,ABA,TP]+(sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1])+extant_thermal[pds[0]+'.'+ABA+'.'+TP])*ramp_rate_percent[TP]*time_diff[H]
model.ramp_up=Constraint(pds,h2,aba,tplants, rule=ramp_up)

#down ramp limit
def ramp_down (model,PD,H,ABA,TP):
    ind=pds.index(PD)
    return model.supply[PD,H+time_diff[H],ABA,TP]>=model.supply[PD,H,ABA,TP]-(sum(model.capacity_therm[PDD,ABA,TP]-model.retire_therm[PDD,ABA,TP] for PDD in pds[:ind+1])+extant_thermal[pds[0]+'.'+ABA+'.'+TP])*ramp_rate_percent[TP]*time_diff[H]
model.ramp_down=Constraint(pds,h2,aba,tplants, rule=ramp_down)

#capacity limit for wind plants
def windcaplimit (model,GL):
    return sum(model.capacity_wind[PD,GL]+model.capacity_wind_recon[PD,GL] for PD in pds) <= maxwind[GL]
model.windcaplimit=Constraint(gl, rule=windcaplimit)

#capacity limit for solar plants
def solarcaplimit (model,GL):
    return sum(model.capacity_solar[PD,GL]+model.capacity_solar_recon[PD,GL] for PD in pds) <= maxsolar[GL]
model.solarcaplimit=Constraint(gl, rule=solarcaplimit)

###################carbon limit constraint can be on or off##################
if provincial_emission_limit:
    def provincialcarbonlimit(model,AP):
        return sum(model.supply[pds[-1],H,ABA,TP]*carbondioxide[TP]/1000000 for H in h for TP in tplants for ABA in aba if AP in ABA)<=carbon_limit[AP]/(365/len(rundays))
    model.provincialcarbonlimit=Constraint(ap, rule=provincialcarbonlimit)
    
if national_emission_limit:
    def nationalcarbonlimit(model,PD):
        return sum(model.supply[PD,H,ABA,TP]*carbondioxide[TP]/1000000 for H in h for TP in tplants for ABA in aba)<=nat_em_limit[PD]/(365/len(rundays))
    model.nationalcarbonlimit=Constraint(pds, rule=nationalcarbonlimit)

if min_installed_LB_PHP:
    def installedgasPHP(model,ABA):
        return sum(model.capacity_storage[PD,ST,ABA] for PD in pds for ST in st)+ba_storage_capacity[ABA+'.PHS']>=LB_PHP_installed_limit[ABA]
    model.installedgasPHP=Constraint(aba, rule=installedgasPHP)

if new_thermal_limit:
    def newthermallimit(model,PD,LTP):
        return sum(model.capacity_therm[PD,ABA,LTP] for ABA in aba)<=0
    model.newthermallimit=Constraint(pds,limited_tplants, rule=newthermallimit)


if thermal_phase_out:
    
    def phaseout(model,PHOT,ABA):
        ind=pds.index(phase_out_type_year[PHOT])
        return extant_thermal[pds[0]+'.'+ABA+'.'+PHOT]+sum(model.capacity_therm[PDD,ABA,PHOT]-model.retire_therm[PDD,ABA,PHOT] for PDD in pds[:ind+1]) <= 0
    model.phaseout=Constraint(ph_out_t,aba, rule=phaseout)    

###New Constraint - Maximum biomass capacity (By Balancing Area)

if max_bio_aba:

    def maxbioaba(model, PD, ABA):
        ind=pds.index(PD)
        return extant_thermal[pds[0]+'.'+ABA+'.'+'biomass']+sum(model.capacity_therm[PDD,ABA,'biomass']-model.retire_therm[PDD,ABA,'biomass'] for PDD in pds[:ind+1]) <= max_bio[ABA]
    model.maxbioaba=Constraint(pds, aba, rule=maxbioaba)

###New Constraint - Maximum nuclear capacity (By Balancing Area)

if max_nuclear_aba:

    def maxnuclearaba(model, PD, ABA):
        ind=pds.index(PD)
        return extant_thermal[pds[0]+'.'+ABA+'.'+'nuclear']+sum(model.capacity_therm[PDD,ABA,'nuclear']-model.retire_therm[PDD,ABA,'nuclear'] for PDD in pds[:ind+1]) <= max_nuclear[ABA]
    model.maxnuclearaba=Constraint(pds, aba,rule=maxnuclearaba)
    
###New Constraint - Maximum SMR capacity (By Balancing Area)

if max_smr_aba:

    def maxsmraba(model, PD, ABA):
        ind=pds.index(PD)
        return extant_thermal[pds[0]+'.'+ABA+'.'+'SMR']+sum(model.capacity_therm[PDD,ABA,'SMR']-model.retire_therm[PDD,ABA,'SMR'] for PDD in pds[:ind+1]) <= max_smr[ABA]
    model.maxsmraba=Constraint(pds, aba,rule=maxsmraba)

if non_emitting_limit:

    def nonemittinglimit(model,PD):
        return sum(model.supply[PD,H,ABA,NEP] for H in h for ABA in aba for NEP in non_emitting_tplants)+sum(model.windout[PD,H,ABA]+model.solarout[PD,H,ABA]+ror_hydroout[PD+'.'+str(H)+'.'+ABA]+model.daystoragehydroout[PD,H,ABA]+model.monthstoragehydroout[PD,H,ABA] for H in h for ABA in aba)>=\
                                                                                                                                                                                                                                    nonemitting_limit[PD]*(sum(model.supply[PD,H,ABA,TP] for H in h for ABA in aba for TP in tplants)+sum(model.windout[PD,H,ABA]+model.solarout[PD,H,ABA]+ror_hydroout[PD+'.'+str(H)+'.'+ABA]+model.daystoragehydroout[PD,H,ABA]+model.monthstoragehydroout[PD,H,ABA] for H in h for ABA in aba))
    model.nonemittinglimit=Constraint(pds, rule=nonemittinglimit) 
    
if limited_new_thermal_gen:
    def limited_new_thermal(model,GT):
        return sum(model.capacity_therm[PD,ABA,GT] for PD in pds for ABA in aba)<=limited_new_thermal_gen[GT]
    model.limited_new_thermal=Constraint(list(limited_new_thermal_gen.keys()),rule=limited_new_thermal)
    
if renewable_portfolio_standards['wind.max']:
    def rpf_wind_max (model,PD,AP):
        ind=pds.index(PD)
        return sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL]) for PDD in pds[:ind+1] for GL in gl if AP==map_gl_to_pr[int(GL)]) + sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'wind'] for GL in gl if str(GL)+'.'+'wind' in extant_wind_solar[0] and AP==map_gl_to_pr[int(GL)])<=wind_max_ap_limit[PD][AP]
    model.rpf_wind_max=Constraint(pds,ap,rule=rpf_wind_max)

if renewable_portfolio_standards['wind.min']:
    def rpf_wind_min (model,PD,AP):
        ind=pds.index(PD)
        return sum((model.capacity_wind[PDD,GL]+model.capacity_wind_recon[PDD,GL]) for PDD in pds[:ind+1] for GL in gl if AP==map_gl_to_pr[int(GL)]) + sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'wind'] for GL in gl if str(GL)+'.'+'wind' in extant_wind_solar[0] and AP==map_gl_to_pr[int(GL)])>=wind_min_ap_limit[PD][AP]
    model.rpf_wind_min=Constraint(pds,ap,rule=rpf_wind_min)

if renewable_portfolio_standards['solar.max']:
    def rpf_solar_max (model,PD,AP):
        ind=pds.index(PD)
        return sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL]) for PDD in pds[:ind+1] for GL in gl if AP==map_gl_to_pr[int(GL)]) + sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'solar'] for GL in gl if str(GL)+'.'+'solar' in extant_wind_solar[0] and AP==map_gl_to_pr[int(GL)]) <=solar_max_ap_limit[PD][AP]
    model.rpf_solar_max=Constraint(pds,ap,rule=rpf_solar_max)

if renewable_portfolio_standards['solar.min']:
    def rpf_solar_min (model,PD,AP):
        ind=pds.index(PD)
        return sum((model.capacity_solar[PDD,GL]+model.capacity_solar_recon[PDD,GL]) for PDD in pds[:ind+1] for GL in gl if AP==map_gl_to_pr[int(GL)]) + sum(extant_wind_solar[pds.index(PD)][str(GL)+'.'+'solar'] for GL in gl if str(GL)+'.'+'solar' in extant_wind_solar[0] and AP==map_gl_to_pr[int(GL)]) >=solar_min_ap_limit[PD][AP]
    model.rpf_solar_min=Constraint(pds,ap,rule=rpf_solar_min)

end=time.time()

print(f'\n==================================================\n\
Creating model time: {round((end-start)/60)} Min and {round((end-start)%60)} Sec \
\n==================================================')

model.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
print('Before Solve Start')
#model.write('problem.lp')
start=time.time()
#import cplex as CPX
#c = CPX.Cplex('problem.lp')
#c.parameters.mip.interval.set(1)
#result_obj=c.solve()
#solving the Optimization Problem
opt = SolverFactory('cplex')
#opt.parameters.mip.interval.set(1)
##optional settings
#opt.options['mipgap'] = 0.05
#opt.options['optimality'] = 0.01
result_obj = opt.solve(model, tee=True)
result_obj.write()
print('Before Result Object Write')
end=time.time()

print(f'\n==================================================\n\
Solving process time: {round((end-start)/60)} Min and {round((end-start)%60)} Sec\
\n==================================================')


############Saving resaults in .csv files################
coordinate = pd.read_excel (r'coordinate.xlsx')

folder_name='outputs'+str(randint(0,10000))+'_ct'+str(ctax['2030'])+'_'+str(ctax['2040'])+'_'+str(ctax['2050'])+'_rd'+str(len(rundays))+'_pds'+str(len(pds))

if test:
   folder_name+='_Test' 
if hierarchical:
    folder_name+='_Hr'
elif downsampling:
    folder_name+='_DS'
    
if not OBPS_on:
    folder_name+='_NoOBPS'
else:
    folder_name+='_OBPS'
if local_gas_price:
    folder_name+='_LGP'    
if not hydro_development:
    folder_name+='_NoHydro'
else:
    folder_name+='_Hydro'
if not national_emission_limit:
    folder_name+='_NoNCL'
else:
    folder_name+='_NCL'
if not provincial_emission_limit:
    folder_name+='_NoPCL'
else:
    folder_name+='_PCL'+str(int(carbon_reduction*100))

if storage_continous:
    folder_name+='_CPHy'
else:
    folder_name+='_NoCPHy'
if not autrarky:
   folder_name+='_NoAr' 
else:
    folder_name+='_Ar'
if SMR_CCS:
   folder_name+='_SMR_CCS' 
else:
    folder_name+='_NoSMR_CCS'   

if CPO:
   folder_name+='_CPO'

if GPS:
   folder_name+='_GPS' 

if technology_evolution_on:
    folder_name+='_TE' 
       
cwd = os.getcwd()
if not os.path.exists(folder_name): 
    os.mkdir(folder_name)
outputdir=cwd+'/'+folder_name
os.chdir(outputdir)



ind=list(model.capacity_therm)
val=list(model.capacity_therm[:,:,:].value)
capacity_thermal = [ i + tuple([j]) for i,j in zip(ind, val)]
#retire_thermal = np.asarray(resultP)
np.savetxt('capacity_thermal.csv', capacity_thermal,fmt='%s',delimiter=',')

ind=list(model.retire_therm)
val=list(model.retire_therm[:,:,:].value)
retire_thermal = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('retire_thermal.csv', retire_thermal,fmt='%s',delimiter=',')

ind=list(model.capacity_wind)
val=list(model.capacity_wind[:,:].value)
capacity_wind = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('capacity_wind.csv', capacity_wind,fmt='%s',delimiter=',')

ind=list(model.capacity_solar)
val=list(model.capacity_solar[:,:].value)
capacity_solar = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('capacity_solar.csv', capacity_solar,fmt='%s',delimiter=',')


ind=list(model.supply)
val=list(model.supply[:,:,:,:].value)
supply = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('supply.csv', supply,fmt='%s',delimiter=',')

ind=list(model.windout)
val=list(model.windout[:,:,:].value)
windout = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('windout.csv', windout,fmt='%s',delimiter=',')

ind=list(model.solarout)
val=list(model.solarout[:,:,:].value)
solarout = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('solarout.csv', solarout,fmt='%s',delimiter=',')

ind=list(model.storageout)
val=list(model.storageout[:,:,:,:].value)
storageout = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('storageout.csv', storageout,fmt='%s',delimiter=',')


ind=list(model.storagein)
val=list(model.storagein[:,:,:,:].value)
storagein = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('storagein.csv', storagein,fmt='%s',delimiter=',')

ind=list(model.storageenergy)
val=list(model.storageenergy[:,:,:,:].value)
storageenergy = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('storageenergy.csv', storageenergy,fmt='%s',delimiter=',')

ind=list(model.daystoragehydroout)
val=list(model.daystoragehydroout[:,:,:].value)
daystoragehydroout = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('daystoragehydroout.csv', daystoragehydroout,fmt='%s',delimiter=',')

ind=list(model.monthstoragehydroout)
val=list(model.monthstoragehydroout[:,:,:].value)
monthstoragehydroout = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('monthstoragehydroout.csv', monthstoragehydroout,fmt='%s',delimiter=',')

ind=list(model.transmission)
val=list(model.transmission[:,:,:,:].value)
transmission = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('transmission.csv', transmission,fmt='%s',delimiter=',')

ind=list(model.capacity_transmission)
val=list(model.capacity_transmission[:,:,:].value)
capacity_transmission = [ i + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('capacity_transmission.csv', capacity_transmission,fmt='%s',delimiter=',')

if hydro_development:
    ind=list(model.ror_renewal_binary)
    val=list(model.ror_renewal_binary[:,:].value)
    ror_renewal_binary = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('ror_renewal_binary.csv', ror_renewal_binary,fmt='%s',delimiter=',')
    
    ind=list(model.day_renewal_binary)
    val=list(model.day_renewal_binary[:,:].value)
    day_renewal_binary = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('day_renewal_binary.csv', day_renewal_binary,fmt='%s',delimiter=',')
    
    ind=list(model.month_renewal_binary)
    val=list(model.month_renewal_binary[:,:].value)
    month_renewal_binary = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('month_renewal_binary.csv', month_renewal_binary,fmt='%s',delimiter=',')
    if not storage_continous:
        ind=list(model.pumphydro)
        val=list(model.pumphydro[:,:].value)
        pumphydro = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
        np.savetxt('pumphydro.csv', pumphydro,fmt='%s',delimiter=',')
    
    ind=list(model.dayrenewalout)
    val=list(model.dayrenewalout[:,:,:].value)
    dayrenewalout = [ i + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('dayrenewalout.csv', dayrenewalout,fmt='%s',delimiter=',')
    
    ind=list(model.monthrenewalout)
    val=list(model.monthrenewalout[:,:,:].value)
    monthrenewalout = [ i + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('monthrenewalout.csv', monthrenewalout,fmt='%s',delimiter=',')

if storage_continous:
    ind=list(model.capacity_storage)
    val=list(model.capacity_storage[:,:,:].value)
    capacity_storage = [ i + tuple([j]) for i,j in zip(ind, val)]
    np.savetxt('capacity_storage.csv', capacity_storage,fmt='%s',delimiter=',')
    
ind=['Objective_function_value']
val=list([model.obj()])
obj = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('obj_value.csv', obj,fmt='%s',delimiter=',')

ind=list(model.capacity_wind_recon)
val=list(model.capacity_wind_recon[:,:].value)
capacity_wind_recon = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('capacity_wind_recon.csv', capacity_wind_recon,fmt='%s',delimiter=',')

ind=list(model.capacity_solar_recon)
val=list(model.capacity_solar_recon[:,:].value)
capacity_solar_recon = [ tuple([i]) + tuple([j]) for i,j in zip(ind, val)]
np.savetxt('capacity_solar_recon.csv', capacity_solar_recon,fmt='%s',delimiter=',')


### Saving config data
number_run_days=len(rundays)
original_stdout = sys.stdout # Save a reference to the original standard output
with open('COPPER config.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(f'Planning for the target year {foryear} considering {refyear} as the reference year')
    print(f'modeled {pds} palnning periods and ran {number_run_days} representative days in each period')
    print(f'Carbon price = {ctax}')
    print(f'reserve margin = {reserve_margin}')
    print(f'pumped hydro retrofit limit = {pump_ret_limit}')
    print(f'down sampling clustering ? {downsampling}')
    print(f'hierarchical clustering ? {hierarchical}')
    print(f'test run ? {test}')
    print(f'hydro development on ? {hydro_development}')
    print(f'autrarky on ? {autrarky}')
    print(f'pump as continous variable ? {storage_continous}')
    print(f'provincial_emission_limit on?  {provincial_emission_limit}   {carbon_reduction} carbon reduction compared to refrence year {emission_limit_ref_year}')
    print(f'national_emission_limit on?  {national_emission_limit}  {nat_em_limit}')
    print(f'local gas price on? {local_gas_price}')
    print(f'OBPS on? {OBPS_on}')
    print(f'SMR and CCS technologies ? {SMR_CCS}')
    print(f'thermal phase out on? {thermal_phase_out} >>>>> if true these types will be phased out by the specified year {phase_out_type_year}')
    print(f'min installed gas PHP requirement on ? {min_installed_LB_PHP}')
    print(f'is developping new thermal banned ? {new_thermal_limit} >>>>> if true development of these types are banned {limited_tplants}')
    print(f'just small hydro on ? {just_small_hydro} >>>>> if true, the model just consider development of small hydro projects (under 100 MW)')
    print(f'DSF tranmission on ? {DSF_tranmission}')
    print(f'is  technology evolution on? {technology_evolution_on}')
    print(f'is  GPS on? {GPS}')
    print(f'is  CPO on? {CPO}')
    print(f'is  tranmission expansion constrained? CTE extant = {CTE_extant} CTE coefficient =>> {CTE_coef}, CTE custom = {CTE_custom}')
    print(f'is  non_emitting limit on? {non_emitting_limit}  {nonemitting_limit}')
    print(f'limited new thermal generation expansion? {limited_new_thermal_gen} ')
    print(f'The discount rate is {discount} and the inflation is {inflation} ')
    sys.stdout = original_stdout # Reset the standard output to its original value


############### Creating the results summary excel sheet that includes generation mix, new generation capacity addition, transmission, carbon, costs, storage, renewable data
#os.chdir('C:/Users\Reza\Desktop\Google Drive\PhD\Thesis\DSF\meeting10\Final results\DSF_BAU_outputs4376_ct170_170_170_rd38_pds6_Hr_NoOBPS_LGP_NoHydro_NCL_NoPCL_CPHy_NoAr_SMR_CCS_CPO_GPS_TE')
capacity_thermal=pd.read_csv(r'capacity_thermal.csv', header=None)

capacity_storage=pd.read_csv(r'capacity_storage.csv', header=None)

retire_thermal =pd.read_csv(r'retire_thermal.csv', header=None)

capacity_wind =pd.read_csv(r'capacity_wind.csv', header=None)

capacity_solar = pd.read_csv(r'capacity_solar.csv', header=None)

supply =pd.read_csv(r'supply.csv', header=None)

windout =pd.read_csv(r'windout.csv', header=None)

solarout =pd.read_csv(r'solarout.csv', header=None)

storageout =pd.read_csv(r'storageout.csv', header=None)

storagein =pd.read_csv(r'storagein.csv', header=None)

storageenergy =pd.read_csv(r'storageenergy.csv', header=None)

monthhydroout=pd.read_csv(r'monthstoragehydroout.csv', header=None)

dayhydroout=pd.read_csv(r'monthstoragehydroout.csv', header=None)

transmission = pd.read_csv(r'transmission.csv', header=None)

capacity_transmission =pd.read_csv(r'capacity_transmission.csv', header=None)
if hydro_development:

    ror_renewal_binary =pd.read_csv(r'ror_renewal_binary.csv', header=None)
    
    day_renewal_binary =pd.read_csv(r'day_renewal_binary.csv', header=None)
    
    month_renewal_binary = pd.read_csv(r'month_renewal_binary.csv', header=None)
    
    dayrenewalout = pd.read_csv(r'dayrenewalout.csv', header=None)
    
    monthrenewalout = pd.read_csv(r'monthrenewalout.csv', header=None)

obj = pd.read_csv(r'obj_value.csv', header=None)


capacity_wind_recon = pd.read_csv(r'capacity_wind_recon.csv',header=None)
capacity_solar_recon =pd.read_csv(r'capacity_solar_recon.csv', header=None)
    
    
######### Generation outline ##########
gen_stor_types=allplants+st
tp_num=len(tplants)
## generation mix in each period
Canada_gen_outline=np.zeros((len(pds),len(gen_stor_types)))
capcitytherm=list(capacity_thermal.iloc[:,3])
retiretherm=list(retire_thermal.iloc[:,3])
## new capacity addition in each period
Total_installed=np.zeros((len(pds),len(gen_stor_types)))
Total_retired=np.zeros((len(pds),len(tplants)))
Total_installed_hydro_aba=dict()
## provincial (BA) breakdown
ABA_generation_mix=np.zeros((len(pds)*len(aba),len(gen_stor_types)))
Total_installed_ABA=np.zeros((len(pds)*len(aba),len(gen_stor_types)))
ABA_investment_costs=np.zeros((len(pds)*len(aba),len(gen_stor_types)))
Canada_investment_costs=np.zeros((len(pds),len(gen_stor_types)))
for PD in pds:
    index_pd=pds.index(PD)
    npv_coef=sum( 1/((1+discount)**(int(PDD)-int(PD))) for PDD in pds[index_pd:])
    for ALP in gen_stor_types:
        index_p=gen_stor_types.index(ALP)
        for ABA in aba:
            index_aba=aba.index(ABA)
                            
            if ALP !='wind' and ALP!='solar' and ALP!='hydro' and ALP not in st:
                index_tp=tplants.index(ALP)
                Canada_gen_outline[index_pd,index_p]+=extant_thermal[pds[0]+'.'+ABA+'.'+ALP]+sum(capcitytherm[ii*len(tplants)*len(aba)+index_aba*tp_num+index_tp] -retiretherm[ii*len(tplants)*len(aba)+index_aba*tp_num+index_tp] for ii in range(index_pd+1))
                Total_installed[index_pd,index_p]+=capcitytherm[index_pd*len(tplants)*len(aba)+index_aba*tp_num+index_tp]
                Total_retired[index_pd,index_tp]+=retiretherm[index_pd*len(tplants)*len(aba)+index_aba*tp_num+index_tp]
                ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=extant_thermal[pds[0]+'.'+ABA+'.'+ALP]+sum(capcitytherm[ii*len(tplants)*len(aba)+index_aba*tp_num+index_tp] -retiretherm[ii*len(tplants)*len(aba)+index_aba*tp_num+index_tp] for ii in range(index_pd+1))
                Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=capcitytherm[index_pd*len(tplants)*len(aba)+index_aba*tp_num+index_tp]
                Canada_investment_costs[index_pd,index_p]+=capcitytherm[index_pd*len(tplants)*len(aba)+index_aba*tp_num+index_tp]*capitalcost[ABA+'.'+ALP]*cap_cost_alter*(npv_coef)
                ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+=capcitytherm[index_pd*len(tplants)*len(aba)+index_aba*tp_num+index_tp]*capitalcost[ABA+'.'+ALP]*cap_cost_alter*(npv_coef)
            elif ALP =='hydro':

                Canada_gen_outline[index_pd,index_p]+=ror_hydro_capacity[PD+'.'+ABA]+day_hydro_capacity[PD+'.'+ABA]+month_hydro_capacity[PD+'.'+ABA]
                ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=ror_hydro_capacity[PD+'.'+ABA]+day_hydro_capacity[PD+'.'+ABA]+month_hydro_capacity[PD+'.'+ABA]
                if hydro_development:
                    Total_installed_hydro_aba[PD+'.'+ABA+'.ror']=0
                    Total_installed_hydro_aba[PD+'.'+ABA+'.day']=0
                    Total_installed_hydro_aba[PD+'.'+ABA+'.month']=0
    
                        
                    for HR_ROR in hr_ror:
                        if ABA==location_renewal[HR_ROR]:
                            index_rn=hr_ror.index(HR_ROR)
                            Canada_gen_outline[index_pd,index_p]+=sum(ror_renewal_binary.iloc[index_rn+ii*len(hr_ror)][2]*capacity_renewal[HR_ROR] for ii in range(index_pd+1))
                            Total_installed[index_pd,index_p]+=ror_renewal_binary.iloc[index_rn+index_pd*len(hr_ror)][2]*capacity_renewal[HR_ROR]
                            Total_installed_hydro_aba[PD+'.'+ABA+'.ror']+=ror_renewal_binary.iloc[index_rn+index_pd*len(hr_ror)][2]*capacity_renewal[HR_ROR]
                            ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=sum(ror_renewal_binary.iloc[index_rn+ii*len(hr_ror)][2]*capacity_renewal[HR_ROR] for ii in range(index_pd+1))
                            Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=ror_renewal_binary.iloc[index_rn+index_pd*len(hr_ror)][2]*capacity_renewal[HR_ROR]
                            Canada_investment_costs[index_pd,index_p]+=ror_renewal_binary.iloc[index_rn+index_pd*len(hr_ror)][2]*cost_ror_renewal[HR_ROR]*cap_cost_alter*(npv_coef)
                            ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+=ror_renewal_binary.iloc[index_rn+index_pd*len(hr_ror)][2]*cost_ror_renewal[HR_ROR]*cap_cost_alter*(npv_coef)

                    for HR_DAY in hr_day:
                        if ABA==location_renewal[HR_DAY]:
                            index_rn=hr_day.index(HR_DAY)
                            Canada_gen_outline[index_pd,index_p]+=sum(day_renewal_binary.iloc[index_rn+ii*len(hr_day)][2]*capacity_renewal[HR_DAY] for ii in range(index_pd+1))
                            Total_installed[index_pd,index_p]+=day_renewal_binary.iloc[index_rn+index_pd*len(hr_day)][2]*capacity_renewal[HR_DAY]
                            Total_installed_hydro_aba[PD+'.'+ABA+'.day']+=day_renewal_binary.iloc[index_rn+index_pd*len(hr_day)][2]*capacity_renewal[HR_DAY]
                            ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=sum(day_renewal_binary.iloc[index_rn+ii*len(hr_day)][2]*capacity_renewal[HR_DAY] for ii in range(index_pd+1))
                            Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=day_renewal_binary.iloc[index_rn+index_pd*len(hr_day)][2]*capacity_renewal[HR_DAY]
                            Canada_investment_costs[index_pd,index_p]+=day_renewal_binary.iloc[index_rn+index_pd*len(hr_day)][2]*cost_day_renewal[HR_DAY]*cap_cost_alter*(npv_coef)
                            ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+=day_renewal_binary.iloc[index_rn+index_pd*len(hr_day)][2]*cost_day_renewal[HR_DAY]*cap_cost_alter*(npv_coef)
    
                    for HR_MO in hr_mo:
                        if ABA==location_renewal[HR_MO]:
                            index_rn=hr_mo.index(HR_MO)
                            Canada_gen_outline[index_pd,index_p]+=sum(month_renewal_binary.iloc[index_rn+ii*len(hr_mo)][2]*capacity_renewal[HR_MO] for ii in range(index_pd+1))
                            Total_installed[index_pd,index_p]+=month_renewal_binary.iloc[index_rn+index_pd*len(hr_mo)][2]*capacity_renewal[HR_MO]
                            Total_installed_hydro_aba[PD+'.'+ABA+'.month']+=month_renewal_binary.iloc[index_rn+index_pd*len(hr_mo)][2]*capacity_renewal[HR_MO]
                            ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=sum(month_renewal_binary.iloc[index_rn+ii*len(hr_mo)][2]*capacity_renewal[HR_MO] for ii in range(index_pd+1))
                            Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=month_renewal_binary.iloc[index_rn+index_pd*len(hr_mo)][2]*capacity_renewal[HR_MO]
                            Canada_investment_costs[index_pd,index_p]+=month_renewal_binary.iloc[index_rn+index_pd*len(hr_mo)][2]*cost_month_renewal[HR_MO]*cap_cost_alter*(npv_coef)
                            ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+=month_renewal_binary.iloc[index_rn+index_pd*len(hr_mo)][2]*cost_month_renewal[HR_MO]*cap_cost_alter*(npv_coef)
                            
                        
            elif ALP=='wind' or ALP=='solar':
                    
                for GL in gl:
                    if map_gl_to_ba[int(GL)]==ABA and str(GL)+'.'+ALP in extant_wind_solar[0]:
                        Canada_gen_outline[index_pd,index_p]+=extant_wind_solar[pds.index(PD)][str(GL)+'.'+ALP]
                        ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=extant_wind_solar[pds.index(PD)][str(GL)+'.'+ALP]               
                    if map_gl_to_ba[int(GL)]==ABA and ALP=='wind':
                        Canada_gen_outline[index_pd,index_p]+=sum(capacity_wind.iloc[len(gl)*ii+int(GL)-1][2]+capacity_wind_recon.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))
                        #Canada_gen_outline[index_pd,index_p]+=sum(capacity_wind.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))+capacity_wind_recon.iloc[len(gl)*0+int(GL)-1][2]
                        Total_installed[index_pd,index_p]+=capacity_wind.iloc[len(gl)*index_pd+int(GL)-1][2]+capacity_wind_recon.iloc[len(gl)*index_pd+int(GL)-1][2]
                        ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=sum(capacity_wind.iloc[len(gl)*ii+int(GL)-1][2]+capacity_wind_recon.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))
                        Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=capacity_wind.iloc[len(gl)*index_pd+int(GL)-1][2]+capacity_wind_recon.iloc[len(gl)*index_pd+int(GL)-1][2]
                        Canada_investment_costs[index_pd,index_p]+=capacity_wind.iloc[len(gl)*index_pd+int(GL)-1][2]*windcost[PD+'.'+GL]*cap_cost_alter*npv_coef+capacity_wind_recon.iloc[len(gl)*index_pd+int(GL)-1][2]*windcost_recon[PD]*cap_cost_alter*npv_coef
                        ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+= capacity_wind.iloc[len(gl)*index_pd+int(GL)-1][2]*windcost[PD+'.'+GL]*cap_cost_alter*npv_coef+capacity_wind_recon.iloc[len(gl)*index_pd+int(GL)-1][2]*windcost_recon[PD]*cap_cost_alter*npv_coef
                    if map_gl_to_ba[int(GL)]==ABA and ALP=='solar':
                        Canada_gen_outline[index_pd,index_p]+=sum(capacity_solar.iloc[len(gl)*ii+int(GL)-1][2]+capacity_solar_recon.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))
                        #Canada_gen_outline[index_pd,index_p]+=sum(capacity_solar.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))+capacity_solar_recon.iloc[len(gl)*0+int(GL)-1][2]
                        Total_installed[index_pd,index_p]+=capacity_solar.iloc[len(gl)*index_pd+int(GL)-1][2]+capacity_solar_recon.iloc[len(gl)*index_pd+int(GL)-1][2]
                        ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=sum(capacity_solar.iloc[len(gl)*ii+int(GL)-1][2]+capacity_solar_recon.iloc[len(gl)*ii+int(GL)-1][2] for ii in range(index_pd+1))
                        Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=capacity_solar.iloc[len(gl)*index_pd+int(GL)-1][2]+capacity_solar_recon.iloc[len(gl)*index_pd+int(GL)-1][2]
                        Canada_investment_costs[index_pd,index_p]+= capacity_solar.iloc[len(gl)*index_pd+int(GL)-1][2]*solarcost[PD+'.'+GL]*cap_cost_alter*npv_coef+capacity_solar_recon.iloc[len(gl)*index_pd+int(GL)-1][2]*solarcost_recon[PD]*cap_cost_alter*npv_coef
                        ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+= capacity_solar.iloc[len(gl)*index_pd+int(GL)-1][2]*solarcost[PD+'.'+GL]*cap_cost_alter*npv_coef+capacity_solar_recon.iloc[len(gl)*index_pd+int(GL)-1][2]*solarcost_recon[PD]*cap_cost_alter*npv_coef
            elif ALP in st:
                if storage_continous:
                    if ST=='LB' and int(PD)>=2040:
                        Canada_gen_outline[index_pd,index_p]+=ba_storage_capacity[ABA+'.'+ALP]+sum(capacity_storage.iloc[ii*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3] for ii in range(index_pd-2,index_pd+1))
                        ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=ba_storage_capacity[ABA+'.'+ALP]+sum(capacity_storage.iloc[ii*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3] for ii in range(index_pd-2,index_pd+1))
                    else:
                        Canada_gen_outline[index_pd,index_p]+=ba_storage_capacity[ABA+'.'+ALP]+sum(capacity_storage.iloc[ii*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3] for ii in range(index_pd+1)) 
                        ABA_generation_mix[len(aba)*index_pd+index_aba,index_p]+=ba_storage_capacity[ABA+'.'+ALP]+sum(capacity_storage.iloc[ii*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3] for ii in range(index_pd+1))
                Total_installed[index_pd,index_p]+=capacity_storage.iloc[index_pd*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3]
                Total_installed_ABA[len(aba)*index_pd+index_aba,index_p]+=capacity_storage.iloc[index_pd*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3]
                Canada_investment_costs[index_pd,index_p]+=capacity_storage.iloc[index_pd*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3]*storage_cost[ST+'.'+PD]*cap_cost_alter*npv_coef
                ABA_investment_costs[len(aba)*index_pd+index_aba,index_p]+=capacity_storage.iloc[index_pd*len(st)*len(aba)+st.index(ALP)*len(aba)+index_aba][3]*storage_cost[ST+'.'+PD]*cap_cost_alter*npv_coef


##### installed transmission ##########3
Installed_transmission=np.zeros((len(pds),len(transmap)))
Transmission_investment=np.zeros((len(pds),len(transmap)))
tr_list=list(capacity_transmission.iloc[:][1])
iter_index=-1  
for PD in pds:
    index_pd=pds.index(PD)
    npv_coef=sum( 1/((1+discount)**(int(PDD)-int(PD))) for PDD in pds[index_pd:])
    for ABA in aba:
        for ABBA in aba:
            if ABA+'.'+ABBA in transmap:
                Installed_transmission[pds.index(PD),transmap.index(ABA+'.'+ABBA)]=capacity_transmission.iloc[pds.index(PD)*len(aba)*len(aba)+aba.index(ABA)*len(aba)+aba.index(ABBA)][3]
                Transmission_investment[pds.index(PD),transmap.index(ABA+'.'+ABBA)]=capacity_transmission.iloc[pds.index(PD)*len(aba)*len(aba)+aba.index(ABA)*len(aba)+aba.index(ABBA)][3]*transcost*distance[ABA+'.'+ABBA]*cap_cost_alter*npv_coef
        
###### Calculating Carbon Emission by BA ######

carbon_ABA=np.zeros((len(pds),len(aba)))
carbon_national_tp=np.zeros((len(pds),len(tplants)))
hours_list=list(supply.iloc[:][0])
ba_list=list(supply.iloc[:][1])
tp_type=list(supply.iloc[:][2])
prod_power=list(supply.iloc[:][3])
carbon_national=dict()
for PD in pds:
    carbon_national[PD]=0
    for H in h:
        for ABA in aba:
            for TP in tplants:
                carbon_ABA[pds.index(PD),aba.index(ABA)]+=supply.iloc[pds.index(PD)*len(h)*len(aba)*len(tplants)+h.index(H)*len(aba)*len(tplants)+aba.index(ABA)*len(tplants)+tplants.index(TP)][4] *carbondioxide[TP]*365/len(run_days)/1000000  
                carbon_national[PD]+=supply.iloc[pds.index(PD)*len(h)*len(aba)*len(tplants)+h.index(H)*len(aba)*len(tplants)+aba.index(ABA)*len(tplants)+tplants.index(TP)][4] *carbondioxide[TP]*365/len(run_days)/1000000
                carbon_national_tp[pds.index(PD),tplants.index(TP)]+=supply.iloc[pds.index(PD)*len(h)*len(aba)*len(tplants)+h.index(H)*len(aba)*len(tplants)+aba.index(ABA)*len(tplants)+tplants.index(TP)][4] *carbondioxide[TP]*365/len(run_days)/1000000
carbon_ap=np.zeros((len(pds),len(ap)))
for AP in ap:
    for PD in pds:
        for ABA in aba:
            if AP in ABA:
                carbon_ap[pds.index(PD),ap.index(AP)]+=carbon_ABA[pds.index(PD),aba.index(ABA)]

Obj=obj.iloc[0][1]*365/len(run_days)/1000000        

#generation by source TWH

j=0
Canada_supply_source=np.zeros((len(pds),len(allplants)))
for i in list(supply[0]):
    pds_index=pds.index(str(i))
    allplants_index=allplants.index(supply[3][j])
    Canada_supply_source[pds_index,allplants_index]+=supply[4][j]
    j+=1
j=0
for i in list(windout[0]):
    pds_index=pds.index(str(i))
    allplants_index=allplants.index('wind')
    Canada_supply_source[pds_index,allplants_index]+=windout[3][j]
    allplants_index=allplants.index('solar')
    Canada_supply_source[pds_index,allplants_index]+=solarout[3][j]
    allplants_index=allplants.index('hydro')
    Canada_supply_source[pds_index,allplants_index]+=monthhydroout[3][j]
    Canada_supply_source[pds_index,allplants_index]+=dayhydroout[3][j]
    Canada_supply_source[pds_index,allplants_index]+=ror_hydroout[str(i)+'.'+str(windout[1][j])+'.'+windout[2][j]]
    j+=1

Canada_supply_source=Canada_supply_source.transpose() 
carbon_ap=carbon_ap.transpose() 
carbon_ABA=carbon_ABA.transpose()
carbon_national_tp=carbon_national_tp.transpose()                            
Canada_gen_outline=Canada_gen_outline.transpose()
Total_installed=Total_installed.transpose()
ABA_generation_mix=ABA_generation_mix.transpose() 
Total_installed_ABA=Total_installed_ABA.transpose()  
Installed_transmission=Installed_transmission.transpose()
ABA_investment_costs=ABA_investment_costs.transpose()
Transmission_investment=Transmission_investment.transpose()  
Canada_investment_costs=Canada_investment_costs.transpose()

can_gen_outline_df= pd.DataFrame(Canada_gen_outline, columns = pds, index=gen_stor_types)
can_gen_outline_df.to_excel("Results_summary.xlsx",sheet_name='Canada_generation_mix') 

Canada_supply_source_df= pd.DataFrame(Canada_supply_source,index=allplants,columns = pds)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Canada_supply_source_df.to_excel(writer, sheet_name='supply_by_source')

Total_installed_df= pd.DataFrame(Total_installed, columns = pds, index=gen_stor_types)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Total_installed_df.to_excel(writer, sheet_name='New_installed_capacity')

ABA_generation_mix_df= pd.DataFrame(ABA_generation_mix, index=gen_stor_types,columns=aba*len(pds))
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    ABA_generation_mix_df.to_excel(writer, sheet_name='ABA_generation_mix')
    
Total_installed_ABA_df= pd.DataFrame(Total_installed_ABA, index=gen_stor_types,columns=aba*len(pds))
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Total_installed_ABA_df.to_excel(writer, sheet_name='New_installed_ABA')
    
   
Installed_transmission_df= pd.DataFrame(Installed_transmission, index=transmap,columns = pds)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Installed_transmission_df.to_excel(writer, sheet_name='New_installed_transmission')

carbon_national_df= pd.DataFrame(list(carbon_national.values()), index=pds,columns=['MT CO2'])
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    carbon_national_df.to_excel(writer, sheet_name='carbon_national')


carbon_ap_df= pd.DataFrame(carbon_ap, index=ap,columns = pds)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    carbon_ap_df.to_excel(writer, sheet_name='carbon_AP')

carbon_ABA_df= pd.DataFrame(carbon_ABA, index=aba,columns = pds)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    carbon_ABA_df.to_excel(writer, sheet_name='carbon_ABA')

carbon_national_tp_df= pd.DataFrame(carbon_national_tp, index=tplants,columns = pds)
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    carbon_national_tp_df.to_excel(writer, sheet_name='carbon_national_tp')

Canada_investment_costs_df= pd.DataFrame(Canada_investment_costs, index=gen_stor_types,columns = pds)    
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Canada_investment_costs_df.to_excel(writer, sheet_name='Canada_investment_costs')
    
ABA_investment_costs_df= pd.DataFrame(ABA_investment_costs, index=gen_stor_types,columns = aba*len(pds))    
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    ABA_investment_costs_df.to_excel(writer, sheet_name='ABA_investment_costs')


Transmission_investment_df= pd.DataFrame(Transmission_investment, index=transmap,columns = pds)    
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    Transmission_investment_df.to_excel(writer, sheet_name='Transmission_investment')

obj_df= pd.DataFrame([Obj],index=['Total cost'],columns = ['M$'])    
with pd.ExcelWriter('Results_summary.xlsx',
                    mode='a',engine="openpyxl") as writer:  
    obj_df.to_excel(writer, sheet_name='OBJ')


### Creates a report when the model is infeasible 
original_stdout = sys.stdout # Save a reference to the original standard output
with open('solver_report.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.    
    #model.load(result_obj) # Loading solution into results object    
    if (result_obj.solver.status == SolverStatus.ok) and (result_obj.solver.termination_condition == TerminationCondition.optimal):
        print('Terminated optimal and noraml ;)\n\n')
    elif (result_obj.solver.termination_condition == TerminationCondition.infeasible):
        print('WARNING:    -------- INFEASIBLE :(  -------\n\n')
        from pyomo.util.infeasible import log_infeasible_constraints
        import logging
        log_infeasible_constraints(model)
        # log_infeasible_constraints(model, log_expression=True, log_variables=True)
        logging.basicConfig(filename='infeasible.log', level=logging.INFO)
    else:
        # Something else is wrong
        print('Solver Status:',  result_obj.solver.status,'\n\n')
    result_obj.write()
    sys.stdout = original_stdout # Reset the standard output to its original value

"""
# Model lp
#my_path = "C:\\Users\\crist\\OneDrive\\Documentos\\UVic\\Courses\\UVIC-2021\\Fall-2021\\Co-op-Work-Term-IESVic\\IESVic-Pyomo-Test"
my_path = "/project/rrg-mcpher16/SHARED2/copper/COPPER6_JGM"
filename = os.path.join(my_path, 'Log-files/pyomo_model_NoCTE_2_5_S1.lp') # LP file of the model instance
model.write(filename, io_options={'symbolic_solver_labels': True})

# Model log
logfile = open("/project/rrg-mcpher16/SHARED2/copper/COPPER6_JGM/Log-files/pyomomodel_NoCTE_2_5_S1.log", "w")
model.pprint(ostream=logfile)
model.solutions.load_from(result_obj)

# Variables log
variable_log = open("/project/rrg-mcpher16/SHARED2/copper/COPPER6_JGM/Log-files/Variables_NoCTE_2_5_S1.log", "w")
for v in model.component_data_objects(Var):
	print(str(v), v.value, file=variable_log)
    
print("Before Duals Log!!!")
    
# Duals log
dual_logfile = open("/project/rrg-mcpher16/SHARED2/copper/COPPER6_JGM/Log-files/Duals_NoCTE_2_5_S1.log", "w") # Create file to print instance duals
print("Model duals:\n", file=dual_logfile)
###Test Consumption
print("Before Duals PPrint!!!")
model.dual.pprint(dual_logfile)
print("After Duals PPrint!!!")
###End Test Consumption
print("\nConstraint slacks:\n", file=dual_logfile)
#print("Constraint || Lower Slack || Upper Slack || Dual Price", file=dual_logfile)
print("Constraint || Lower Slack || Upper Slack", file=dual_logfile)
for c in model.component_objects(Constraint, active=True):
	for index in c:
		print(c[index], " || ", c[index].lslack(), " || ", c[index].uslack(), file=dual_logfile)


try:
    # display all duals
    print("Duals\n", file=dual_logfile)
    for c in model.component_objects(Constraint, active=True):
        print("   Constraint",c, file=dual_logfile)
        for index in c:
            print("      ", index, model.dual[c], file=dual_logfile)
except:
    print("An exception occurred - Error printing the Dual Prices!")
    print("An exception occurred - Error printing the Dual Prices!", file=dual_logfile)

print("After Duals Log!!!")

try:
    # Display the objective
    print("Objective\n", file=dual_logfile)
    for o in model.component_objects(Objective, active=True):
        print("   Objective = ",o, file=dual_logfile)
except:
    print("An exception occurred - Error printing the Objective Value!")
    print("An exception occurred - Error printing the Objective Value!", file=dual_logfile)
"""

os.chdir(cwd)
