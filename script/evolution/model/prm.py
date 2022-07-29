
prm_list = {
    
    "feedback": False, ## ab renewable or not
    "output": True, ### write out info and warn
    "randomAb":True,
    "useWeek":False,
    "useSim":True,
    "debug":False,
    "useBell":False,
    "goodDiff":False,
    "useBinomial":True,
    "dump_eta": False,
    "no_selection": False,
    "dynamical_force": False, ## use dynamical force to extract antigen
    
    "death": "random",
    
    "update_rule": "topK", #### "all, ", "random", "topK"
    
    "potential": "linear-cubic", ## "cusp" or "linear-cubic" or "bell"
    
    "cag":100, ## antigen concentration
    
    "Eb": 14, ## naive B cell affinity
    "Ea": 14, ## tether affinity
    
    "xb1": 1.5,   ## APC-Ag bond stiffness
    "xb2": 2.0,   ## BCR-Ag bond stiffness
    
    ### plasma cell property
    "Td": 0, ### delay betwin diff to Ab feedback
    "Tp": 1000, ## plasma cell lifetime, 0 do not die
    
    "cutoff": 200, ## cycle cutoff to calculate the average maturation rate 
    
    "w0": 8,  ### max number of divisions in each cycle
    "pd": 0.05, ### B cell differentiation rate
    "pa": 0.7,  ### 1 - B cell  death rate
    "pm": 0.5, #### B cell mutation rate
    
    "pm_xb": 0.5, ### probability of bond length mutation
    
    "Nc": 2000,  ### GC carrying capacity
    "Npc": 5000, ### plasma cell capacity
    
    "N0": 1000, ## initial B cell number
    "Nab": 100, ### Ab pool size
    
    "dE": 0.1, ### kT, mutation step
    "dxb": 0.1, ### nm, mutation in bond length
    
    "f":0,  # const force
    "df": 0, ## force increasing step size per GC cycle
    "eta0": 0.5, 
    "tm": 300,
    "r": 0.0001, ### ramping rate when extracting antigens
    
    "Eb_min": 2.0, ### minimal barrier height that analytical expression apply, the higher the more accurate,
    "eta_file": "fixedEb",
}