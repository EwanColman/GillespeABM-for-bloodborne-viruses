# all non zero transitions, rates are per_year e.g. 1/rate = expected number of events in one year
progression_rates={
        # HCV reversion 
        ('SI','SS'):0.2, # HCV reversion is no HIV
        ('AI','AS'):0.1, # HCV reversion if HIV acute
        ('CI','CS'):0.06, # HCV reversion if HIV chronic
        ('RI','RS'):0.06, # HCV reversion if HIV ART enrolled
                 
        # HIV progression
        ('AS','CS'):4, # HIV progression
        ('AI','CI'):4, # HIV progression if HCV+
        
        # ART enrolment
        ('AS','RS'):0.0, # ART enrolment (acute and HepC-)
        ('AI','RI'):0.0, # ART enrolment (acute and HepC+)
        ('CS','RS'):0.0, # ART enrolment (chronic and HepC-)
        ('CI','RI'):0.0, # ART enrolment (chronic and HepC+)

        # cessation
        ('SS','XX'):0.15, # susceptible cessation or death
        ('SI','XX'):0.15, # HCV only cessation or death
        
        ('AS','XX'):0.15, # acute HIV only cessation or death
        ('AI','XX'):0.15, # acute HIV and HCV coinfection cessation or death
        
        ('CS','XX'):0.15, # cronic HIV only cessation or death
        ('CI','XX'):0.15, # cronic HIV and HCV coinfection cesation or death
        
        ('RS','XX'):0.15, # on ART only cessation or death
        ('RI','XX'):0.15 # on ART and HCV coinfection cesation or death

        }

transmission_rates={}
    
transmission_rates[('Baseline','Baseline')]={
        # (state of node, transition state, state of neighbour) : infection rate per infector
    
        # HCV infection
        # if HIV negative
        ('SS','SI','SI'):0.06, # from HCV+ contacts that are HIV negative
        ('SS','SI','AI'):0.06, # from HCV+ contacts that are HIV acute
        ('SS','SI','CI'):0.06, # from HCV+ contacts that are HIV chronic
        ('SS','SI','RI'):0.06, # from HCV+ contacts that are ART enrolled

        ('SS','SI','SS'):0, # from HCV- contacts that are HIV negative
        ('SS','SI','AS'):0, # from HCV- contacts that are HIV acute
        ('SS','SI','CS'):0, # from HCV- contacts that are HIV chronic
        ('SS','SI','RS'):0, # from HCV- contacts that are ART enrolled        

        # if HIV acute
        ('AS','AI','SI'):0.3, # from HCV+ contacts that are HIV negative
        ('AS','AI','AI'):0.3, # from HCV+ contacts that are HIV acute
        ('AS','AI','CI'):0.3, # from HCV+ contacts that are HIV chronic
        ('AS','AI','RI'):0.3, # from HCV+ contacts that are ART enrolled

        ('AS','AI','SS'):0, # from HCV- contacts that are HIV negative
        ('AS','AI','AS'):0, # from HCV- contacts that are HIV acute
        ('AS','AI','CS'):0, # from HCV- contacts that are HIV chronic
        ('AS','AI','RS'):0, # from HCV- contacts that are ART enrolled
        
        # if HIV chronic            
        ('CS','CI','SI'):0.3, # from HCV+ contacts that are HIV negative
        ('CS','CI','AI'):0.3, # from HCV+ contacts that are HIV acute
        ('CS','CI','CI'):0.3, # from HCV+ contacts that are HIV chronic
        ('CS','CI','RI'):0.3, # from HCV+ contacts that are ART enrolled
        
        ('CS','CI','SS'):0, # from HCV- contacts that are HIV negative
        ('CS','CI','AS'):0, # from HCV- contacts that are HIV acute
        ('CS','CI','CS'):0, # from HCV- contacts that are HIV chronic
        ('CS','CI','RS'):0, # from HCV- contacts that are ART enrolled
        
        # if on ART
        ('RS','RI','SI'):0.3, # from HCV+ contacts that are HIV negative
        ('RS','RI','AI'):0.3, # from HCV+ contacts that are HIV acute
        ('RS','RI','CI'):0.3, # from HCV+ contacts that are HIV chronic
        ('RS','RI','RI'):0.3, # from HCV+ contacts that are ART enrolled
        
        ('RS','RI','SS'):0, # from HCV- contacts that are HIV negative
        ('RS','RI','AS'):0, # from HCV- contacts that are HIV acute
        ('RS','RI','CS'):0, # from HCV- contacts that are HIV chronic
        ('RS','RI','RS'):0, # from HCV- contacts that are ART enrolled
        
        
        # HIV infection
        
        # if HCV negative
        ('SS','AS','SS'):0, # from HIV susceptible contacts that are HCV negative
        ('SS','AS','SI'):0, # from HIV susceptible contacts that are HCV positive
        
        ('SS','AS','AS'):0.3, # from HIV acute contacts that are HCV negative
        ('SS','AS','AI'):0.3, # from HIV acute contacts that are HCV positive

        ('SS','AS','CS'):0.03, # from HIV chronic contacts that are HCV negative
        ('SS','AS','CI'):0.03, # from HIV chronic contacts that are HCV positive

        ('SS','AS','RS'):0.01, # from ART enrolled contacts that are HCV negative
        ('SS','AS','RI'):0.01, # from ART enrolled contacts that are HCV positive
        
        # if HCV positive
        ('SI','AI','SS'):0, # from HIV susceptible contacts that are HCV negative
        ('SI','AI','SI'):0, # from HIV susceptible contacts that are HCV positive
        
        ('SI','AI','AS'):0.3, # from HIV acute contacts that are HCV negative
        ('SI','AI','AI'):0.3, # from HIV acute contacts that are HCV positive

        ('SI','AI','CS'):0.03, # from HIV chronic contacts that are HCV negative
        ('SI','AI','CI'):0.03, # from HIV chronic contacts that are HCV positive

        ('SI','AI','RS'):0.01, # from ART enrolled contacts that are HCV negative
        ('SI','AI','RI'):0.01, # from ART enrolled contacts that are HCV positive  
    }

transmission_rates[('Baseline','Intervention')]=transmission_rates[('Baseline','Baseline')].copy()

for rate in transmission_rates[('Baseline','Intervention')]:
    transmission_rates[('Baseline','Intervention')][rate]=transmission_rates[('Baseline','Intervention')][rate]*0.1


transmission_rates[('Intervention','Baseline')]=transmission_rates[('Baseline','Baseline')].copy()

for rate in transmission_rates[('Intervention','Baseline')]:
    transmission_rates[('Intervention','Baseline')][rate]=transmission_rates[('Intervention','Baseline')][rate]*0.1


transmission_rates[('Intervention','Intervention')]=transmission_rates[('Baseline','Baseline')].copy()

for rate in transmission_rates[('Baseline','Intervention')]:
    transmission_rates[('Baseline','Intervention')][rate]=0




