# code to generate data for the EDAM assessment
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import choice
import pandas as pd
import random
import pickle as pk
from time import time as timer
# SIRS model

def initialise_states():
    
    global state_of
    global initial_degree_of
    global DoB_of
    global DoI_of
    global time_of_HIV_diagnosis
    global number_of_dens
    
    #HIV_seeds=choice(nodes,int(len(nodes)*HIV_prevalence))
    HCV_seeds=choice(nodes,int(len(nodes)*HCV_prevalence))    
    
    
    state_of={}
    initial_degree_of={}
    DoB_of={}
    DoI_of={}
    time_of_HIV_diagnosis={}
    number_of_dens={}
    
    for node in range(N):
        #state_str=''
        
        state_str=random.choices(['S','A','C','R'],[S_fraction,A_fraction,C_fraction,R_fraction])[0]
        '''
        if node in HIV_seeds:
            if random.random()<chronic_fraction:
                state_str=state_str+'C'
            else:
                state_str=state_str+'A'
        else:
            state_str=state_str+'S'
        '''
        
        if node in HCV_seeds:
            state_str=state_str+'I'
        else:
            state_str=state_str+'S'
             
        state_of[node]=state_str

        initial_degree_of[node]=4
        DoB_of[node]=-params['mean_age']
        DoI_of[node]=0        
        number_of_dens[node]=0

        
def initialise_variables():
    global time
    global old_time
    global total_rate
    global P
    global state_of

    time=0
    old_time=0
    # the live transition probabilities (the main engine of the code)
    P={}
    
    
    # initialise with 0s
    for hiv_state in HIV_status:
        for hcv_state in HCV_status:
            for node in nodes:
                
                P[(node,hiv_state+hcv_state)]=0

            
    for node in nodes:

        # DISEASE PROGRESSION TRANSITIONS
        # use the keys from params to loop over all possible transitions
        for (current_state,trans_state) in progression_rates:
            # if transition is relevant (if the state of "node" is in the "from" part of the transition)
            if state_of[node]==current_state:
                # then some probability of transitioning (to trans state)
                P[(node,trans_state)]=progression_rates[(current_state,trans_state)]

        # INFECTION TRANSITIONS
        # get the current state
        current_state=state_of[node]
        # if its susceptible to HIV infection
        if current_state[0]=='S':
            # get the state that the noe would transition to if infected
            trans_state='A'+state_of[node][1]
            # calulate the total infection pressure
            P[(node,trans_state)]=sum(transmission_rates[(intervention_status[contact],intervention_status[node])][(current_state,trans_state,state_of[contact])] for contact in social_neighbours[node])

        # if its susceptible to HCV infection
        if current_state[1]=='S':
            # get the state that the noe would transition to if infected
            trans_state=state_of[node][0]+'I'
            # calulate the total infection pressure
            P[(node,trans_state)]=sum(transmission_rates[(intervention_status[contact],intervention_status[node])][(current_state,trans_state,state_of[contact])] for contact in social_neighbours[node])
 
    total_rate=sum(P.values())

def initialise_output():
    global state_list
    
    global HIV_population
    global HCV_population
    # keep a record of the state of every node
    
    HIV_population={'S':[],'A':[],'C':[],'R':[]}
    HCV_population={'S':[],'I':[]}
    
    for HIV_status in HIV_population:
        HIV_population[HIV_status].append(len([node for node in nodes if state_of[node][0]==HIV_status]))
    
    for HCV_status in HCV_population:
        HCV_population[HCV_status].append(len([node for node in nodes if state_of[node][1]==HCV_status]))


def update_variables(node,transition):
    
    global P
    global state_of
    global total_rate
    # change the transition rate from whatever it is to 0

    # action the effects the transition has on other nodes
    leave(state_of[node])

    # change the state (move compartment)
    state_of[node]=transition
   
    # consequences of event
    enter(state_of[node])


def update_output():
    # add to the state list for the whole range
    for t in range(int(old_time*time_steps_per_year),min(int(time*time_steps_per_year),end_time*time_steps_per_year)):
        # this can be a non-loop if time has not crossed an integer

        for HIV_status in HIV_population:
            HIV_population[HIV_status].append(len([node for node in nodes if state_of[node][0]==HIV_status]))
        
        for HCV_status in HCV_population:
            HCV_population[HCV_status].append(len([node for node in nodes if state_of[node][1]==HCV_status]))
        


def plot_output(disease):

    global state_list

    if disease=='HIV':
        x=0
        col={'S':'lightgrey','A':'r','C':'purple','R':'c'}
        population=HIV_population
        
    else:
        x=1
        col={'S':'lightgrey','I':'r'}
        population=HCV_population

    fs=12
    plt.figure(figsize=(5,3))
    plt.title(disease)
    for s in population:

        # % presentation
        # percent_pop=[100*population[s][t]/sum(population[i][t] for i in population) for t in range(end_time*time_steps_per_year)]
        # plt.plot(percent_pop,c=col[s],label=s)
        
        # raw numbers presentation
        plt.plot(population[s],c=col[s],label=s)
    
    plt.legend(loc=1,ncol=4,prop={'size':fs})#ncol=2
    
    
    plt.xlim([0,end_time*time_steps_per_year])
    #plt.ylim([0,100])
    plt.ylim([0,1.3*max([max(population[i]) for i in population])])
    #plt.xticks([])
    #plt.yticks([])
    #plt.xlabel('Time',size=fs)
    #plt.ylabel('Percentage of IDUs',size=fs)
    plt.ylabel('Number of IDUs',size=fs)
    plt.xticks([i*10*time_steps_per_year for i in range(5)],['1990','2000','2010','2020','2030'],size=fs)
    #plt.yticks([0,50,100,150],size=fs)
    
    plt.savefig('figs/'+disease+'.png',format='png',dpi=300,bbox_inches='tight')


##########################
def leave(state):
    global total_rate
    
    # set ALL transistions for the node to 0
    for hiv_state in HIV_status:
        for hcv_state in HCV_status:
            if (node,hiv_state+hcv_state) in P:
                total_rate=total_rate-P[(node,hiv_state+hcv_state)]
            P[(node,hiv_state+hcv_state)]=0

    # and the cessation rate to 0
    if (node,'XX') in P:
        total_rate=total_rate-P[(node,'XX')]
    P[(node,'XX')]=0
    
    
    
    # reduce HIV pressure on neighbours
    if state[0]=='A':
        
        # decrease the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:#+sexual_neighbours[node]: # this will need splitting into two loops
            # only if they are susceptible to HIV
            if state_of[neighbour][0]=='S':
                # trans_state is the state the neighbour can transition to
                trans_state='A'+state_of[neighbour][1]
                # decrease by the amount that the node was contributing to the infection pressure of neighbour
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]
                # apply same change to the total rate
                total_rate=total_rate-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]

    # reduce HIV pressure on neighbours
    if state[0]=='C':
        
        # decrease the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:#+sexual_neighbours[node]: # this will need splitting into two loops
            # only if they are susceptible to HIV
            if state_of[neighbour][0]=='S':
                # trans_state is the state the neighbour can transition to
                trans_state='A'+state_of[neighbour][1]
                # 90% reduction as going from acute to chronic   
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]
                # apply same change to the total rate
                total_rate=total_rate-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]

    # reduce HCV pressure on neighbours    
    if state[1]=='I':
        
        # decrease the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:
            # only if they are susceptible to HCV
            if state_of[neighbour][1]=='S':
                # trans_state is the state the neighbour can transition to
                trans_state=state_of[neighbour][0]+'I'
                # this bit of code is same as above. Can it be written more efficiently?
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]
    
                total_rate=total_rate-transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state)]
        
#~~~~~~~~~~~~~~~~~~#
def enter(state):
    global total_rate
    
    # Entering HIV susceptible
    if state[0]=='S':
        
        trans_state='A'+state[1]
        
        infection_pressure=sum(transmission_rates[(intervention_status[contact],intervention_status[node])][(state,trans_state,state_of[contact])] for contact in social_neighbours[node])
        
        P[(node,trans_state)]=infection_pressure
        # update total too       
        total_rate=total_rate+infection_pressure

    
    # Entering HIV acute
    if state[0]=='A':
        # trans state is the state the node could potentially transistion into
        trans_state='C'+state[1]
        # change the recovery rate of the node from 0 to whatever it is 
        
        P[(node,trans_state)]=progression_rates[(state,trans_state)]
    
        total_rate=total_rate+progression_rates[(state,trans_state)]
    
        # increase the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:#+sexual_neighbours[node]:
            # only if they are susceptible to HIV
            if state_of[neighbour][0]=='S':
                
                trans_state='A'+state_of[neighbour][1]
                # increase by the amount that the node contributes to the infection pressure of neighbour
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]+transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]
            
                total_rate = total_rate + transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]
            
    # Entering HIV chronic
    elif state[0]=='C': 
        # trans state is the state the node could potentially transistion into
        trans_state='R'+state[1]
        # change the recovery rate of the node from 0 to whatever it is 
        
        P[(node,trans_state)]=progression_rates[(state,trans_state)]
    
        total_rate=total_rate+progression_rates[(state,trans_state)]
        

        # increase the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:#+sexual_neighbours[node]:
            # only if they are susceptible to HIV
            if state_of[neighbour][0]=='S':
                
                trans_state='A'+state_of[neighbour][1]
                # increase by the amount that the node contributes to the infection pressure of neighbour
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]+transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]

                total_rate = total_rate + transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]
                
    elif state[0]=='R':
        time_of_HIV_diagnosis[node]=time
    #~~~~#
    
    # Entering HCV susceptible
    if state[1]=='S':
        
        trans_state=state[0]+'I'
        
        infection_pressure=sum(transmission_rates[(intervention_status[contact],intervention_status[node])][(state,trans_state,state_of[contact])] for contact in social_neighbours[node])
        
        # change the infection rate of the event node from 0 to whatever it is
        P[(node,trans_state)]=infection_pressure
    
        total_rate=total_rate+infection_pressure
    
    # Entering HCV infected
    elif state[1]=='I':
    
        trans_state=state[0]+'S'
        # change the recovery rate of the node from 0 to whatever it is 
        P[(node,trans_state)]=progression_rates[(state,trans_state)]
    
        total_rate=total_rate+progression_rates[(state,trans_state)]
    
        # increase the infection rate of the neighbours of the event node
        for neighbour in social_neighbours[node]:
            # only if they are susceptible to HCV
            if state_of[neighbour][1]=='S':
                
                trans_state=state_of[neighbour][0]+'I'
                P[(neighbour,trans_state)]=P[(neighbour,trans_state)]+transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]

                total_rate=total_rate+transmission_rates[(intervention_status[node],intervention_status[neighbour])][(state_of[neighbour],trans_state,state_of[node])]

    # whatever state they enter (except the "ceased" state), give them a cessation rate
    if state!='XX':
        # change the cessation rate of the node from 0 to whatever it is 
        P[(node,'XX')]=progression_rates[(state,'XX')]
        total_rate=total_rate+progression_rates[(state,'XX')]

    # leaving the network
    else:
    
        # destroy the node and its links
        nodes.remove(node)
        # remove from other nodes' neighbour lists
        for i in social_neighbours[node]:
            social_neighbours[i].remove(node)
        # remove from dictionary
        social_neighbours.pop(node)
        DoB_of.pop(node)
        DoI_of.pop(node)
        state_of.pop(node)
        initial_degree_of.pop(node)
        
        # remove from dens
        for den in users_of:
            if node in users_of[den]:
                users_of[den].remove(node)
        
        
        # loop over all trasitions relevant to the node 
        for hiv_state in HIV_status:
            for hcv_state in HCV_status:
                total_rate=total_rate-P[(node,hiv_state+hcv_state)]
                P.pop((node,hiv_state+hcv_state))
        
        total_rate=total_rate-P[(node,'XX')]
        P.pop((node,'XX'))

# creates a node and connects it to others in the network
def create_node():
    global N
    global event
    
    node=N+1
    N=N+1
    
    # decide what state it is in when it enters
    
    trans_state=random.choices(['S','A','C','R'],[S_fraction,A_fraction,C_fraction,R_fraction])[0]
  
    if random.random()<HCV_prevalence:
        trans_state=trans_state+'I'
    else:
        trans_state=trans_state+'S'

    # and create links
    initial_degree_of[node]=int(np.random.lognormal(np.log(params['mean_degree']),np.log(params['degree_disp']))) # change to a random variable in future
    
    age_of_initiation=np.random.normal(params['mean_age'],params['age_sd'])
    DoB_of[node]=time-age_of_initiation 
    DoI_of[node]=time
    
    # to speed things up, take a sample of 100 candidate nodes
    #sample_nodes=random.sample(nodes,min(len(nodes),100))

    # weight the sample nodes with probability proportional to initial degree and similarity in age
    #weights=[initial_degree_of[i]*np.exp(-((DoB_of[i]-DoB_of[node])/params['age_divergence'])**2) for i in sample_nodes]
    
    #social_neighbours[node]=list(set(random.choices(sample_nodes,weights,k=min(len(nodes),initial_degree_of[node]))))
 
    # den-based connections
    social_neighbours[node]=[]
    # choose the number of dens
    number_of_dens[node]=np.random.geometric(1/1.05)
    
    # select the specific dens
    dens=random.sample(list(users_of.keys()),number_of_dens[node])

    
    # choose the number in each den
    connections_per_den=int(initial_degree_of[node]/number_of_dens[node])
    # choose the den
    for den in dens:
        # select the people in the den
        social_neighbours[node]=social_neighbours[node]+list(set(random.sample(users_of[den],min(len(users_of[den]),connections_per_den))))
        # add the new node to the den
        users_of[den].append(node)
        

    # the older simpler method
    #social_neighbours[node]=list(set(random.sample(nodes,min(len(nodes),initial_degree_of[node]))))
    
    for neighbour in social_neighbours[node]:

        social_neighbours[neighbour].append(node)
    # add now (not earlier as we don't want self-loops)
    nodes.append(node)
    

    state_of[node]='XX'

    intervention_status[node]='Baseline'

    event=(node,trans_state)

# this function changes all the rates affected by a rate change. "rate" is the state-to-state transition, "new_value" is the new value
def global_rate_change(rate,new_value):
    global rates
    global total_rate
    
    progression_rates[rate]=new_value

    # find transitions that the change applies to
    for trans in P:
        
        node=trans[0]
        
        # if the "from" and "to" states are the same
        if rate==(state_of[node],trans[1]):
            
            old_rate=P[trans]
            
            # if its not a network change then its a simple change of value
                
            new_rate=new_value

            P[trans]=new_rate
            # total rate needs to change too
           
            total_rate=total_rate+new_rate-old_rate

 

#######################################
seed=1
random.seed(seed)
np.random.seed(seed)

N=100

#HIV_prevalence=0.1
S_fraction,A_fraction,C_fraction,R_fraction=[0.95,0.01,0.04,0]
HCV_prevalence=0.05
end_time=40 # duration in years
# separately set the initiation rate (initiations per time unit)
#initiation_rate=200
time_steps_per_year=12 # for polotting
chronic_fraction=0.8
art_fraction=0
total_number_of_dens=10
initial_den_size=10

number_of_interventions=100

users_of={}
for i in range(total_number_of_dens):
    users_of['D'+str(i)]=random.sample(range(N),initial_den_size)


# the rate dictionaries were getting a bit big so I moved them to a separate .py file
from rate_lists import progression_rates, transmission_rates


params={'age_divergence':11,
              'mean_duration_injecting':7,
              'mean_age':25,
              'age_sd':10,
              'mean_degree':4,
              'degree_disp':1.5}



# list the possible states
HIV_status=['S','A','C','R'] # susceptible , acute, chronic, on ART
HCV_status=['S','I'] # susceptible, infectious



# create the initial nodes
nodes=[i for i in range(N)]

intervention_status={}


# initial nides have no neighbours
social_neighbours={}
for node in nodes:
    social_neighbours[node]=[]
    intervention_status[node]='Baseline'

initialise_states()

initialise_variables()

initialise_output()

start=timer()
section_time=0


# choose which event is happening
while time<end_time:

    
    # update the state record
    old_time=time
    
    
    initiation_rate=0.5*len(nodes)*(1-len(nodes)/2000)

    #print(sum(P.values()),total_rate)
    # get time of next event
    time=time+np.random.exponential(1/(total_rate+initiation_rate))
    
    # if it crosses a time point
    art_rates={15:0.05,23:0.2}
    HIV_fractions={15:[0.95,0.01,0.04,0],23:[0.95,0.005,0.015,0.03]}
    for rate_change_time in art_rates:
        if time>rate_change_time and old_time<rate_change_time:
            global_rate_change(('AS','RS'),art_rates[rate_change_time])
            global_rate_change(('AI','RI'),art_rates[rate_change_time])
            global_rate_change(('CS','RS'),art_rates[rate_change_time])
            global_rate_change(('CI','RI'),art_rates[rate_change_time])
            
            # change background rates
            S_fraction,A_fraction,C_fraction,R_fraction=HIV_fractions[rate_change_time]

    # do an intervention to remove some nodes from the network
    # if it crosses a time point
    '''
    for intervention_time in [30+i for i in range(10)]:
        if time>intervention_time and old_time<intervention_time:
            
            # choose a bunch of (random) nodes
            available_nodes=[node for node in nodes if intervention_status[node]=='Baseline']
            intervention_nodes=random.sample(available_nodes,min(number_of_interventions,len(available_nodes)))
            
            # alternatively
            # sort nodes by their number of dens
            sorted_nodes=sorted([(node,number_of_dens[node]) for node in available_nodes], key = lambda item: item[1],reverse=True)
            print(sorted_nodes)
            sorted_nodes=[d[0] for d in sorted_nodes]
            intervention_nodes=sorted_nodes[:number_of_interventions]
            
            
            
            # find the nodes that have transitioned to ART in last year
            newly_diagnosed=[node for node in time_of_HIV_diagnosis if time_of_HIV_diagnosis[node]>intervention_time-1]
            # get the dens
            
            den_list=[(den,len([user for user in users_of[den] if user in newly_diagnosed])) for den in users_of]
            den_list=sorted(den_list, key = lambda item: item[1])
            
            for den in den_list:
                print(den)
            den_list=[d[0] for d in den_list]          
            print(den_list)   
            
            # for the targeted version choose from list of dens
            nodes_by_den=sum([users_of[d] for d in den_list],[])
            
            # intervention nodes (use set to avoid repeats)
            print(len(set(nodes_by_den)))
            
            intervention_nodes=[]
            while nodes_by_den and len(intervention_nodes)<number_of_interventions:
                n=nodes_by_den.pop()
                if n not in intervention_nodes and intervention_status[node]=='Baseline':
                    intervention_nodes.append(n)
            
    
            # for each one
            for node in intervention_nodes:
                
                # remove them
                leave(state_of[node])
                # change the status of the node
                intervention_status[node]='Intervention'
                # put them back in the same  
                enter(state_of[node])
    '''
    
    # choose the event
    # is it a birth? 
    if random.random()<random.choices([True, False],[initiation_rate,total_rate])[0]:
        section_start=timer()  
        # add a new node to the network. This function also creates an event
        create_node()
        
    # if its a transition between states
    else:
        
        section_start=timer()
        # try speeding it up with subsample
        #sample_events=random.sample(list(P.keys()),min(len(P),100))
        sample_events=random.choices(list(P.keys()),k=min(len(P),100))
        
        event=random.choices(sample_events,[P[event] for event in sample_events])[0]
        
        #event=random.choices(list(P.keys()),list(P.values()))[0]
        section_time=section_time+timer()-section_start
    
    
    node=event[0]
    transition=event[1]
   
    #print('At time '+str(round(time,2))+' node '+str(node)+' transitions from '+state_of[node]+' to '+transition)
       
    update_variables(node,transition)
    
    update_output()
        
    


print(timer()-start,section_time)
    
########### PLOTTING ############
print(N)
plot_output('HIV')
plot_output('HCV')



