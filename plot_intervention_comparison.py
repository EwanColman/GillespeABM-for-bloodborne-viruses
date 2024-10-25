import pickle as pk
import matplotlib.pyplot as plt
########### PLOTTING ############
output_list={}
output_list['Random']=pk.load(open('pickles/model_outputs_random_interventions.p','rb'))
output_list['Targeted']=pk.load(open('pickles/model_outputs_targeted_interventions.p','rb'))
output_list['No intervention']=pk.load(open('pickles/model_outputs_no_interventions.p','rb'))

runs=100
################
# make the plot and set the gridspec
fig = plt.figure(figsize=(13,6))
gs = fig.add_gridspec(2,3)
plt.subplots_adjust(hspace=0,wspace=0.3)



#################

j=0
#intervention='Targeted'
for intervention in ['No intervention','Random','Targeted']:

    # HIV    
    HIV_prevalence_list=[[] for t in range(40*12)]
    for output in output_list[intervention]['HIV']:
        
        for t in range(40*12):
            HIV_prevalence_list[t].append(sum(output[state][t] for state in ['A','C','R']))
 
    HIV_prevalence_median=[]
    HIV_prevalence_upper=[]
    HIV_prevalence_lower=[]
    
    for t in range(40*12):
        #print(HIV_prevalence_list[t])
        
        sorted_outputs=sorted(HIV_prevalence_list[t])
        
        
        HIV_prevalence_median.append(sorted_outputs[int(runs*0.5)])
        HIV_prevalence_upper.append(sorted_outputs[int(runs*0.975)])
        HIV_prevalence_lower.append(sorted_outputs[int(runs*0.025)])
    
    # HCV
    HCV_prevalence_list=[[] for t in range(40*12)]
    for output in output_list[intervention]['HCV']:
        
        for t in range(40*12):
            HCV_prevalence_list[t].append(sum(output[state][t] for state in ['I']))
 
    HCV_prevalence_median=[]
    HCV_prevalence_upper=[]
    HCV_prevalence_lower=[]
    
    for t in range(40*12):
        #print(HIV_prevalence_list[t])
        
        sorted_outputs=sorted(HCV_prevalence_list[t])
        
        
        HCV_prevalence_median.append(sorted_outputs[int(runs*0.5)])
        HCV_prevalence_upper.append(sorted_outputs[int(runs*0.975)])
        HCV_prevalence_lower.append(sorted_outputs[int(runs*0.025)])


    
        
    ax=fig.add_subplot(gs[0,j])
    
    plt.title(intervention)
    plt.text(20,290,'Hepatitis C')
    plt.plot(HCV_prevalence_median)
    plt.fill_between(range(40*12),HCV_prevalence_upper,HCV_prevalence_lower,alpha=0.5)
    
    plt.xlim([0,40*12])
    plt.ylim([0,350])
    plt.xticks([])
    plt.yticks([0,100,200,300])
    #plt.xticks([0,120,240,360,480],[1990,2000,2010,2020,2030])

    plt.ylabel('HCV Positive IDUs')
        
    
    if intervention!='No intervention':
        for i in range(10):
            plt.axvline((30+i)*12,linestyle='-',linewidth=0.5,color='r')
    #if intervention=='Targeted':
    #    plt.text(52,90,'Top 10% most\nconnected\nnodes blocked',color='r',verticalalignment='top')
    #else:
    #    plt.text(52,90,'Random 10% of\nnodes blocked',color='r',verticalalignment='top')
    
    
    ax=fig.add_subplot(gs[1,j])
    if intervention!='No intervention':
        for i in range(10):
            plt.axvline((30+i)*12,linestyle='-',linewidth=0.5,color='r')
    
    
    plt.text(20,290,'HIV')
    plt.plot(HIV_prevalence_median)
    plt.fill_between(range(40*12),HIV_prevalence_upper,HIV_prevalence_lower,alpha=0.5)
    
    plt.xlim([0,40*12])
    plt.ylim([0,350])
    #plt.xticks([])
    plt.xticks([0,120,240,360,480],[1990,2000,2010,2020,2030])
    plt.yticks([0,50,100,150,200,250,300])
    plt.xlabel('Time')
    plt.ylabel('HIV Positive IDUs')  
    j=j+1
    
plt.savefig('figs/interventions.png',format='png',dpi=300,bbox_inches='tight')

