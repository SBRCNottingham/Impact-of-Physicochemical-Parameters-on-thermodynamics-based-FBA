from scipy.optimize import linprog
import numpy as np

R=0.008314 #universal gas constant, KJ/(mol x K)
T=298.15 #temperature, K

ignore_these = ['h_c','h_p','h2o_c']

#eQuilibrator uses the following concentration constraints (in mM):
#lower_bound = 0.001 (based on enzyme affinity)
#upper_bound = 10 (based on toxicity/osmotic pressure/solubility)

#ARGUMENTS:
#m = E. coli core model
#v = flux distribution computed with m
#g = standard Gibbs free energies for reactions in v
#objective_metabolites = metabolites whose concentration should be maximised/minimised
#objective_direction = optimisation direction ('Max' or 'Min')
#metabolite_bounds = a tuple (l,u), with lower and upper metabolite concentration bounds
#upper_dG_limit = all reaction dGs must be lower than this value (default = -1)
def run_mdf_LP(m,v,g,objective_metabolites=False,objective_direction='Min',metabolite_bounds=(1e-6,0.01),upper_dG_limit=0,mdf=True):  

    """
    m = E coli core COBRA model,
    v = flux dict,
    g = standard Gibbs energy dict (keys = reaction IDs, vals = Gibbs energy [J/mol]),
    metabolite_bounds=(l,u)
    """

    if metabolite_bounds:
        c_lb=-np.log(metabolite_bounds[0])
        c_ub=np.log(metabolite_bounds[1])
    else:
        print('WARNING: NO CONCENTRATION BOUNDS SET!')
        c_lb=-np.inf
        c_ub=np.inf
    
    #DEFINE METABOLITES INVOLVED IN SOLUTION
    #FOR E. COLI CORE
    internal_reactions = list(v.keys())
    if 'ATPM' in internal_reactions:
        internal_reactions.remove('ATPM')

    #FOR E. COLI CORE
    metabolites = set()
    for reaction_id in internal_reactions: 
        reaction_metabolites = m.reactions.get_by_id(reaction_id).metabolites
        metabolite_ids = [x.id for x in reaction_metabolites.keys()]
        metabolites = metabolites.union(metabolite_ids)

    #FOR E. COLI CORE
    metabolites = list(metabolites.difference(ignore_these))       

    #OBJECTIVE FUNCTION
    ###MDF STUFF: 
    if mdf:
        metabolites.append('B')
        objective_metabolites=['B']
        objective_direction='Max'
    
    c = np.zeros(len(metabolites)+len(internal_reactions)) 

    if objective_direction == 'Min':
        objective_coefficient = 1
    elif objective_direction == 'Max':
        objective_coefficient = -1
    else:
        print("INVALID OBJECTIVE DIRECTION")
        print("MUST BE 'Max' OR 'Min'")
        return
    
    if objective_metabolites:
        if type(objective_metabolites) == list:
            for obj_met in objective_metabolites: 
                if obj_met in metabolites:
                    metabolite_index = metabolites.index(obj_met)
                    c[metabolite_index] = objective_coefficient
                else:
                    print("CANNOT MAXIMISE/MINIMISE "+obj_met+" IN OBJECTIVE, NOT IN SOLUTION VECTOR.")
                    return
        elif type(objective_metabolites) == dict:
            for obj_met in objective_metabolites.keys():
                if obj_met in metabolites:
                    metabolite_index = metabolites.index(obj_met)
                    c[metabolite_index] = objective_metabolites[obj_met]*objective_coefficient
                    print(obj_met,objective_metabolites[obj_met]*objective_coefficient)
                else:
                    print("CANNOT MAXIMISE/MINIMISE "+obj_met+" IN OBJECTIVE, NOT IN SOLUTION VECTOR.")
                    return
    
    upper_bound_vector=np.full(len(metabolites),c_ub)
    lower_bound_vector=np.full(len(metabolites),c_lb)

    ###MDF STUFF###
    if mdf:
        upper_bound_vector=np.full(len(metabolites)-1,c_ub)
        lower_bound_vector=np.full(len(metabolites)-1,c_lb)

    #DEFINE GIBBS ENERGY OF REACTION INEQUALITY CONSTRAINTS
    if type(v) == dict:
        NT_v = np.empty((len(internal_reactions),len(metabolites))) #NT_v = transpose of stoichiometry matrix for reactions in flux distribution, v.
        b_ub = np.concatenate((np.full(len(internal_reactions)*3,upper_dG_limit),upper_bound_vector,lower_bound_vector)) #right hand side, actual Gibbs energy must be below -1. 
        #b_ub = np.concatenate((np.zeros(len(internal_reactions)*3),np.full(len(metabolites),c_ub),np.full(len(metabolites),c_lb))) #right hand side

        for i in range(len(internal_reactions)):
            reaction_id = internal_reactions[i]
            standard_gibbs_energy = g[reaction_id][0] #must be energy change for relevant direction
            standard_gibbs_energy_error = g[reaction_id][1]
            b_ub[len(internal_reactions)+i] = ((standard_gibbs_energy+standard_gibbs_energy_error)/(R*T)) 
            b_ub[len(internal_reactions)*2+i] = -1*((standard_gibbs_energy-standard_gibbs_energy_error)/(R*T))

            #FOR E. COLI CORE
            metabolite_dict_i = m.reactions.get_by_id(reaction_id).metabolites 
            stoich_dict_i = dict([(x.id,metabolite_dict_i[x]) for x in metabolite_dict_i.keys() if x.id not in ignore_these])
            #print stoich_dict_i
            
            for j in range(len(metabolites)):
                metabolite_id = metabolites[j] 

                ###MDF STUFF###
                if mdf and metabolite_id == 'B':
                    NT_v[i,j] = 1
                    continue
                ###MDF STUFF###
                
                #FOR E. COLI CORE
                if metabolite_id in stoich_dict_i.keys():
                    stoichiometric_coefficient = stoich_dict_i[metabolite_id]
                    flux = v[reaction_id]
                    reaction_direction = flux/abs(flux) #1 if forward, -1 if reverse
                    NT_v[i,j] = reaction_direction*stoichiometric_coefficient
                else:
                    NT_v[i,j] = 0

        #CONSTRUCT FULL LP PROBLEM
        I_n_by_n = np.identity(len(internal_reactions)) #n-by-n identity matrix (n = number of reactions, diagonal 1s are 'stoichiometry' for dG values)
        I_m_by_m = np.identity(len(metabolites))
        Z_2n_by_m = np.zeros((2*len(internal_reactions),len(metabolites)))
        Z_2m_by_n = np.zeros((2*len(metabolites),len(internal_reactions)))

        ###MDF STUFF
        if mdf:
            I_m_by_m = I_m_by_m[0:-1] #i.e. remove last row so that driving force is unbounded by conc constraints
            Z_2m_by_n = Z_2m_by_n[0:-2]
        ###MDF STUFF
        
        row_1 = np.c_[NT_v,I_n_by_n]
        row_2 = np.c_[Z_2n_by_m,np.concatenate((I_n_by_n,-I_n_by_n))]
        row_3 = np.c_[np.concatenate((I_m_by_m,-I_m_by_m)),Z_2m_by_n]
        A_ub = np.concatenate((row_1,row_2,row_3))
    
        #print '# internal reactions:',len(internal_reactions),internal_reactions
        #print '# internal metabolites:',len(metabolites),metabolites
        #return c,A_ub,b_ub #return matrices to check everything's in order.

        #bounds=[(np.log(0.0001),np.log(100))]*len(metabolites)+[(None,None)]*len(internal_reactions)
        bounds=(-np.inf,np.inf)
            
        #SOLVE
        sol = linprog(c,A_ub=A_ub,b_ub=b_ub,bounds=bounds,method='interior-point') #'bounds' for all variables, 'metabolite_bounds' restrict metabolite concentrations. options={"bland":True}
#        sol = linprog(c,A_ub=A_ub,b_ub=b_ub,bounds=bounds) #'bounds' for all variables, 'metabolite_bounds' restrict metabolite concentrations. options={"bland":True}
        
        if sol.success: #sort data into dictionary structure

            sol_dict={}
            
            ###MDF STUFF
            if mdf:
                #print('mdf:',sol.x[len(metabolites)+1])
                #sol_dict['mdf']=sol.x[len(metabolites)+1]
                kjmol=-sol.x[len(metabolites)-1]*R*T/1000.0
                sol_dict['mdf']=kjmol
                print('mdf: %.4f' % kjmol,'kJ/mol')
                print('B:',sol.x[len(metabolites)-1])
	
            #CONVERT LOG METABOLITE CONCENTRATIONS TO mM
            mM_conc = 1000*np.exp(sol.x[:len(metabolites)])
            metabolites_dict = dict(zip(metabolites,mM_conc))

            #GET dG_STANDARD VALUES
            dG_standard = (R*T)*sol.x[len(metabolites):]
            dG_standard_dict = dict(zip(internal_reactions,dG_standard))

            #GET dG VALUES
            b_sol = np.dot(A_ub,sol.x)
            dG = b_sol[:len(internal_reactions)]
            dG_dict = dict(zip(internal_reactions,dG))
            
            sol_dict['metabolites']=metabolites_dict
            sol_dict['dG_standard']=dG_standard_dict
            sol_dict['dG']=dG_dict
            
            return sol_dict
        else:
            print('NO FEASIBLE SOLUTION')
            return sol.success
