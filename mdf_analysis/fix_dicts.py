def correct_reversibilities(run,dG_dict,flines,analysis_type): #analysis_type='FBA'/'TFA'
    
    #header_list=flines[0].strip().split('\t')
    #reaction_id_index=header_list.index('Reaction')
    header_list=flines[1].strip().split('\t') #for updated FBA/TFA tables
    reaction_id_index=header_list.index(analysis_type+'_fluxes1') #for updated FBA/TFA tables

    corrected_dG_dict={} #dictionary of corrected dG values for reverse reactions
    flux_dict={} #dictionary of flux values as floats
    flux_index=header_list.index(run)
    for fline in flines[1:]:
        fline_list=fline.split('\t')
        #fline_list=fline.split(' ') #for updated FBA/TFA tables
        reaction_id=fline_list[reaction_id_index]
        flux=fline_list[flux_index].strip()
        if flux and reaction_id in dG_dict.keys() and 'E-' not in flux:
            dG=dG_dict[reaction_id][0]
            if float(flux) < 0 and dG != 0:
                reverse_dG=dG*-1
                dG_error=dG_dict[reaction_id][1]
                corrected_dG_tuple=(reverse_dG,dG_error)
                corrected_dG_dict[reaction_id]=corrected_dG_tuple
                flux_dict[reaction_id]=float(flux)
            elif float(flux) > 0 and dG != 0:
                corrected_dG_dict[reaction_id]=dG_dict[reaction_id]
                flux_dict[reaction_id]=float(flux)

    corrected_dG_dict=dict([(x.replace('NF_',''),corrected_dG_dict[x]) for x in corrected_dG_dict.keys()])
    flux_dict=dict([(x.replace('NF_',''),flux_dict[x]) for x in flux_dict.keys()])
    
    return corrected_dG_dict,flux_dict
