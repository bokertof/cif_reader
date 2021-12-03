import glob, os, re, datetime, random


def reader(file_name):
    
    '''reads cif-file

    args:
    file_name - name of cif file

    output:
    chemical_formula_sum,cell_volume,cell_length_a,cell_length_b,cell_length_c,cell_angle_alpha,cell_angle_beta,cell_angle_gamma,cell_formula_units_Z,
    chemical_formula_weight,crystal_density_diffrn,absorpt_coefficient_mu,diffrn_radiation_wavelength,diffrn_radiation_type,
    
    '''
    lines = file_name.readlines()
    
    chemical_formula_sum = '?'
    crystal_system = '?'
    group = '?'
    cell_volume = '?'
    cell_length_a = '?'
    cell_length_b = '?'
    cell_length_c = '?'
    cell_angle_alpha = '?'
    cell_angle_beta = '?'
    cell_angle_gamma = '?'
    cell_formula_units_Z = '?'
    chemical_formula_weight = '?'
    crystal_density_diffrn = '?'
    F000 = '?'
    absorpt_coefficient_mu = '?'
    diffrn_radiation_wavelength = '?'
    diffrn_radiation_type = '?'
    diffrn_ambient_temperature = '?'
    diffrn_reflns_number = '?'
    reflns_number_total = '?'
    reflns_number_gt = '?'
    goodness_of_fit = '?'
    peak = '?'
    hole = '?'
    h_max = '?'
    h_min = '?'
    k_max = '?'
    k_min = '?'
    l_max = '?'
    l_min = '?'
    theta_min = '?' 
    theta_max = '?'
    R_int = '?'
    R_factor_gt = '?'
    R_factor_all = '?'
    wR_factor_gt = '?'
    wR_factor_all = '?'
    params = '?'
    restraints = '?'


    
    formula = True
    param0 = 0
    for line in lines:

        if re.findall("\s-*\d\s-*\d\s-*\d\s-*\d*\.*\d*\s*\d", line):
            break

        if formula == 1:
             if not re.findall("_chemical_formula_weight", line):
                 chemical_formula_sum = re.findall(r".*", line)[0].replace("'",'')
                 formula = False
             continue
                 
        if re.findall("_chemical_formula_sum\s*('.*?'){0,1}", line):
            chemical_formula_sum = re.findall("_chemical_formula_sum\s*('.*?'){0,1}", line)[0].replace("'",'').replace(' ','')
            if chemical_formula_sum:
                formula = 0
            else:
                formula = 1
            continue
    

        if re.findall("_chemical_formula_weight", line):
            formula = 2
            chemical_formula_weight = re.findall("_chemical_formula_weight\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_volume", line):
            "STOP if the CIF file has more than one structure"
            param0 += 1
            if param0 > 1:
                break
            
            cell_volume = re.findall("_cell_volume\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_length_a", line): 
            cell_length_a = re.findall("_cell_length_a\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue
            
        if re.findall("_cell_length_b", line): 
            cell_length_b = re.findall("_cell_length_b\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue
            
        if re.findall("_cell_length_c", line): 
            cell_length_c = re.findall("_cell_length_c\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_angle_alpha", line):
            cell_angle_alpha = re.findall("_cell_angle_alpha\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_angle_beta", line): 
            cell_angle_beta = re.findall("_cell_angle_beta\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_angle_gamma", line): 
            cell_angle_gamma = re.findall("_cell_angle_gamma\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_cell_formula_units_Z", line): 
            cell_formula_units_Z = re.findall("_cell_formula_units_Z\s*(\d*)", line)[0]
            continue

        if re.findall("_exptl_crystal_density_diffrn", line): 
            crystal_density_diffrn = re.findall("_exptl_crystal_density_diffrn\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_exptl_crystal_F_000", line): 
            F000 = re.findall("_exptl_crystal_F_000\s*(\d*)", line)[0]
            continue

        if re.findall("_exptl_absorpt_coefficient_mu", line): 
            absorpt_coefficient_mu = re.findall("_exptl_absorpt_coefficient_mu\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_diffrn_radiation_wavelength", line): 
            diffrn_radiation_wavelength = re.findall("_diffrn_radiation_wavelength\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_diffrn_radiation_type", line):
            diffrn_radiation_type = re.findall(r"_diffrn_radiation_type\s*(.*)", line)[0].replace("'",'')
            continue

        if re.findall("_diffrn_ambient_temperature", line): 
            diffrn_ambient_temperature = re.findall("_diffrn_ambient_temperature\s*(\d*\.*\d*\(*\d*\.*\d*\)*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_number", line): 
            diffrn_reflns_number = re.findall("_diffrn_reflns_number\s*(\d*)", line)[0]
            continue

        if re.findall("_reflns_number_total", line): 
            reflns_number_total = re.findall("_reflns_number_total\s*(\d*)", line)[0]
            continue

        if re.findall("_reflns_number_gt", line): 
            reflns_number_gt = re.findall("_reflns_number_gt\s*(\d*)", line)[0]
            continue

        if re.findall("_refine_ls_goodness_of_fit_ref", line): 
            goodness_of_fit = re.findall("_refine_ls_goodness_of_fit_ref\s*(\d*\.*\d*)", line)[0]
            continue

        if re.findall("_refine_diff_density_max", line): 
            peak = re.findall("_refine_diff_density_max\s*(-*\d*\.*\d*)", line)[0]
            continue

        if re.findall("_refine_diff_density_min", line): 
            hole = re.findall("_refine_diff_density_min\s*(-*\d*\.*\d*)", line)[0]
            continue

        if re.findall("_refine_ls_number_parameters", line): 
            params = re.findall("_refine_ls_number_parameters\s*(\d*)", line)[0]
            continue

        if re.findall("_refine_ls_number_restraints", line): 
            restraints = re.findall("_refine_ls_number_restraints\s*(\d*)", line)[0]
            continue
        
        if re.findall("_diffrn_reflns_theta_min", line): 
            theta_min = re.findall("_diffrn_reflns_theta_min\s*(\d*\.*\d*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_theta_max", line): 
            theta_max = re.findall("_diffrn_reflns_theta_max\s*(\d*\.*\d*)", line)[0]
            continue
        
        if re.findall("_symmetry_cell_setting", line):
            if re.findall("_symmetry_cell_setting\s*(.*)", line): 
                crystal_system = re.findall("_symmetry_cell_setting\s*(.*)", line)[0].replace("'",'')        
        else:
            if re.findall("_space_group_crystal_system\s*(.*)", line):
                crystal_system = re.findall("_space_group_crystal_system\s*(.*)", line)[0].replace("'",'')  
                continue

        if re.findall("_space_group_name_H-M_alt", line):
            if re.findall("_space_group_name_H-M_alt\s*(.*)", line): 
                group = re.findall("_space_group_name_H-M_alt\s*(.*)", line)[0].replace("'",'')        
        else:
            if re.findall("_symmetry_space_group_name_H-M\s*(.*)", line):
                group = re.findall("_symmetry_space_group_name_H-M\s*(.*)", line)[0].replace("'",'')  
                continue

        if re.findall("_diffrn_reflns_av_R_equivalents", line): 
            R_int = re.findall("_diffrn_reflns_av_R_equivalents\s*(\d*\.*\d*)", line)[0]
            continue
        
        if re.findall("_diffrn_reflns_limit_h_max", line): 
            h_max = re.findall("_diffrn_reflns_limit_h_max\s*(\d*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_limit_h_min", line): 
            h_min = re.findall("_diffrn_reflns_limit_h_min\s*(-*\d*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_limit_k_max", line): 
            k_max = re.findall("_diffrn_reflns_limit_k_max\s*(\d*)", line)[0]
            continue
        
        if re.findall("_diffrn_reflns_limit_k_min", line): 
            k_min = re.findall("_diffrn_reflns_limit_k_min\s*(-*\d*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_limit_l_max", line): 
            l_max = re.findall("_diffrn_reflns_limit_l_max\s*(\d*)", line)[0]
            continue

        if re.findall("_diffrn_reflns_limit_l_min", line): 
            l_min = re.findall("_diffrn_reflns_limit_l_min\s*(-*\d*)", line)[0]
            continue
        
        if re.findall("_refine_ls_wR_factor_ref", line): 
            wR_factor_all = re.findall("_refine_ls_wR_factor_ref\s*(\d*\.*\d*)", line)[0]
            continue
        
        if re.findall("_refine_ls_wR_factor_gt", line): 
            wR_factor_gt = re.findall("_refine_ls_wR_factor_gt\s*(\d*\.*\d*)", line)[0]
            continue

        if re.findall("_refine_ls_R_factor_gt", line): 
            R_factor_gt = re.findall("_refine_ls_R_factor_gt\s*(\d*\.*\d*)", line)[0]
            continue

        if re.findall("_refine_ls_R_factor_all", line): 
            R_factor_all = re.findall("_refine_ls_R_factor_all\s*(\d*\.*\d*)", line)[0]
            continue
        
        if re.findall("_refine_ls_R_factor_gt", line): 
            R_factor_gt = re.findall("_refine_ls_R_factor_gt\s*(\d*\.*\d*)", line)[0]
            continue
        

            
    peak_hole = str(peak)+' and '+str(hole)
    hkl = str(h_min)+','+str(h_max)+'/'+str(k_min)+','+str(k_max)+'/'+str(l_min)+','+str(l_max)
    params_restr = str(params)+'/'+str(restraints)
    
    return [chemical_formula_sum, chemical_formula_weight, crystal_system, group, cell_length_a, cell_length_b, cell_length_c, cell_angle_alpha, cell_angle_beta,\
           cell_angle_gamma, cell_volume,\
           cell_formula_units_Z, crystal_density_diffrn, absorpt_coefficient_mu, F000, diffrn_radiation_wavelength,\
           diffrn_radiation_type, diffrn_ambient_temperature, hkl, diffrn_reflns_number,\
           reflns_number_total, reflns_number_gt, params_restr, R_int, R_factor_all, wR_factor_all, R_factor_gt, wR_factor_gt, goodness_of_fit, theta_min, theta_max, peak_hole]
           



number = 0

labels = ['File name', 'Molecular formula', 'Molecular weight', 'Crystal system', 'Space group', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'Cell volume',\
           'Z', 'Density', 'Absorption coefficient (mu)', 'F000', 'Radiation wavelength', 'Radiation type', 'Temperature', 'h/k/l (min, max)', 'Reflections measured',\
           'Independent reflections', 'Observed reflections [I > 2s(I)]', 'Parameters/restraints', 'R int', 'R1 (all data)', 'wR2 (all data)', 'R1 [I > 2s(I)]',  'wR2 [I > 2s(I)]','Goodness-of-fit',\
          'Theta min', 'Theta max', 'Largest diff. peak and hole']



outname = str(datetime.datetime.now().time()).replace(':','_') + '.csv'

general_info = '#'*5+' CIF READER v. 2.0. (28.10.2021)  '+'#'*5+'\n'+'In case of an unexpected error: stas_melnikov_1@mail.ru'+2*'\n'


quotes = ['Feci, quod potui, faciant meliora potentes', 'Festina lente', 'Carum est, quod rarum est', 'Clavus clavo pellitur', 'Dictum est factum',
          'Divide et impera', 'Docendo discimus', 'Facile dictu, difficile factu', 'Festina lente!', 'Fiat lux!', 'Hic mortui vivunt, hic muti loquuntur',
          'Hodie mihi, cras tibi', 'Homo doctus in se semper divitias habet', 'Homo quisque fortunae faber', 'Id agas, ut sis felix, non ut videaris',
          'In hoc signo vinces', 'In optima forma', 'In vino veritas', 'Invenit et perfecit', 'Is fecit, cui prodest', 'Latrante uno, latrat statim et alter canis',
          'Legem brevem esse oportet', 'Mens sana in corpore sano', 'Nec sibi, nec alteri', 'Nil adsuetudine majus', 'Terra incognita', 'Una hirundo non facit ver',
          'Usus est optimus magister', 'Ut ameris, amabilis esto', 'Ut salutas, ita salutaberis', 'Ut vivas, igitur vigila', 'Verba volant, scripta manent',
          'Vivat Academia! Vivant professores!', 'Vivere est cogitare', 'Volens nolens']


res = [labels]
for file_name in glob.glob("*.cif"):
    amount = len(glob.glob("*.cif"))
    with open(file_name, 'r') as file:
        
        features = reader(file)
        features.insert(0, file_name)
        
        res.append(features)
        
    number += 1

res = list(map(list, zip(*res)))



with open(outname, 'w') as out:
    out.write(general_info)
    for i in res:
        out.write('; '.join(i)+'\n')

    out.write('\n')
    out.write(random.choice(quotes))
    

if number == 0:
    with open(outname, 'w') as out:
        out.write(general_info)
        out.write('No CIF files were found!')
        out.write('\n')
        out.write('Errare humanum est!')
