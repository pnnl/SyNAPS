import os

def case_files(base_path):
    '''
    This function is used to get the file names based on the case specificed by the user.
    '''
    print('User provided Case Name:',base_path.split(os.path.sep)[-1])
    if base_path.split(os.path.sep)[-1] == 'Case_14bus':
        dc_case_name = os.path.join(base_path, "mtdc_net_14AC.xlsx")
        ac_case_name_base = os.path.join(base_path, "ieee14AC.xlsx")
        ac_case_name_pf = os.path.join(base_path, "ieee14AC_with_dc_mod.xlsx")
        ac_case_name = os.path.join(base_path, "ieee14AC_ss.xlsx")
        eq_data_file = os.path.join(base_path, "pf_info_"+ base_path.split(os.path.sep)[-1]+".xlsx")
        eq_data_file_extended = os.path.join(base_path, "pf_info_extended"+ base_path.split(os.path.sep)[-1] + ".xlsx")

    if base_path.split(os.path.sep)[-1] == 'Case_5bus':
        dc_case_name = os.path.join(base_path, "mtdc_net_5AC.xlsx")
        ac_case_name_base = os.path.join(base_path, "5BusAC.xlsx")
        ac_case_name_pf = os.path.join(base_path, "5BusAC_with_dc_mod.xlsx")
        ac_case_name = os.path.join(base_path, "5BusAC_ss.xlsx")
        eq_data_file = os.path.join(base_path, "pf_info_"+ base_path.split(os.path.sep)[-1]+".xlsx")
        eq_data_file_extended = os.path.join(base_path, "pf_info_extended"+ base_path.split(os.path.sep)[-1] + ".xlsx")

    if base_path.split(os.path.sep)[-1] == 'Case_5bus_v2':
        dc_case_name = os.path.join(base_path, "mtdc_net_5AC.xlsx")
        ac_case_name_base = os.path.join(base_path, "5BusAC.xlsx")
        ac_case_name_pf = os.path.join(base_path, "5BusAC_with_dc_mod.xlsx")
        ac_case_name = os.path.join(base_path, "5BusAC_ss.xlsx")
        eq_data_file = os.path.join(base_path, "pf_info_"+ base_path.split(os.path.sep)[-1]+".xlsx")
        eq_data_file_extended = os.path.join(base_path, "pf_info_extended"+ base_path.split(os.path.sep)[-1] + ".xlsx")
    
    return dc_case_name,ac_case_name_base,ac_case_name_pf,ac_case_name,eq_data_file,eq_data_file_extended
    