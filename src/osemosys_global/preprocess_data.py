"""Pre-process OSeMOSYS data file to reduce matrix generation time

This script pre-processes an OSeMOSYS input data file by adding lines that list 
commodity-technology-mode combinations that data is provided for. Pre-processing 
a data file before starting a model run significantly reduces the time taken 
for matrix generation. 

 Pre-processing consists of the following steps:

1. Reading the ``InputActivityRatio`` and ``OutputActivityRatio`` sections of the 
data file to identify commodity-technology-mode combinations that data has been 
explicitly provided for.
2. Adding a set entry for each commodity that lists all technology-mode combinations 
that are associated with it.  
3. Values from the ``InputActivityRatios`` and ``OutputActivityRatios`` sections are 
added to the sets ``MODExTECHNOLOGYperFUELin`` and ``MODExTECHNOLOGYperFUELout`` respectively.
4. Values from the ``TechnologyToStorage`` and ``TechnologyFromStorage`` sections 
are added to the sets ``MODExTECHNOLOGYperSTORAGEto`` and ``MODExTECHNOLOGYperSTORAGEfrom`` respectively.
5. All values for technology-mode combinations are added to the sets 
``MODEperTECHNOLOGY``.

 In order to start a model run with a pre-processed data file, the following sets 
need to be introduced to its associated OSeMOSYS model file::

    set MODEperTECHNOLOGY{TECHNOLOGY} within MODE_OF_OPERATION;
    set MODExTECHNOLOGYperFUELout{COMMODITY} within MODE_OF_OPERATION cross TECHNOLOGY;
    set MODExTECHNOLOGYperFUELin{COMMODITY} within MODE_OF_OPERATION cross TECHNOLOGY;
    set MODExTECHNOLOGYperSTORAGEto{STORAGE} within MODE_OF_OPERATION cross TECHNOLOGY;
    set MODExTECHNOLOGYperSTORAGEfrom{STORAGE} within MODE_OF_OPERATION cross TECHNOLOGY;

"""

import pandas as pd
import os, sys
from collections import defaultdict


def main(data_format, data_infile, data_outfile):

    lines = []

    with open(data_infile, 'r') as f1:
        for line in f1:
            if not line.startswith(('set MODEper','set MODEx', 'end;')):
                lines.append(line)

    with open(data_outfile, 'w') as f2:
        f2.writelines(lines)

    parsing = False
    parsing_year = False
    parsing_tech = False
    parsing_fuel = False
    parsing_mode = False
    parsing_storage = False
    parsing_emission = False

    year_list = []
    fuel_list = []
    tech_list = []
    storage_list = []
    mode_list = []
    emission_list = []

    data_all = []
    data_out = []
    data_inp = []
    output_table = []
    storage_to = []
    storage_from = []
    emission_table = []

    param_check = ['OutputActivityRatio', 'InputActivityRatio', 'TechnologyToStorage', 'TechnologyFromStorage', 'EmissionActivityRatio']

    with open(data_infile, 'r') as f:
        for line in f:
            if parsing_year:
                year_list += [line.strip()] if line.strip() not in ['', ';'] else []
            if parsing_fuel:
                fuel_list += [line.strip()] if line.strip() not in ['', ';'] else [] 
            if parsing_tech:
                tech_list += [line.strip()] if line.strip() not in ['', ';'] else [] 
            if parsing_storage:
                storage_list += [line.strip()] if line.strip() not in ['', ';'] else [] 
            if parsing_mode:
                mode_list += [line.strip()] if line.strip() not in ['', ';'] else [] 
            if parsing_emission:
                emission_list += [line.strip()] if line.strip() not in ['', ';'] else []

            if line.startswith('set YEAR'):
                if len(line.split('=')[1]) > 1:
                    year_list = line.split(' ')[3:-1]
                else:
                    parsing_year = True
            if line.startswith('set COMMODITY'):  # Extracts list of COMMODITIES from data file. Some models use FUEL instead. 
                if len(line.split('=')[1]) > 1:
                    fuel_list = line.split(' ')[3:-1]
                else:
                    parsing_fuel = True
            if line.startswith('set FUEL'):  # Extracts list of FUELS from data file. Some models use COMMODITIES instead. 
                if len(line.split('=')[1]) > 1:
                    fuel_list = line.split(' ')[3:-1]
                else:
                    parsing_fuel = True
            if line.startswith('set TECHNOLOGY'):
                if len(line.split('=')[1]) > 1:
                    tech_list = line.split(' ')[3:-1]
                else:
                    parsing_tech = True
            if line.startswith('set STORAGE'):
                if len(line.split('=')[1]) > 1:
                    storage_list = line.split(' ')[3:-1]
                else:
                    parsing_storage = True
            if line.startswith('set MODE_OF_OPERATION'):
                if len(line.split('=')[1]) > 1:
                    mode_list = line.split(' ')[3:-1]
                else:
                    parsing_mode = True
            if line.startswith('set EMISSION'):
                if len(line.split('=')[1]) > 1:
                    emission_list = line.split(' ')[3:-1]
                else:
                    parsing_emission = True

            if line.startswith(";"):
                parsing_year = False
                parsing_tech = False
                parsing_fuel = False
                parsing_mode = False
                parsing_storage = False
                parsing_emission = False

    start_year = year_list[0]

    if data_format == 'momani':
        with open(data_infile, 'r') as f:
            for line in f:
                if line.startswith(";"):
                    parsing = False
                if parsing:
                    if line.startswith('['):
                        fuel = line.split(',')[2]
                        tech = line.split(',')[1]
                        emission = line.split(',')[2]
                    elif line.startswith(start_year):
                        years = line.rstrip(':= ;\n').split(' ')[0:]
                        years = [i.strip(':=') for i in years]
                    else:
                        values = line.rstrip().split(' ')[1:]
                        mode = line.split(' ')[0]

                        if param_current == 'OutputActivityRatio':    
                            data_out.append(tuple([fuel, tech, mode]))
                            for i in range(0, len(years)):
                                output_table.append(tuple([tech, fuel, mode, years[i], values[i]]))

                        if param_current == 'InputActivityRatio':
                            data_inp.append(tuple([fuel, tech, mode]))   

                        data_all.append(tuple([tech, mode]))

                        if param_current == 'TechnologyToStorage':
                            if not line.startswith(mode_list[0]):
                                storage = line.split(' ')[0]
                                values = line.rstrip().split(' ')[1:]
                                for i in range(0, len(mode_list)):
                                    if values[i] != '0':
                                        storage_to.append(tuple([storage, tech, mode_list[i]]))

                        if param_current == 'TechnologyFromStorage':
                            if not line.startswith(mode_list[0]):
                                storage = line.split(' ')[0]
                                values = line.rstrip().split(' ')[1:]
                                for i in range(0, len(mode_list)):
                                    if values[i] != '0':
                                        storage_from.append(tuple([storage, tech, mode_list[i]]))

                        if param_current == 'EmissionActivityRatio':
                            emission_table.append(tuple([emission, tech, mode]))

                if line.startswith(('param OutputActivityRatio', 'param InputActivityRatio', 'param TechnologyToStorage', 'param TechnologyFromStorage', 'param EmissionActivityRatio')):
                    param_current = line.split(' ')[1]
                    parsing = True

    if data_format == 'otoole':
        with open(data_infile, 'r') as f:
            for line in f:
                if line.startswith(";"):
                    parsing = False
                if parsing:
                    if len(line.split(' ')) > 1:
                        if param_current == 'OutputActivityRatio':
                            tech = line.split(' ')[1].strip()
                            fuel = line.split(' ')[2].strip()
                            mode = line.split(' ')[3].strip()
                            year = line.split(' ')[4].strip()
                            value = line.split(' ')[5].strip()

                            data_out.append(tuple([fuel, tech, mode]))
                            output_table.append(tuple([tech, fuel, mode, year, value]))

                        if param_current == 'InputActivityRatio':
                            tech = line.split(' ')[1].strip()
                            fuel = line.split(' ')[2].strip()
                            mode = line.split(' ')[3].strip()

                            data_inp.append(tuple([fuel, tech, mode]))

                        data_all.append(tuple([tech, mode]))

                        if param_current == 'TechnologyToStorage':
                            tech = line.split(' ')[1].strip()
                            storage = line.split(' ')[2].strip()
                            mode = line.split(' ')[3].strip()

                            storage_to.append(tuple([storage, tech, mode]))

                        if param_current == 'TechnologyFromStorage':
                            tech = line.split(' ')[1].strip()
                            storage = line.split(' ')[2].strip()
                            mode = line.split(' ')[3].strip()

                            storage_from.append(tuple([storage, tech, mode]))

                        if param_current == 'EmissionActivityRatio':
                            tech = line.split(' ')[1].strip()
                            emission = line.split(' ')[2].strip()
                            mode = line.split(' ')[3].strip()

                            emission_table.append(tuple([emission, tech, mode]))

                if any(param in line for param in param_check):
                    param_current = line.split(' ')[-2]
                    parsing = True

    data_out = list(set(data_out))
    data_inp = list(set(data_inp))
    data_all = list(set(data_all))
    storage_to = list(set(storage_to))
    storage_from = list(set(storage_from))
    emission_table = list(set(emission_table))

    dict_out = defaultdict(list)
    dict_inp = defaultdict(list)
    dict_all = defaultdict(list)
    dict_stt = defaultdict(list)
    dict_stf = defaultdict(list)
    dict_emi = defaultdict(list)

    for fuel, tech, mode in data_out:
        dict_out[fuel].append((mode, tech))

    for fuel, tech, mode in data_inp:
        dict_inp[fuel].append((mode, tech))

    for tech, mode in data_all:
        if mode not in dict_all[tech]:
            dict_all[tech].append(mode)

    for storage, tech, mode in storage_to:
        dict_stt[storage].append((mode, tech))

    for storage, tech, mode in storage_from:
        dict_stf[storage].append((mode, tech))

    for emission, tech, mode in emission_table:
        dict_emi[emission].append((mode, tech))

    def file_output_function(if_dict, str_dict, set_list, set_name, extra_char):
        for each in set_list:
            if each in if_dict.keys():
                line = set_name + str(each) + ']:=' + str(str_dict[each]) + extra_char
                if set_list == tech_list:
                    line = line.replace(',','').replace(':=[',':= ').replace(']*','').replace("'","")
                else:
                    line = line.replace('),',')').replace('[(',' (').replace(')]',')').replace("'","")
            else:
                line = set_name + str(each) + ']:='
            file_out.write(line + ';' + '\n')

    # Append lines at the end of the data file
    with open(data_outfile, 'w') as file_out:  # 'a' to open in 'append' mode

        file_out.writelines(lines)

        file_output_function(dict_out, dict_out, fuel_list, 'set MODExTECHNOLOGYperFUELout[', '')
        file_output_function(dict_inp, dict_inp, fuel_list, 'set MODExTECHNOLOGYperFUELin[', '')
        file_output_function(dict_all, dict_all, tech_list, 'set MODEperTECHNOLOGY[', '*')

        if len(storage_list) > 1:
            file_output_function(dict_stt, dict_out, storage_list, 'set MODExTECHNOLOGYperSTORAGEto[', '')
            file_output_function(dict_stf, dict_out, storage_list, 'set MODExTECHNOLOGYperSTORAGEfrom[', '*')

        if len(emission_list) > 1:
            file_output_function(dict_emi, dict_emi, emission_list, 'set MODExTECHNOLOGYperEMISSION[', '')

        file_out.write('end;')


if __name__ == '__main__':

    if len(sys.argv) != 4:
        msg = "Usage: python {} <otoole/momani> <infile> <outfile>"
        print(msg.format(sys.argv[0]))
        sys.exit(1)
    else:
        data_format = sys.argv[1]
        data_infile = sys.argv[2]
        data_outfile = sys.argv[3]
        main(data_format, data_infile, data_outfile)
