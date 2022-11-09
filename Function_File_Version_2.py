#Used packages:

import matplotlib.pyplot as plt
import pandas as pd
import csv 
import numpy as np
import scipy.optimize as opt
import array
from itertools import combinations
import ipynb

# _____________________________________________________________________________________________________________________________________________________________
# 1) Importing Data
# _____________________________________________________________________________________________________________________________________________________________



# _____________________________________________________________________________________________________________________________________________________________
# 2) Import and creation of Header and Data Dataframe
# _____________________________________________________________________________________________________________________________________________________________


def create_dataframe_new(file_path):
    """
    Reads in a csv file and creates a header and a data section

    Args:
        file_path (str): a relative or abolute pathway to the file, which contains the whole information

    Returns:
        header_section (dataframe): pandas dataframe of the header section
        data_section_raw (dataframe): pandas dataframe of the data section

    """

    dataframe = pd.read_csv(file_path)

    header_dataframe = dataframe[0:12][dataframe.columns[0:2]]
    print("______________________________________________________")
    print("Header Section:")
    print("______________________________________________________")
    print(header_dataframe)
    print("______________________________________________________")

    column_names = dataframe.loc[dataframe.index == 14].values.flatten().tolist()
    #print(column_names)

    data_dataframe = dataframe[15:]
    data_dataframe.columns = column_names
    data_dataframe = data_dataframe.dropna(axis=0,how="all")
    #drop rows that are completely empty
    data_dataframe = data_dataframe.fillna(0)
    #replacing remaing NaN with 0 as they cause trouble later on
    #wirklich notwendig?

    return header_dataframe, data_dataframe





# _____________________________________________________________________________________________________________________________________________________________
# 3) Adjusting global parameters of data dataframe
#   + creating local parameters csv
# _____________________________________________________________________________________________________________________________________________________________

def adjust_dataframe_global(data_section, header_section, local_parameters_template, local_parameters_path):
    """
    Turns the raw data into plottable FRET data regarding global parameters

    Args:
        data_section (dataframe): a raw dataframe containing the to be adjusted data
        header_section (dataframe): contains information about the dataframe like the global adjusting parameter
        local_parameters_template (dict): template entry for creating the Local Parameter CSV File
        local_parameters_path (str): path where Local Parameter CSV File gets saved

    Returns:
        data_section (dataframe): the adjusted dataframe
    """   

    #GLOBAL PARAMETERS ----------------------------------------------------------------------------------------------------------------------------------------

    # X-axis: Frames per second to seconds
    # Y-axis: Photon counts to frequency (Counts/s)

    data_section = data_section.rename(columns={'Frame index': 'Time', 'Donor Photons': 'Frequency Donor', 'Acceptor Photons': 'Frequency Acceptor'})

    header_series = header_section.set_index("Settings:")["Unnamed: 1"]
    frames_per_second = header_series["Frames per second"]
    #create series to access values


    frames_per_second = float(frames_per_second)
    frames_per_second = round(frames_per_second, 2)


    for row in range(len(data_section)):

        data_section.iloc[row]["Time"] = float(data_section.iloc[row]["Time"]) / frames_per_second

        try:
            data_section.iloc[row]["Frequency Donor"] = float(data_section.iloc[row]["Frequency Donor"]) / frames_per_second
        except ZeroDivisionError:
            data_section.iloc[row]["Frequency Donor"] = 0

        try:
            data_section.iloc[row]["Frequency Acceptor"] = float(data_section.iloc[row]["Frequency Acceptor"]) / frames_per_second
        except ZeroDivisionError:
            data_section.iloc[row]["Frequency Acceptor"] = 0


        convert_dict = {"ID": int,
                        "Time": int,
                        "Groups": int,
                        "Frequency Donor": float,
                        "Frequency Acceptor": float,
                        "FRET efficiency": float,
                        "X (nm)": float,
                        "Y (nm)": float
                             }
 
        data_section = data_section.astype(convert_dict)

    # Creating the Local Parameter CSV File ---------------------------------------------------------------------------------

    ID_set = set()

    for row in range(len(data_section)):
        ID_set.add(int(data_section.iloc[row]["ID"]))

    ID_dict = {}

    for ID in ID_set:
        ID_dict[ID] = local_parameters_template
        
    df_ID_dict = pd.DataFrame.from_dict(ID_dict, orient='index')

    df_ID_dict.to_csv (local_parameters_path, header=True)

    return data_section










# _____________________________________________________________________________________________________________________________________________________________
# 4) Adjusting local parameters of data dataframe
# _____________________________________________________________________________________________________________________________________________________________


def adjust_dataframe_local(data_section, local_parameters_path):
    """
    Turns the raw data into plottable FRET data regarding local parameters

    Args:
        data_section (dataframe): a raw dataframe containing the to be adjusted data
        local_parameters_path (str): csv file, found at this path, contains the adjusting parameteres for every ID {IDX:[Group, Bin, Bleach, Blink]}

    Returns:
        data_section_new (dataframe): the adjusted, ready to print dataframe
    """

    # LOCAL PARAMETERS: Preparations ---------------------------------------------------------------------------------------------------------------------------

    column_names = list(data_section.columns)
    data_section_new = pd.DataFrame(columns = column_names)
    #This empty dataframe will be filled at the very end at "binning"
    
    ID_set = set()

    for row in range(len(data_section)):
        ID_set.add(int(data_section.iloc[row]["ID"]))
        #The ID Set is requirred for blinking and bleaching to deal with the IDs seperately

    local_parameters_dict = dict()
    local_parameters_csv = open(local_parameters_path)
    #The lp csv will be used to fill the lp dict
    #The lp dict will be needed to access the tresholds for blinking/bleaching/binning for each ID

    for counter, line in enumerate(local_parameters_csv):

        if counter == 0:
            continue

        line = line. strip('\n')
        ID, Group,  Bin,  Bleach_range,  Bleach_ratio,  Blink,  Additional_cutoff = line.split(",")
        local_parameters_dict[int(ID)] = {"Group": int(Group), "Bin" : int(Bin), "Bleach range" : int(Bleach_range), "Bleach ratio" : int(Bleach_ratio), "Blink" : int(Blink), "Additional cutoff" : int(Additional_cutoff)}

    # BLEACHING ---------------------------------------------------------------------------------------------------------------

    output_data_section = data_section.copy()  

    for current_ID in ID_set:

        iter_data_section = data_section.copy()

        # Der übergeordnete Loop ist zwar umständlich (alle Reihen für jede ID),
        # aber sonst ist es schwierig den row Loop zu beenden wenn Zeilen gelöscht wurden
        
        bleaching_range = local_parameters_dict[current_ID]["Bleach range"]
        bleaching_ratio = local_parameters_dict[current_ID]["Bleach ratio"]
        additional_cutoff = local_parameters_dict[current_ID]["Additional cutoff"]

        bleaching_counter = 0

        for row in range(len(iter_data_section)):

            if iter_data_section.iloc[row]["ID"] != current_ID:
                continue            

            if iter_data_section.iloc[row]["Frequency Donor"] > iter_data_section.iloc[row]["Frequency Acceptor"] * bleaching_ratio:
                bleaching_counter = bleaching_counter + 1
            
            if bleaching_counter >= bleaching_range:

                print("__________________________________________________________________________________")
                print(f"current row: {row}")
                print(f"current ID: {current_ID}")


                to_be_removed_time_start = iter_data_section.iloc[row - bleaching_range - additional_cutoff]["Time"]
                print(f"Time to be removed start: {to_be_removed_time_start}")
                to_be_removed_time_end = iter_data_section[iter_data_section.ID == current_ID]["Time"].max()
                print(f"Time to be removed end: {to_be_removed_time_end}")

                data_section_ID = data_section[data_section.ID == current_ID]
                first_row_org_df = data_section_ID[data_section_ID["Time"] == to_be_removed_time_start].index.tolist()
                print(f"first row index = {first_row_org_df[0]}")
                last_row_org_df = data_section_ID[data_section_ID["Time"] == to_be_removed_time_end].index.tolist()
                print(f"last row index = {last_row_org_df[0]}")


                output_data_section["Frequency Donor"][first_row_org_df[0] : last_row_org_df[0]] = 0
                output_data_section["Frequency Acceptor"][first_row_org_df[0] : last_row_org_df[0]] = 0
                #Indexing problem gets resolved by not removing lines, but replacing them with 0
                #Zeros won't be plotted and/or can be removed afterwards 

                print(f"Bleaching happended at ID = {current_ID}, time start = {to_be_removed_time_start} and time end = {to_be_removed_time_end}")
                
                break


    # Blinking ----------------------------------------------------------------------------------------------------------------------------------------------


    for current_ID in ID_set:

        iter_data_section = data_section.copy()
        #hier nicht Kopie von output data_section?

        blinking_value = local_parameters_dict[current_ID]["Blink"]
        additional_cutoff = local_parameters_dict[current_ID]["Additional cutoff"]

        for row in range(len(iter_data_section)):

            if iter_data_section.iloc[row]["ID"] != current_ID:
                continue

            if iter_data_section.iloc[row]["Frequency Acceptor"] != 0:
                continue

            list_of_values = []

            for line in range(blinking_value):
                if iter_data_section.iloc[row + line].ID == current_ID:
                    list_of_values.append(iter_data_section.iloc[row+line]["Frequency Acceptor"])
                else:
                    continue

            
            if all([ values == 0.0 for values in list_of_values]) and list_of_values != [] and len(list_of_values) == blinking_value:

                print("__________________________________________________________________________________")
                print(f"current row: {row}")
                print(f"current ID: {current_ID}")
                print(f"list of values: {list_of_values}")

                to_be_removed_time_start = iter_data_section.iloc[row-additional_cutoff]["Time"]
                print(f"Time to be removed start: {to_be_removed_time_start}")
                to_be_removed_time_end = iter_data_section[iter_data_section.ID == current_ID]["Time"].max()
                print(f"Time to be removed end: {to_be_removed_time_end}") 

                data_section_ID = output_data_section[output_data_section['ID']==current_ID]
                first_row_org_df = data_section_ID[data_section_ID["Time"] == to_be_removed_time_start].index.tolist()
                print(f"first row index = {first_row_org_df[0]}")
                last_row_org_df = data_section_ID[data_section_ID["Time"] == to_be_removed_time_end].index.tolist()
                print(f"last row index = {last_row_org_df[0]}")


                output_data_section["Frequency Donor"][first_row_org_df[0] : last_row_org_df[0]] = 0
                output_data_section["Frequency Acceptor"][first_row_org_df[0] : last_row_org_df[0]] = 0


                print(f"Blinking happended at ID = {current_ID}, time start = {to_be_removed_time_start} and time end = {to_be_removed_time_end}")

                break
            

    # BINNING -----------------------------------------------------------------------------------------------------------------------------------------------
    
    row_number = 0
    counter = 0

    for row in range(len(output_data_section)):
        
        current_ID = int(output_data_section.iloc[row]["ID"])
        bin = local_parameters_dict[current_ID]["Bin"]

        if counter % bin == 0 and current_ID == int(output_data_section.iloc[row]["ID"]):

            data_bin = output_data_section.iloc[counter : counter + bin]
            data_bin = data_bin.astype(float)
            #um .mean() ausführen zu können brauchen wir floats

            data_means = data_bin.mean(axis=0)


            if int(data_means["ID"]) != data_bin["ID"].iloc[0]:
                counter = counter +1
                continue
            #Inhomogener Bin wird hier verworfen


            data_section_new.loc[row_number] = ["","","","","","","",""]

            data_section_new.iloc[row_number]["ID"] = int(data_means["ID"])
            data_section_new.iloc[row_number]["Time"] = int(data_means["Time"])
            data_section_new.iloc[row_number]["Groups"] = int(data_means["Groups"])
            data_section_new.iloc[row_number]["Frequency Donor"] = data_means["Frequency Donor"]
            data_section_new.iloc[row_number]["Frequency Acceptor"] = data_means["Frequency Acceptor"]
            data_section_new.iloc[row_number]["FRET efficiency"] = data_means["FRET efficiency"]
            data_section_new.iloc[row_number]["X (nm)"] = data_means["X (nm)"]
            data_section_new.iloc[row_number]["Y (nm)"] = data_means["Y (nm)"]

            row_number = row_number + 1
            counter = counter +1

        else: counter = counter +1


    return data_section_new


# _____________________________________________________________________________________________________________________________________________________________
# 5) Visualisation of Data with subplots
# _____________________________________________________________________________________________________________________________________________________________


def plot_dataframe(data_section_global_and_local, data_section_global, local_parameters_path ,group_dict):
    """
    Prints the dataframe, which one subplot per ID

    Args:
        data_section_global_and_local (dataframe): Contains the globally and locally adjusted, to be plotted data
        data_section_global (dataframe): Contains the to globally adjusted, to be plotted data
        local_parameters_dict (dict): Binning, Bleaching and Blinking Parameters for plot description
        group_dict (dict): Contains the references for groups for the plot description

    Returns:
        None
    """

    # Collecting occuring ID Values

    ID_set = set()

    for row in range(len(data_section_global)):
        ID_set.add(int(data_section_global.iloc[row]["ID"]))

    print(f'ID_set: {ID_set}')

    # How many are there?

    amount_IDs = len(ID_set)
    print(f'There are {amount_IDs} different IDs')


    # Building the Parameters Dict

    local_parameters_dict = dict()
    local_parameters_csv = open(local_parameters_path)

    for counter, line in enumerate(local_parameters_csv):

        if counter == 0:
            continue

        line = line. strip('\n')
        ID, Group,  Bin,  Bleach_range,  Bleach_ratio,  Blink,  Additional_cutoff = line.split(",")
        local_parameters_dict[int(ID)] = {"Group": int(Group), "Bin" : int(Bin), "Bleach range" : int(Bleach_range), "Bleach ratio" : int(Bleach_ratio), "Blink" : int(Blink), "Additional cutoff" : int(Additional_cutoff)}

    # For each ID make one subplot!

    fig, axes = plt.subplots(amount_IDs, 2, figsize=(25, amount_IDs*7))
    fig.suptitle("smFRET", fontsize=30)
    fig.subplots_adjust(hspace = 0.5, wspace = 0.2)
    #fig.tight_layout(h_pad = 2)


    # axes.spines["right"].set_visible(False)
    # axes.spines["top"].set_visible(False)

    for number, element in enumerate(ID_set):

        # Print datasection global on the left side -------------------------------------------------------------------------
        
        x_values_g = []
        yD_values_g = []
        yA_values_g = []

        filtered_data_g = data_section_global.loc[data_section_global.ID == int(element)]
    
        for line in range(len(filtered_data_g)):
            x_values_g.append(filtered_data_g.iloc[line]["Time"])
            yD_values_g.append(filtered_data_g.iloc[line]["Frequency Donor"])
            yA_values_g.append(filtered_data_g.iloc[line]["Frequency Acceptor"])


        x_limit_g = data_section_global.Time.loc[data_section_global.ID == int(element)].max()
        x_limit_g = x_limit_g * 1.05
        y_limit_Acceptor_g = data_section_global.loc[:,"Frequency Acceptor"].loc[data_section_global.ID == int(element)].max()
        y_limit_Donor_g = data_section_global.loc[:,"Frequency Donor"].loc[data_section_global.ID == int(element)].max()
        y_limit_g = max([y_limit_Acceptor_g,y_limit_Donor_g])*1.2


        row_g = number
        coloum_g = 0

        axes[row_g, coloum_g].plot(x_values_g, yD_values_g, label = "Donor")
        axes[row_g, coloum_g].plot(x_values_g, yA_values_g, label = "Acceptor")
        axes[row_g, coloum_g].set_title(f"ID: {element}\nData only adjusted to fit axes")
        #axes[row_g, coloum_g].set_xticks([0,x_tick_1,x_tick_2,x_tick_3])
        axes[row_g, coloum_g].set_xlim(0,x_limit_g)
        #axes[row_g, coloum_g].set_yticks([0,y_tick_1,y_tick_2,y_tick_3])
        axes[row_g, coloum_g].set_ylim(0,y_limit_g)
        axes[row_g, coloum_g].legend()
        #axes[row, coloum].xticks(rotation=45)
        axes[row_g, coloum_g].set_axisbelow(True)
        #axes[row_g, coloum_g].grid()


        # Print datasection global and local on the right side  -------------------------------------------------------------

        group_of_ID = local_parameters_dict[element]["Group"]
        group_name = group_dict[group_of_ID]

        binsize = local_parameters_dict[element]["Bin"]
        bleach_ratio = local_parameters_dict[element]["Bleach ratio"]
        bleach_range = local_parameters_dict[element]["Bleach range"]
        blink = local_parameters_dict[element]["Blink"]


        x_values_g_l = []
        yD_values_g_l = []
        yA_values_g_l = []

        filtered_data_g_l = data_section_global_and_local.loc[data_section_global_and_local.ID == element]

        for line in range(len(filtered_data_g_l)):

            x_values_g_l.append(filtered_data_g_l.iloc[line]["Time"])
            yD_values_g_l.append(filtered_data_g_l.iloc[line]["Frequency Donor"])
            yA_values_g_l.append(filtered_data_g_l.iloc[line]["Frequency Acceptor"])

        x_limit_g_l = data_section_global_and_local.Time.loc[data_section_global_and_local.ID == element][data_section_global_and_local["Frequency Donor"] != 0].max()
        x_limit_g_l = x_limit_g_l *1.2
        y_limit_Acceptor_g_l = data_section_global_and_local.loc[:,"Frequency Acceptor"].loc[data_section_global_and_local.ID == element].max()
        y_limit_Donor_g_l = data_section_global_and_local.loc[:,"Frequency Donor"].loc[data_section_global_and_local.ID == element].max()
        y_limit_g_l = max([y_limit_Acceptor_g_l,y_limit_Donor_g_l])*1.2

        row_g_l = number
        coloum_g_l = 1

        axes[row_g_l, coloum_g_l].plot(x_values_g_l, yD_values_g_l, label = "Donor")
        axes[row_g_l, coloum_g_l].plot(x_values_g_l, yA_values_g_l, label = "Aceptor")
        axes[row_g_l, coloum_g_l].set_title(f"ID: {element}\nGroup: {group_name}, Binsize: {binsize} frames\nCutoffs - Bleaching range:{bleach_range} frames, Bleaching ratio: {bleach_ratio}, Blinking:{blink} frames")
        #axes[row_g_l, coloum_g_l].set_xticks([0,x_tick_1,x_tick_2,x_tick_3])
        axes[row_g_l, coloum_g_l].set_xlim(0,x_limit_g_l)
        #axes[row_g_l, coloum_g_l].set_yticks([0,y_tick_1,y_tick_2,y_tick_3])
        axes[row_g_l, coloum_g_l].set_ylim(0,y_limit_g_l)
        axes[row_g_l, coloum_g_l].legend()
        #axes[row, coloum].xticks(rotation=45)
        axes[row_g_l, coloum_g_l].set_axisbelow(True)
        #axes[row_g_l, coloum_g_l].grid()

    for axes in axes.flat:
        axes.set(xlabel='Time [s]', ylabel='Frequency [Photons/s]')
        #axes.yticks([0,200,400,600,800])
        #axes.xticks([0,200,400,600])
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        #axes.spines["left"].set_visible(False)
        #axes.spines["bottom"].set_visible(False)

    #fig.legend(ID_set, ncol = len(ID_set))



# _____________________________________________________________________________________________________________________________________________________________
# 6) Haar Cycles
# _____________________________________________________________________________________________________________________________________________________________


def haar_cyles (dataframe, number_cycles):

    # FRET Werte aus Dataframe in eine lange Liste geben 

    FRET_list = []
    for element in dataframe["FRET efficiency"]:
        FRET_list.append(float(element))

    # Je ein Dict erstellen - Mit je zwei Werten und Mittelwert - 1 lange liste mit Dicts
        # Wenn einer der Werte 0 - überspringen
        # Wenn einer der Werte aus neuer ID - überspringen (noch nicht implementiert)

    new_list = []

    for number, element in enumerate(FRET_list):
        if number % 2 == 0:
            continue
        
        new_dict = {"First": FRET_list[number-1] ,"Second": FRET_list[number], "Mean": (FRET_list[number]+FRET_list[number-1])/2}
        
        if new_dict["First"] != 0 and new_dict["Second"] != 0:
            new_list.append(new_dict)

    #print(new_list)

    # Pro Zyklus: (0 Zyklen = überspringen)
        # Von zwei Dicts hintereinander die Mittelwerte als Werte für neues Dict
        # Mittelwert dieser berechnen und ebenfalls in Dict
        # Belieblig oft wiederholen

    while number_cycles > 0:

        iter_list = new_list.copy()
        new_list = []

        #print("____________________________________________")
        #print("One cycle happend")
        #print("____________________________________________")

        for number, element in enumerate(iter_list):

            try:

                if number % 2 == 0:
                    continue
        
                new_dict = {"First": iter_list[number-1]["Mean"] ,"Second": iter_list[number]["Mean"], "Mean": (iter_list[number]["Mean"]+iter_list[number-1]["Mean"])/2}
                new_list.append(new_dict)
                
            except IndexError:
                continue

        number_cycles = number_cycles - 1

        #print (new_list)

    # Abweichung von Werten zu  Mittelwert berechnen
    # Diese Abweichung in Liste geben

    noises = []

    for element in new_list:
        noise = element["Mean"] - element["First"]
        noises.append(abs(noise))
        #noises.append(noise)

    # print("\n-------------------------------\n")
    # print("List of noises")
    # print(noises)

    return noises



# _____________________________________________________________________________________________________________________________________________________________
# 7) Histogram fit
# _____________________________________________________________________________________________________________________________________________________________

def histogram_fit (noises, binsize, cutoff = 68.2):

    y_values, bin_edges, third_variable = plt.hist(noises, binsize)
    
    bin_edges_list = []

    for element in bin_edges:
        bin_edges_list.append(element)


    bin_width = (bin_edges_list[1] - bin_edges_list [0])

    x_values = []

    for element in bin_edges_list:
        current_value = element + bin_width/2
        x_values.append(current_value)

    x_values.pop()
    # Letzter Wert muss weg weil letzte Bin-Grenze + 1/2 Binsize keinem Datenpunkt entspricht!

    plt.scatter(x_values, y_values, c = "red")
    #plt.figure(figsize=(50, 100), dpi = 80)
    plt.rcParams['figure.figsize'] = [10, 10]
    plt.show

    y_values_cumu = []
    current_value = 0

    for element in y_values:
        current_value += element
        y_values_cumu.append(current_value)

    #print(y_values_cumu)

    plt.scatter(x_values, y_values_cumu, c = "green")
    plt.show

    hundred_percent = y_values_cumu[-1]
    value_crit = hundred_percent * cutoff/100

    print(f"The value at {cutoff}% is {value_crit}")

    plt.axhline(y=value_crit, color="yellow", linestyle='-')
    plt.show

    y_values_cumu = array.array("f", y_values_cumu)
    #?????????????????????????????????????????

    return y_values, y_values_cumu, x_values, value_crit




# _____________________________________________________________________________________________________________________________________________________________
# 8) Gaussian Curve fit
# _____________________________________________________________________________________________________________________________________________________________


def gaussian_curve_fit (x_values, y_values):

    # das unterhalb funktioniert nur wenn die y-werte arrays sind!
    # y values ist array, aber y values cumu ist liste!

    #y_values = array.array("f", y_values)

    print(y_values)
    print(type(y_values))

    random_list = [4,2,6,2,7,7,5,3,2,5,6,3]

    random_array = array.array("f", random_list)
    print(random_array)
    print(type(random_array))

    # random_numpy = np.ndarray(shape = (),dtype = float, buffer = random_list)
    # print(random_numpy)


    #Dieser Teil derzeit nur Copy Paste!

    n = len(x_values)                         
    mean = sum(x_values*y_values)/n                  
    sigma = sum(y_values*(x_values-mean)**2)/n 

    def gauss_function(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
        # a = height of the curves peak
        # x0 = position of center of peak
        # sigma = controlls the width of the "bell"
        
    popt,pcov = opt.curve_fit(f = gauss_function, xdata = x_values, ydata = y_values, p0=[1,mean,sigma], maxfev = 10000)
    #maxfev = 5000 um error zu vermeiden: "RuntimeError: Optimal parameters not found: Number of calls to function has reached maxfev = 800"
    #Gauß fittet nicht nur sigma sondern auch müh etc!
    #sigma zuvor zwar "richtiger" aber "blind" für alles andere

    plt.plot(x_values,y_values,'b+:',label='data')
    plt.plot(x_values,gauss_function(x_values,*popt),'ro:',label='fit')
    plt.legend()
    plt.title('Titel')
    plt.xlabel('x label')
    plt.ylabel('y label')
    plt.show()

    
    print("popt: letzter der drei Werte ist Sigma")
    print(popt)
    print("simga 2:")
    print(sigma)




# _____________________________________________________________________________________________________________________________________________________________
# 9) Find transition points
# _____________________________________________________________________________________________________________________________________________________________

def find_transition_points(liste_y, sigma, transition_search_loops, transition_treshold):

    print("___________________________________________________")
    print("FIND TRANSITION POINTS")
    print("___________________________________________________")

    liste_x = []
    number = 1
    while number < len(liste_y)+1:
        liste_x.append(number)
        number += 1

    transition_points = [{"index": 0, "value": 0}]
    sigma = 2

    number = 0

    print(len(liste_y))

    #while len(transition_points) < 4
    while number < transition_search_loops:

        print("NEW ROUND---------------------------------------------")

        candidates_for_transition = []

        for i,transition_point in enumerate(transition_points):
            start = transition_points[i]["index"]
            print("START:")
            print(start)

            try:
                end = transition_points[i+1]["index"]
            except IndexError:
                end = -1

            print("END")
            print(end)

            current_list = liste_y[start:end]
            #print(current_list)

            values = []

            for index, element in enumerate (current_list):

                try:
                    one_by_N_minus_i = 1/(len(current_list)-index)
                except ZeroDivisionError:
                    one_by_N_minus_i = 0

                try:
                    one_by_i = 1/index
                except ZeroDivisionError:
                    one_by_i = 0

                try:
                    average_before = sum(current_list[0:index])/len(current_list[0:index])
                except ZeroDivisionError:
                    average_before = 0
                try:
                    average_after = sum(current_list[index:-1])/len(current_list[index:-1])
                except ZeroDivisionError:
                    average_after = 0

                numerator = average_after - average_before
                denominator = sigma*((one_by_N_minus_i+ one_by_i)**0.5)
                value = numerator/denominator
                values.append(value)

            #print(values)

            plt.plot(liste_x[start:end], values)

            try:
                max_value = max(values[1:], key=abs)
            except ValueError:
                continue
            print("MAX VALUE:")
            print(max_value)
            max_index_local = values.index(max_value)
            max_index_global = max_index_local + start
            print("MAX INDEX GLOBAL:")
            print(max_index_global)
            max_dict = {"index": max_index_global , "value": abs(max_value)}
            candidates_for_transition.append(max_dict)

        print("candidates_for_transition:")
        print(candidates_for_transition)

        max_candidate = max(candidates_for_transition,  key=lambda x:x["value"])
        max_transitions = max(transition_points, key=lambda x:x["value"])

        # max_candidate = max(item["value"] for item in candidates_for_transition)
        # max_transitions = max(item["value"] for item in transition_points)

        if max_candidate["value"] > max_transitions["value"]*transition_treshold:
            transition_points.append(max_candidate)
        else:
            print("Value rejectedXXXXXXXXXXXXXXXXXXXXXXX")

        number = number + 1

    print(transition_points)


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    plt.plot(liste_x[1:], liste_y[1:])
    #plt.plot(liste_x[1:], values[1:])

    for element in transition_points:
        plt.axvline(x = element["index"], color="black", linestyle='-')

    plt.axhline(y = 0, color="black", linestyle='-')

    plt.show()

    return transition_points

#_____________________________________________________________________________________________


# x_array = np.array(liste_x,dtype = np.int)
# y_array = np.array(values, dtype = np.float)
# slope = np.diff(y_array)/np.diff(x_array)

# print (slope)
# print(np.diff(slope))

# plt.plot(liste_x[2:-1], slope[2:])
# plt.plot(liste_x[2:-2], np.diff(slope[2:]))

#_____________________________________________________________________________________________



# von values höchsten absolut wert finden
# den dazugehörigen x wert in transition point liste speichern 
# dann die liste der y werte teilen an diesem punkt
# für die zwei listen nun wieder values berechnen
# den höchsten der beiden werte (absolut) nun wieder in transition points abspeichern
# diese liste wieder spalten
# mal ein while loop mit bestimmten wiederholungen, sonst vl mit statsistischem tool?






# _____________________________________________________________________________________________________________________________________________________________
# 10) Find levels
# _____________________________________________________________________________________________________________________________________________________________

def find_levels(liste_y, transition_points, blending_cycles, max_diff):

    print("___________________________________________________")
    print("FIND LEVELS")
    print("___________________________________________________")

    
    means = []

    liste_x = []
    number = 1
    while number < len(liste_y)+1:
        liste_x.append(number)
        number += 1

    plt.plot(liste_x[1:], liste_y[1:])

    for i,transition_point in enumerate(transition_points):
        start = transition_points[i]["index"]
        print("START:")
        print(start)

        try:
            end = transition_points[i+1]["index"]
        except IndexError:
            end = -1

        print("END")
        print(end)

        current_list = liste_y[start:end]

        try:
            current_mean = sum(current_list)/len(current_list)
        except ZeroDivisionError:
            continue
        
        means.append({"mean": current_mean, "start": start, "end": end})

        if end != -1:
            plt.plot((start, end), (current_mean, current_mean), 1)
        else:
            plt.plot((start, max(liste_x)), (current_mean, current_mean), 1)

    print ("MEANS:")
    print(means)

    plt.show()

    binsize = 20

    means_for_histo = []
    for element in means:
        means_for_histo.append(element["mean"])

    lol = histogram_fit (means_for_histo, binsize, cutoff = 68.2)
    plt.show()

    number_of_cycles = 0

    while number_of_cycles < blending_cycles:

        print ("___________________________________________________________________________________")
        print(f"Blending Cycle Nr. {number_of_cycles+1}")
        print ("___________________________________________________________________________________")
        
        pairings = []
        temp = combinations(means, 2)
        for element in list(temp):
            current_diff = element[0]["mean"]-element[1]["mean"]
            pairings.append({"elements": element, "diff": abs(current_diff)})
            
        print("pairings:")
        print(pairings)

        # min_diff = min(pairings, key=lambda x:x["diff"])

        #[{'elements': ({'mean': 2.5609756097560976, 'start': 0, 'end': 41}, {'mean': 7.476190476190476, 'start': 41, 'end': 62}), 
        #'diff': 4.915214866434379}, {'elements': ({'mean': 2.5609756097560976, 'start': 0, 'end': 41}, {'mean': 2.0, 'start': 62, 'end': 80}), 'diff': 0.5609756097560976}, 

        current_min = 0

        for index, _ in enumerate(pairings):
            if pairings[index - 1]["diff"] > 0:
                current_min = pairings[index - 1]["diff"]
                min_diff = pairings[index - 1]
                print(f"New curren min is {current_min}")

        if current_min == 0:
            print ("THERE IS ONLY ONE MEAN VALUE LEFT")
            break

        for element in pairings:
            if element["diff"] > 0 and element["diff"] < current_min:
                current_min = element["diff"]
                min_diff = element

        print("min_diff")
        print(min_diff)

        if min_diff["diff"] > max_diff:
            print(f"The current difference of {min_diff['diff']} is higher than the maximal difference parameter of {max_diff}")
            break

        combine_values_1 = min_diff["elements"][0]["start"]
        combine_values_2 = min_diff["elements"][1]["start"]
        new_value = (min_diff["elements"][0]["mean"] + min_diff["elements"][1]["mean"])/2

        to_be_counified_values_1 = min_diff["elements"][0]["mean"]
        to_be_counified_values_2 = min_diff["elements"][1]["mean"]

        #{'elements': ({'mean': 7.476190476190476, 'start': 41}, {'mean': 7.375, 'start': 80}), 'diff': 0.10119047619047628}

        print("combine_values_1:")    
        print(combine_values_1)

        print("combine_values_2:")    
        print(combine_values_2)

        #MEANS:
        #[{'mean': 2.5609756097560976, 'start': 0}, {'mean': 7.476190476190476, 'start': 41}, {'mean': 2.0, 'start': 62}, {'mean': 7.375, 'start': 80}, {'mean': 1.8, 'start': 88}]

        print("MEAN BEFORE:")
        print(means)

        for element in means:
            if element["start"] == combine_values_1 or element["start"] == combine_values_2:
                element["mean"] = new_value
            if element["mean"] == to_be_counified_values_1 or element["mean"] == to_be_counified_values_2:
                element["mean"] = new_value
                
        print("MEANS AFTER")
        print(means)

        plt.plot(liste_x[1:], liste_y[1:])
        
        for element in means:
            if element["end"] != -1:
                plt.plot((element["start"], element["end"]), (element["mean"], element["mean"]), 1)
            else:
                plt.plot((element["start"], max(liste_x)), (element["mean"], element["mean"]), 1)


        #[{'mean': 2.5609756097560976, 'start': 0, 'end': 41}, {'mean': 7.425595238095238, 'start': 41, 'end': 62}, 
        #{'mean': 2.0, 'start': 62, 'end': 80}, {'mean': 7.425595238095238, 'start': 80, 'end': 88}, {'mean': 1.8, 'start': 88, 'end': -1}]        

        vertical_lines = []
        for element in transition_points:
            vertical_lines.append({"x": element["index"], "y_1": 0, "y_2":0 })

        for element in vertical_lines:
            for item in means:
                if item["start"] == element ["x"]:
                    element["y_2"] = item["mean"]
                if item["end"] == element ["x"]:
                    element["y_1"] = item["mean"]


        print("vertical_lines:")
        print(vertical_lines)

        for element in vertical_lines:
            if element ["x"] != 0:
                plt.plot((element["x"], element["x"]), (element["y_1"], element["y_2"]), 1)

        plt.show()

        number_of_cycles = number_of_cycles + 1