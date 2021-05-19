import re
import sys
import itertools


sys.tracebacklimit = None

def Singlecutter(input_sequence, cutting_pattern, desired_length=''):
    if ' ' in input_sequence == True:
        raise  Exception("Please provide a valid sequence. Spaces are not allowed in the input sequence")
    elif bool(re.search(r'[BJOUXZ]', input_sequence)) == True:
        raise  Exception("Please provide a valid sequence. BJOUXZ are not allowed in the input sequence")
    elif input_sequence.isalpha() == False:
        raise  Exception("Please provide a valid sequence. Digits are not allowed in the input sequence")
    else:
        p = re.compile(cutting_pattern)
        pointer = 0
        array = []
        list_of_input_fragments = []
        temp = ""
        for m in p.finditer(input_sequence):
            part = input_sequence[pointer:m.start()]
            array.append(part)
            pointer = m.start()
        array.append(input_sequence[pointer:len(input_sequence)])

        desired_length = int(desired_length or 1)
        fragmentcaller = FragmentCaller(array)
        fragmentcallercheck = FragmentCallerCheck(fragmentcaller, input_sequence)
        filteredFragments = Length_filter_sequential(fragmentcallercheck, desired_length)
        return filteredFragments
        # p = re.compile(cutting_pattern)
        # pointer = 0
        # array = []
        # list_of_input_fragments = []
        # temp = ""
        # for m in p.finditer(input_sequence):
        #     part = input_sequence[pointer:m.start()+len(m.group())]
        #     array.append(part)
        #     pointer = m.start()+len(m.group())
        #     split = re.split(cutting_pattern, part)
        #     list_of_input_fragments.append(temp + split[0] + split[1])
        #     temp = split[2]+split[3]
        # temp += input_sequence[pointer:len(input_sequence)]
        # if len(temp) > 0:
        #     list_of_input_fragments.append(temp)
        # desired_length = int(desired_length or 1)
        # filteredFragments = Length_filter(list_of_input_fragments, desired_length)
        # fragmentcaller = FragmentCaller(filteredFragments)
        # fragmentcallercheck = FragmentCallerCheck(fragmentcaller, input_sequence)
        # return fragmentcallercheck

def FragmentCaller(fragments):
    filteredFragments = []
    for idx, element in enumerate(fragments):
        index_within_list = idx < len(fragments)
        if index_within_list:
            start_joining = True
            counter = idx
            temp_string = element
            while counter < len(fragments):
                counter += 1
                try:
                    temp_string += fragments[counter]
                    filteredFragments.append(element) if element not in filteredFragments else None
                    filteredFragments.append(temp_string) if temp_string not in filteredFragments else None
                except IndexError as e:
                    pass
    else:
        filteredFragments.append(element) if element not in filteredFragments else None
    return filteredFragments

def FragmentCallerCheck(fragments, sequence):
    filteredFragmentsCheck = []
    for x in fragments:
        if len(sequence) > len(x):
            filteredFragmentsCheck.append(x) if x not in filteredFragmentsCheck else None
    return filteredFragmentsCheck

def Length_filter(list_of_input_fragments, desired_length):
    desired_length = int(desired_length)
    filteredFragments=[]
    for idx, element in enumerate(list_of_input_fragments):
        index_within_list = idx < len(list_of_input_fragments)
        last_index = idx == len(list_of_input_fragments) -1
        desired_length_not_reached = len(element) < desired_length
        if desired_length_not_reached and index_within_list:
            start_joining = True
            counter = idx
            temp_string = element
            while start_joining and counter < len(list_of_input_fragments):
                # print(len(list_of_input_fragments))
                counter += 1
                # print(counter)
                try:
                    temp_string += list_of_input_fragments[counter]
                    if len(temp_string) >= desired_length:
                        filteredFragments.append(temp_string)
                        filteredFragments.append(temp_string) if temp_string not in filteredFragments else filteredFragments
                        start_joining = False
                
                except IndexError as e:
                    start_joining = False
                    if desired_length_not_reached and last_index:
                        reverse_joining = True
                        reverse_counter = idx
                        reverse_string = element
                        while reverse_joining:
                            reverse_counter -= 1
                            reverse_string = list_of_input_fragments[reverse_counter]+reverse_string
                            if len(reverse_string) >= desired_length:
                                filteredFragments.append(reverse_string) if reverse_string not in filteredFragments else filteredFragments
                                reverse_joining = False
        else:
            filteredFragments.append(element) 
    return filteredFragments

def Get_pattern(input_sequence, list_of_proteases, desired_length=''):
    protease_fragments_pair = {}
    """
    Calling the dictionary from the Peptidases script for single cutter?
    """
    patterndict={
    "BACE1" : r'(?<=[EG][VIL].[LF])(.)(?=[AV].[VF])',
    "CPQ" : r'(?<=.{4})(F)(?=.{3})',
    "CTSB" : r'(?<=.{3}[RG])(.)(?=.{3})',
    "CTSC" : r'(?<=.{2}(?<!R|K)[SE])(?<!P)(?=.{3})',
    "CTSD" : r'(?<=.{3}[LFAVIPMW])(.)(?=.{3})',
    "CTSF" : r'(?<=.{2}[FLVK].)(.)(?=.{3})',
    "CTSH" : r'(?<=QEVR)(S)(?=PS.)',
    "CTSK" : r'(?<=.{2}[LIVPMF].)(.)(?=.{3})',
    "CTSL" : r'(?<=.{2}[LVFI].)(.)(?=.{3})',
    "CTSS" : r'(?<=.{2}[LV].)(.)(?=.{3})',
    "CTSV" : r'(?<=.{2}[LVI].)(.)(?=.{3})',
    "LGMN" : r'(?<=.{3}[ND])(.)(?=.{3})',
    "PRCP" : r'(?<=.{3}P)(.)(?=.{3})',
    "SPPL2A" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    "SPPL2B" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    "TPP1" : r'(?<=.{2}[GP][FMG])(F)(?=[RL].P)'
    }
    # patterndict={
    # "BACE1" : r'([EG][VIL].[LF])(.[AV].[VF])',
    # "CPQ" : r'(....)(F...)',
    # "CTSB" : r'(...[RG])(....)',
    # "CTSC" : r'(..(?<!R|K)[SE])((?<!P)...)',
    # "CTSD" : r'(...[LFAVIPMW])(....)',
    # "CTSF" : r'(..[FLVK].)(....)',
    # "CTSH" : r'(QEVR)(SPS.)',
    # "CTSK" : r'(..[LIVPMF].)(....)',
    # "CTSL" : r'(..[LVFI].)(....)',
    # "CTSS" : r'(..[LV].)(....)',
    # "CTSV" : r'(..[LVI].)(....)',
    # "LGMN" : r'(...[ND])(....)',
    # "PRCP" : r'(...P)(....)',
    # "SPPL2A" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    # "SPPL2B" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    # "TPP1" : r'(..[GP][FMG])(F[RL].P)'
    # }
    dictfilt_function = lambda x, y: dict([ (i,x[i]) for i in x if i in set(y) ])
    Selected_proteases = dictfilt_function(patterndict, list_of_proteases)
    for proteases, cutting_pattern in Selected_proteases.items():
        protease_fragments_pair[proteases] = Singlecutter(input_sequence, cutting_pattern, desired_length)
    return protease_fragments_pair

#######################
def Get_pattern_output(input_sequence, list_of_proteases, desired_length=''):
    protease_fragments_pair = {}
    """
    Calling the dictionary from the Peptidases script for single cutter?
    """
    patterndict={
    "BACE1" : r'(?<=[EG][VIL].[LF])(.)(?=[AV].[VF])',
    "CPQ" : r'(?<=.{4})(F)(?=.{3})',
    "CTSB" : r'(?<=.{3}[RG])(.)(?=.{3})',
    "CTSC" : r'(?<=.{2}(?<!R|K)[SE])(?<!P)(?=.{3})',
    "CTSD" : r'(?<=.{3}[LFAVIPMW])(.)(?=.{3})',
    "CTSF" : r'(?<=.{2}[FLVK].)(.)(?=.{3})',
    "CTSH" : r'(?<=QEVR)(S)(?=PS.)',
    "CTSK" : r'(?<=.{2}[LIVPMF].)(.)(?=.{3})',
    "CTSL" : r'(?<=.{2}[LVFI].)(.)(?=.{3})',
    "CTSS" : r'(?<=.{2}[LV].)(.)(?=.{3})',
    "CTSV" : r'(?<=.{2}[LVI].)(.)(?=.{3})',
    "LGMN" : r'(?<=.{3}[ND])(.)(?=.{3})',
    "PRCP" : r'(?<=.{3}P)(.)(?=.{3})',
    "SPPL2A" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    "SPPL2B" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    "TPP1" : r'(?<=.{2}[GP][FMG])(F)(?=[RL].P)'
    }
    # patterndict={
    # "BACE1" : r'([EG][VIL].[LF])(.[AV].[VF])',
    # "CPQ" : r'(....)(F...)',
    # "CTSB" : r'(...[RG])(....)',
    # "CTSC" : r'(..(?<!R|K)[SE])((?<!P)...)',
    # "CTSD" : r'(...[LFAVIPMW])(....)',
    # "CTSF" : r'(..[FLVK].)(....)',
    # "CTSH" : r'(QEVR)(SPS.)',
    # "CTSK" : r'(..[LIVPMF].)(....)',
    # "CTSL" : r'(..[LVFI].)(....)',
    # "CTSS" : r'(..[LV].)(....)',
    # "CTSV" : r'(..[LVI].)(....)',
    # "LGMN" : r'(...[ND])(....)',
    # "PRCP" : r'(...P)(....)',
    # "SPPL2A" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    # "SPPL2B" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    # "TPP1" : r'(..[GP][FMG])(F[RL].P)'
    # }
    dictfilt_function = lambda x, y: dict([ (i,x[i]) for i in x if i in set(y) ])
    Selected_proteases = dictfilt_function(patterndict, list_of_proteases)
    for proteases, cutting_pattern in Selected_proteases.items():
        protease_fragments_pair[proteases] = Sequence_output_cut(input_sequence, cutting_pattern, desired_length)
    return protease_fragments_pair


def Sequence_output_cut(input_sequence, cutting_pattern, desired_length=''):
    if ' ' in input_sequence == True:
        raise  Exception("Please provide a valid sequence. Spaces are not allowed in the input sequence")
    elif bool(re.search(r'[BJOUXZ]', input_sequence)) == True:
        raise  Exception("Please provide a valid sequence. BJOUXZ are not allowed in the input sequence")
    elif input_sequence.isalpha() == False:
        raise  Exception("Please provide a valid sequence. Digits are not allowed in the input sequence")
    else:
        p = re.compile(cutting_pattern)
        pointer = 0
        array = []
        list_of_input_fragments = []
        temp = ""
        for m in p.finditer(input_sequence):
            part = input_sequence[pointer:m.start()]
            array.append(part)
            pointer = m.start()
        array.append(input_sequence[pointer:len(input_sequence)])

        fragmentcallercheck = FragmentCallerCheck(array, input_sequence)

        return fragmentcallercheck
        # p = re.compile(cutting_pattern)
        # pointer = 0
        # array = []
        # list_of_input_fragments = []
        # temp = ""
        # for m in p.finditer(input_sequence):
        #     part = input_sequence[pointer:m.start()+len(m.group())]
        #     array.append(part)
        #     pointer = m.start()+len(m.group())
        #     split = re.split(cutting_pattern, part)
        #     list_of_input_fragments.append(temp + split[0] + split[1])
        #     temp = split[2]+split[3]
        # temp += input_sequence[pointer:len(input_sequence)]
        # if len(temp) > 0:
        #     list_of_input_fragments.append(temp)
        # # print(list_of_input_fragments)

        # # fragmentcallercheck = FragmentCallerCheck(array, input_sequence)
        # return list_of_input_fragments

############################################################################################################################################

def Array_pattern_sequential(input, cutting_pattern, desired_length=''):
    if isinstance(input, str):
        input = [input]
    output = []
    for input_sequence in input:
        output = output + Fragment_forming_sequential(input_sequence, cutting_pattern, desired_length)
    return output

def Get_pattern_sequential(input_sequence, list_of_proteases, desired_length=''):
    output = input_sequence
    # patterndict={
    # "BACE1" : r'(?<=[EG][VIL].[LF])(.)(?=[AV].[VF])',
    # "CPQ" : r'(?<=.{4})(F)(?=.{3})',
    # "CTSB" : r'(?<=.{3}[RG])(.)(?=.{3})',
    # "CTSC" : r'(?<=.{2}(?<!R|K)[SE])(?<!P)(?=.{3})',
    # "CTSD" : r'(?<=.{3}[LFAVIPMW])(.)(?=.{3})',
    # "CTSF" : r'(?<=.{2}[FLVK].)(.)(?=.{3})',
    # "CTSH" : r'(?<=QEVR)(S)(?=PS.)',
    # "CTSK" : r'(?<=.{2}[LIVPMF].)(.)(?=.{3})',
    # "CTSL" : r'(?<=.{2}[LVFI].)(.)(?=.{3})',
    # "CTSS" : r'(?<=.{2}[LV].)(.)(?=.{3})',
    # "CTSV" : r'(?<=.{2}[LVI].)(.)(?=.{3})',
    # "LGMN" : r'(?<=.{3}[ND])(.)(?=.{3})',
    # "PRCP" : r'(?<=.{3}P)(.)(?=.{3})',
    # "SPPL2A" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    # "SPPL2B" : r'(?<=[LFT][SFLC][LF][FSLC])([SLIH])(?=[FVL][LSAHG][IFGV])',
    # "TPP1" : r'(?<=.{2}[GP][FMG])(F)(?=[RL].P)'
    # }
    patterndict={
    "BACE1" : r'([EG][VIL].[LF])(.[AV].[VF])',
    "CPQ" : r'(....)(F...)',
    "CTSB" : r'(...[RG])(....)',
    "CTSC" : r'(..(?<!R|K)[SE])((?<!P)...)',
    "CTSD" : r'(...[LFAVIPMW])(....)',
    "CTSF" : r'(..[FLVK].)(....)',
    "CTSH" : r'(QEVR)(SPS.)',
    "CTSK" : r'(..[LIVPMF].)(....)',
    "CTSL" : r'(..[LVFI].)(....)',
    "CTSS" : r'(..[LV].)(....)',
    "CTSV" : r'(..[LVI].)(....)',
    "LGMN" : r'(...[ND])(....)',
    "PRCP" : r'(...P)(....)',
    "SPPL2A" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    "SPPL2B" : r'([LFT][SFLC][LF][FSLC])([SLIH][FVL][LSAHG][IFGV])',
    "TPP1" : r'(..[GP][FMG])(F[RL].P)'
    }
    pattern = []
    emptylist = []
    notempty = []
    #This is the order of entering
    for x in list_of_proteases:
        pattern.append(patterndict.get(x))
    not_none_values = filter(None.__ne__, pattern)
    list_not_none = list(not_none_values)
    for x in list_not_none:
        output = Array_pattern_sequential(output, x, desired_length)
        if not output:
            emptylist.append(output) 
        else:
            notempty.append(output)
    for x in notempty:
        fragmentcallercheck = FragmentCallerCheck(x, input_sequence)
    lengthfilter = Length_filter_sequential(fragmentcallercheck, desired_length)
    return lengthfilter

def Fragment_forming_sequential(input_sequence, cutting_pattern, desired_length=''):
    # p = re.compile(cutting_pattern)
    # pointer = 0
    # array = []
    # list_of_input_fragments = []
    # temp = ""
    # for m in p.finditer(input_sequence):
    #     part = input_sequence[pointer:m.start()]
    #     array.append(part)
    #     pointer = m.start()
    # array.append(input_sequence[pointer:len(input_sequence)])
    # desired_length = int(desired_length or 1)
    # fragmentcaller = FragmentCaller(array)
    # fragmentcallercheck = FragmentCallerCheck(fragmentcaller, input_sequence)
    # filteredFragments = Length_filter_sequential(fragmentcallercheck, desired_length)
    # return filteredFragments

    p = re.compile(cutting_pattern)
    pointer = 0
    array = []
    list_of_input_fragments = []
    temp = ""
    for m in p.finditer(input_sequence):
        part = input_sequence[pointer:m.start()+len(m.group())]
        array.append(part)
        pointer = m.start()+len(m.group())
        split = re.split(cutting_pattern, part)
        list_of_input_fragments.append(temp + split[0] + split[1])
        temp = split[2]+split[3]
    temp += input_sequence[pointer:len(input_sequence)]
    if len(temp) > 0:
        list_of_input_fragments.append(temp)
    desired_length = int(desired_length or 1)
    filteredFragments = Length_filter(list_of_input_fragments, desired_length)
    fragmentcaller = FragmentCaller(filteredFragments)
    fragmentcallercheck = FragmentCallerCheck(fragmentcaller, input_sequence)
    return fragmentcallercheck

def Length_filter_sequential(list_of_input_fragments, desired_length):
    """
    To remove the fragments after all the output is formed for the chosen protease(s). A list of fragments > desired length is returned.
    """
    desired_length = int(desired_length)
    length_list = []
    for x in list_of_input_fragments:
        if len(x) > desired_length:
            length_list.append(x)
    return length_list

