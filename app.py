from flask import Flask, render_template, request, redirect, session, url_for
import Peptidases
import requests as r
from Bio import SeqIO
from io import StringIO
import pandas as pd  
import textwrap


app = Flask(__name__)
app.secret_key = "hello"

@app.route('/')
def Home():
    return render_template("Home.html")

@app.route('/Documentation')
def About():
    return render_template("Documentation.html")

@app.route('/Table.html')
def Table():
    return render_template("Table.html")

@app.route('/Tool', methods=['GET', 'POST'])
def Tool():    
    if request.method == "POST":
        result = request.form["sequence"]
        length = request.form["length"]
        radio_options = request.form["radio_options"]
        checkbox = request.form.getlist("enzymes")

        session["result"] = result
        session["length"] = length
        session["radio_options"] = radio_options
        session["checkbox"] = checkbox
        return redirect(url_for("Results_cutter"))
    else:
        return render_template('Tool.html')

@app.route('/Sequentialcutter', methods=['GET', 'POST'])
def Sequentialcutter():
    if request.method == "POST":
        result = request.form["sequence"]
        length = request.form["length"]
        checkbox = request.form.getlist("enzymes")

        session["result"] = result
        session["length"] = length
        session["checkbox"] = checkbox
        return redirect(url_for("Results_sequential"))
    else:
        return render_template("Sequentialcutter.html")

@app.route('/Results')
def Results_cutter():
    """ 

    """
    sequence = session["result"]
    sequence = Seq_input(sequence)
    length_fragment = session["length"]
    length_fragment = Check_length(length_fragment)
    peptidase = session["checkbox"]
    table = Html_table(sequence, peptidase, length_fragment)

    wraptext=[]
    wrapper = textwrap.TextWrapper(width=100)
    words = wrapper.wrap(text=sequence)
    for element in words:
        wraptext.append(element)
    hardcut = Peptidases.Get_pattern_output(sequence, peptidase, length_fragment)
    # print(hardcut)
    listcut = []
    listempty = [] 
    Cutters_worked = []
    Cutters_notworked = []
    for key,value in hardcut.items():
        if not value:
            listempty.append(value)
            Cutters_notworked.append(key)
        else:
            listcut.append(value)
            Cutters_worked.append(key) 

    return render_template("Output_Simple_Cutter.html", sequence=wraptext, peptidase=peptidase, cuts=listcut, sequence_zip=zip(Cutters_worked, listcut), 
                            tables=[table.to_html(classes='data' , justify='center')], titles=table.columns.values,listcut=hardcut, not_cutlist=Cutters_notworked)


@app.route('/Results_sequential')
def Results_sequential():
    """ 

    """
    sequence = session["result"]
    sequence = Seq_input(sequence)
    length_fragment = session["length"]
    length_fragment = Check_length(length_fragment)
    peptidase = session["checkbox"]
    table = Html_table_sequential_cutter(sequence, peptidase, length_fragment)

    wraptext=[]
    wrapper = textwrap.TextWrapper(width=100)
    words = wrapper.wrap(text=sequence)
    for element in words:
        wraptext.append(element)
    hardcut = Peptidases.Get_pattern_output(sequence, peptidase, length_fragment)
    listcut = []
    listempty = [] 
    Cutters_worked = []
    Cutters_notworked = []
    
    for key,value in hardcut.items():
        if not value:
            listempty.append(value)
            Cutters_notworked.append(key)
        else:
            listcut.append(value)
            Cutters_worked.append(key) 
    return render_template("Output_Sequential_Cutter.html", sequence=wraptext, peptidase=peptidase, cuts=listcut, sequence_zip=zip(Cutters_worked, listcut), 
                            tables=[table.to_html(classes='data' , justify='center')], titles=table.columns.values,listcut=hardcut, not_cutlist=Cutters_notworked)

@app.route('/InputError')
def Uniprot_id(uID):
    """

    """
    uniprotUrl="http://www.uniprot.org/uniprot/"
    currentUrl=uniprotUrl+uID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))
    sequence = str(pSeq[0].seq)
    return sequence
        # RENDER TEMPLATE FOR INPUT ERRORS HANDLING

################################################################################################################################
def Check_length(length_option):
    if length_option == "":
        length_option = 1
    else:
        length_option = length_option
    return length_option

def Fasta_sequence(Fasta):
    for record in SeqIO.parse(StringIO(Fasta), "fasta"):
        sequence = record.seq
    return str(sequence)

def Seq_input(sequence):
    """

    """
    sequence=sequence.strip()
    if '>' in sequence:
        sequence = Fasta_sequence(sequence)
        return sequence
    elif sequence.isalpha() == False: #For now it may work but not always. We will have to add multiple checks
        sequence = Uniprot_id(sequence)
        return sequence
    else:
        sequence = sequence.upper()
    return sequence
        ## CHECK IF PROTEIN SEQUENCE OR NOT, IF NOT DO NOT RETURN THE SEQUENCE
##############################################################################
"""
REGULAR CUTTER PART!
"""
def Html_table(sequence, peptidase, length_option):
    """
    The function will create a dataframe with the data to render in the .to_html 
    """
    inputdict = {}
    inputdict = Peptidases.Get_pattern(sequence,peptidase,length_option)
    outputdf = pd.DataFrame(inputdict.items())
    outputdf.columns = ['Protease', 'Fragment']
    outputdf1 = outputdf.set_index(['Protease'])['Fragment'].apply(pd.Series).stack()
    outputdf1 = outputdf1.reset_index()

    #Fragment needs to be removed
    outputdf1.columns = ['Protease', 'Fragment', 'Fragment Sequence']
    outputdf1['Length of fragment'] = outputdf1['Fragment Sequence'].apply(len)

    return outputdf1
##############################################################################
"""
SEQUENTIAL CUTTER PART!
"""
def Html_table_sequential_cutter(sequence, peptidase, length_option):
    """
    The function will create a dataframe with the data to render in the .to_html 
    """
    inputdict = {}
    inputdict = Peptidases.Get_pattern_sequential(sequence,peptidase,length_option)
    outputdf = pd.DataFrame(inputdict)
    outputdf1 = outputdf
    outputdf1.columns = ['Fragment Sequence']
    outputdf1['Length of fragment'] = outputdf1['Fragment Sequence'].apply(len)
    return outputdf1
    

if __name__ == '__main__':
    from waitress import serve
    serve(app, host="0.0.0.0", port=8080)