from flask import Flask, render_template, request, redirect
from subprocess import Popen, PIPE
import json 
import pandas as pd
import numpy as np

app=Flask(__name__)
app.config['SECRET_KEY'] = '#DtfrL98G5t1dC*4'

@app.route("/")
def root():
    return render_template('index.html')

@app.route("/eqtl")
def eqtl():
    return render_template('eqtl.html')

@app.route("/scatter", methods=['GET', 'POST'])
def plots():
    if request.method == 'POST':
        term = request.form
        chr, bp=term["loc"].split(":")
        process=Popen(['bash', 'export_eqtl_data.sh', chr, bp], stdout=PIPE, stderr=PIPE)
        stdout, stderr=process.communicate()
        file=str(chr)+"_"+(bp)+"_"+"win"+str(2000000)+".tab"
        df=pd.read_csv(file, delimiter="\t")
        dfplot=df.iloc[:,3:5]
        dfplot.columns=['bp', 'p']
        dfplot['logp'] = -np.log10(dfplot['p'])
        chart_data=dfplot.to_dict(orient='records')
        chart_data=json.dumps(chart_data,indent=2)
        data={'chart_data':chart_data}
        return render_template("d3.html", data=data, loc=term)
if __name__ == '__main__':
    app.run(debug=True)

