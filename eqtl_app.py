from flask import Flask, render_template, request, redirect
from subprocess import Popen, PIPE
import sys

app=Flask(__name__)
app.config['SECRET_KEY'] = '#DtfrL98G5t1dC*4'

@app.route("/")
def root():
    return render_template('index.html')

@app.route("/eqtl")
def eqtl():
    return render_template('eqtl.html')

@app.route("/plots", methods=['GET', 'POST'])
def plots():
    if request.method == 'POST':
        term = request.form
        chr, bp=term["loc"].split(":")
        process=Popen(['bash', 'export_eqtl_data.sh', chr, bp], stdout=PIPE, stderr=PIPE)
        stdout, stderr=process.communicate()
        #print (bp, file=sys.stderr)
        return render_template("plots.html", loc=term, data=stdout.split())

if __name__ == '__main__':
    app.run(debug=True)

