import os


def getCutDict ( file_name ) :
    
    file = open ( file_name, "r" )
    file_contents = file.read()
    file.close() 
    
    raw_table = file_contents.split("DATA")[1].split("\n\n")[0].strip()
    column_names = raw_table.split("\n")[0].split()
    rows = raw_table.split("\n")
    
    d_cutName_cutData = {} 

    for row in rows[1:]:
        row = row.strip()
        row_fields = row.split() 
        cut_data = dict ( zip ( column_names, row_fields ) )
        cut_name = row_fields[0]
        d_cutName_cutData [ cut_name ] = cut_data 

    return d_cutName_cutData

dat_file_name     = os.environ["LQDATA"]+"/enujj_analysis/enujj/WZSherpa_scaled_output_cutTable_lq_enujj/analysisClass_lq_enujj_tables.dat"

dat_file = open ( dat_file_name, "r" ) 

d_cutName_cutTitle = { 
    "ST":"S_{T}^{e\\nu} > (GeV)",
    "MET":"E_{T}^{miss} > (GeV)",
    "Mej":"M(e,jet) > (GeV)" 
    }

cutNames = [ 
    "ST",
    "MET",
    "Mej" 
    ]


LQ_masses = [ 
    250, 
    350, 
    400, 
    450, 
    500, 
    550, 
    600, 
    650, 
    750, 
    850 
]
 
d_cutName_cutData = getCutDict ( dat_file_name ) 

latex_file_name = "table_selection_enujj.tex"
latex_file = open ( latex_file_name, "w" )  

latex_file.write("\documentclass{article}\n")
latex_file.write("\usepackage{multirow} \n")
latex_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
latex_file.write("\pagestyle{empty} \n")
latex_file.write("\\begin{document}\n")
latex_file.write("\\begin{table} \n")
latex_file.write("\\small \n")
line = "\\begin{tabular}{| l | "
for LQ_mass in LQ_masses:
    line = line + "c | "
line = line + "} \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
line = "$M_{LQ} (GeV)$ "
for LQ_mass in LQ_masses:
    line = line + " & " + str ( LQ_mass ) 
line = line + " \\\ \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
latex_file.write("  \hline \n")
for cutName in cutNames:
    line = "$ " + d_cutName_cutTitle [ cutName ] + " $ "
    for LQ_mass in LQ_masses:
        this_cutName  = cutName + "_LQ" + str ( LQ_mass ) 
        this_cutValue = d_cutName_cutData [ this_cutName ]["min1"]
        line = line + " & " + str ( int ( float ( ( this_cutValue ) ) ) )
    line = line + " \\\ \n" 
    latex_file.write ( line ) 
    latex_file.write("  \hline \n")
latex_file.write("\end{tabular}\n")
latex_file.write("\end{table}\n")
latex_file.write("\end{document}\n")
latex_file.close()

os.system ( "pdflatex " + latex_file_name ) 
os.system ( "rm *.aux *.log" )
                    
pdf_file_name = "/mnt/lxplus/scratch0/rootNtupleAnalyzer/CMSSW_4_2_3/src/rootNtupleAnalyzerV2/" + latex_file_name
pdf_file_name = pdf_file_name.replace (".tex", ".pdf" )

print "open " + pdf_file_name
