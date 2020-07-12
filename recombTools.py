## file recombTools.py
## maps of recombination frequencies obtained from http://www.stanford.edu/~lipatov/recombination/methods.html

import os, string, MySQLdb, cPickle, tkFont
from Tkinter import *

#import psyco ## comment out these lines
#psyco.full() ## if you don't want to bother with psyco

def main():
    app = recombulator()
    app.master.title("Recombulator")
    app.mainloop()

class recombulator(Frame):              
    def __init__(self, master=None):
        Frame.__init__(self, master)   
        self.grid()                    
        self.createWidgets()
        self.answercoord = 0
    def createWidgets(self):
##         self.menubutton = Button(self)
##         menu = Menu(self.menubutton)
##         menu.add_cascade(label="quit")
        self.gene1box = Entry(self, bg="white")
        self.gene1box.insert(0,"jeb")
        self.gene1box.grid()
        self.gene2box = Entry(self,bg="white")
        self.gene2box.grid()
        self.gene2box.insert(0,"2L:100000..200000")
        self.fireButton = Button(self, text = "calculate",
                                 command=self.calculate)
        self.fireButton.grid()
        self.answerbox = Canvas(self,width = 150, height=80)
#        self.answerbox.create_bitmap(0,0,bitmap=@"recomLogo.xbm")
        self.answerbox.grid()
        self.quitButton = Button(self, text="Quit",
            command=self.quit)        
        self.quitButton.grid()
    def calculate(self):
        gene1string = self.gene1box.get()
        gene2string = self.gene2box.get()
        if ('..' and ':') in gene1string:
            gene1 = interval(gene1string)
            gene1string = 'intvl'
        else: gene1 = locus(gene1string)
        if ('..' and ':') in gene2string:
            gene2 = interval(gene2string)
            gene2string = 'intvl'
        else: gene2 = locus(gene2string)
        answer = predictRF(gene1, gene2)*float(100)
        self.answerbox.create_text(5,self.answercoord+10, anchor= W,text="%0s-%0s: %%%.5f" % (gene1string, gene2string, answer), font=tkFont.Font(weight="bold"))
        self.setAnswerCoord()
    def setAnswerCoord(self):
        self.answercoord = self.answercoord + 20

def calculateX(physical_location):
    "takes nt position in units of Megabases as argument"
    x = physical_location
    genetic_map = -0.0097*x**3 + 0.2996*x**2 + 1.1626*x - 1.8904
    return genetic_map

def calculate2L(physical_location):
    "takes nt position in units of Megabases as argument"
    x = physical_location
    genetic_map = -0.0099*x**3 + 0.2087*x**2 + 2.5183*x - 0.8057
    return genetic_map

def calculate2R(physical_location):
    "takes nt position in units of Megabases as argument"
    x = physical_location
    genetic_map = -0.0083*x**3 + 0.3627*x**2 - 1.5224*x + 58.252
    return genetic_map

def calculate3L(physical_location):
    "takes nt position in units of Megabases as argument"
    x = physical_location
    genetic_map = -0.006*x**3 + 0.1091*x**2 + 2.6663*x -1.6899
    return genetic_map

def calculate3R(physical_location):
    "takes nt position in units of Megabases as argument"
    x = physical_location
    genetic_map = -0.0038*x**3 + 0.233*x**2 - 1.5567*x + 50.127
    return genetic_map

def mapPosition(location_nt, chromosome_arm):
    "This function returns a genetic map distance coordinate for any specific nucleotide location in the genome of *Drosophila melanogaster*. It does not handle chromosome IV."
    calculator_chooser = {'X':calculateX,'2L':calculate2L,'2R':calculate2R,'3L':calculate3L,'3R':calculate3R,'x':calculateX,'2l':calculate2L,'2r':calculate2R,'3l':calculate3L,'3r':calculate3R}
    calculate_position = calculator_chooser[chromosome_arm]
    map_position = calculate_position(location_nt/float(1e6))
    return map_position

class interval:
    "Interval class handles basic info about any genetic element that can be described as an 'interval', i.e. the sequence contained between specified nucleotide posistions in a genome. It was designed for use with *Drosophila* but is generic enough for use with any genome."
    def __init__(self,location,mapType='physical'):
        "location should be formatted like: '3R:20857..39123', or '2L:complement(50213..55480)'."
        genomic_coordinates = string.split(location,':')
        self.type = mapType
        self.linkage_group = genomic_coordinates[0]
        coord = genomic_coordinates[1]
#        print(location)
        if 'complement' in coord:
            coord = coord[11:-1]
        self.breakpoints = string.split(coord,'..')
        self.name = 'none'
        
    def proximalBreakpoint(self):
        return eval(self.breakpoints[0])

    def distalBreakpoint(self):
        return eval(self.breakpoints[1])

    def getLinkageGroup(self):
        return self.linkage_group

    def mapProximal(self):
        return mapPosition(self.proximalBreakpoint(),self.getLinkageGroup())

    def mapDistal(self):
        return mapPosition(self.distalBreakpoint(),self.getLinkageGroup())

    def setName(self,name):
        self.name = name

    def Name(self):
        return self.name


def mapDistance(locus1,arm1,locus2,arm2):
    distance = abs(mapPosition(locus1,arm1)-mapPosition(locus2,arm2))
    return distance

def estimateRF(distance):
    "estimate the recombinant fraction given two genetic map loci"
    if distance < 20:
        rf = distance/100
    else:
        ## adapted from Ashburner "Drosophila A Laboratory Handbook", pp 458
        ## see Kosambi's formula
        y = estimateRF(distance/2)
        rf = 2*y/(1+4*y**2)
    return rf


def locus(stable_id,filelocation='geneheader.txt'):
    "Creates an instance of class *interval* corresponding to the given gene name."
    location = os.popen2("grep -e 'name=%0s;' %0s | grep -o -e 'loc[=][0-9]*[A-Z]*[:][complement(]*[0-9]*[.][.][0-9]*[)]*[;]'" % (stable_id, filelocation))[1].readline()[4:-2]
    location = string.replace(location,'complement\(', '')
    location = string.replace(location,'\)','')
    linkageGroup = string.split(location,':')[0]
    new_locus = interval(location)
    new_locus.setName(stable_id)
    del(location)
    return new_locus
        
def predictRF(locus1,locus2):
    "Estimates recombination frequency between two intervals (class interval) associated with genes or deficiencies."
    if locus1.getLinkageGroup()[0]==locus2.getLinkageGroup()[0]: #indexing preserves calculations for different arms of the same chrom.
        locus1p = locus1.mapProximal()
        locus1d = locus1.mapDistal()
        locus2p = locus2.mapProximal()
        locus2d = locus2.mapDistal()
        a = abs(locus1p-locus2d)
        b = abs(locus1d-locus2p)
        if a < b:
            return estimateRF(a)
        else:
            return estimateRF(b)
    else: return 0.5

def alkrecomb(defcyFile = '/home/hazmat/Desktop/drosdeldfweiss.csv'):
    outfile = open('recombmap.csv','w')
    defcys = open(defcyFile,'r').readlines()
    defcys = [string.split(df,',') for df in defcys]
    chrom2df = []
    alk = locus('Alk')
    ## >>> defcys[0] ##example
    ## ['#1', '9052', '"Df(1)ED6396', ' P{w[+mW.Scer\\FRT.hs3]=3\'.RS5+3.3\'}ED6396 w[1118]/FM7h"', '"01B05;1B8 (R4 estimated cytology)', ' X:350351..380452 (R4)"\n']
    outfile.write("num,stock,genotype,cytol,physloc,rf\n")
    for df in defcys:
        num,stock = df[0][1:],df[1]
        genotype = 'none'
        for i in range(2,len(df)):
            if 'Df(' in df[i]:
                genotype = df[i][1:]
##             if 'Scer' in df[i]:
##                 full_genotype = df[i][1:-1]
            if 'cytology' in df[i]:
                cytol = string.split(df[i])[0][1:]
            if '..' in df[i]:
                physloc = string.split(df[i])[0]
        if '2' == physloc[0]:
            outfile.write("%0s,%0s,%0s,%0s,%0s,%0f\n" % (num,stock,genotype,cytol,physloc,predictRF(alk,interval(physloc))))
    outfile.close()

def makeAllGenes():
    genenamedata = os.popen2("grep -e '[>]' ~/data/DrosophilaGenome/dmel-all-gene-r4.3.fasta | grep -o -e 'Name[=][^;]*'")[1].readlines()
    genenamedata = [string.replace(i, 'Name=', '')[:-1] for i in genenamedata]
    print("%d gene entries, " % len(genenamedata)),
    non_redundant = []
    for label in genenamedata:
        if not label in non_redundant:
            non_redundant.append(label)
    del(genenamedata)
    for i in range(len(non_redundant)-1,-1,-1):
        if ('[' or ']') in non_redundant[i]:
            non_redundant.pop(i)    
    genes = [locus(i) for i in non_redundant]
    del(non_redundant)
    return genes

## allgenes = makeAllGenes()
## tmpoutfile = open('allgene.pickle','w')
## cPickle.dump(allgenes,tmpoutfile)
## tmpoutfile.close()

    
def mapGeneOverlap(df,list_of_genes):
    "Computes a list of genes on or near any deficiency (df). Requires defcy from defmap.py"
    dfpm = df.physical_location()
    dfi = interval("%s:%s..%s" % (df.chromosome(),dfpm[0],dfpm[1]))
    positiveOverlaps = []
    for gene in list_of_genes:
        if dfi.getLinkageGroup()==gene.getLinkageGroup():
            if (dfi.proximalBreakpoint()<=gene.proximalBreakpoint())&(dfi.distalBreakpoint()>=gene.proximalBreakpoint()):
                assoc_type='direct overlap'
                positiveOverlaps.append((gene.Name(),assoc_type))
            elif (dfi.proximalBreakpoint()<=gene.distalBreakpoint())&(dfi.distalBreakpoint()>=gene.distalBreakpoint()):
                assoc_type='direct overlap'
                positiveOverlaps.append((gene.Name(),assoc_type))
            elif abs(dfi.proximalBreakpoint()-gene.distalBreakpoint())<=10000:
                assoc_type='within 10kb'
                positiveOverlaps.append((gene.Name(),assoc_type))
            elif abs(gene.proximalBreakpoint()-dfi.distalBreakpoint())<=10000:
                assoc_type='within 10kb'
                positiveOverlaps.append((gene.Name(),assoc_type))
            else:
                continue
        else:
            continue
    return positiveOverlaps

def buildOverlapDB(list_of_dfs, list_of_genes, start_from=0,end_at=295,drop_table=False):
    "list_of_dfs==list of elements class defcy (from defmap.py)"
    db_cnx = MySQLdb.connect(host='localhost',passwd='Och0Bush',db='alkscreen')
    curs = db_cnx.cursor()
    if drop_table==True:
        curs.execute("DROP TABLE IF EXISTS dfcy_genes;")
        curs.execute("CREATE TABLE dfcy_genes (defcy_strain VARCHAR(25), gene_stable_id VARCHAR(15), association VARCHAR(25));")
    for df in list_of_dfs[start_from:end_at]:
        associated_genes = mapGeneOverlap(df,list_of_genes)
        print(df.label()),
        temporary_geneinfo = []
        for gene in associated_genes:
            genename,overlapclass = gene[0],gene[1]
            if "'" in genename:
                genename = string.replace(genename,"'","")
            entry = "\'%s\', \'%s\', \'%s\'"%(df.label(),genename,overlapclass)
            if not (entry in temporary_geneinfo):
                temporary_geneinfo.append(entry)
        for eachline in temporary_geneinfo:
            sql_cmd = "INSERT INTO dfcy_genes VALUES(%s);" % eachline
            curs.execute(sql_cmd)
        del(temporary_geneinfo,associated_genes)
        print('..finished')
    db_cnx.close()
    return None

if __name__=="__main__": main()
