
"""
Created on Wed May 21 17:34:14 2014
@author: Wei He, Zhi Liang
@Email: hwkobe.1027@gmail.com
mFASD: A structure-based algorithm for discriminating different types of metal binding sites
"""


from __future__ import division
from Bio.PDB.PDBParser import PDBParser
import os, argparse, logging, sys



def GetLigand(PDBDir,PDBCode,Ligand):
    try:
       p = PDBParser(PERMISSIVE=1)
       PDBFileName = os.path.join(PDBDir,PDBCode)
       s = p.get_structure(PDBCode,PDBFileName)
    
       chain = s[0]['A']
       for residue in chain.get_residues():
           residue_id = residue.get_id()
           hetfield = residue_id[0]
           if hetfield == 'H_'+ Ligand:
               return residue
    except Exception,e:
         print e
     


def GetFunctionatomsFromResidue(PDBDir,PDBCode,ResidueList):
    atom_dict = {}
  
    try:
       p = PDBParser(PERMISSIVE=1)
       PDBFileName = os.path.join(PDBDir,PDBCode)
       s = p.get_structure(PDBCode,PDBFileName)
    except Exception:
        pass
    
    for j in ResidueList:
        re = s[0]['A'][j]
        for atom in re.get_list():
            if atom.get_name() not in ['C','O','N','CA']:
                
               atom_dict[atom] = atom_type[re.get_resname()][atom.get_name()]
    return atom_dict
    


def GetAtomInteractDictFromResidue(PDBDir,PDBCode,ResidueList):
    atom_dicts = GetFunctionatomsFromResidue(PDBDir,PDBCode,ResidueList)
    interaction_dict_all = {}
    value_dict =  {'1~1':0,'1~2':0,'1~3':0,'1~4':0,'1~5':0,'1~6':0,'2~1':0,'2~2':0,'2~3':0,'2~4':0,'2~5':0,
                   '2~6':0,'3~1':0,'3~2':0,'3~3':0,'3~4':0,'3~5':0,'3~6':0,'4~1':0,'4~2':0,'4~3':0,'4~4':0,
                   '4~5':0,'4~6':0,'5~1':0,'5~2':0,'5~3':0,'5~4':0,'5~5':0,'5~6':0,'6~1':0,'6~2':0,'6~3':0,
                   '6~4':0,'6~5':0,'6~6':0}
   
    
    newdict = atom_dicts
    
    for i in sorted(atom_dicts):
        for j in sorted(newdict):
            if 0 < j-i < 5:
               value_dict[str(atom_dicts[i])+'~'+str(newdict[j])] = value_dict[str(atom_dicts[i])+'~'+str(newdict[j])] + 1
        interaction_dict_all[i] = value_dict
        value_dict = {'1~1':0,'1~2':0,'1~3':0,'1~4':0,'1~5':0,'1~6':0,'2~1':0,'2~2':0,'2~3':0,'2~4':0,'2~5':0,
                      '2~6':0,'3~1':0,'3~2':0,'3~3':0,'3~4':0,'3~5':0,'3~6':0,'4~1':0,'4~2':0,'4~3':0,'4~4':0,
                      '4~5':0,'4~6':0,'5~1':0,'5~2':0,'5~3':0,'5~4':0,'5~5':0,'5~6':0,'6~1':0,'6~2':0,'6~3':0,
                      '6~4':0,'6~5':0,'6~6':0}
   
    return interaction_dict_all
    
         
    
def Get_Functionatoms_Dict(PDBDir,PDBCode,Ligand):
    atom_dict = {}
    re = GetLigand(PDBDir,PDBCode,Ligand)
    
    try:
       atoms1 = re.get_parent().get_atoms()
    
       for atom1 in atoms1:
          for atom2 in re:
             if atom1-atom2 < 5 :
               try:
                  atom_dict[atom1] = atom_type[atom1.get_parent().get_resname()][atom1.get_name()]
               except KeyError:
                  pass
    except Exception:
           pass              
    return atom_dict    
    


def GetAtomInteractDict(PDBDir,PDBCode,Ligand):    
    interaction_dict_all = {}
    value_dict =  {'1~1':0,'1~2':0,'1~3':0,'1~4':0,'1~5':0,'1~6':0,'2~1':0,'2~2':0,'2~3':0,'2~4':0,'2~5':0,'2~6':0,
                   '3~1':0,'3~2':0,'3~3':0,'3~4':0,'3~5':0,'3~6':0,'4~1':0,'4~2':0,'4~3':0,'4~4':0,'4~5':0,'4~6':0,
                   '5~1':0,'5~2':0,'5~3':0,'5~4':0,'5~5':0,'5~6':0,'6~1':0,'6~2':0,'6~3':0,'6~4':0,'6~5':0,'6~6':0}
   
    atom_dicts = Get_Functionatoms_Dict(PDBDir,PDBCode,Ligand)
    newdict = atom_dicts
    
    for i in sorted(atom_dicts):
        for j in sorted(newdict):
            if 0 < j-i < 5:
               value_dict[str(atom_dicts[i])+'~'+str(newdict[j])] = value_dict[str(atom_dicts[i])+'~'+str(newdict[j])] + 1
        interaction_dict_all[i] = value_dict
        value_dict = {'1~1':0,'1~2':0,'1~3':0,'1~4':0,'1~5':0,'1~6':0,'2~1':0,'2~2':0,'2~3':0,'2~4':0,'2~5':0,'2~6':0,
                      '3~1':0,'3~2':0,'3~3':0,'3~4':0,'3~5':0,'3~6':0,'4~1':0,'4~2':0,'4~3':0,'4~4':0,'4~5':0,'4~6':0,
                      '5~1':0,'5~2':0,'5~3':0,'5~4':0,'5~5':0,'5~6':0,'6~1':0,'6~2':0,'6~3':0,'6~4':0,'6~5':0,'6~6':0}
   
    return interaction_dict_all
    


def CalculateDistance(PDBDir,Query,ResidueList,Reference,ReferenceLigand):  
    ReferenceDict = GetAtomInteractDict(PDBDir,Reference,ReferenceLigand)
    QueryDict = GetAtomInteractDictFromResidue(PDBDir,Query,ResidueList)
    f = len(QueryDict)
    dist_list = []
    Distance = 0
    
    for AAq in sorted(QueryDict):
       for AAr in sorted(ReferenceDict):
            a = 0
            if AAq.get_name() == AAr.get_name():
               a = 0
            else: a = 1.0
            c = 0
            b = 0
            for AAIq in sorted(QueryDict[AAq]):
                for AAIr in sorted(ReferenceDict[AAr]):
                    if AAIq == AAIr:
                      
                       b = b + abs(QueryDict[AAq][AAIq]-ReferenceDict[AAr][AAIr])
                       if QueryDict[AAq][AAIq] >= ReferenceDict[AAr][AAIr]:
                           c = c + QueryDict[AAq][AAIq]
                       else: c = c + ReferenceDict[AAr][AAIr]
                           
            try:           
              d = 0.8*a + 0.2*(b/c)
              dist_list.append(d)
            except Exception:
               pass
       #print dist_list
       try:        
          Distance = Distance + min(dist_list)
       except Exception:
           logging.warning('No metal binding sites is found in reference PDB! Please check whether PDB is in the PDBDir!')
       dist_list = []
          
    if f!= 0:   
       return Distance/f
    else:
       return 1.0
             

def ChangePDBCodeToList(LigandFileName):
    PDBList = []
    for ln in file(LigandFileName,'rt'):
        PDBList.extend(ln.strip().split())
    return PDBList           


if __name__ == '__main__': 
      
     atom_type = {'ALA':{'O':2, 'C':6, 'CA':6, 'N':3,'CB':4},
                  'GLY':{'O':2, 'C':6, 'CA':6, 'N':3 },
                  'PRO':{'O':2, 'C':6, 'CA':6, 'CB':4,'N':3, 'CG':4, 'CD':4},
                  'ASN':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':6, 'OD1':2, 'ND2':3},
                  'ASP':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':6, 'OD1':2, 'OD2':1},
                  'PHE':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':5, 'CD1':5, 'CD2':5, 'CE1':5, 'CE2':5, 'CZ':5},
                  'LYS':{'O':2, 'C':6, 'CA':6, 'N':3, 'CG':4, 'CD':4, 'CE':6, 'NZ':3},
                  'ILE':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG1':4, 'CG2':4, 'CD1':4},
                  'LEU':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':4, 'CD1': 4, 'CD2':4},
                  'ARG':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':4, 'CD':6, 'NE':3, 'CZ':4, 'NH1':3, 'NH2':1},
                  'CYS':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':6, 'SG':6},
                  'MET':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':6, 'SG':6, 'CE':6},
                  'THR':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':6, 'CG2':4, 'OG1':1},
                  'TYR':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':5, 'CD1':5, 'CD2':5,'CE1':5, 'CE2':5, 'CZ':5, 'OH':1},
                  'HIS':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':5, 'ND1':3, 'CD2':5, 'CE1':5, 'NE2':2},
                  'VAL':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG1':4, 'CG2':4},
                  'SER':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'OG':1},
                  'GLU':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':3, 'CD':6, 'OE1':2,'OE2':1},
                  'GLN':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':3, 'CD':6, 'OE1':2,'NE2':3},
                  'TRP':{'O':2, 'C':6, 'CA':6, 'N':3, 'CB':4, 'CG':6, 'CD1':4, 'CD2':5, 'NE1':3, 'CE2':5, 'CE3':5,'CZ2':5, 'CZ3':5, 'CH2':5}}
     
     
     ## Set logging format
     logging.basicConfig(level=logging.DEBUG,  
                    format='%(levelname)s:%(asctime)s @%(message)s',  
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filemode='a')
     
     ## Add arguments for user input with command line
     parser = argparse.ArgumentParser(description='mFASD:Structure-based algorithm for discriminating different types of metal binding sites.')
     
     parser.add_argument('-i','--input',type=str,help='PDB codes of inqurey protein structures for metal binding sites prediction, if a structure has no PDB code yet, user can assign a temporal one',required=True)
     parser.add_argument('-m','--metal_type',type=str,help='The metal type for which the user want to predict for given structures, mFASD now support prediction for CA,CU,FE,MN,MG,ZN',required=True)
     parser.add_argument('-r','--residue_list',type=str,help='The residue list of certain protein region for metal prediction,eg:20,47,214',required=True)
     parser.add_argument('-t1','--threshold1',type=float,help='The threshold to determine whether two metal binding sites are similar',default=0.3)
     parser.add_argument('-t2','--threshold2',type=float,help='The threshold to determine the majority of votes',default=0.5)
     parser.add_argument('-d','--pdb_directory',type=str,help='The directory where the input protein structure are stored',default='PDBFiles')
     parser.add_argument('-o','--outputdir',type=str,help='The output directory to save the result files',default='Output')
     
     args = parser.parse_args()
     metal = args.metal_type.upper()
     QueryPDBCode = args.input
     ReferenceSet = ChangePDBCodeToList(os.path.join(os.getcwd(),metal+'_reference'))
     ResidueList = [int(i) for i in args.residue_list.split(',')]
     t1 = args.threshold1
     t2 = args.threshold2
     PDBDir = os.path.join(os.getcwd(),args.pdb_directory)
      
     ## To generate the outputdir for saving result files
     try:
        os.mkdir(args.outputdir)
        logging.info('Creat the outputdir {} to place result files'.format(args.outputdir))
     except OSError:
        logging.warning('outputdir {} already exist'.format(args.outputdir))
     
     if QueryPDBCode not in os.listdir(PDBDir):
         logging.error('No protein structure named %s in the PDBDir, please save your qurey PDB in the PDBDir!'%QueryPDBCode)
         sys.exit(-1)
         
     for PDB in ReferenceSet:
         if PDB not in os.listdir(PDBDir):
            logging.warning('Reference protein structure %s is not available, Download and save it to PDBDir!'%PDB)
            os.system('curl https://files.rcsb.org/download/'+PDB+' -o '+PDBDir+'/'+PDB)
       
     logging.info('Predicting '+metal+' binding sites for '+QueryPDBCode+' at region:'+str(ResidueList))
     outputdir = os.path.join(os.getcwd(),args.outputdir) ## Set the path to outputdir
     
     d = 0
     f = open(os.path.join(outputdir,'Prediction_'+QueryPDBCode+'_'+metal+'.txt'),'w')
     f.write('  '.join(['ReferencePDB','FASD','Judge'])+'\n')
     vote_all = 0
     for ReferencePDBCode in ReferenceSet:
         vote = 0
         logging.info('Compare the qurey binding sites to metal binding sites of %s'%ReferencePDBCode)
         d = CalculateDistance(PDBDir,QueryPDBCode,ResidueList,ReferencePDBCode,' '+metal)
         if d < t1: #if the distance is less than t, two binding sites are thought to be similar, here t = 0.3 
            vote = 1
            vote_all += 1
            
         f.write('  '.join([ReferencePDBCode,str(d),str(vote)])+'\n')
     
     Judge = 'FALSE'
     PassRate = vote_all/len(ReferenceSet)    
     if PassRate > t2:  # if more than half of structures in reference set say yes, output yes else no
         logging.info('Yes! %s is predicted to bind %s!'%(args.input,args.metal_type))
         Judge = 'TRUE'
     else: logging.info('No! %s is predicted do not bind %s!'%(args.input,args.metal_type))
     
     f.write('VoteAll: '+str(vote_all) + '  PassRate: '+str(PassRate) + '  Judge: '+Judge)
     
     f.close()






