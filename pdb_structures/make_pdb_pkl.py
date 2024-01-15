# from utfils import *
import pandas as pd,numpy as np
from colorama import Fore
import glob,os,multiprocessing
from subprocess import check_output
import gemmi
from functools import partial
# import logging
# logging.basicConfig(filename='pdb.log',filemode='w', level=logging.DEBUG)

if not os.path.exists('dat'):os.mkdir('dat')

pkl_file='dat/pdb.pkl'
col_names={
    'sp':'space_group_hm','params':'cell_parameters','res':'resolution',
    'nb_chains':'number_of_chains','ligands':'ligands','nb_models':'nb_models',
    'nb_atoms' : 'Number of atoms','V':'Volume(A)',
    'exp':'exp_method','title':'title',
}

def process_entry(folder,i,nf):
    df = pd.DataFrame(columns=list(col_names.keys()))
    files=glob.glob('%s/*' %folder)

    folder_base=os.path.basename(folder)
    f=open('dat/pdb_%s.log' %folder_base,'w')
    msg='%d/%d : %s with %d files\n' %(i+1,nf,folder,len(files))    #;print(msg)
    f.write(msg)
    for filename in files:
        try:
            p = gemmi.read_structure(filename)
            pdb_code = p.info['_entry.id']
            df.loc[pdb_code] = '?'
            df.loc[pdb_code,'exp']      = p.info['_exptl.method']
            df.loc[pdb_code,'title']    = p.info['_struct.title']
            df.loc[pdb_code,'V']        = int(p.cell.volume)
            df.loc[pdb_code,'params']   = str(p.cell.parameters)
            df.loc[pdb_code,'sp']  = p.spacegroup_hm
            df.loc[pdb_code,'res'] = p.resolution
            df.loc[pdb_code,'nb_models'] = len(p)
            df.loc[pdb_code,'nb_chains'] = len(p[0])
            df.loc[pdb_code,'nb_atoms']  = p[0].count_atom_sites()
            df.loc[pdb_code,'ligands']   = ' '.join(np.unique([l.name for model in p[0] for l in model.get_ligands()]).tolist())

        except Exception as e:
            f.write('filename:%s\n' %filename)
            f.write(e.__str__()+'\n')

    if len(df):
        pkl_file='dat/pdb_%s.pkl' %folder_base
        df.to_pickle(pkl_file )
        f.write('file saved : %s' %pkl_file)
        print(Fore.GREEN+'file saved :\n\t'+Fore.YELLOW+pkl_file+Fore.RESET)
    f.close()

def make_pickle(parallel=True):
    all_folders=[f.split('/')[-1] for f in glob.glob('resources/*')]
    cmp_folders=[f.split('/')[-1].split('pdb_')[1].split('.')[0] for f in glob.glob('dat/pdb_*.pkl')]
    folders=np.setdiff1d(all_folders,cmp_folders)
    folders=['resources/%s' %s for s in folders]

    nf=len(folders)
    if parallel:
        p = multiprocessing.Pool()
        p.starmap(process_entry,zip(folders,list(range(nf)),[nf]*nf))

    else:
        for i,folder in enumerate(folders):
            process_entry(folder,i,nf)


    # check_output('cd dat;cat pdb_*.log>pdb.log;rm  pdb_*.log',shell=True)
    df=concat_df()
    return df

def concat_df():
    print('...concatenating all dataframes...')
    df=pd.concat([pd.read_pickle(pkl) for pkl in glob.glob('dat/*.pkl')])
    df=df.drop_duplicates()

    print_info = lambda df,col: print(df[col])
    df.quick_info = lambda :print_info(df,['sp','V','nb_atoms','nb_chains','ligands','res'])

    print('saving pickle')
    df.to_pickle(pkl_file)
    print(Fore.GREEN+'file saved :\n\t'+Fore.YELLOW+pkl_file+Fore.RESET)
    return df



if __name__=="__main__":
    # df = pd.read_pickle(pdbf.pkl_file)
    df = make_pickle(parallel=True)
    # df = concat_df()
    df.quick_info()
