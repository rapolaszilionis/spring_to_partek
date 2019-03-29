#########################
### import statements ###
#########################


from sys import argv
import os
import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import scipy.io
import io
import gzip


#################
### functions ###
#################


def text_to_sparse_in_chunks(
    path,
    sep = ',',
    chunksize = 100,
    skiprows = 1,
    skipcols = 1,
    compressed = True,
    save_skipped_rows = True,
    save_skipped_columns = True,
    comment = '#',
    verbose = True,
    ):
    
    """ for reading and simultaneously sparsifying giant csv/tsv of count data.
    input:
        path - path to a counts table.
        sep - separator
        chunksize - how many lines to load at a time before sparsifying
        skiprows - nr of rows to skip
        skipcols - nr of columns to skip
        compressed - whether gzip or not
        save_skipped_rows - if True, will return skipped rows as a dictionary
                            of the form {rows_number: line_as_string}
        
        save_skipped_columns - if True, will return skipped columns as a dictionary
                            of the form {rows_number: [col_1,col_2...col_skipcols]}
                            only for those rows that were not skipped
        
    
    """
    
    if compressed:
        
        f = io.TextIOWrapper(io.BufferedReader(gzip.open(path)))
        
    else:
        f = open(path,"r")
    
    skipped_rows = {}
    skipped_columns = {}
    
    
    counter = 0
    chunks = []
    frame = []

    for line in f:
        counter += 1
        if (counter <= skiprows)|(line.startswith(comment)):
            if verbose:
                print("skipping row starting with:",line[:25])
            if save_skipped_rows:
                #add line to dictionary
                skipped_rows[counter-1] = line
            continue

        l = line.strip('\n').split(sep)

        # save skipped columns, but only for rows that are not skipped.
        skipped_columns[counter-1] = l[:skipcols]

        frame.append(l[skipcols:])
        if float(counter/chunksize) == int(counter/chunksize):
            if verbose:
                print(counter)
            frame = np.array(frame).astype(np.float)
            frame = scipy.sparse.csc_matrix(frame)
            chunks.append(frame)

            # reset frame
            del frame
            frame = []
    
    # in case the total number of lines is a multiple of 
    if not (float(counter/chunksize) == int(counter/chunksize)):
        print(counter)
        frame = np.array(frame).astype(np.float)
        frame = scipy.sparse.csc_matrix(frame)
        chunks.append(frame)
        
        # reset frame
        del frame
        frame = []
    
    f.close()
    
    print("concatenating chunks...")
    E = scipy.sparse.vstack(chunks)
    print("turning into a csc matrix...")
    E = E.tocsc()
    print("done")

    return {'E':E,'skipped_rows':skipped_rows,'skipped_columns':skipped_columns}

###########################
### things happen below ###
###########################


# directory with files downloaded from spring
sdl = argv[1]

# add slash if necessary
if sdl[-1]!='/':
    sdl = sdl+'/'
    
# make a directory for saving
outdir = 'for_partek/'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
# load expression data
epath = sdl+'expr.csv.gz'
Edict = text_to_sparse_in_chunks(epath,chunksize=5e3,skipcols=1,skiprows=0)



# sparse expression data matrix
E = Edict['E']


# get gene names
genes = Edict['skipped_columns']
genes = np.array([genes[i] for i in np.arange(len(genes))]).T[0]


# load original cell index
oix = np.loadtxt(sdl+'original_cell_indices.txt',dtype=int)

# fake some barcode names, these are not saved during spring
barcodes = ['cell_%d'%i for i in oix]


# save 10x-like barcodes:
np.savetxt(outdir+'barcodes.tsv',barcodes,fmt='%s')


# save 10x-like genes:
feat = pd.DataFrame([
    np.repeat('its_a_hack',len(genes)), #hack!
    genes,
    np.repeat('Gene Expression',len(genes)) #hack!
])
feat.T.to_csv(outdir+'features.tsv',sep='\t',header=None,index=None)


# save mtx
scipy.io.mmwrite(outdir+'matrix.mtx',E)

# Save attributes for partek
# load categorcial coloring:
cg = pd.read_csv(sdl+'cell_groupings.csv',header=None,index_col=0).T
cg.index = barcodes


# load numeric coloring
num = pd.read_csv(sdl+'custom_colors.csv',header=None,index_col=0).T
num.index = barcodes
num.columns = ['numeric_%s'%i for i in num.columns]

# load coordinates
coo = pd.read_csv(sdl+'coordinates.csv',header=None,index_col=0)
coo.index = barcodes
coo.columns = ['x','y']
coo['y'] = -coo['y']
coo.columns = ['numeric_%s'%i for i in coo.columns]


# concatenate:
attr = pd.concat([num,cg,coo],axis=1)

# add index name
attr.index.name = 'Cell name'


# save
attr.to_csv(outdir+'attributes.txt',sep='\t')
