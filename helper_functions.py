import pandas as pd
import os
import glob
import matplotlib.patches as Rectangle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_48_hd_files(file_dir, ignore, trim=True):
    """
    Loads 48hd files as one large dataframe
    Args:
        file_dir (str): path to 48hd files
        ignore (bool): specifies if nuc sequences larger than 24 should be excluded

    Returns:

    """
    all_files = glob.glob(os.path.join(file_dir, "*.txt"))
    li_files = []

    for file_ in all_files:
        df = pd.read_csv(file_, index_col=0, header=0, sep=' ')
        df.drop(['mindex', 'Primer', 'Mod', 'AA'], axis=1, inplace=True)
        df.set_index('Nuc', inplace=True)
        print(f'Length of DF {file_} before: {len(df)}')
        print(len(df[df.index.duplicated(keep='first')]))
        df = df.groupby(by=df.index).sum()
        print(f'Len after {len(df)}')
        cols = df.columns[~df.columns.isin(['Nuc'])]
        df.rename(columns=dict(
            zip(cols, str(os.path.basename(file_)[:-4]) + '__' + cols))
            , inplace=True)
        li_files.append(df)

    df_full = pd.concat(li_files, axis=1)

    # Delete first row containing totals:
    first_row = df_full.loc[['XX']]
    df_full.drop(index=df_full.loc[['XX']].index.tolist(), inplace=True)

    if ignore:
        # Drop rows that have sequence length longer than 24:
        print(f'Length of DF before dropping rows: {len(df_full)}')
        seq_to_drop_idx = [x for x in df_full.index if len(x) > 24]
        df_full.drop(index=seq_to_drop_idx, inplace=True)
        print(f'Length of DF after dropping rows: {len(df_full)}')

    if trim:
        df_full.index = df_full.index.str[-24:-3]
        print(f'Final trimming for 48hd files')
        df_full = df_full.groupby(by=df_full.index).sum()
        print(f'Len after {len(df_full)}')

    df_full.fillna(0, inplace=True)
    return df_full, first_row


def load_vt_unfiltered_df(file_path, trim, keep):
    """Load vt_unfiltered_df files as one DF

    Args:
        file_path (str): Path to df_vt_unfiltered
        trim (bool): Trim nuc sequence to last 24 characters.
        keep (str): 'sum','first','last' Specifies how duplicates should be
        reconciled.

    Returns:
        [df]: Dataframe for VT_unfiltered
    """
    df = pd.read_fwf(file_path, colspecs='infer')
    # Drops the first row containing the totals
    first_row = df.iloc[:1, :]
    df = df[1:]
    df.rename(columns={'Tag': 'Nuc'}, inplace=True)
    df.set_index('Nuc', inplace=True)
    print(f'Length of DF_VT before: {len(df)}')
    df.drop(columns=['X'], inplace=True)
    df = df.apply(pd.to_numeric)
    df.fillna(0, inplace=True)
    df = df.groupby(by=df.index).sum()
    print(f'Len after {len(df)}')

    if trim:
        df.index = df.index.str[-24:-3]

        if keep == 'sum':
            print(f'Length of DF {file_path} before: {len(df)}')
            df = df.groupby(by=df.index).sum()
            print(f'Len after {len(df)}')

        else:
            print(f'Length of DF {file_path} before: {len(df)}')
            df = df.loc[~df.index.duplicated(keep=keep)]
            print(f'Len after {len(df)}')

    return df, first_row


def clean_48_zip_files(file_path):
    """

    Args:
        file_path:

    Returns:

    """
    all_files = glob.glob(os.path.join(file_path, "*.txt"))
    for file_ in all_files:
        with open(file_, 'r') as txt_file:
            lines = txt_file.readlines()

        with open(file_, 'w') as f:
            f.writelines(lines[16:-1])


def load_48_aa(file_dir, trim, reliable=3):
    all_Files = glob.glob(os.path.join(file_dir, "*.txt"))
    li = []

    for file_ in all_Files:
        df = pd.read_csv(file_, index_col=0, header=0, sep=' ')
        df.drop(['mindex','Primer','Mod', 'Nuc'], axis=1, inplace=True)
        df.set_index('AA', inplace=True)
        df = df.iloc[1:] # removes the first row that contains total number of AA for each file
        print(f'File Name: {os.path.basename(file_)}')

        # Trims the AA to the last 4 in sequence.
        if trim:
            df.index = df.index.str[-4:]

        # Sums up duplicates.
        print(f'Length of file before trimming and(or) combining: {len(df)}')
        df = df.groupby(by=df.index).sum()
        print(f'Length of file after trimming and(or) combining: {len(df)}')

        cols = df.columns[~df.columns.isin(['Nuc'])]

        # Only include a row if the sum of that row is > 3 -> remove rows that have very few values.
        # Based on matlab script!
        if reliable:
            len_b = len(df)
            df = df.loc[(df.sum(axis=1) > reliable)]
            print(f'Reliability check of {reliable} removed {len_b - len(df)} sequences')

        df.rename(columns = dict(
            zip(cols, str(os.path.basename(file_)[:-4])+ '__' + cols))
        , inplace = True)
        li.append(df)

    df_full = pd.concat(li, axis=1)

    df_full.fillna(0, inplace=True)
    df_full = df_full.drop(index='blank')
    return df_full

def load_vt_aa(file_path, trim, reliable=3):
    df = pd.read_fwf(file_path, colspecs='infer')
    df = df[1:] # remove first row containing totals

    df.drop(columns=['Tag'], inplace=True) # remove nucleotide column.
    df.rename(columns={'X':'AA'}, inplace=True)
    df.set_index('AA', inplace=True)

    df = df.apply(pd.to_numeric)
    df.fillna(0, inplace=True)

    if trim:
        df.index = df.index.str[18:22]

    print(f'Length before trimming and(or) combining {len(df)}')
    df = df.groupby(by=df.index).sum()
    print(f'Length after trimming and(or) combining {len(df)}')

    if reliable:
        len_b = len(df)
        df = df.loc[(df.sum(axis=1) > reliable)]
        print(f'Reliability check of {reliable} removed {len_b - len(df)} sequences')

    return df


def ttlocate(seqx):
    aa = ['R', 'K', 'Q', 'E', 'D', 'N', 'Y', 'P', 'T', 'S', 'H', 'A', 'G', 'W', 'M', 'F', 'L', 'V', 'I', 'C']
    row_i1 = aa.index(seqx[0])
    row_i2 = aa.index(seqx[2])
    col_i1 = aa.index(seqx[1])
    col_i2 = aa.index(seqx[3])
    row_i = (row_i1 * 20) + row_i2
    col_i = (col_i1 * 20) + col_i2
    return (row_i, col_i)


# function to get 20x20 matrix given array of aa sequences and DC label
def ttmatrix(seqs, dc):
    # define empty matrix to fill in
    ttmat = np.full((400, 400), np.nan)
    for i, s in enumerate(seqs):
        # ttmat[ttlocate(s)[0],ttlocate(s)[1]] = dc[i]
        ttmat[ttlocate(s)] = dc[i]
#     mask = np.isnan(ttmat)
    mask = np.where(ttmat < 1, True, False)
    return ttmat, mask


# plot the ttmatrix, and highlight the query if present
def ttplot(seqs, 
           dc, 
           query=None, 
           title='Position of queried sequence/s in 20x20 plot (log (deep conversion))',
           vmin=None,
           vmax=None):
    """Plots a 20x20 Matrix based on provided Sequences and LogDC values. The query is optional and draws a yellow
    square around the queried sequence
    Args:
        vmin:
        vmax:
        title:
        seqs ([String]): array of sequences to build the ttplot with
        dc ([Float]): array of log DC values corresponding description to the sequences.
        query (list, optional): List of query sequences to see where they might show up on the 20x20 Plot. Defaults to [].
    Returns:
        [type]: [description]
    """
    if query is None:
        query = []
    ttmat, mask = ttmatrix(seqs, dc)

    aa = ['R', 'K', 'Q', 'E', 'D', 'N', 'Y', 'P', 'T', 'S', 'H', 'A', 'G', 'W', 'M', 'F', 'L', 'V', 'I', 'C']
    ticks = [[a[0]] + [''] * 19 for a in aa]
    ticks = [j for i in ticks for j in i]
    sns.set(rc={'figure.figsize': (9, 6)})
    cm = 'jet'
    ax = sns.heatmap(ttmat, cmap=cm, xticklabels=ticks, yticklabels=ticks, mask=mask, vmin=0, vmax=200)
    ax.set_title(title, pad=20)
    ax.set(facecolor='#F5F5F5')

    if len(query) > 0:
        for q in query:
            query_pos = ttlocate(q)[1], ttlocate(q)[0]
            ax.add_patch(Rectangle(query_pos, 1, 1, fill=True, edgecolor='red', lw=8))
            

    return ax.figure