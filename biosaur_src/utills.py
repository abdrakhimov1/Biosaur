import logging
import warnings
warnings.formatwarning = lambda msg, *args, **kw: str(msg) + '\n'
logger = logging.getLogger(__name__)
from pyteomics import pepxml, achrom, auxiliary as aux, mass, fasta, mzid, parser
import pandas as pd


def add_mod_info(df_raw, mod):
    sequence, mods_counter = df_raw['peptide'], df_raw['mods_counter']
    mod_aa = mod.split(' at ')[1]
    if 'term' not in mod_aa and mod_aa not in sequence:
        return -1
    else:
        return (1 if mods_counter.get(mod, 0) >= 1 else 0)

def is_decoy(proteins, decoy_prefix, decoy_infix=False):
    if not decoy_infix:
        return all(z.startswith(decoy_prefix) for z in proteins)
    else:
        return all(decoy_infix in z for z in proteins)
    
def is_decoy_2(proteins, decoy_set):
    return all(z in decoy_set for z in proteins)

def prepare_mods(df):
    all_mods = set()
    for cnt in df['mods_counter'].values:
        for k in cnt.keys():
            all_mods.add(k)
    for mod in all_mods:
        df[mod] = df.apply(add_mod_info, axis=1, mod=mod)
        
def calc_RT(seq, RC):
    try:
        return achrom.calculate_RT(seq, RC)
    except Exception:
        return 0
    
def remove_column_hit_rank(df):
    if 'hit_rank' in df.columns:
        return df[df['hit_rank'] == 1]
    else:
        return df
        
def parse_mods(df_raw):
    mods_counter = {}
    sequence, mods = df_raw['peptide'], df_raw['modifications']
    if isinstance(mods, list):
        for mod in mods:
            mod_mass, aa_ind = mod.split('@')
            mod_mass = float(mod_mass)
            aa_ind = int(aa_ind)
            if aa_ind == 0:
                aa = 'N_term'
                mod_mass = round(mod_mass - 1.007825, 3)
            elif aa_ind == len(sequence) + 1:
                aa = 'C_term'
                mod_mass = round(mod_mass - 17.002735, 3)
            else:
                aa = sequence[aa_ind - 1]
                mod_mass = round(mod_mass - mass.std_aa_mass[aa], 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] = mods_counter.get(mod_name, 0) + 1
    return mods_counter

def parse_mods_msgf(df_raw):
    mods_counter = {}
    sequence, mods = df_raw['peptide'], df_raw['Modification']
    if isinstance(mods, list):
        for mod in mods:
            mod_mass, aa_ind = mod['monoisotopicMassDelta'], mod['location']
            if aa_ind == 0:
                aa = 'N_term'
                mod_mass = round(mod_mass, 3)
            elif aa_ind == len(sequence) + 1:
                aa = 'C_term'
                mod_mass = round(mod_mass, 3)
            else:
                aa = sequence[aa_ind - 1]
                mod_mass = round(mod_mass, 3)
            mod_name = 'mass shift %.3f at %s' % (mod_mass, aa)
            mods_counter[mod_name] = mods_counter.get(mod_name, 0) + 1
    return mods_counter

def prepare_dataframe(infile_path, decoy_prefix=None, decoy_infix=False, cleavage_rule=False, fdr=0.01, decoy2set=None):
    if not cleavage_rule:
        cleavage_rule = parser.expasy_rules['trypsin']
    if infile_path.lower().endswith('.pep.xml') or infile_path.lower().endswith('.pepxml'):
        df1 = pepxml.DataFrame(infile_path)
        ftype = 'pepxml'
    elif infile_path.lower().endswith('.mzid'):
        df1 = mzid.DataFrame(infile_path)
    else:
        raise WrongInputError()
    if not df1.shape[0]:
        raise EmptyFileError()

    if 'Morpheus Score' in df1.columns:
        df1 = df1[df1['Morpheus Score'] != 0]
        df1['expect'] = 1 / df1['Morpheus Score']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))

    if 'MS-GF:EValue' in df1.columns:
        # MSGF search engine
        ftype = 'msgf'
        df1['peptide'] = df1['PeptideSequence']
        df1['num_missed_cleavages'] = df1['peptide'].apply(lambda x: parser.num_sites(x, rule=cleavage_rule))
        df1['assumed_charge'] = df1['chargeState']
        df1['spectrum'] = df1['spectrumID']
        df1['massdiff'] = (df1['experimentalMassToCharge'] - df1['calculatedMassToCharge']) * df1['assumed_charge']
        df1['calc_neutral_pep_mass'] = df1['calculatedMassToCharge'] * df1['chargeState'] - df1['chargeState'] * 1.00727649
        df1['protein'] = df1['accession']
        df1['protein_descr'] = df1['protein description']
        df1['expect'] = df1['MS-GF:EValue']

    if set(df1['protein_descr'].str[0]) == {None}:
        # MSFragger
        logger.debug('Adapting MSFragger DataFrame.')
        logger.debug('Proteins before: %s', df1.loc[1, 'protein'])
        protein = df1['protein'].apply(lambda row: [x.split(None, 1) for x in row])
        df1['protein'] = protein.apply(lambda row: [x[0] for x in row])
        try:
            df1['protein_descr'] = protein.apply(lambda row: [x[1] for x in row])
        except IndexError:
            df1['protein_descr'] = protein.apply(lambda row: ['' for x in row])
        logger.debug('Proteins after: %s', df1.loc[1, 'protein'])

    df1.loc[pd.isna(df1['protein_descr']), 'protein_descr'] = df1.loc[pd.isna(df1['protein_descr']), 'protein']

    df1 = df1[~pd.isna(df1['peptide'])]
    if 'MS1Intensity' not in df1:
        df1['MS1Intensity'] = 0.0
    df1['length'] = df1['peptide'].apply(len)
    df1 = df1[df1['length'] >= 6]
    df1['spectrum'] = df1['spectrum'].apply(lambda x: x.split(' RTINS')[0])
    if 'retention_time_sec' not in df1.columns:
        if 'scan start time' in df1.columns:
            df1['RT exp'] = df1['scan start time']
            df1 = df1.drop(['scan start time', ], axis=1)
        else:
            df1['RT exp'] = 0
    else:
        df1['RT exp'] = df1['retention_time_sec'] / 60
        df1 = df1.drop(['retention_time_sec', ], axis=1)

    df1['massdiff_int'] = df1['massdiff'].apply(lambda x: int(round(x, 0)))
    df1['massdiff_ppm'] = 1e6 * (df1['massdiff'] - df1['massdiff_int'] * 1.003354) / df1['calc_neutral_pep_mass']

    df1 = remove_column_hit_rank(df1)

    if ftype == 'pepxml':
        df1['mods_counter'] = df1.apply(parse_mods, axis=1)
    elif ftype == 'msgf':
        df1['mods_counter'] = df1.apply(parse_mods_msgf, axis=1)
    prepare_mods(df1)




    try:
        logger.info('Calibrating retention model...')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

        
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, retention_coefficients))

        logger.info('RT model training results: R^2 = %f , std = %f', r_value**2, std_value)
        df1['RT diff'] = df1['RT pred'] - df1['RT exp']
        logger.info('Retention model calibrated successfully.')
    except Exception:
        logger.warning('Retention times are probably missing in input file.')
        df1['RT pred'] = df1['peptide'].apply(lambda x: calc_RT(x, achrom.RCs_krokhin_100A_tfa))
        df1['RT diff'] = df1['RT exp']
    return df1, decoy2set