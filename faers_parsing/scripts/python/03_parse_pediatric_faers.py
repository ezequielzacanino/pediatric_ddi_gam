#!/usr/bin/env python
# coding: utf-8

# In[1]:


import glob
import numpy as np
import pandas as pd
import pickle
import os
import gzip
import gc

script_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))
data_root = os.path.join(project_dir, "data", "raw")
data_dir = os.path.join(data_root, "openFDA_drug_event") + os.sep
er_dir = os.path.join(data_dir, 'er_tables') + os.sep

vocab_path = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
vocab_dir = os.path.join(
    vocab_path,
    'data',
    'vocabulary',
    'vocabulary_SNOMED_MEDDRA_RxNorm_ATC'
)


# In[2]:


primarykey='safetyreportid'


def load_atc_5th_concepts(vocab_dir, chunksize=500000):
    # 1- Reads vocabulary by chunks.
    atc_chunks = []
    concept_path = os.path.join(vocab_dir, 'CONCEPT.csv')
    concept_cols = ['concept_id', 'concept_code', 'concept_name',
                    'concept_class_id', 'vocabulary_id']

    for chunk in pd.read_csv(
        concept_path,
        sep='\t',
        usecols=concept_cols,
        dtype={
            'concept_id': 'Int64',
            'concept_code': 'str',
            'concept_name': 'str',
            'concept_class_id': 'str',
            'vocabulary_id': 'str'
        },
        chunksize=chunksize
    ):
        filtered = chunk[
            (chunk['vocabulary_id'] == 'ATC') &
            (chunk['concept_class_id'] == 'ATC 5th')
        ].copy()

        if not filtered.empty:
            atc_chunks.append(
                filtered.loc[:, ['concept_id', 'concept_code',
                                 'concept_name', 'concept_class_id']]
            )

    atc_concepts = (
        pd.concat(atc_chunks, ignore_index=True)
        .drop_duplicates()
        .rename(columns={
            'concept_id': 'ATC_concept_id',
            'concept_code': 'ATC_concept_code',
            'concept_name': 'ATC_concept_name',
            'concept_class_id': 'ATC_concept_class_id'
        })
    )
    atc_concepts['ATC_concept_id'] = atc_concepts['ATC_concept_id'].astype(int)
    return atc_concepts


def build_complete_rxnorm_to_atc_map(standard_drugs_path, vocab_dir,
                                     chunksize=500000):
    # 1- Identify only RxNorm used by openFDA.
    openfda_rxnorm_ids = set()

    for chunk in pd.read_csv(
        standard_drugs_path,
        compression='gzip',
        usecols=['RxNorm_concept_id'],
        dtype={'RxNorm_concept_id': 'Int64'},
        chunksize=chunksize
    ):
        ids = chunk['RxNorm_concept_id'].dropna().astype(int).unique()
        openfda_rxnorm_ids.update(ids)

    print('RxNorm openFDA únicos:', len(openfda_rxnorm_ids))

    # 2- Loads ATC 5th concepts.
    atc_concepts = load_atc_5th_concepts(vocab_dir=vocab_dir,
                                         chunksize=chunksize)
    atc_concept_ids = set(atc_concepts['ATC_concept_id'].astype(int).unique())

    # 3- Extracts RxNorm -> ATC relationships in OMOP.
    direct_relation_chunks = []
    relationship_path = os.path.join(vocab_dir, 'CONCEPT_RELATIONSHIP.csv')

    for chunk in pd.read_csv(
        relationship_path,
        sep='\t',
        usecols=['concept_id_1', 'concept_id_2', 'relationship_id'],
        dtype={
            'concept_id_1': 'Int64',
            'concept_id_2': 'Int64',
            'relationship_id': 'str'
        },
        chunksize=chunksize
    ):
        filtered = chunk[
            (chunk['relationship_id'] == 'RxNorm - ATC') &
            (chunk['concept_id_2'].isin(atc_concept_ids))
        ].copy()

        if not filtered.empty:
            direct_relation_chunks.append(
                filtered.loc[:, ['concept_id_1', 'concept_id_2']]
            )

    direct_relations = (
        pd.concat(direct_relation_chunks, ignore_index=True)
        .drop_duplicates()
    )
    direct_relations['concept_id_1'] = direct_relations['concept_id_1'].astype(int)
    direct_relations['concept_id_2'] = direct_relations['concept_id_2'].astype(int)

    print('Relaciones directas RxNorm -> ATC:', direct_relations.shape[0])

    # 4- Keeps openFDA concepts with direct relationship.
    descendant_map_chunks = []
    direct_openfda = direct_relations[
        direct_relations['concept_id_1'].isin(openfda_rxnorm_ids)
    ].copy()

    if not direct_openfda.empty:
        descendant_map_chunks.append(
            direct_openfda.rename(columns={
                'concept_id_1': 'RxNorm_concept_id',
                'concept_id_2': 'ATC_concept_id'
            }).loc[:, ['RxNorm_concept_id', 'ATC_concept_id']]
        )

    # 5- Expands by CONCEPT_ANCESTOR
    direct_rxnorm_ids = set(direct_relations['concept_id_1'].astype(int).unique())
    ancestor_path = os.path.join(vocab_dir, 'CONCEPT_ANCESTOR.csv')

    for chunk in pd.read_csv(
        ancestor_path,
        sep='\t',
        usecols=['ancestor_concept_id', 'descendant_concept_id'],
        dtype={
            'ancestor_concept_id': 'Int64',
            'descendant_concept_id': 'Int64'
        },
        chunksize=chunksize
    ):
        filtered = chunk[
            (chunk['descendant_concept_id'].isin(openfda_rxnorm_ids)) &
            (chunk['ancestor_concept_id'].isin(direct_rxnorm_ids))
        ].copy()

        if filtered.empty:
            continue

        expanded = (
            filtered.merge(
                direct_relations,
                left_on='ancestor_concept_id',
                right_on='concept_id_1',
                how='inner'
            )
            .loc[:, ['descendant_concept_id', 'concept_id_2']]
            .rename(columns={
                'descendant_concept_id': 'RxNorm_concept_id',
                'concept_id_2': 'ATC_concept_id'
            })
            .drop_duplicates()
        )

        if not expanded.empty:
            descendant_map_chunks.append(expanded)

    rxnorm_to_atc_map = (
        pd.concat(descendant_map_chunks, ignore_index=True)
        .drop_duplicates()
        .merge(atc_concepts, on='ATC_concept_id', how='left')
        .dropna(subset=['ATC_concept_id'])
        .drop_duplicates()
    )

    rxnorm_to_atc_map['RxNorm_concept_id'] = (
        rxnorm_to_atc_map['RxNorm_concept_id'].astype(int)
    )
    rxnorm_to_atc_map['ATC_concept_id'] = (
        rxnorm_to_atc_map['ATC_concept_id'].astype(int)
    )

    print('RxNorm openFDA con mapeo ATC final:',
          rxnorm_to_atc_map['RxNorm_concept_id'].nunique())

    return rxnorm_to_atc_map


def write_complete_pediatric_drugs_reactions(
    pediatric_reporter_df,
    ped_reports,
    standard_drugs_path,
    standard_reactions_path,
    rxnorm_to_atc_map,
    output_path,
    chunksize=250000
):
    # Preparing the reaction table once maintains the join by chunks, avoiding re-reading full file
    reactions = (
        pd.read_csv(
            standard_reactions_path,
            compression='gzip',
            dtype={'safetyreportid': 'str'}
        )
        .query('safetyreportid in @ped_reports')
        .dropna(subset=['safetyreportid', 'MedDRA_concept_id'])
        .drop_duplicates()
    )
    reactions['safetyreportid'] = reactions['safetyreportid'].astype(str)
    reactions['MedDRA_concept_id'] = reactions['MedDRA_concept_id'].astype(int)
    reactions['MedDRA_concept_code'] = reactions['MedDRA_concept_code'].astype(int)

    reporter = pediatric_reporter_df.copy()
    reporter['safetyreportid'] = reporter['safetyreportid'].astype(str)
    reporter = reporter.drop_duplicates(subset=['safetyreportid'])

    rxnorm_to_atc_map = rxnorm_to_atc_map.drop_duplicates().copy()
    rxnorm_to_atc_map['RxNorm_concept_id'] = (
        rxnorm_to_atc_map['RxNorm_concept_id'].astype(int)
    )

    written_rows = 0
    written_reports = set()
    chunk_number = 0
    first_write = True

    with gzip.open(output_path, mode='wt', newline='') as gz_file:
        for chunk in pd.read_csv(
            standard_drugs_path,
            compression='gzip',
            usecols=['safetyreportid', 'RxNorm_concept_id'],
            dtype={
                'safetyreportid': 'str',
                'RxNorm_concept_id': 'Int64'
            },
            chunksize=chunksize
        ):
            chunk_number += 1

            filtered_drugs = (
                chunk
                .dropna(subset=['safetyreportid', 'RxNorm_concept_id'])
                .query('safetyreportid in @ped_reports')
                .copy()
            )

            if filtered_drugs.empty:
                continue

            filtered_drugs['RxNorm_concept_id'] = (
                filtered_drugs['RxNorm_concept_id'].astype(int)
            )

            pediatric_drugs_atc = (
                filtered_drugs
                .drop_duplicates()
                .merge(rxnorm_to_atc_map, on='RxNorm_concept_id', how='inner')
                .drop_duplicates()
            )

            if pediatric_drugs_atc.empty:
                continue

            merged = (
                reporter
                .merge(pediatric_drugs_atc, on='safetyreportid', how='inner')
                .merge(reactions, on='safetyreportid', how='inner')
                .drop_duplicates()
            )

            if merged.empty:
                continue

            merged = merged.reindex(np.sort(merged.columns), axis=1)
            merged.to_csv(gz_file, index=False, header=first_write)
            first_write = False

            written_rows += merged.shape[0]
            written_reports.update(merged['safetyreportid'].astype(str).unique())

            print(
                'Chunk', chunk_number,
                '- filas acumuladas:', written_rows,
                '- reportes acumulados:', len(written_reports)
            )

            del filtered_drugs
            del pediatric_drugs_atc
            del merged
            gc.collect()

    print('Filas finales escritas:', written_rows)
    print('Reportes finales escritos:', len(written_reports))


# In[3]:


patients = pd.read_csv(er_dir+'patient.csv.gz',
                       compression='gzip',
                       dtype={
                           'safetyreportid' : 'str',
                           'patient_custom_master_age' : 'float'
                       })


# In[4]:


age_col='patient_custom_master_age'
aged = patients[patients[age_col].notnull()].reset_index(drop=True).copy()


# In[5]:


col = 'nichd'

neonate = aged[age_col].apply(lambda x : float(x)>0 and float(x)<=(1/12))
infant = aged[age_col].apply(lambda x : float(x)>(1/12) and float(x)<=1)
toddler = aged[age_col].apply(lambda x : float(x)>1 and float(x)<=2)
echildhood = aged[age_col].apply(lambda x : float(x)>2 and float(x)<=5)
mchildhood = aged[age_col].apply(lambda x : float(x)>5 and float(x)<=11)
eadolescence = aged[age_col].apply(lambda x : float(x)>11 and float(x)<=18)
ladolescence = aged[age_col].apply(lambda x : float(x)>18 and float(x)<=21)

aged[col] = None

aged.loc[neonate,col] = 'term_neonatal'
aged.loc[infant,col] = 'infancy'
aged.loc[toddler,col] = 'toddler'
aged.loc[echildhood,col] = 'early_childhood'
aged.loc[mchildhood,col] = 'middle_childhood'
aged.loc[eadolescence,col] = 'early_adolescence'
aged.loc[ladolescence,col] = 'late_adolescence'


# In[6]:


col = 'ich_ema'

term_newborn_infants = (aged[age_col].
                        apply(lambda x : float(x)>0 and float(x)<=(1/12)))
infants_and_toddlers = (aged[age_col].
                       apply(lambda x : float(x)>(1/12) and float(x)<=2))
children = aged[age_col].apply(lambda x : float(x)>2 and float(x)<=11)
adolescents = aged[age_col].apply(lambda x : float(x)>11 and float(x)<=17)

aged[col] = np.nan

aged.loc[term_newborn_infants,col] = 'term_newborn_infants'
aged.loc[infants_and_toddlers,col] = 'infants_and_toddlers'
aged.loc[children,col] = 'children'
aged.loc[adolescents,col] = 'adolescents'


# In[7]:


col = 'fda'

neonates = (aged[age_col].
                        apply(lambda x : float(x)>0 and float(x)<(1/12)))
infants = (aged[age_col].
                       apply(lambda x : float(x)>=(1/12) and float(x)<2))
children = aged[age_col].apply(lambda x : float(x)>=2 and float(x)<11)
adolescents = aged[age_col].apply(lambda x : float(x)>=11 and float(x)<16)

aged[col] = np.nan

aged.loc[neonates,col] = 'neonates'
aged.loc[infants,col] = 'infants'
aged.loc[children,col] = 'children'
aged.loc[adolescents,col] = 'adolescents'


# In[8]:


pediatric_patients = (aged.
                      dropna(subset=['nichd']).
                      reset_index(drop=True))
print(pediatric_patients.shape)
print(pediatric_patients.head())


# In[9]:


del patients
del aged


# In[10]:


pediatric_patients.head()


# In[11]:


report = (pd.read_csv(er_dir+'report.csv.gz',
                      compression='gzip',
                     dtype={
                         'safetyreportid' : 'str'
                     }))
report.head()


# In[12]:


df1 = pediatric_patients.copy()
ped_reports = df1.safetyreportid.unique()
df2 = report.copy()
print(df1.shape)
print(df2.shape)
df1[primarykey] = df1[primarykey].astype(str)
df2[primarykey] = df2[primarykey].astype(str)
pediatric_patients_report = pd.merge(df1,
         df2,
         on=primarykey,
         how='inner').query('safetyreportid in @ped_reports')
print(pediatric_patients_report.shape)


# In[13]:


del pediatric_patients
del report


# In[14]:


report_serious = pd.read_csv(er_dir+'report_serious.csv.gz',compression='gzip')
report_serious.head()


# In[15]:


df1 = pediatric_patients_report.copy()
df2 = report_serious.copy()
print(df1.shape)
print(df2.shape)
df1[primarykey] = df1[primarykey].astype(str)
df2[primarykey] = df2[primarykey].astype(str)
pediatric_patients_report_serious = pd.merge(df1,
         df2,
         on=primarykey,
         how='inner')
print(pediatric_patients_report_serious.shape)


# In[16]:


pediatric_patients_report_serious.head()


# In[17]:


del report_serious
del pediatric_patients_report


# In[18]:


reporter = pd.read_csv(er_dir+'reporter.csv.gz',compression='gzip')
reporter.head()


# In[19]:


df1 = pediatric_patients_report_serious.copy()
df2 = reporter.copy()
print(df1.shape)
print(df2.shape)
df1[primarykey] = df1[primarykey].astype(str)
df2[primarykey] = df2[primarykey].astype(str)
pediatric_patients_report_serious_reporter = pd.merge(df1,
         df2,
         on=primarykey,
         how='inner')
print(pediatric_patients_report_serious_reporter.shape)


# In[20]:


pediatric_patients_report_serious_reporter.head()


# In[21]:


pediatric_patients_report_serious_reporter.info()


# In[22]:


del reporter


# In[23]:


del pediatric_patients_report_serious


# In[24]:


(pediatric_patients_report_serious_reporter.
 to_csv(os.path.join(data_root, 'pediatric_patients_report_serious_reporter.csv.gz'),
       compression='gzip')
)


# In[25]:


ped_reports = pediatric_patients_report_serious_reporter.safetyreportid.astype(str).unique()
len(ped_reports)


# In[26]:


pediatric_patients_report_serious_reporter = (pd.
 read_csv(os.path.join(data_root, 'pediatric_patients_report_serious_reporter.csv.gz'),
       compression='gzip',
         index_col=0)
)
pediatric_patients_report_serious_reporter.head()


# In[27]:

standard_drugs_path = os.path.join(er_dir, 'standard_drugs.csv.gz')
standard_reactions_path = os.path.join(er_dir, 'standard_reactions.csv.gz')
final_output_path = os.path.join(
    data_root,
    'pediatric_patients_report_serious_reporter_drugs_reactions.csv.gz'
)

rxnorm_to_atc_map = build_complete_rxnorm_to_atc_map(
    standard_drugs_path=standard_drugs_path,
    vocab_dir=vocab_dir
)
print(rxnorm_to_atc_map.head())

write_complete_pediatric_drugs_reactions(
    pediatric_reporter_df=pediatric_patients_report_serious_reporter,
    ped_reports=ped_reports,
    standard_drugs_path=standard_drugs_path,
    standard_reactions_path=standard_reactions_path,
    rxnorm_to_atc_map=rxnorm_to_atc_map,
    output_path=final_output_path
)


# In[28]:


del rxnorm_to_atc_map
gc.collect()


# In[33]:


del pediatric_patients_report_serious_reporter


# In[34]:


pediatric_standard_drugs = (pd.
                            read_csv(os.path.join(er_dir, 'standard_drugs.csv.gz'),
                                     compression='gzip',
                                    dtype={
                                        'safetyreportid' : 'str'
                                    }).
                            query('safetyreportid in @ped_reports')
                           )
pediatric_standard_drugs.safetyreportid = pediatric_standard_drugs.safetyreportid.astype(str) 
pediatric_standard_drugs.RxNorm_concept_id = pediatric_standard_drugs.RxNorm_concept_id.astype(int)
pediatric_standard_drugs.head()


# In[35]:


concept = pd.read_csv(
    os.path.join(vocab_dir, 'CONCEPT.csv'),
    sep='\t',
    dtype={
        'concept_id': 'Int64'
    }
)
concept_relationship = pd.read_csv(
    os.path.join(vocab_dir, 'CONCEPT_RELATIONSHIP.csv'),
    sep='\t',
    dtype={
        'concept_id_1': 'Int64',
        'concept_id_2': 'Int64'
    }
)

brand_concepts = (concept.
                  query('vocabulary_id=="RxNorm" & concept_class_id=="Brand Name"').
                  loc[:, ['concept_id', 'concept_code', 'concept_name', 'concept_class_id']].
                  drop_duplicates().
                  rename(columns={
                      'concept_id': 'brand_concept_id',
                      'concept_code': 'brand_concept_code',
                      'concept_name': 'brand_concept_name',
                      'concept_class_id': 'brand_concept_class_id'
                  }))
tobrand = (concept_relationship.
           loc[:, ['concept_id_1', 'concept_id_2']].
           drop_duplicates().
           merge(
               brand_concepts,
               left_on='concept_id_2',
               right_on='brand_concept_id',
               how='inner'
           ))


# In[37]:


a = pediatric_standard_drugs.copy()
print(a[primarykey].nunique())
m = (pd.merge(
    a,
    tobrand,
    left_on='RxNorm_concept_id',
    right_on='concept_id_1'
)
)
m[primarykey].nunique()


# In[38]:


m_renamed = (m.
 loc[:,
     [primarykey, 'brand_concept_class_id', 'brand_concept_code',
      'brand_concept_name', 'brand_concept_id']
    ].
 rename(columns={
     'brand_concept_class_id' : 'RxNorm_concept_class_id',
     'brand_concept_code' : 'RxNorm_concept_code',
     'brand_concept_name' : 'RxNorm_concept_name',
     'brand_concept_id' : 'RxNorm_concept_id'})
)


# In[39]:


(m_renamed.
 to_csv(os.path.join(data_root, 'pediatric_patients_report_drug_brands.csv.gz'),
       compression='gzip')
)


# In[ ]:
