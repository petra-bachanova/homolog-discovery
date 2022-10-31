from datetime import datetime
import pandas as pd
import xmltodict


def read_xml_to_df(path):
    """parse the xml database file, save as a dictionary"""
    with open(path, encoding="utf-8") as xml_file:
        data_dict = xmltodict.parse(xml_file.read())
    data = pd.DataFrame.from_dict(data_dict["conoserver"]["entry"])
    return data


today = datetime.now().strftime('%y%m%d')
# imports
protein = read_xml_to_df("data/conoserver/conoserver_protein.xml")
structure = read_xml_to_df("data/conoserver/conoserver_structure.xml")
# convert database dictionaries into "PDB" and "BMRB" columns
database_list = []
for row in structure.linkOut:
    if type(row) is list:
        concise_dict = {}
        for index in range(len(row)):
            concise_dict[row[index]['@database']] = row[index]['#text']
        database_list.append(concise_dict)
    elif type(row) is float:
        database_list.append({})
    else:
        database_list.append({row['@database']: row['#text']})
structure = structure.join(pd.DataFrame(database_list))
structure.PDB = structure.PDB.astype('string')
structure = structure.dropna(subset=['PDB'])
structure = structure.drop_duplicates(subset=['proteinId'])

# check how many PDB entries there are
len([row for row in structure.PDB if type(row) == str])  # 198
len([row for row in structure.BMRB if type(row) == str])  # 110

cono_data = pd.merge(protein, structure[["proteinId", "PDB"]], left_on="id", right_on="proteinId", how="left")

# get PDB and UniProt identifiers from the linkOut dictionary column
database_list = []
for row in cono_data['linkOut']:
    if type(row) is list:
        concise_dict = {}
        for index in range(len(row)):
            concise_dict[row[index]['@database']] = row[index]['#text']
        database_list.append(concise_dict)
    elif type(row) is dict:
        database_list.append({row['@database']: row['#text']})
    else:
        database_list.append({})
cono_data = cono_data.join(pd.DataFrame(database_list))
cono_data.UniProt = cono_data['UniProt'].astype('string')

# conopeptide P05971 had its sequence in the 'note' column
cono_data.loc[cono_data['id'] == 'P05971', 'sequence'] = cono_data.loc[cono_data['id'] == 'P05971']['note']
# add a sequence length column
cono_data["sequenceLength"] = [len(row) for row in cono_data['sequence']]
# export
cono_data.to_csv(f"data/{today}_conoserver_data.csv", index=False)
