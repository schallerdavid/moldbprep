from moldbprep.io import update_progress
import pandas as pd
from molvs.standardize import Standardizer, LargestFragmentChooser, Uncharger
from molvs.tautomer import TautomerCanonicalizer
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import time


def protonate_mol(mol):
    """
    This function protonates molecules based on substructure patterns.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        An RDKit molecule.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        A protonated RDKit molecule.

    """
    transformations = {
        'acids C=C-OH': '[H][O;$(OC=CC(=O))!$([-1])!$(OC[-1]):1][C:2]>>[C:2][O-;$(OC=CC(=O)):1]',
        'acids NH': '[H][NH1;$(N(C=O)S(=O)(=O)[#6])!$([-1])!$(N*[-1]):1]([C:3])[S:2]>>[C:3][N-;$(N(C(=O))S(=O)(=O)):1][S:2]',
        'acids OH': '[H][O;$(O[C,S,P]=O)!$([-1])!$(O[C,S][-1]):1][A;C,S,P:2]>>[O-;$(O[C,S,P]=O):1][A;C,S,P:2]',
        'aliphatic amines': '[#1,#6:4][N;!$([+1])!$(N=*)!$(N[+1])!$(N*([+1,#8,#15,#16,F,Cl]))!$(N*#*)!$(N**[+1])!$(N[#7,#8,#15,#16,c])!$(NC=[#7,#8,#16])!$(N#*)!$(N@*@*@O)!$(NC(F)F)!$(NCC(F)F)!$(NC=C[C,S,N]~[NH0,OH0]):1]([#1,#6:2])[#1,#6:3]>>[H][N+:1]([#1,#6:4])([#1,#6:2])[#1,#6:3]',
        'amidines': '[N:3][C:2]=[N;!$([+1])!$(N[#7,#8])!$(N=CN[#7,#8])!$(N=C([N,O])N)!$(NS(=O)(=O))!$(N=CNS(=O)(=O))!$(N(c)=CS)!$(N=C(S)Nc)!$(N=CN=*)!$(NC#N)!$(N=CNC#N)!$(N=CNC=[#8,#16])!$(N(C=*)=CNC=*)!$([N](=[CX3])[CX3]):1][#1,#6:4]>>[H][N+:1]([#1,#6:4])=[C:2][N:3]',
        'guanidines': '[N:4][C:2]([N:3])=[N;!$([+1])!$(NC[+1])!$(N=CN=*)!$(N(C(=O))=CNC(=O))!$(N=C(NC=O)NC=O)!$(N=CN*(~O)~O)!$(NC(=O)[CX3])!$(NC#N):1][#1,#6:5]>>[H][N+:1]([#1,#6:5])=[C:2]([N:4])[N:3]',
        'tetrazole1': '[c:1]1[n:2][n:3]([H])[n:4][n:5]1>>[c:1]1[n:2][n-:3][n:4][n:5]1',
        'tetrazole2': '[c:1]1[n:2][n:3][n:4][n:5]1([H])>>[c:1]1[n:2][n:3][n:4][n-:5]1'
    }
    mol = Uncharger().uncharge(mol)
    mol = Chem.AddHs(mol)
    for smarts in transformations.values():
        products = [0]
        rxn = ReactionFromSmarts(smarts)
        while len(products) > 0:
            products = rxn.RunReactants((mol,))
            if len(products) > 0:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(products[0][0]))
                mol = Chem.AddHs(mol)
    mol = Chem.RemoveHs(mol)
    return mol


def enumerate_stereo_isomers(mol, max_stereo_isomers):
    """
    This function emumerates stereo isomers.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        An RDKit molecule.

    max_stereo_isomers: int
        Maximal number of stereo isomers to generate.

    Returns
    -------
    isomers : tuple
        A tuple of enumerated RDKit molecules.

    """
    options = StereoEnumerationOptions(tryEmbedding=False, unique=True, maxIsomers=max_stereo_isomers)
    isomers = tuple(EnumerateStereoisomers(mol, options=options))
    return isomers


def standardize_mols(jobs, mol_counter, num_mols, results, start_time, vendors, max_stereo_isomers, failures,
                     verbose=False):
    """
    This function passes molecules to the standardization functions.

    Parameters
    ----------
    jobs: multiprocessing.manager.list
        A list containing job information as dictionaries.

    mol_counter: multiprocessing.manager.value
        A counter keeping track of processed molecules.

    num_mols: int
        Total number of molecules to be processed.

    results: multiprocessing.manager.list
        A list containing lists describing the processed molecules.

    start_time: float
        Starting time of molecule processing.

    vendors: list
        List of vendors.

    max_stereo_isomers: int
        Maximal number of stereo isomers to generater per molecule.

    verbose : bool
        If RDKit warning should be displayed.

    """
    if not verbose:
        RDLogger.DisableLog('rdApp.*')
    job = 'initiate'
    processed_mols = []
    while job is not None:
        try:
            job = jobs.pop(0)
            vendor_position = vendors.index(job['vendor'])
            supplier = Chem.SDMolSupplier(job['sdf_path'])
            for mol_id in range(job['mol_start'], job['mol_end'] + 1):
                mol = supplier[mol_id]
                identifier = mol.GetProp(job['identifier_field'])
                try:
                    # default standardization from molvs
                    mol = Standardizer().standardize(mol)
                    # choose largest fragment
                    mol = LargestFragmentChooser().choose(mol)
                    # canonicalize tautomer
                    mol = TautomerCanonicalizer().canonicalize(mol)
                    # protonate mol
                    mol = protonate_mol(mol)
                    # enumerate stereo isomers and append mols
                    for mol in enumerate_stereo_isomers(mol, max_stereo_isomers):
                        mol_as_list = [Chem.MolToSmiles(mol)] + [''] * len(vendors)
                        mol_as_list[1 + vendor_position] = identifier
                        processed_mols.append(mol_as_list)
                except:
                    failures.append(Chem.MolToSmiles(mol))
                with mol_counter.get_lock():
                    mol_counter.value += 1
                update_progress(mol_counter.value / num_mols, 'Progress of standardization',
                                ((time.time() - start_time) / mol_counter.value) * (num_mols - mol_counter.value))
        except IndexError:
            job = None
    results += processed_mols
    return


def merge_ids(results, vendors):
    """
    This function merges identifiers from vendors for the same molecule.

    Parameters
    ----------
    results: pandas.DataFrame
        A dataframe with columns smiles, vendor1, ..., vendorx.

    vendors: list
        A list containing vendor names matching the ones in results.

    Returns
    -------
    merged_results: pandas.DataFrame
        A dataframe with columns smiles, vendor1, ..., vendorx. Duplicates of smiles are removed by merging vendor
        identifiers with a comma.

    """
    grouped_per_vendor = []
    print('Merging molecules and identifiers for each vendor...')
    for vendor in vendors:
        print('Merging {} database...'.format(vendor))
        vendor_results = results[results[vendor] != '']
        joined_ids = vendor_results.groupby(['smiles'])[vendor].apply(','.join).to_frame().reset_index()
        other_vendors = pd.DataFrame([[''] * (len(vendors) - 1)] * joined_ids.shape[0], columns=[x for x in vendors
                                                                                                 if x != vendor])
        grouped_per_vendor.append(pd.concat([joined_ids, other_vendors], axis=1))
    grouped_per_vendor = pd.concat(grouped_per_vendor).reset_index(drop=True)
    grouped = []
    print('Merging molecules and identifiers into main database...')
    for vendor in vendors:
        print('Merging {} database...'.format(vendor))
        if len(grouped) == 0:
            grouped.append(grouped_per_vendor.groupby(['smiles'])[vendor].apply(','.join).str.strip(',').reset_index())
        else:
            grouped.append(grouped_per_vendor.groupby(['smiles'])[vendor].apply(','.join
                                                                                ).str.strip(',').reset_index(drop=True))
    merged_results = pd.concat(grouped, axis=1)
    return merged_results


