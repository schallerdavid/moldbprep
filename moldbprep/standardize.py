from moldbprep.io import update_progress
from molvs.charge import Uncharger
from molvs.fragment import LargestFragmentChooser
from molvs import Standardizer
from rdkit import Chem
from rdkit.Chem.AllChem import ReactionFromSmarts
import time


def largest_fragment(mol):
    """
    This function standardizes molecules.

    Parameters
    ----------
    mol - rdkit.Chem.rdchem.Mol
        An RDKit molecule

    Returns
    -------
    smiles - str
        A list containing SMILES, vendor and vendor ID.

    """
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(LargestFragmentChooser().choose(mol)))
    return mol


def protonate_mol(mol):
    """

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
        'aliphatic amines': '[#1,#6:4][N;!$([+1])!$(N=*)!$(N[+1])!$(N*([+1,#7,#8,#15,#16,F,Cl]))!$(N*#*)!$(N**[+1])!$(N[#7,#8,#15,#16,c])!$(NC=[#7,#8,#16])!$(N#*)!$(N@*@*@O)!$(NC(F)F)!$(NCC(F)(F)F)!$(NC=C[C,S,N]~[NH0,OH0]):1]([#1,#6:2])[#1,#6:3]>>[H][N+:1]([#1,#6:4])([#1,#6:2])[#1,#6:3]',
        'amidines': '[N:3][C:2]=[N;!$([+1])!$(N[#7,#8])!$(N=CN[#7,#8])!$(N=C([N,O])N)!$(NS(=O)(=O))!$(N=CNS(=O)(=O))!$(N(c)=CS)!$(N=C(S)Nc)!$(N=CN=*)!$(NC#N)!$(N=CNC#N)!$(N=CNC=[#8,#16])!$(N(C=*)=CNC=*)!$([N](=[CX3])[CX3]):1][#1,#6:4]>>[H][N+:1]([#1,#6:4])=[C:2][N:3]',
        'guanidines': '[N:4][C:2]([N:3])=[N;!$([+1])!$(NC[+1])!$(N=CN=*)!$(N(C(=O))=CNC(=O))!$(N=C(NC=O)NC=O)!$(N=CNN(~O)~O)!$(NC(=O)[CX3]):1][#1,#6:5]>>[H][N+:1]([#1,#6:5])=[C:2]([N:4])[N:3]',
        'tetrazoles': '[#7,#6:2][#7;r5;$(*@n@n)!$(*[#6][#6])!$([#7][-1])!$([#7][#6]=O)!$([#7]*[#6](=O)):1]([#7,#6:3])[H]>>[#7,#6:2][*-:1][#7,#6:3]'
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


def standardize_mols(jobs, mol_counter, num_mols, results, start_time):
    """
    This function passes molecules to the standardize_mol function.

    Parameters
    ----------
    file_path : str
        Full path to sdf file.

    first_mol : int
        Position of the first molecule in sdf-file to standardize.

    last_mol : int
        Position of the last molecule in sdf-file to standardize.

    """
    #
    #for mol_id in range(first_mol, last_mol):

    #    mol_standardized = standardize_mol(suppl[mol_id])
    #    if mol_counter is not None:
    #        with mol_counter.get_lock():
    #            mol_counter.value += 1
    job = 'initiate'
    processed_mols = []
    while job is not None:
        try:
            job = jobs.pop(0)
            supplier = Chem.SDMolSupplier(job['sdf_path'])
            for mol_id in range(job['mol_start'], job['mol_end'] + 1):
                mol = supplier[mol_id]
                identifier = mol.GetProp(job['identifier_field'])
                mol = Standardizer().standardize(mol)
                mol = largest_fragment(mol)
                mol = protonate_mol(mol)
                with mol_counter.get_lock():
                    mol_counter.value += 1
                    update_progress(mol_counter.value / num_mols, 'Progress of standardization',
                                    ((time.time() - start_time) / mol_counter.value) * (num_mols - mol_counter.value))
                processed_mols.append([Chem.MolToSmiles(mol), Chem.MolToInchiKey(mol), job['vendor'],
                                       identifier])
        except IndexError:
            job = None
    results += processed_mols
    return
