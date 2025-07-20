# tetrahydra_molecules.py
# Part of the Metalalchemist Neurochemical Series — AngelNET Exclusive

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction
import os


def generate_molecule(smiles: str, name: str):
    """
    Creates an RDKit molecule object from a SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string for {name}: {smiles}")
    AllChem.Compute2DCoords(mol)
    return mol


def draw_molecules(molecule_data):
    """
    Takes a list of tuples (smiles, name) and renders the molecules side by side.
    """
    mols = [generate_molecule(smiles, name) for smiles, name in molecule_data]
    legends = [name for _, name in molecule_data]
    img = Draw.MolsToImage(mols, legends=legends)
    return img


def define_reaction(smarts: str):
    """
    Defines a chemical reaction from a SMARTS string.
    """
    reaction = ChemicalReaction()
    reaction = AllChem.ReactionFromSmarts(smarts)
    return reaction


def run_reaction(reactants, reaction_smarts):
    """
    Runs a chemical reaction on given reactants and returns the products.
    """
    reaction = define_reaction(reaction_smarts)
    ps = reaction.RunReactants(reactants)
    return ps


if __name__ == "__main__":
    # Molecule data: (SMILES, Name)
    molecule_data = [
        ("C1CCC2(CC1)CCCC2O", "Tetrahydrahypanol"),
        ("C1=CC2=C(C=C1)CC(C2)O", "Tetrahydrahypenol"),
    ]

    # Draw original conceptual molecules
    image = draw_molecules(molecule_data)
    image.save("tetrahydra_molecules.png")

    # === Tetrahydrahypanol Synthesis ===
    hypanone = generate_molecule("C1CCC2=CCCCC2C1=O", "Hypanone")
    hydrogenation_reaction = "[C:1]=[C:2]>>[C:1][C:2]"  # Mock hydrogenation
    reduction_reaction = "[C:1](=O)[C:2]>>[C:1](O)[C:2]"  # Mock ketone to alcohol

    hydrogenated = run_reaction((hypanone,), hydrogenation_reaction)
    reduced = run_reaction((hydrogenated[0][0],), reduction_reaction)

    intermediate_mols = [
        ("C1CCC2=CCCCC2C1=O", "Hypanone (Start)"),
        (Chem.MolToSmiles(hydrogenated[0][0]), "Hydrogenated Intermediate"),
        (Chem.MolToSmiles(reduced[0][0]), "Tetrahydrahypanol (Product)"),
    ]
    intermediate_img = draw_molecules(intermediate_mols)
    intermediate_img.save("tetrahydrahypanol_synthesis.png")

    # === Tetrahydrahypenol Synthesis ===
    hypenone = generate_molecule("C1=CC2=CCCCC2C1=O", "Hypenone")
    partial_hydrogenation = "[c:1]=[c:2]>>[C:1][C:2]"  # Mock partial hydrogenation (preserving aromaticity)
    enolization = "[C:1](=O)[C:2]>>[C:1](O)[C:2]"        # Enol tautomerization

    partially_hydrogenated = run_reaction((hypenone,), partial_hydrogenation)
    enolized = run_reaction((partially_hydrogenated[0][0],), enolization)

    hypenol_mols = [
        ("C1=CC2=CCCCC2C1=O", "Hypenone (Start)"),
        (Chem.MolToSmiles(partially_hydrogenated[0][0]), "Partially Hydrogenated"),
        (Chem.MolToSmiles(enolized[0][0]), "Tetrahydrahypenol (Product)"),
    ]
    hypenol_img = draw_molecules(hypenol_mols)
    hypenol_img.save("tetrahydrahypenol_synthesis.png")

    print("Full synthesis complete. All visualizations saved.")
