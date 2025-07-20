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


def draw_molecules(molecule_data, filename):
    """
    Takes a list of tuples (smiles or mol, name) and renders the molecules side by side.
    """
    mols = [mol if isinstance(mol, Chem.Mol) else generate_molecule(mol, name) for mol, name in molecule_data]
    legends = [name for _, name in molecule_data]
    img = Draw.MolsToImage(mols, legends=legends)
    img.save(filename)
    print(f"Saved visualization: {filename}")


def define_reaction(smarts: str):
    """
    Defines a chemical reaction from a SMARTS string.
    """
    return AllChem.ReactionFromSmarts(smarts)


def run_reaction(reactants, reaction_smarts):
    """
    Runs a chemical reaction on given reactants and returns the products.
    """
    reaction = define_reaction(reaction_smarts)
    ps = reaction.RunReactants(reactants)
    return ps


def run_forward_synthesis():
    """
    Runs the synthesis steps for Tetrahydrahypanol and Tetrahydrahypenol.
    """
    # Tetrahydrahypanol synthesis
    hypanone = generate_molecule("C1CCC2=CCCCC2C1=O", "Hypanone")
    hydrogenated = run_reaction((hypanone,), "[C:1]=[C:2]>>[C:1][C:2]")
    reduced = run_reaction((hydrogenated[0][0],), "[C:1](=O)[C:2]>>[C:1](O)[C:2]")

    draw_molecules([
        ("C1CCC2=CCCCC2C1=O", "Hypanone (Start)"),
        (hydrogenated[0][0], "Hydrogenated Intermediate"),
        (reduced[0][0], "Tetrahydrahypanol (Product)")
    ], "tetrahydrahypanol_synthesis.png")

    # Tetrahydrahypenol synthesis
    hypenone = generate_molecule("C1=CC2=CCCCC2C1=O", "Hypenone")
    partial_hydrogenated = run_reaction((hypenone,), "[c:1]=[c:2]>>[C:1][C:2]")
    enolized = run_reaction((partial_hydrogenated[0][0],), "[C:1](=O)[C:2]>>[C:1](O)[C:2]")

    draw_molecules([
        ("C1=CC2=CCCCC2C1=O", "Hypenone (Start)"),
        (partial_hydrogenated[0][0], "Partially Hydrogenated"),
        (enolized[0][0], "Tetrahydrahypenol (Product)")
    ], "tetrahydrahypenol_synthesis.png")


def run_retro_synthesis():
    """
    Simulates retrosynthetic steps for both molecules.
    """
    # Retro for Tetrahydrahypanol
    tetrahydrahypanol = generate_molecule("C1CCC2(CC1)CCCC2O", "Tetrahydrahypanol")
    oxidized = run_reaction((tetrahydrahypanol,), "[C:1](O)[C:2]>>[C:1](=O)[C:2]")
    desaturated = run_reaction((oxidized[0][0],), "[C:1][C:2]>>[C:1]=[C:2]")

    draw_molecules([
        ("C1CCC2(CC1)CCCC2O", "Tetrahydrahypanol (Start)"),
        (oxidized[0][0], "Oxidized Intermediate"),
        (desaturated[0][0], "Hypanone (Recovered)")
    ], "retro_tetrahydrahypanol.png")

    # Retro for Tetrahydrahypenol
    tetrahydrahypenol = generate_molecule("C1=CC2=C(C=C1)CC(C2)O", "Tetrahydrahypenol")
    keto_form = run_reaction((tetrahydrahypenol,), "[C:1](O)[C:2]>>[C:1](=O)[C:2]")
    re_aromatized = run_reaction((keto_form[0][0],), "[C:1][C:2]>>[c:1]=[c:2]")

    draw_molecules([
        ("C1=CC2=C(C=C1)CC(C2)O", "Tetrahydrahypenol (Start)"),
        (keto_form[0][0], "Keto Tautomer"),
        (re_aromatized[0][0], "Hypenone (Recovered)")
    ], "retro_tetrahydrahypenol.png")


if __name__ == "__main__":
    # Visualize original target molecules
    base_molecules = [
        ("C1CCC2(CC1)CCCC2O", "Tetrahydrahypanol"),
        ("C1=CC2=C(C=C1)CC(C2)O", "Tetrahydrahypenol")
    ]
    draw_molecules(base_molecules, "tetrahydra_molecules.png")

    # Synthesis flow
    run_forward_synthesis()

    # Retrosynthesis flow
    run_retro_synthesis()

    print("✅ All synthesis and retrosynthesis operations complete.")
